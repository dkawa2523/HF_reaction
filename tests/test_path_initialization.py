from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest

from hfauto.backends.ts.nwchem_neb import NWChemNEBEngine
from hfauto.backends.ts.nwchem_path_support import prepare_initial_path
from hfauto.chemistry.geometry import kabsch_rmsd
from hfauto.chemistry.path_initialization import initialize_reaction_path
from hfauto.chemistry.reactions import coordinate_term_value
from hfauto.chemistry.xyz import XYZ, read_xyz
from hfauto.chemistry.xyz_trajectory import read_xyz_trajectory
from hfauto.core.schemas.artifact import Artifact
from hfauto.core.schemas.chemistry import ReactionCoordinateTermRecord

ROOT = Path(__file__).resolve().parents[1]
HONO_TRANS = ROOT / "examples" / "m3_hono_isomerization" / "trans.xyz"
HONO_CIS = ROOT / "examples" / "m3_hono_isomerization" / "cis.xyz"
HONO_COORDINATE = {
    "reaction_coordinate": {
        "min_change": 1.0,
        "terms": [
            {
                "kind": "dihedral",
                "atoms": [0, 1, 2, 3],
                "coefficient": 1.0,
                "label": "H-O-N-O torsion",
            }
        ],
    }
}


def _species(path: Path, artifact_id: str) -> Artifact:
    return Artifact(
        artifact_id=artifact_id,
        artifact_type="species_optimized",
        paths={"xyz": str(path)},
        data={
            "species_id": artifact_id,
            "xyz_path": str(path),
            "charge": 0,
            "multiplicity": 1,
        },
    )


def test_hono_dihedral_path_is_smooth_collision_free_and_exact() -> None:
    start, end = read_xyz(HONO_TRANS), read_xyz(HONO_CIS)

    initialized = initialize_reaction_path(
        HONO_COORDINATE, start, end, image_count=7
    )

    assert initialized is not None
    assert initialized.metadata["strategy"] == "acyclic_dihedral_rotation"
    assert initialized.metadata["minimum_nonbonded_distance_A"] > 0.55
    assert initialized.metadata["minimum_covalent_bond_compression_ratio"] > 0.70
    assert initialized.metadata["maximum_aligned_segment_rmsd_A"] < 0.5
    assert kabsch_rmsd(initialized.images[0].coords, start.coords) < 1.0e-10
    assert kabsch_rmsd(initialized.images[-1].coords, end.coords) < 1.0e-8
    term = ReactionCoordinateTermRecord(
        kind="dihedral", atoms=(0, 1, 2, 3)
    )
    values = [coordinate_term_value(term, image.coords) for image in initialized.images]
    increments = [
        float(np.arctan2(np.sin(right - left), np.cos(right - left)))
        for left, right in pairwise(values)
    ]
    assert abs(sum(increments)) > 3.0
    assert max(abs(value) for value in increments) < 0.7
    oh_lengths = [
        np.linalg.norm(image.coords[0] - image.coords[1])
        for image in initialized.images
    ]
    assert min(oh_lengths) > 0.8


def test_unsupported_distance_coordinate_does_not_claim_chemical_path() -> None:
    start, end = read_xyz(HONO_TRANS), read_xyz(HONO_CIS)
    reaction = {
        "reaction_coordinate": {
            "terms": [{"kind": "distance", "atoms": [0, 3]}]
        }
    }

    assert initialize_reaction_path(reaction, start, end, image_count=7) is None


def test_bonded_center_atom_transfer_uses_external_arc() -> None:
    symbols = ["C", "N", "H"]
    start = XYZ(
        symbols,
        np.asarray([[0.0, 0.0, 0.0], [1.16, 0.0, 0.0], [-1.06, 0.0, 0.0]]),
        "HCN",
    )
    end = XYZ(
        symbols,
        np.asarray([[0.0, 0.0, 0.0], [1.16, 0.0, 0.0], [2.16, 0.0, 0.0]]),
        "HNC",
    )
    reaction = {
        "bond_changes": [
            {"kind": "form", "atoms": [1, 2]},
            {"kind": "break", "atoms": [0, 2]},
        ]
    }

    initialized = initialize_reaction_path(reaction, start, end, image_count=7)

    assert initialized is not None
    assert initialized.metadata["strategy"] == "bonded_center_transfer_arc"
    middle = initialized.images[3].coords
    assert np.linalg.norm(middle[2] - middle[0]) > 0.8
    assert np.linalg.norm(middle[2] - middle[1]) > 0.8
    assert np.linalg.norm(np.cross(middle[1] - middle[0], middle[2] - middle[0])) > 0.5
    donor_distances = [
        np.linalg.norm(image.coords[2] - image.coords[0])
        for image in initialized.images[:3]
    ]
    acceptor_distances = [
        np.linalg.norm(image.coords[2] - image.coords[1])
        for image in initialized.images[-3:]
    ]
    assert min(donor_distances) > 0.9
    assert min(acceptor_distances) > 0.9
    assert kabsch_rmsd(initialized.images[0].coords, start.coords) < 1.0e-10
    assert kabsch_rmsd(initialized.images[-1].coords, end.coords) < 1.0e-8


def test_cyclic_dihedral_fails_closed() -> None:
    symbols = ["C", "C", "C", "C"]
    start = XYZ(
        symbols,
        np.asarray(
            [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [1.5, 1.5, 0.0], [0.0, 1.5, 0.0]]
        ),
        "planar ring",
    )
    end = XYZ(symbols, start.coords.copy(), "puckered ring")
    end.coords[0, 2] = 0.5
    reaction = {
        "reaction_coordinate": {
            "terms": [
                {"kind": "dihedral", "atoms": [0, 1, 2, 3]}
            ]
        }
    }

    with pytest.raises(ValueError, match="absent or cyclic"):
        initialize_reaction_path(reaction, start, end, image_count=7)


def test_neb_renderer_receives_generated_xyz_path(tmp_path: Path) -> None:
    reactant = _species(HONO_TRANS, "hono_trans")
    product = _species(HONO_CIS, "hono_cis")
    reaction = Artifact(
        artifact_id="rxn_hono",
        artifact_type="reaction",
        data={"reaction_id": "rxn_hono", **HONO_COORDINATE},
    )
    start, end = read_xyz(HONO_TRANS), read_xyz(HONO_CIS)
    generated, metadata = prepare_initial_path(
        reaction,
        start,
        end,
        {},
        tmp_path / "generated.xyz",
        image_count=7,
    )
    assert generated is not None
    method = {
        "functional": "pbe0",
        "basis": "def2-svpd",
        "nbeads": 7,
        "path_initial_trajectory": str(generated),
        "path_initialization": metadata,
    }

    input_path = NWChemNEBEngine().render_neb_input(
        reactant, product, method, tmp_path / "neb" / "path.nw"
    )

    rendered = input_path.read_text(encoding="utf-8")
    assert "xyz_path initial_path.xyz" in rendered
    assert "nbeads 7" in rendered
    assert len(read_xyz_trajectory(input_path.parent / "initial_path.xyz")) == 7
