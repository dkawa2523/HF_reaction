from __future__ import annotations

from hfauto.backends.conformer.crest import CRESTConformerBackend
from hfauto.backends.conformer.rdkit import RDKitConformerBackend
from hfauto.backends.db.atct import ATcTProvider
from hfauto.backends.db.cas_common_chemistry import CASCommonChemistryProvider
from hfauto.backends.db.cccbdb import CCCBDBProvider
from hfauto.backends.db.comptox import CompToxProvider
from hfauto.backends.db.dummy import DummyDBProvider
from hfauto.backends.db.niosh import NIOSHProvider
from hfauto.backends.db.nist_webbook import NISTWebBookProvider
from hfauto.backends.db.pubchem import PubChemProvider
from hfauto.backends.kinetics.cantera import CanteraEngine
from hfauto.backends.kinetics.tst import TSTKineticsEngine
from hfauto.backends.qm.dummy import DummyQMEngine
from hfauto.backends.qm.nwchem import NWChemEngine
from hfauto.backends.qm.orca import OrcaEngine
from hfauto.backends.qm.xtb import XTBEngine
from hfauto.backends.reaction_discovery.readuct import ReaDuctDiscoveryBackend
from hfauto.backends.thermo.arkane import ArkaneEngine
from hfauto.backends.thermo.goodvibes import GoodVibesEngine
from hfauto.backends.ts.dummy import DummyTSEngine
from hfauto.backends.ts.nwchem_neb import NWChemNEBEngine
from hfauto.backends.ts.nwchem_saddle import NWChemSaddleEngine
from hfauto.backends.ts.nwchem_string import NWChemStringEngine
from hfauto.backends.ts.orca_nebts import ORCANEBTSEngine
from hfauto.backends.ts.pysisyphus import PysisyphusEngine
from hfauto.backends.ts.pysisyphus_saddle import PysisyphusSaddleEngine


def _normalize_spec(spec, explicit_kwargs: dict | None = None) -> tuple[str, dict]:
    explicit_kwargs = explicit_kwargs or {}
    if isinstance(spec, str):
        return spec, explicit_kwargs
    if isinstance(spec, dict):
        name = spec.get("name") or spec.get("provider") or spec.get("engine")
        kwargs = {k: v for k, v in spec.items() if k not in {"name", "provider", "engine"}}
        kwargs.update(explicit_kwargs)
        return name, kwargs
    raise TypeError(f"Invalid backend spec: {spec!r}")


def get_qm_engine(spec: str | dict, **kwargs):
    engines = {"dummy": DummyQMEngine, "nwchem": NWChemEngine, "orca": OrcaEngine, "xtb": XTBEngine}
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown QM engine: {name}")
    return engines[name](**init_kwargs)


def get_conformer_engine(spec: str | dict, **kwargs):
    engines = {"rdkit": RDKitConformerBackend, "crest": CRESTConformerBackend}
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown conformer engine: {name}")
    return engines[name](**init_kwargs)


def get_reaction_discovery_engine(spec: str | dict, **kwargs):
    engines = {"readuct": ReaDuctDiscoveryBackend}
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown reaction-discovery engine: {name}")
    return engines[name](**init_kwargs)


def get_ts_engine(spec: str | dict, **kwargs):
    engines = {
        "dummy": DummyTSEngine,
        "nwchem_neb": NWChemNEBEngine,
        "nwchem_saddle": NWChemSaddleEngine,
        "nwchem_string": NWChemStringEngine,
        "orca_nebts": ORCANEBTSEngine,
        "pysisyphus": PysisyphusEngine,
        "pysisyphus_saddle": PysisyphusSaddleEngine,
    }
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown TS/IRC engine: {name}")
    return engines[name](**init_kwargs)


def get_db_provider(spec: str | dict, **kwargs):
    providers = {
        "dummy": DummyDBProvider,
        "pubchem": PubChemProvider,
        "nist": NISTWebBookProvider,
        "nist_webbook": NISTWebBookProvider,
        "comptox": CompToxProvider,
        "niosh": NIOSHProvider,
        "atct": ATcTProvider,
        "cccbdb": CCCBDBProvider,
        "cas_common_chemistry": CASCommonChemistryProvider,
        "cas": CASCommonChemistryProvider,
    }
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in providers:
        raise KeyError(f"Unknown DB provider: {name}")
    return providers[name](**init_kwargs)


def get_thermo_engine(spec: str | dict, **kwargs):
    engines = {
        "goodvibes": GoodVibesEngine,
        "internal": GoodVibesEngine,
        "internal_quasirrho": GoodVibesEngine,
        "quasi_rrho": GoodVibesEngine,
        "dummy": GoodVibesEngine,
        "arkane": ArkaneEngine,
    }
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown thermo engine: {name}")
    return engines[name](**init_kwargs)


def get_kinetics_engine(spec: str | dict, **kwargs):
    engines = {
        "tst": TSTKineticsEngine,
        "simple": TSTKineticsEngine,
        "internal": TSTKineticsEngine,
        "internal_tst": TSTKineticsEngine,
        "dummy": TSTKineticsEngine,
        "cantera": CanteraEngine,
        "arkane": ArkaneEngine,
    }
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown kinetics engine: {name}")
    return engines[name](**init_kwargs)


def get_arkane_engine(spec: str | dict = "arkane", **kwargs):
    engines = {"arkane": ArkaneEngine, "dummy": ArkaneEngine}
    name, init_kwargs = _normalize_spec(spec, kwargs)
    if name not in engines:
        raise KeyError(f"Unknown Arkane engine: {name}")
    return engines[name](**init_kwargs)
