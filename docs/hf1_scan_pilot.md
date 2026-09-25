# HF1 restrained-scan pilot

This pilot is the small, auditable replacement for sending collapsed RC/IP
geometries directly to TS search. It evaluates four HF1 cases from the existing
production preoptimization manifest.

## Run

From the repository root in the production WSL environment:

```bash
export PATH="$HOME/.local/bin:$PATH"
~/.venvs/hfauto-prod/bin/python scripts/run_hf1_scan_pilot.py
~/.venvs/hfauto-prod/bin/python scripts/report_hf1_scan_pilot.py
```

The run is resumable. Each xTB scan and each NWChem single point is written to
`runs/hf1_scan_validation_001/cases/<case>/case_result.json` before the next job
starts.

## Scientific contract

- Reaction coordinate: `q = r(H-F) - r(B-H)`.
- Neutral complex: `q < -0.15 A`; ion pair: `q > +0.15 A`.
- GFN2-xTB: 11-point concerted restrained relaxed scan.
- Scan acceptance: exact point count, finite energies, correct endpoint classes,
  and maximum distance-constraint deviation no greater than 0.12 A.
- DFT validation: PBE0/def2-SVPD NWChem single points at scan points 0, 5, and 10.
- DFT acceptance: real QM execution, SCF convergence, normal termination, and no
  dummy fallback.
- An endpoint maximum is never labeled as a TS. TS/IRC work requires a resolved
  internal maximum and distinct validated minima on both sides.

## Current result

All 44 xTB points and all 12 NWChem points passed their quality gates. Every
profile has a shallow neutral-side minimum followed by a rise to the ion-pair
endpoint, with no internal maximum. The current deliverable is therefore an
electronic-energy endpoint assessment, not an activation barrier or free-energy
profile.

Open `runs/hf1_scan_validation_001/report/report.html` for the self-contained
technical report. PNG/SVG figures, four 11-frame GIF animations, multi-XYZ
trajectories, and reviewed CSV files are in the same report directory.
