# Phase B offline tests

Run: 2026-11-12 15:35 ET

Commands: `/workspace/lab-dasch/.conda/bin/python phaseB/test_fit.py`; `python phaseB/02_fit.py --dry-run --case injected|null`. No network; prereg.md unchanged.

## T1 — end-to-end, both serializations

- injected, gaia_array: **CLAIM SUPPORTED**; S=[2.4187, 2.3928, 23.2036, 24.8022, 24.9308]; gates={'G1': True, 'G2': True, 'G3': True, 'G4': True}.

- injected, standard: **CLAIM SUPPORTED**; S=[2.4187, 2.3928, 23.2036, 24.8022, 24.9308]; gates={'G1': True, 'G2': True, 'G3': True, 'G4': True}.

- null, gaia_array: **face value (9.4 MJup) excluded**; S=[1.7337, 1.6959, 2.8346, 2.7135, 2.6204]; gates={'G1': True, 'G2': True, 'G3': True, 'G4': True}.

- null, standard: **face value (9.4 MJup) excluded**; S=[1.7337, 1.6959, 2.8346, 2.7135, 2.6204]; gates={'G1': True, 'G2': True, 'G3': True, 'G4': True}.

## T2 — 40 independent seeds per truth

|Fixed truth|Any DETECT|W4 DETECT|Median W4 UL MJup|Median s2 W4|All gates fraction|
|---|---:|---:|---:|---:|---:|---:|
|null|0.000|0.000|3.902|1.006|1.000|
|injected|0.750|0.700|11.830|1.045|1.000|

## Claim-true posterior — 60 seeds

Seeds 2000–2059; independent RNG 20261112 selects uniform W1–W5 and K=162+z*(69 if z<0 else 76), rejecting K<0. Mass maps K at the window midpoint, e=0.09, sin i=1; generator draws the remaining truth parameters. Central noise=0.75; 100 controls and relaxed G3.
|Truth|n|Any DETECT|Own DETECT|Own raw >=13.82|Central raw forecast|
|---|---:|---:|---:|---:|---:|
|Overall|60|0.550|0.550|0.833|—|
|W1|11|0.455|0.455|0.818|~0.70|
|W2|19|0.632|0.632|0.895|~0.70|
|W3|11|0.636|0.636|0.727|0.80|
|W4|11|0.545|0.545|0.909|0.81|
|W5|8|0.375|0.375|0.750|0.84|

Forecast uses a raw statistic and nominal central assumptions; full DETECT adds amplitude, opposite-phase and 19.1 criteria. Small per-window samples limit comparison precision.

## T3 — empirical DR3 scales and gate failures

- DR3-derived per-star scales: FPR=[0.0, 0.0, 0.0, 0.0, 0.0], G3 PASS; control s2 min/median/max=1.0000/1.0048/1.4228. Errors track the per-star scale, so s2 is near 1, clipped below at 1.

- RUWE 1.6: G1 failed and headline INCONCLUSIVE (data quality), PASS.

- 40 FoV: G1 failed and quality veto applied, PASS.

- Catalogue parallax offset +1 mas: G2 failed, both barycentric signs diagnosed, quality veto applied, PASS.

- IPD harmonic amplitude 10: G4 failed and quality veto applied, PASS.

- Only 199 controls: G3 failed, PASS.

- Errors divided by 4: s2=15.2275; S uses each Mc residual rescale (M0 diagnostic shown), PASS.

## T4 — real-mode safety

- prereg.md: exit 2; no result written, PASS.

- prereg_v2.md: exit 2; no result written, PASS.

- phaseB/02_fit.py: exit 2; no result written, PASS.

- SYNTHETIC_MARKER: exit 2; no result written, PASS.

- missing-v2: exit 2; no result written, PASS.

- missing-fit: exit 2; no result written, PASS.

## Additional checks

- Analytic flat-amplitude Gaussian 95% limit, exact mass inversion, zero displacement at transit for both conventions, null/NaN parsing, and AST equality of copied generator orbital functions: PASS.

## Fixed-mass comparison

- Null W4 median UL 3.902 MJup versus forecast ~4.1. Injected DETECT fraction 0.750 versus ~0.98 Phase A raw-statistic power. The posterior comparison above uses a different estimand.

Elapsed 225.4 seconds.

## v1 comparison
|Rule|Injected any|Injected W4|Null any FPR|Null median W4 UL|
|---|---:|---:|---:|---:|
|v1|0.675|0.600|0/40|3.89|
|v2|0.750|0.700|0/40|3.902|

Reduced runs explicitly relax G3 to 100 controls for both truths; production requires 200. DETECT fractions precede the quality veto. No network. Real DR4 performance remains unverified.

## G3 raise and independent regression checks

```text
Contaminated controls: FPR [0.   0.   0.05 0.05 0.05] thresholds [ 13.82        13.82       341.84285978 371.3540144  365.97291875] G3 passes without quality veto PASS
Brace arrays, t/f booleans, NaN/null padding through collapse PASS
Restricted-result headline preservation PASS
Full-model residual, phase-specific rescale, statistic and covariance UL PASS
```

## CLI and frozen files

- Injected/null CLI dry runs: exit 0, all gates pass; hash statuses v1=True, v2=False (absent), fitter=True. Both serializations match in T1.
- prereg.md, phaseB/02_fit_v1.py and results/phaseB_tests.md were left untouched. prereg_v2.md was neither created nor modified.
- Latest fitter SHA256 appended to results/hashes.txt; stale v1 fitter entry retained. Real mode requires the operator’s prereg_v2.md and matching hash entry. Dry-run memo marks its hash MISSING.
- Shared-memory append is outside writable roots; local dated note is results/progress-cli-phaseb-v2.md.
