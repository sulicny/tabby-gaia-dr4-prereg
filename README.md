# Pre-registered Gaia DR4 test of a proposed companion to Tabby's star

Author: Steve Ulicny, Independent Researcher

This repository holds a pre-registered analysis plan and the code that will
test a proposed ~9 Jupiter-mass companion to KIC 8462852 (Tabby's star;
Gaia DR3 source_id 2081900940499099136) using Gaia Data Release 4 epoch
astrometry. The claim under test is that of Madurga-Favieres et al. (2026),
MNRAS, 550, stag1273, [doi:10.1093/mnras/stag1273](https://doi.org/10.1093/mnras/stag1273),
[arXiv:2607.20643](https://arxiv.org/abs/2607.20643): a companion of mass
9.4 +4.9/−4.1 Jupiter masses, radial-velocity semi-amplitude
K = 162 +76/−69 m/s, in one of five photometric period windows between
1029.7 and 1368.3 days.

The analysis rules were written and frozen before any Gaia DR4 data were
accessed. Gaia DR4 is scheduled for 2026-12-02. Until that date, the code
can be run only on simulated data (`--dry-run`).

## Pre-registration versions

Two versions of the rules are included, both as frozen, byte-identical
copies of the files locked on 2026-11-12 (America/Detroit, UTC-5). They
must not be edited. Version 2 governs the planned 2026-12-02 run; version 1
is kept unchanged for the record.

| File | Frozen (America/Detroit) | SHA-256 |
| --- | --- | --- |
| `prereg.md` (v1) | 2026-11-12 14:30:54 | `c5374c6e41de69304483b3c1faba2804e66e57cf33240e7d0bf290e53e2eb5d6` |
| `prereg_v2.md` (v2) | 2026-11-12 15:42:01 | `432664ad096b266eb67546215240ea2cbf0fe5a5e8af8c6131d25197d21f4a33` |

The v2 fitter, `phaseB/02_fit.py`, was frozen at 2026-11-12 15:42:09
(SHA-256 `5e79d81632f8500472c4a57295ec2e59c8046582bc6aa2307e424adbab45e783`).
The superseded v1 fitter is `phaseB/02_fit_v1.py`.

The documents still contain the wording used at freeze time ("LOCAL ONLY",
"not posted anywhere"). That records their status when they were locked,
before this public posting. A path mentioned in those files
(`/workspace/memory/progress-lab.md`) is a local lab note and is not part
of this package.

### Why version 2 exists

Offline tests of the v1 fitter showed that its noise-rescaling step,
s² = max(1, χ²_M0 / (N−5)), partly absorbed a real companion signal
because χ²_M0 includes that signal. With a 9.4 Jupiter-mass companion
injected in window W4, the v1 detection rate was 67.5% (60% in W4
itself), with 0/40 false detections on null data.

Version 2 makes two changes, and only these:

1. Noise rescale: s²_W = max(1, χ²_Mc(W) / (N−7)), the reduced χ² of the
   fit that includes that window's orbit term.
2. Gate G3: a high false-positive rate among comparison stars raises that
   window's detection threshold to the controls' 99th percentile of S_W.
   It does not declare the result inconclusive. G3 fails only if fewer
   than 200 valid controls are available. In v1 the raise and the
   inconclusive outcome were written in a way that could not both apply.

Under the same offline tests, v2 detected the injected W4 companion in
75% of runs (70% in W4 itself), still with 0/40 false detections and a
median W4 95% mass upper limit of 3.9 Jupiter masses. When the injected
truth is drawn from the paper's posterior (random window and K; 60
seeds), the detection rate is about 55%. Details are in
`results/phaseB_tests_v2.md`. Phase A forecasts from the earlier
simulator (`scripts/sim.py`) are in `results/forecast.csv` and
`results/forecast.png`. Those Phase A figures use an unrescaled
statistic and are not a v2 detection-rate forecast.

## How to run

Requires Python 3.10 or later. The offline tests were run with Python
3.13 and the pins in `requirements.txt` / `environment.yml`.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run all commands from the repository root.

### Dry run (simulated data; no network)

```bash
python phaseB/02_fit.py --dry-run --case null
python phaseB/02_fit.py --dry-run --case injected
```

This generates a temporary simulated DR4-format data set, applies the v2
rule, and writes `memo/phaseB_dryrun_<case>.md` plus
`results/phaseB_dryrun_<case>.png`. Use `--skip-plots` to skip the figure.
The default `--hashes-file` is `results/hashes.txt`. A dry run does not
abort on a hash mismatch; it reports the match status and continues.

`--case null` injects no companion. `--case injected` injects a 9.4
Jupiter-mass companion in window W4 (override with `--seed` and the
generator's window options if calling `phaseB/gen_sim_dr4.py` directly).

The dry-run path and the Phase A simulator need the three CSV files in
`data/` (see `data/README.md`).

### After Gaia DR4 (2026-12-02 or later)

1. Download Tabby's catalog rows, flags, NSS tables, and epoch
   astrometry, plus the pre-registered control sample:

   ```bash
   python phaseB/01_fetch_dr4.py
   ```

   This talks to the public Gaia TAP and DataLink services (no login).
   Table names follow the draft DR4 data model and should be checked
   against TAP_SCHEMA once the archive is live. A schema-only probe is
   `python phaseB/01_fetch_dr4.py --dry` (prints the HEALPix source_id
   range; this is not the same flag as `--dry-run` on the fitter).

2. Run the frozen v2 fitter in real mode:

   ```bash
   python phaseB/02_fit.py --data-dir data_dr4 --out memo/phaseB_result.md --plot results/phaseB_fit.png
   ```

   Real mode exits with status 2 and writes no result if `prereg.md`,
   `prereg_v2.md`, or `phaseB/02_fit.py` does not match the corresponding
   last entry in `results/hashes.txt`, or if a `SYNTHETIC_MARKER` file is
   present in the data directory.

Do not edit `prereg.md`, `prereg_v2.md`, or `phaseB/02_fit.py`. Any
departure from the locked rule has to be labelled post hoc.

### Tests and the Phase A simulator

```bash
python phaseB/test_extra.py
python phaseB/test_fit.py          # long: 40+60 simulation seeds
python scripts/sim.py              # long: Phase A forecast, 1500 sims/case
```

`test_fit.py` as published here uses `sys.executable` for the
subprocess checks in section T4. The historical test log in
`results/phaseB_tests_v2.md` records the local interpreter path used
when that log was written.

## Repository contents

- `prereg.md`, `prereg_v2.md` — locked analysis rules
- `phaseB/01_fetch_dr4.py` — Gaia Archive download (after 2026-12-02)
- `phaseB/02_fit.py` — v2 fitter (governs the Dec 2 run)
- `phaseB/02_fit_v1.py` — v1 fitter, retained for comparison
- `phaseB/gen_sim_dr4.py` — simulated DR4-format data generator
- `phaseB/test_fit.py`, `phaseB/test_extra.py` — offline tests of the v2 rule
- `scripts/sim.py` — Phase A forecast simulator
- `results/hashes.txt` — frozen hashes in the format the fitter checks
- `SHA256SUMS` — SHA-256 of every file in this package except itself
- `results/phaseB_tests_v2.md`, `results/forecast.csv`, `results/forecast.png`
- `data/` — public Gaia DR3 extracts and GOST predicted transits

## Hashes and integrity

`results/hashes.txt` lists the hashes that were frozen on 2026-11-12 and
is what `phaseB/02_fit.py` consults in real mode. `SHA256SUMS` covers the
full package (including README, licenses, and data files) and can be
checked with `sha256sum -c SHA256SUMS`.

## License

Code is released under the MIT License (see `LICENSE`). Documentation,
including the pre-registration files and the results write-ups, is
released under CC BY 4.0 (see `NOTICE`).

## AI-assistance disclosure

Analysis code and text drafted with AI coding assistants under the
author's direction.

## Acknowledgements

This work will use, and the `data/` extracts already use, data from the
European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia),
processed by the Gaia Data Processing and Analysis Consortium (DPAC).
The official Gaia DR4 credit line and DOI should be added when they are
posted at https://www.cosmos.esa.int/web/gaia/data-release-4.
