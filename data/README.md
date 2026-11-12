# Input tables for the offline simulators

These files are public catalog extracts and predicted-transit times. They are
required by `scripts/sim.py` (Phase A forecast) and by `phaseB/gen_sim_dr4.py`
(the Phase B dry-run generator used by `phaseB/02_fit.py --dry-run`).

| File | Source | Role |
| --- | --- | --- |
| `tabby_dr3.csv` | Gaia DR3 `gaia_source` row for source_id 2081900940499099136 | Tabby's DR3 astrometric noise scale |
| `controls_raw.csv` | Gaia DR3 `gaia_source` rows near Tabby, selected as in the pre-registration | Per-star noise scales for simulated controls |
| `gost_tabby_dr4window.csv` | ESA Gaia Observation Schedule Tool (GOST), DR4 time window | Predicted field-of-view transit times, scan angles, and parallax factors |

Gaia data: ESA/Gaia/DPAC. Cite the Gaia mission and DR3 release papers when
reusing these extracts.
