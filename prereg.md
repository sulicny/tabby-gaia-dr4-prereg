# PRE-REGISTRATION (LOCAL ONLY – not posted anywhere) — Gaia DR4 test of the Tabby's-star companion claim
Written before any Gaia DR4 data were accessed (DR4 release scheduled 2026-12-02; draft data model only). Local file; its sha256 is recorded in /workspace/memory/progress-lab.md. Not published; no OSF/Zenodo/etc.

## 1. Claim under test
arXiv:2607.20643: a ~9.4 +4.9/−4.1 MJup companion to KIC 8462852 (Gaia source_id 2081900940499099136; DR3 ϖ = 2.2545 mas), K = 162 +76/−69 m/s, one of five period windows (W1 1029.7–1031.5, W2 1034.5–1035.3, W3 1301.0–1303.6, W4 1345.648±0.008, W5 1367.3–1368.3 d), transit anchor BJD 2458730.164. Nominal astrometric amplitudes α_nom (K=162, e=0.09, ϖ=2.2545 mas): 34.5, 34.6, 43.5, 45.0, 45.7 µas (W1–W5). See memo/paper_extract.md.

## 2. Data (only these)
(a) Tabby's DR4 epoch_astrometry (DataLink) + gaia_source/all_source_flags/NSS rows; (b) controls: DR4 sources within the HEALPix(level-6) cells covering a 2.5° disc around Tabby, 11.2<G<12.4, |bp_rp−0.776|<0.25, astrometric_params_solved=31, RUWE<1.4, with epoch astrometry; the 500 (minimum 200) with the lowest source_id are used. No other data/selection changes after data access.

## 3. Processing (frozen)
Use only CCD samples with used_by_agis_al=True; weight 1/(centroid_pos_error_al²+agis_source_excess_noise²); collapse the CCD samples of each FoV transit to one weighted abscissa (time = weighted mean of obs_time_tcb + obs_time_bary_corr; scan angle and parallax_factor_al from the table). Model: w = dα*·sinψ + dδ·cosψ + ϖ·pf + (μα* sinψ + μδ cosψ)·t + orbit. M0 = 5 parameters. Mc(W) = M0 + orbit term (cx sinψ + cy cosψ)·X′_W(t), where X′_W is the unit edge-on stellar reflex template for P = window centre, e=0.09, ω_s=141°, phase anchored to transit BJD 2458730.164 with ω_p=ω_s+π (opposite convention reported as a sensitivity). Rescaling: s² = max(1, χ²_M0/(N−5)); all Δχ² below are divided by s². Statistic S_W = (χ²_M0 − χ²_Mc)/s²; amplitude a_W = hypot(cx,cy).

## 4. Calibration gates (all must pass, else the outcome is INCONCLUSIVE (data quality) and the cause is reported)
G1 Tabby has epoch astrometry with ≥60 FoV transits and DR4 RUWE<1.4. G2 The 5-parameter fit reproduces the catalogue parallax within 3σ. G3 Over the controls, S_W (each W) has fraction above 13.82 (χ²₂, p=0.001) ≤ 1.5% (expected 0.1%); if higher, the detection threshold for that window is raised to the controls' 99th percentile of S_W. G4 Tabby's ipd_frac_multi_peak and ipd_gof_harmonic_amplitude ≤ the 95th percentile of the controls (blend/neighbour-contamination screen). Phase A sims: a scan-angle-windowed leak needs ≥10% of the G=17.6 neighbour's flux inside Tabby's windows to fake a detection; static offsets are absorbed.

## 5. Decision rule (applied per window, then overall)
DETECT(W): S_W ≥ max(T_W, 19.1) [T_W = 13.82 unless raised by G3; 19.1 = ΔBIC≥10 for N≈93] AND 0.4 ≤ a_W/α_nom(W) ≤ 2.0 AND the same-window fit with the opposite-phase template does not give S larger by more than 4 (phase consistency). Overall DETECTED if any window satisfies DETECT; windows ranked by S; DR4 NSS orbit/acceleration rows for Tabby (period within ±5% of any window) are reported as supporting evidence only.
EXCLUDE(W): not DETECT(W) and the 95% credible upper limit on amplitude (flat prior on a ≥ 0, uniform orientation, Gaussian likelihood in (cx,cy), other parameters profiled) converts to a mass limit M95(W) < 5.3 MJup (= 9.4 − 4.1, the paper's 1σ lower bound; conversion K→M uses M*=1.19 Msun, e=0.09, sin i=1). WEAK-EXCLUDE(W): M95(W) < 9.4 MJup (rules out the face-value claim but not its low-mass tail).
Overall: CLAIM SUPPORTED if DETECTED; CLAIM EXCLUDED (strong) if all five windows are EXCLUDE; "face value (9.4 MJup) excluded" if all five are at least WEAK-EXCLUDE; otherwise INCONCLUSIVE (report per-window S_W, a_W, M95). Windows may be dropped as inactive only on the basis of published, citable transit non-detections dated before 2026-12-02 and listed explicitly; the unrestricted five-window result is always also reported. Every deviation from this document is labelled post hoc in the memo.

## 6. Phase B analysis outline (code to be written under this rule before data access, not altered afterwards)
01_fetch_dr4.py (ready) → 02_fit.py (M0, M1 acceleration, Mc, Mfree grid 500–2500 d, circular, BIC) → control null and real-noise injection (Mc orbits of 9.4 MJup into controls' abscissae, reporting recovery) → memo. Free-period search (Mfree) is exploratory and not used in the rule.

## 7. Phase A forecast (recorded for later comparison; sims frozen as scripts/sim.py, results/forecast.csv; per-FoV σ = 41.7 µas × r, r ∈ {1,0.75,0.5} from DR3)
Claim true: P(DETECT) central ≈ 0.69–0.83 (W1<W5), pessimistic 0.48–0.65, optimistic 0.86–0.90. Claim false: P(strong EXCLUDE) central 0.66–0.86, pessimistic 0.16–0.44, optimistic ≥0.98; median M95 4.1–4.8 MJup (central), 5.5–6.5 (pessimistic). Most likely overall outcome: no detection; strong exclusion or inconclusive depending on DR4 noise.

## 8. Not covered
Statistical significance of any Tabby-specific anomaly outside this rule; any public posting; any DR4 NSS-orbit-only claim. Nothing is posted or sent without Steve's separate yes.
Frozen: 2026-11-12T14:30:54-05:00 (America/Detroit)
