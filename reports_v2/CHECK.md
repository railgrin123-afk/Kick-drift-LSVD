# reports_v2 — build check

Built by `python/kd_present.py`, read-only from `python/out/` (stored run results) and `reports/` (the previous report). Nothing outside `reports_v2/` was written.

## Existing outputs unchanged

- SHA-256 of every file in `reports/` and `python/out/`, before and after: 85 files compared; changed: none; new files: none

## Previous report reproduced

- The stored run files, re-analysed into a scratch folder with the unchanged report code, reproduce `reports/kd_tables.json` (all tables and checks), `reports/kd_tables.md` and `reports/kd_findings.json` exactly.

## Pages and terms

- paper1_stretching.pdf: 6 pages; banned terms and internal codes found: none.
- paper2_reentry.pdf: 5 pages; banned terms and internal codes found: none.
- appendix.pdf (all pages before this check): banned terms found: none; the estimator for censored data is named once, in Methods.
- 2. Figures: 21 built; every text element is ≥ 11 pt at printed size (figures are placed at 100%); no tick labels overlap; no banned terms; at most 3 panels each.
- 3. Numbers: 186 quoted numbers; 164 taken unchanged from the previous report, 11 derived arithmetically from previous-report values (ratios, ranges), 11 new (read directly from the stored run files; listed in CHECK.md), 0 differ.
- 4. The stored run files, re-analysed in memory, reproduce the previous report's tables, checks and summary exactly; no existing output was written (checksums before and after in CHECK.md).
- Terms scanned in the main reports: cosine, LSVD, Eulerian, compressed, Kaplan-Meier / Kaplan–Meier, an_br, S_THR, Pearson r, window, lattice, row0, rows/8, cycleData, "± 0.0000", and part codes such as 1A or 4C. In the appendix the same list except row0, rows/8 and cycleData, and Kaplan–Meier must appear exactly once (Methods).

## Figures

| figure file | printed size (in) | min font (pt) | overlapping tick labels | banned terms |
|---|---|---|---|---|
| p1_fig1_branched_flow | 6.46 × 2.75 | ≥ 11 | none | none |
| appA_lnS_shape | 6.46 × 2.95 | ≥ 11 | none | none |
| p1_fig3_why | 6.46 × 3.00 | ≥ 11 | none | none |
| appA_caustics | 6.46 × 3.35 | ≥ 11 | none | none |
| p1_fig5_lowS_branches | 6.46 × 3.20 | ≥ 11 | none | none |
| p2_fig1_return | 6.46 × 3.50 | ≥ 11 | none | none |
| p2_fig2_hazard | 6.46 × 3.60 | ≥ 11 | none | none |
| p2_fig3_depth | 6.46 × 3.70 | ≥ 11 | none | none |
| appA_branched_low | 6.46 × 2.75 | ≥ 11 | none | none |
| appA_branched_high | 6.46 × 2.75 | ≥ 11 | none | none |
| appA_manifold_annealed | 6.46 × 3.50 | ≥ 11 | none | none |
| appA_manifold_quenched | 6.46 × 3.50 | ≥ 11 | none | none |
| appA_manifold_periodic | 6.46 × 3.50 | ≥ 11 | none | none |
| appA_quantiles | 6.46 × 3.40 | ≥ 11 | none | none |
| appA_tail | 6.46 × 3.60 | ≥ 11 | none | none |
| appA_rate_function | 6.46 × 3.40 | ≥ 11 | none | none |
| appA_lowr_branches | 6.46 × 3.20 | ≥ 11 | none | none |
| appA_r_vs_S | 6.46 × 3.40 | ≥ 11 | none | none |
| appA_metric | 6.46 × 3.60 | ≥ 11 | none | none |
| appA_leave_again | 6.46 × 3.40 | ≥ 11 | none | none |
| appA_cohorts | 6.46 × 3.50 | ≥ 11 | none | none |

## Every figure and table, and the question that references it

| item | where | referenced by |
|---|---|---|
| Paper 1, Figure 1 | paper1_stretching.pdf | Paper 1 Q1 |
| Paper 1, Figure 2 | paper1_stretching.pdf | Paper 1 Q3 |
| Paper 1, Figure 3 | paper1_stretching.pdf | Paper 1 Q5 |
| Paper 2, Figure 1 | paper2_reentry.pdf | Paper 2 Q1 |
| Paper 2, Figure 2 | paper2_reentry.pdf | Paper 2 Q2 |
| Paper 2, Figure 3 | paper2_reentry.pdf | Paper 2 Q5 |
| Paper 1, Table 1 | paper1 main text | Paper 1 Q2 |
| Paper 1, Table 2 | paper1 main text | Paper 1 Q4 |
| Paper 1, Table 3 | paper1 main text | Paper 1 Q6 |
| Paper 2, Table 1 | paper2 main text | Paper 2 Q3 |
| Paper 2, Table 2 | paper2 main text | Paper 2 Q4 |
| Figure A1 | appendix | Paper 1 Q1 |
| Figure A2 | appendix | Paper 1 Q1 |
| Figure A3 | appendix | Paper 1 Q4 |
| Figure A4 | appendix | Paper 1 Q4 |
| Figure A5 | appendix | Paper 1 Q4 |
| Figure A6 | appendix | Paper 1 Q2 |
| Figure A7 | appendix | Paper 1 Q2 |
| Figure A8 | appendix | Paper 1 Q3 |
| Figure A9 | appendix | Paper 1 Q5 |
| Figure A10 | appendix | Paper 1 Q5 |
| Figure A11 | appendix | Paper 1 Q6 |
| Figure A12 | appendix | Paper 2 Q3 |
| Figure A13 | appendix | Paper 2 Q2 |
| Figure A14 | appendix | Paper 1 Q2 |
| Figure A15 | appendix | Paper 1 Q4 |
| Table B1 | appendix | Paper 1 Q1, Paper 1 Q2, Paper 2 Q1 |
| Table B2 | appendix | Paper 1 Q2 |
| Table B2 (seed 0) | appendix | Paper 1 Q2 |
| Table B3 | appendix | Paper 1 Q2 |
| Table B4 | appendix | Paper 1 Q2 |
| Table B5 | appendix | Paper 1 Q3 |
| Table B5 (shares) | appendix | Paper 1 Q3 |
| Table B6 | appendix | Paper 1 Q3 |
| Table B7 | appendix | Paper 1 Q4 |
| Table B8 | appendix | Paper 1 Q5 |
| Table B8 (rarefaction) | appendix | Paper 1 Q5 |
| Table B9 | appendix | Paper 1 Q5 |
| Table B10 | appendix | Paper 1 Q5 |
| Table B11 | appendix | Paper 1 Q5 |
| Table B12 | appendix | Paper 1 Q6 |
| Table B13 | appendix | Paper 2 Q2 |
| Table B14 | appendix | Paper 2 Q1, Paper 2 Q3 |
| Table B15 | appendix | Paper 2 Q5 |

Removed from the previous report because no question needs them: ridge plot and heat map of the distribution of S; stacked shares of S bands; band projections of the tracer parcels onto x; per-mode hazard figures with at-risk bars; pooled age-bin bar chart; phase-plane snapshots at N = 40 and 60; correlation-against-N figure (its numbers are in Table B10).

## Number ledger

186 numbers: 164 unchanged from the previous report, 11 derived from previous-report values, 11 new, 0 differ.

| where | source in the previous report (table, column, row) | previous | shown | status |
|---|---|---|---|---|
| figures | exact caustic count, seed 0, annealed, N = 5 | (not in previous report) | 32 | NEW |
| figures | exact caustic count, seed 0, annealed, N = 10 | (not in previous report) | 82 | NEW |
| figures | exact caustic count, seed 0, annealed, N = 20 | (not in previous report) | 1,114 | NEW |
| figures | exact caustic count, seed 0, quenched, N = 5 | (not in previous report) | 22 | NEW |
| figures | exact caustic count, seed 0, quenched, N = 10 | (not in previous report) | 152 | NEW |
| figures | exact caustic count, seed 0, quenched, N = 20 | (not in previous report) | 1,124 | NEW |
| figures | exact caustic count, seed 0, periodic, N = 5 | (not in previous report) | 2 | NEW |
| figures | exact caustic count, seed 0, periodic, N = 10 | (not in previous report) | 44 | NEW |
| figures | exact caustic count, seed 0, periodic, N = 20 | (not in previous report) | 1,866 | NEW |
| figures | 4C / s=2: rate (ratio) / {'mode': 'annealed'} | 0.1085 ± 0.0018 (1.51) | 0.1085 ± 0.0018 (1.51) | same as previous report |
| figures | 4C / pred a²/2v / {'mode': 'annealed'} | 0.0717 ± 0.0014 | 0.0717 ± 0.0014 | same as previous report |
| figures | 4C / s=2: rate (ratio) / {'mode': 'quenched'} | 0.0332 ± 0.0070 (0.92) | 0.0332 ± 0.0070 (0.92) | same as previous report |
| figures | 4C / pred a²/2v / {'mode': 'quenched'} | 0.0370 ± 0.0061 | 0.0370 ± 0.0061 | same as previous report |
| figures | 4C / s=2: rate (ratio) / {'mode': 'cosine'} | 0.0210 ± 0.0000 (1.68) | 0.0210 ± n/a (deterministic) (1.68) | same as previous report |
| figures | 4C / pred a²/2v / {'mode': 'cosine'} | 0.0125 ± 0.0000 | 0.0125 ± n/a (deterministic) | same as previous report |
| figures | iid / λ-CDF I(0+), N 40→60 → corrected for (ln S)² / {'model': 'i.i.d. model (1e7)'} | 0.0702 → 0.1107 | 0.0702 → 0.1107 | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'annealed', 'N': 5} | 4.18 ± 0.36 | 4.18 ± 0.36 | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'quenched', 'N': 5} | 4.44 ± 0.23 | 4.44 ± 0.23 | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'cosine', 'N': 5} | 3.04 ± 0.00 | 3.04 ± n/a (deterministic) | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'annealed', 'N': 60} | 1.50 ± 0.05 | 1.50 ± 0.05 | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'quenched', 'N': 60} | 2.28 ± 0.24 | 2.28 ± 0.24 | same as previous report |
| Paper 1, Q1 | 1C / contrast, all rays / {'mode': 'cosine', 'N': 60} | 3.08 ± 0.00 | 3.08 ± n/a (deterministic) | same as previous report |
| Paper 1, Q1 | previous summary text | 96.0% within 2% | 96.0% | same as previous report |
| Paper 1, Q1 | previous summary text | 97.8% within 2% | 97.8% | same as previous report |
| Paper 1, Q2 | 4C / alpha / {'mode': 'annealed'} | 0.1526 ± 0.0022 | 0.1526 ± 0.0022 | same as previous report |
| Paper 1, Q2 | 4C / alpha / {'mode': 'quenched'} | 0.1452 ± 0.0094 | 0.1452 ± 0.0094 | same as previous report |
| Paper 1, Q2 | 4C / alpha / {'mode': 'cosine'} | 0.0236 ± 0.0000 | 0.0236 ± n/a (deterministic) | same as previous report |
| Paper 1, Q2 | 4C / v / {'mode': 'annealed'} | 0.1623 ± 0.0028 | 0.1623 ± 0.0028 | same as previous report |
| Paper 1, Q2 | 4C / v / {'mode': 'quenched'} | 0.2891 ± 0.0378 | 0.2891 ± 0.0378 | same as previous report |
| Paper 1, Q2 | 4C / v / {'mode': 'cosine'} | 0.0223 ± 0.0000 | 0.0223 ± n/a (deterministic) | same as previous report |
| Paper 1, Q2 | mom / skew / {'mode': 'annealed', 'N': 60} | 0.132 ± 0.008 | 0.132 ± 0.008 | same as previous report |
| Paper 1, Q2 | mom / skew / {'mode': 'quenched', 'N': 60} | -0.265 ± 0.158 | -0.265 ± 0.158 | same as previous report |
| Paper 1, Q2 | mom / skew / {'mode': 'cosine', 'N': 60} | 1.957 ± 0.000 | 1.957 ± n/a (deterministic) | same as previous report |
| Paper 1, Q2 | mom / ex. kurt / {'mode': 'annealed', 'N': 60} | -0.057 ± 0.011 | -0.057 ± 0.011 | same as previous report |
| Paper 1, Q2 | mom / ex. kurt / {'mode': 'quenched', 'N': 60} | -0.292 ± 0.230 | -0.292 ± 0.230 | same as previous report |
| Paper 1, Q2 | mom / ex. kurt / {'mode': 'cosine', 'N': 60} | 15.366 ± 0.000 | 15.366 ± n/a (deterministic) | same as previous report |
| Paper 1, Q2 | 4C / pred a²/2v / {'mode': 'annealed'} | 0.0717 ± 0.0014 | 0.0717 ± 0.0014 | same as previous report |
| Paper 1, Q2 | 4C / s=2: rate (ratio) / {'mode': 'annealed'} | 0.1085 ± 0.0018 (1.51) | 0.1085 ± 0.0018 (1.51) | same as previous report |
| Paper 1, Q2 | min–max of the s = 1.25…10 means, 4C, annealed | derived from previous-report values | 0.0996–0.1096 | derived from previous report |
| Paper 1, Q2 | 4C / pred a²/2v / {'mode': 'quenched'} | 0.0370 ± 0.0061 | 0.0370 ± 0.0061 | same as previous report |
| Paper 1, Q2 | 4C / s=2: rate (ratio) / {'mode': 'quenched'} | 0.0332 ± 0.0070 (0.92) | 0.0332 ± 0.0070 (0.92) | same as previous report |
| Paper 1, Q2 | min–max of the s = 1.25…10 means, 4C, quenched | derived from previous-report values | 0.0307–0.0381 | derived from previous report |
| Paper 1, Q2 | 4C / pred a²/2v / {'mode': 'cosine'} | 0.0125 ± 0.0000 | 0.0125 ± n/a (deterministic) | same as previous report |
| Paper 1, Q2 | 4C / s=2: rate (ratio) / {'mode': 'cosine'} | 0.0210 ± 0.0000 (1.68) | 0.0210 ± n/a (deterministic) (1.68) | same as previous report |
| Paper 1, Q2 | min–max of the s = 1.25…10 means, 4C, cosine | derived from previous-report values | 0.0183–0.0218 | derived from previous report |
| Paper 1, Q2 | 4C_seed0 / predicted / {'mode__(seed__0)': 'annealed'} | 0.0724 | 0.0724 | same as previous report |
| Paper 1, Q2 | 4C_seed0 / measured s=2 / {'mode__(seed__0)': 'annealed'} | 0.1088 | 0.1088 | same as previous report |
| Paper 1, Q2 | 4C_seed0 / predicted / {'mode__(seed__0)': 'quenched'} | 0.0319 | 0.0319 | same as previous report |
| Paper 1, Q2 | 4C_seed0 / measured s=2 / {'mode__(seed__0)': 'quenched'} | 0.0266 | 0.0266 | same as previous report |
| Paper 1, Q2 | previous summary text | agree to 1e-4 | 10⁻⁴ | same as previous report |
| Paper 1, Q3 | iid / decay of P(S≤2) = I(0) / {'model': 'i.i.d. model (1e7)'} | 0.1090 | 0.1090 | same as previous report |
| Paper 1, Q3 | iid / decay of P(S≤2) = I(0) / {'model': '2-D annealed (10 seeds)'} | 0.1085 ± 0.0018 | 0.1085 ± 0.0018 | same as previous report |
| Paper 1, Q3 | iid / Gaussian a²/2v / {'model': 'i.i.d. model (1e7)'} | 0.0731 | 0.0731 | same as previous report |
| Paper 1, Q3 | iid / −min L(q): N 20→40 / 40→60 / {'model': 'i.i.d. model (1e7)'} | 0.1181 / 0.1011 | 0.1181 / 0.1011 | same as previous report |
| Paper 1, Q3 | iid / λ-CDF I(0+), N 40→60 → corrected for (ln S)² / {'model': 'i.i.d. model (1e7)'} | 0.0702 → 0.1107 | 0.0702 → 0.1107 | same as previous report |
| Paper 1, Q3 | previous summary text | measured P(1.25)/P(2) = 0.097 | 0.097 | same as previous report |
| Paper 1, Q3 | previous summary text | (ln 1.25/ln 2)² = 0.104 | 0.104 | same as previous report |
| Paper 1, Q4 | caustics / caustic growth (fit) / {'mode': 'annealed', 'launch__row': 'row0'} | 0.2417 ± 0.0085 | 0.2417 ± 0.0085 | same as previous report |
| Paper 1, Q4 | caustics / caustic growth (fit) / {'mode': 'quenched', 'launch__row': 'row0'} | 0.1983 ± 0.0166 | 0.1983 ± 0.0166 | same as previous report |
| Paper 1, Q4 | caustics / caustic growth (fit) / {'mode': 'cosine', 'launch__row': 'row0'} | 0.2066 ± 0.0000 | 0.2066 ± n/a (deterministic) | same as previous report |
| Paper 1, Q4 | caustics / alpha+beta/4 / {'mode': 'annealed', 'launch__row': 'row0'} | 0.2344 ± 0.0057 | 0.2344 ± 0.0057 | same as previous report |
| Paper 1, Q4 | caustics / alpha+beta/4 / {'mode': 'quenched', 'launch__row': 'row0'} | 0.1848 ± 0.0376 | 0.1848 ± 0.0376 | same as previous report |
| Paper 1, Q4 | caustics / alpha+beta/4 / {'mode': 'cosine', 'launch__row': 'row0'} | 0.0606 ± 0.0000 | 0.0606 ± n/a (deterministic) | same as previous report |
| Paper 1, Q4 | caustics / line-length growth (exact) / {'mode': 'annealed', 'launch__row': 'row0'} | 0.2423 ± 0.0039 | 0.2423 ± 0.0039 | same as previous report |
| Paper 1, Q4 | caustics / line-length growth (exact) / {'mode': 'quenched', 'launch__row': 'row0'} | 0.1919 ± 0.0174 | 0.1919 ± 0.0174 | same as previous report |
| Paper 1, Q4 | caustics / line-length growth (exact) / {'mode': 'cosine', 'launch__row': 'row0'} | 0.1756 ± 0.0000 | 0.1756 ± n/a (deterministic) | same as previous report |
| Paper 1, Q4 | caustics / fit N / {'mode': 'annealed', 'launch__row': 'row0'} | 10–31 | 10–31 | same as previous report |
| Paper 1, Q4 | caustics / ratio caustic/pred / {'mode': 'annealed', 'launch__row': 'row0'} | 1.03 | 1.03 | same as previous report |
| Paper 1, Q4 | caustics / ratio caustic/pred / {'mode': 'quenched', 'launch__row': 'row0'} | 1.11 | 1.11 | same as previous report |
| Paper 1, Q4 | caustics / ratio caustic/pred / {'mode': 'cosine', 'launch__row': 'row0'} | 3.41 | 3.41 (deterministic) | same as previous report |
| Paper 1, Q5 | 1C / low-S share of branch mass / {'mode': 'annealed', 'N': 5} | 0.3312 ± 0.0653 | 0.3312 ± 0.0653 | same as previous report |
| Paper 1, Q5 | 1C / low-S share of branch mass / {'mode': 'annealed', 'N': 20} | 0.0379 ± 0.0148 | 0.0379 ± 0.0148 | same as previous report |
| Paper 1, Q5 | 1C / low-S share of branch mass / {'mode': 'annealed', 'N': 60} | 0.0006 ± 0.0007 | 0.0006 ± 0.0007 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'annealed', 'N': 5} | +0.197 ± +0.084 | 0.197 ± 0.084 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'annealed', 'N': 20} | +0.040 ± +0.015 | 0.040 ± 0.015 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'annealed', 'N': 60} | +0.001 ± +0.001 | 0.001 ± 0.001 | same as previous report |
| Paper 1, Q5 | 1C / corr(ρ, ρ_noLow) / {'mode': 'annealed', 'N': 5} | 0.858 ± 0.039 | 0.858 ± 0.039 | same as previous report |
| Paper 1, Q5 | 1C / corr(ρ, ρ_noLow) / {'mode': 'annealed', 'N': 20} | 0.989 ± 0.007 | 0.989 ± 0.007 | same as previous report |
| Paper 1, Q5 | 1C / corr(ρ, ρ_noLow) / {'mode': 'annealed', 'N': 60} | 1.000 ± 0.000 | 1.000 ± < 0.001 | same as previous report |
| Paper 1, Q5 | 1C / enrichment / {'mode': 'annealed', 'N': 10} | 1.94 ± 0.31 | 1.94 ± 0.31 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'quenched', 'N': 20} | +0.122 ± +0.022 | 0.122 ± 0.022 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'quenched', 'N': 40} | +0.117 ± +0.040 | 0.117 ± 0.040 | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'quenched', 'N': 60} | +0.109 ± +0.040 | 0.109 ± 0.040 | same as previous report |
| Paper 1, Q5 | 1C / low-S share of branch mass / {'mode': 'cosine', 'N': 20} | 0.4130 ± 0.0000 | 0.4130 ± n/a (deterministic) | same as previous report |
| Paper 1, Q5 | 1C / contrast lost / {'mode': 'cosine', 'N': 20} | +0.414 ± +0.000 | 0.414 ± n/a (deterministic) | same as previous report |
| Paper 1, Q5 | 1C / corr(ρ, ρ_noLow) / {'mode': 'cosine', 'N': 20} | 0.668 ± 0.000 | 0.668 ± n/a (deterministic) | same as previous report |
| Paper 1, Q5 | 1C / null lost (S / r) / {'mode': 'annealed', 'N': 20} | +0.001 ± +0.001 / +0.001 ± +0.003 | 0.001 ± 0.001 / +0.001 ± +0.003 | same as previous report |
| Paper 1, Q5 | 1C / branch overlap (Jaccard) / {'mode': 'annealed', 'N': 20} | 0.90 ± 0.03 | 0.90 ± 0.03 | same as previous report |
| Paper 1, Q6 | 4E / alpha / {'mode': 'annealed', 'y__scale__f': '1'} | 0.1526 ± 0.0022 | 0.1526 ± 0.0022 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=1.5 / {'mode': 'annealed', 'y__scale__f': '1'} | 0.1093 ± 0.0021 | 0.1093 ± 0.0021 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'annealed', 'y__scale__f': '1'} | 0.1085 ± 0.0018 | 0.1085 ± 0.0018 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=3 / {'mode': 'annealed', 'y__scale__f': '1'} | 0.1068 ± 0.0018 | 0.1068 ± 0.0018 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'annealed', 'y__scale__f': '0.5'} | 0.1066 ± 0.0018 | 0.1066 ± 0.0018 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'annealed', 'y__scale__f': '2'} | 0.1079 ± 0.0022 | 0.1079 ± 0.0022 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'annealed', 'y__scale__f': '4'} | 0.1055 ± 0.0019 | 0.1055 ± 0.0019 | same as previous report |
| Paper 1, Q6 | largest seed sd among the decay rates, 4E, annealed | derived from previous-report values | 0.0022 | derived from previous report |
| Paper 1, Q6 | 4E / alpha / {'mode': 'quenched', 'y__scale__f': '1'} | 0.1452 ± 0.0094 | 0.1452 ± 0.0094 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=1.5 / {'mode': 'quenched', 'y__scale__f': '1'} | 0.0316 ± 0.0066 | 0.0316 ± 0.0066 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'quenched', 'y__scale__f': '1'} | 0.0332 ± 0.0070 | 0.0332 ± 0.0070 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=3 / {'mode': 'quenched', 'y__scale__f': '1'} | 0.0348 ± 0.0074 | 0.0348 ± 0.0074 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'quenched', 'y__scale__f': '0.5'} | 0.0317 ± 0.0065 | 0.0317 ± 0.0065 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'quenched', 'y__scale__f': '2'} | 0.0356 ± 0.0078 | 0.0356 ± 0.0078 | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'quenched', 'y__scale__f': '4'} | 0.0373 ± 0.0081 | 0.0373 ± 0.0081 | same as previous report |
| Paper 1, Q6 | largest seed sd among the decay rates, 4E, quenched | derived from previous-report values | 0.0081 | derived from previous report |
| Paper 1, Q6 | 4E / alpha / {'mode': 'cosine', 'y__scale__f': '1'} | 0.0236 ± 0.0000 | 0.0236 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=1.5 / {'mode': 'cosine', 'y__scale__f': '1'} | 0.0218 ± 0.0000 | 0.0218 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'cosine', 'y__scale__f': '1'} | 0.0210 ± 0.0000 | 0.0210 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=3 / {'mode': 'cosine', 'y__scale__f': '1'} | 0.0195 ± 0.0000 | 0.0195 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'cosine', 'y__scale__f': '0.5'} | 0.0225 ± 0.0000 | 0.0225 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'cosine', 'y__scale__f': '2'} | 0.0142 ± 0.0000 | 0.0142 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | 4E / rate S_THR=2 / {'mode': 'cosine', 'y__scale__f': '4'} | 0.0137 ± 0.0000 | 0.0137 ± n/a (deterministic) | same as previous report |
| Paper 1, Q6 | max /rate − rate(S≤2, unscaled)/ / rate, annealed, from 4E | derived from previous-report values | 3% | derived from previous report |
| Paper 1, Q6 | 4E / P(N=60) / {'mode': 'annealed', 'y__scale__f': '1'} | 6.65e-05 | 6.65e-05 | same as previous report |
| Paper 1, Q6 | 4E / P(N=60) / {'mode': 'annealed', 'y__scale__f': '4'} | 5.00e-05 | 5.00e-05 | same as previous report |
| Paper 2, Q1 | KM / F(1) / {'mode': 'annealed', 'direction': 're-entry'} | 0.121 ± 0.010 | 0.121 ± 0.010 | same as previous report |
| Paper 2, Q1 | KM / F(1) / {'mode': 'quenched', 'direction': 're-entry'} | 0.128 ± 0.016 | 0.128 ± 0.016 | same as previous report |
| Paper 2, Q1 | KM / F(1) / {'mode': 'cosine', 'direction': 're-entry'} | 0.085 ± 0.000 | 0.085 ± n/a (deterministic) | same as previous report |
| Paper 2, Q1 | KM / F(5) / {'mode': 'annealed', 'direction': 're-entry'} | 0.226 ± 0.014 | 0.226 ± 0.014 | same as previous report |
| Paper 2, Q1 | KM / F(5) / {'mode': 'quenched', 'direction': 're-entry'} | 0.257 ± 0.027 | 0.257 ± 0.027 | same as previous report |
| Paper 2, Q1 | KM / F(5) / {'mode': 'cosine', 'direction': 're-entry'} | 0.334 ± 0.000 | 0.334 ± n/a (deterministic) | same as previous report |
| Paper 2, Q1 | KM / whole run / {'mode': 'annealed', 'direction': 're-entry'} | 0.283 ± 0.015 (k=59) | 0.283 ± 0.015 (k=59) | same as previous report |
| Paper 2, Q1 | KM / whole run / {'mode': 'quenched', 'direction': 're-entry'} | 0.330 ± 0.026 (k=59) | 0.330 ± 0.026 (k=59) | same as previous report |
| Paper 2, Q1 | KM / whole run / {'mode': 'cosine', 'direction': 're-entry'} | 0.366 ± 0.000 (k=59) | 0.366 ± n/a (deterministic) (k=59) | same as previous report |
| Paper 2, Q1 | KM / naive / {'mode': 'annealed', 'direction': 're-entry'} | 0.283 ± 0.015 | 0.283 ± 0.015 | same as previous report |
| Paper 2, Q1 | KM / naive / {'mode': 'quenched', 'direction': 're-entry'} | 0.328 ± 0.026 | 0.328 ± 0.026 | same as previous report |
| Paper 2, Q1 | KM / naive / {'mode': 'cosine', 'direction': 're-entry'} | 0.359 ± 0.000 | 0.359 ± n/a (deterministic) | same as previous report |
| Paper 2, Q1 | F(5)/whole run, re-entry, annealed | derived from previous-report values | 80% | derived from previous report |
| Paper 2, Q1 | F(5)/whole run, re-entry, quenched | derived from previous-report values | 78% | derived from previous report |
| Paper 2, Q1 | F(5)/whole run, re-entry, cosine | derived from previous-report values | 91% | derived from previous report |
| Paper 2, Q1 | pooled / at risk (mean) / {'mode': 'annealed', 'age': '1'} | 1284674 | 1284674 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'annealed', 'age': '1'} | 0.1215 ± 0.0100 | 0.1215 ± 0.0100 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'annealed', 'age': '2'} | 0.0420 ± 0.0028 | 0.0420 ± 0.0028 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'annealed', 'age': '6–10'} | 0.0093 ± 0.0006 | 0.0093 ± 0.0006 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'annealed', 'age': '21–40'} | 0.0003 ± 0.0000 | 0.0003 ± < 0.0001 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'annealed', 'age': '>40'} | 0.0000 ± 0.0000 | 0.0000 ± < 0.0001 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'quenched', 'age': '1'} | 0.1276 ± 0.0158 | 0.1276 ± 0.0158 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'quenched', 'age': '2'} | 0.0572 ± 0.0047 | 0.0572 ± 0.0047 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'quenched', 'age': '6–10'} | 0.0117 ± 0.0017 | 0.0117 ± 0.0017 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'quenched', 'age': '21–40'} | 0.0006 ± 0.0001 | 0.0006 ± 0.0001 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'quenched', 'age': '>40'} | 0.0002 ± 0.0000 | 0.0002 ± < 0.0001 | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'cosine', 'age': '1'} | 0.0847 ± 0.0000 | 0.0847 ± n/a (deterministic) | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'cosine', 'age': '2'} | 0.1036 ± 0.0000 | 0.1036 ± n/a (deterministic) | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'cosine', 'age': '6–10'} | 0.0035 ± 0.0000 | 0.0035 ± n/a (deterministic) | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'cosine', 'age': '21–40'} | 0.0005 ± 0.0000 | 0.0005 ± n/a (deterministic) | same as previous report |
| Paper 2, Q2 | pooled / h_re / {'mode': 'cosine', 'age': '>40'} | 0.0011 ± 0.0000 | 0.0011 ± n/a (deterministic) | same as previous report |
| Paper 2, Q2 | age-1 return hazard, 10th percentile over exit cycles (risk ≥ 200), annealed, seeds pooled | (not in previous report) | 0.106 | NEW |
| Paper 2, Q2 | age-1 return hazard, 90th percentile over exit cycles (risk ≥ 200), annealed, seeds pooled | (not in previous report) | 0.125 | NEW |
| Paper 2, Q3 | KM / F(1) / {'mode': 'annealed', 'direction': 'exit (entrants)'} | 0.401 ± 0.010 | 0.401 ± 0.010 | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'annealed', 'direction': 'exit (entrants)'} | 0.612 ± 0.009 | 0.612 ± 0.009 | same as previous report |
| Paper 2, Q3 | KM / F(10) / {'mode': 'annealed', 'direction': 'exit (entrants)'} | 0.992 ± 0.001 | 0.992 ± 0.001 | same as previous report |
| Paper 2, Q3 | KM / F(1) / {'mode': 'annealed', 'direction': 'exit (orig. cohort)'} | 0.116 ± 0.019 | 0.116 ± 0.019 | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'annealed', 'direction': 'exit (orig. cohort)'} | 0.395 ± 0.022 | 0.395 ± 0.022 | same as previous report |
| Paper 2, Q3 | KM / F(10) / {'mode': 'annealed', 'direction': 'exit (orig. cohort)'} | 0.989 ± 0.002 | 0.989 ± 0.002 | same as previous report |
| Paper 2, Q3 | KM / F(1) / {'mode': 'quenched', 'direction': 'exit (entrants)'} | 0.388 ± 0.023 | 0.388 ± 0.023 | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'quenched', 'direction': 'exit (entrants)'} | 0.598 ± 0.037 | 0.598 ± 0.037 | same as previous report |
| Paper 2, Q3 | KM / F(10) / {'mode': 'quenched', 'direction': 'exit (entrants)'} | 0.983 ± 0.013 | 0.983 ± 0.013 | same as previous report |
| Paper 2, Q3 | KM / F(1) / {'mode': 'quenched', 'direction': 'exit (orig. cohort)'} | 0.116 ± 0.019 | 0.116 ± 0.019 | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'quenched', 'direction': 'exit (orig. cohort)'} | 0.395 ± 0.019 | 0.395 ± 0.019 | same as previous report |
| Paper 2, Q3 | KM / F(10) / {'mode': 'quenched', 'direction': 'exit (orig. cohort)'} | 0.978 ± 0.011 | 0.978 ± 0.011 | same as previous report |
| Paper 2, Q3 | KM / F(1) / {'mode': 'cosine', 'direction': 'exit (entrants)'} | 0.451 ± 0.000 | 0.451 ± n/a (deterministic) | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'cosine', 'direction': 'exit (entrants)'} | 0.821 ± 0.000 | 0.821 ± n/a (deterministic) | same as previous report |
| Paper 2, Q3 | KM / F(1) / {'mode': 'cosine', 'direction': 'exit (orig. cohort)'} | 0.123 ± 0.000 | 0.123 ± n/a (deterministic) | same as previous report |
| Paper 2, Q3 | KM / F(2) / {'mode': 'cosine', 'direction': 'exit (orig. cohort)'} | 0.355 ± 0.000 | 0.355 ± n/a (deterministic) | same as previous report |
| Paper 2, Q3 | KM / F(10) / {'mode': 'cosine', 'direction': 'exit (orig. cohort)'} | 0.875 ± 0.000 | 0.875 ± n/a (deterministic) | same as previous report |
| Paper 2, Q3 | largest seed sd among the leaving probabilities shown, KM table | derived from previous-report values | 0.037 | derived from previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'annealed', 'stratum': 'all depths (crude)'} | 1.16 ± 0.08 | 1.16 ± 0.08 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'annealed', 'stratum': 'MH, age-adjusted'} | 0.98 ± 0.06 | 0.98 ± 0.06 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'annealed', 'stratum': 'MH, depth-adjusted'} | 0.94 ± 0.04 | 0.94 ± 0.04 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'annealed', 'stratum': 'MH, depth × age'} | 0.94 ± 0.04 | 0.94 ± 0.04 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'quenched', 'stratum': 'all depths (crude)'} | 1.77 ± 0.40 | 1.77 ± 0.40 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'quenched', 'stratum': 'MH, age-adjusted'} | 1.34 ± 0.24 | 1.34 ± 0.24 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'quenched', 'stratum': 'MH, depth-adjusted'} | 1.03 ± 0.17 | 1.03 ± 0.17 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'quenched', 'stratum': 'MH, depth × age'} | 1.03 ± 0.15 | 1.03 ± 0.15 | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'all depths (crude)'} | 2.57 ± 0.00 | 2.57 ± n/a (deterministic) | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'MH, age-adjusted'} | 1.60 ± 0.00 | 1.60 ± n/a (deterministic) | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'MH, depth-adjusted'} | 0.90 ± 0.00 | 0.90 ± n/a (deterministic) | same as previous report |
| Paper 2, Q4 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'MH, depth × age'} | 1.09 ± 0.00 | 1.09 ± n/a (deterministic) | same as previous report |
| Paper 2, Q4 | previous summary text | earlier 1.8× | 1.8× | same as previous report |
| Paper 2, Q4 | previous summary text | earlier 2.6× | 2.6× | same as previous report |
| Paper 2, Q4 | annealed (crude − 1)/sd, depth table | derived from previous-report values | 2 | derived from previous report |
| Paper 2, Q5 | depth / h(returned before) / h(first) / {'mode': 'annealed', 'stratum': 'S_{N-1} (2,2.5]'} | 0.93 ± 0.05 | 0.93 ± 0.05 | same as previous report |
| Paper 2, Q5 | depth / h(returned before) / h(first) / {'mode': 'quenched', 'stratum': 'S_{N-1} (2,2.5]'} | 1.01 ± 0.15 | 1.01 ± 0.15 | same as previous report |
| Paper 2, Q5 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'S_{N-1} (2,2.5]'} | 0.91 ± 0.00 | 0.91 ± n/a (deterministic) | same as previous report |
| Paper 2, Q5 | depth / h(returned before) / h(first) / {'mode': 'cosine', 'stratum': 'S_{N-1} (3,5]'} | 0.43 ± 0.00 | 0.43 ± n/a (deterministic) | same as previous report |
| Paper 2, Q5 | depth / events/at risk returned (per seed) / {'mode': 'annealed', 'stratum': 'S_{N-1} >5'} | 40/11293294 | 40/11293294 | same as previous report |
| Paper 2, Q5 | depth / events/at risk first (per seed) / {'mode': 'annealed', 'stratum': 'S_{N-1} >5'} | 105/34455341 | 105/34455341 | same as previous report |
