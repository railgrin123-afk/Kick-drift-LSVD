# Python companion to `Kick-drift LSVD`

Fast multi-seed versions of Parts 1, 2 and 4 (launched lines and caustics,
crossing ages, distribution of S). Same map, conventions and random stream as
the MATLAB script; seed 0 is the MATLAB `rng(0,'twister')` run.

    python3 kd_run.py                 # 3 modes x seeds 0-9 + i.i.d. model -> out/  (~25 min, 3 workers)
    python3 kd_run.py --only-part1c   # just Part 1C (seconds per run)
    python3 kd_report.py              # tables + figures -> ../reports/
    python3 kd_findings.py            # summary / method / caveats text
    python3 kd_pdf.py                 # ../reports/kd_report.pdf

`kd_engine.run(Params(...), kicks=dp_x_hist)` reruns any MATLAB run on its own
kicks for a direct comparison. Needs numpy, scipy, matplotlib, reportlab.
