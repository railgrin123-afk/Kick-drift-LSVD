Reference scripts for the band / caustic / cusp / branch figures (Python).
They use the same map and kicks as the MATLAB script (seed 0, 1920 x 1920 frame).
kd_engine.py is the engine from the earlier python/ folder; core.py imports it.

core.py               starting picture, kick-drift map, exact pullback picture, S, overlay
compute.py            band amounts per 4-px strip, every cycle, three modes -> stats.npz
figs.py               Figures 1-5 (picture, bands -> branched flow, low-stretch, shares, decay)
peaks.py              "what a peak is": picture + red band amount curve, cycle 3
branch_in_picture.py  band material in a branch coloured cyan, annealed vs quenched
cuspshow.py           exact cusp moment (cycle 2.44, x ~ 765, red band) with in-between moments
branches.py           automatic branch outlines (>= 2x, >= 1 cycle) and lifetimes
Run compute.py first, then e.g.:  python figs.py "fig1()"
