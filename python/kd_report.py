"""
kd_report -- aggregate the multi-seed runs: tables, PASS/FAIL checks, figures, PDF.

    python3 kd_report.py [--out python/out] [--report reports]

Reads <out>/run_<mode>_s<seed>.pkl.gz and <out>/iid.pkl.gz (from kd_run.py) and
<out>/matlab_checks.json if present (written by the Octave verification of the
MATLAB port). Writes <report>/kd_report.pdf, <report>/kd_tables.md and the
figure PNGs in <report>/fig/.
"""
import argparse, glob, json, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import kd_engine as K
from kd_run import load

MODES = ['annealed', 'quenched', 'cosine']
COL = {'annealed': '#2a78d6', 'quenched': '#eb6834', 'cosine': '#1baf7a'}
BANDC = ['#2a78d6', '#eda100', '#e34948']          # S<=2, 2<S<=5, S>5
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
S_FIT = [1.25, 1.5, 2, 3, 5, 10]
MIN_RISK = 200

plt.rcParams.update({
    'font.size': 8.5, 'axes.titlesize': 9, 'axes.labelsize': 8.5, 'legend.fontsize': 7.5,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'lines.linewidth': 1.6, 'figure.dpi': 110, 'savefig.dpi': 160,
    'savefig.bbox': 'tight', 'legend.frameon': False})


# ============================================================ helpers
def ms(x, axis=0):
    x = np.asarray(x, float)
    m = np.nanmean(x, axis=axis)
    s = np.nanstd(x, axis=axis, ddof=1) if x.shape[axis] > 1 else np.zeros_like(m)
    return m, s


def pm(m, s, f='{:.4f}'):
    return (f + ' ± ' + f).format(m, s)


def slope(N, y):
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    return np.polyfit(N[ok], y[ok], 1)[0]


def band(ax, x, Y, color, label=None, ls='-', lw=1.6, logy=False):
    m, s = ms(Y)
    ax.plot(x, m, color=color, ls=ls, lw=lw, label=label)
    if Y.shape[0] > 1 and np.any(s > 0):
        lo, hi = m - s, m + s
        if logy:
            lo = np.where(lo > 0, lo, np.nan)
        ax.fill_between(x, lo, hi, color=color, alpha=0.18, lw=0)


# ============================================================ per-run analysis
def fit_4c(o, Nmin=20, key='P_le', mom=None):
    C = o['cycles']; N = np.arange(C + 1); sel = N >= Nmin
    m = o['mom'] if mom is None else mom
    al = slope(N[sel], m[sel, 0]); v = slope(N[sel], m[sel, 1])
    pred = al ** 2 / (2 * v)
    meas = []
    for s in S_FIT:
        col = int(np.nonzero(K.S_LIST == s)[0][0])
        P = o[key][:, col]
        cnt = P * o['n_win']
        ok = sel & (cnt >= 10)
        meas.append(-slope(N[ok], np.log(P[ok])) if ok.sum() >= 5 else np.nan)
    return al, v, pred, np.array(meas)


def caustic_fit(o, ri=1, Nmin=10):
    C = o['cycles']; N = np.arange(C + 1)
    c = o['caustics_exact'][ri]
    ok = (N >= Nmin) & np.isfinite(c) & (c > 0)
    rate = slope(N[ok], np.log(c[ok])) if ok.sum() >= 5 else np.nan
    sel = N >= 10
    al = slope(N[sel], o['line_mean_r'][ri, sel]); be = 2 * slope(N[sel], o['line_var_r'][ri, sel])
    Nmax = int(N[ok].max()) if ok.any() else 0
    p = o['exact_pieces'][ri]
    okp = ok & (p > 0)
    # pieces ~ total variation of x_N(x0) = integral of |a| dx0: the exact
    # growth of the projected line length, which Kaplan's alpha+beta/4 models
    lrate = slope(N[okp], np.log(p[okp])) if okp.sum() >= 5 else np.nan
    return rate, al, be, al + be / 4, Nmax, lrate


def part2(o):
    C = o['cycles']
    r_re, e_re = o['risk_re'].sum(0), o['ev_re'].sum(0)            # by age
    r_ex, e_ex = o['risk_ex'][1:].sum(0), o['ev_ex'][1:].sum(0)    # entrants only
    r_co, e_co = o['risk_ex'][0], o['ev_ex'][0]                    # original cohort
    h_re = K.hazard(e_re, r_re, MIN_RISK); h_ex = K.hazard(e_ex, r_ex, MIN_RISK)
    h_co = K.hazard(e_co, r_co, MIN_RISK)
    F_re, F_ex, F_co = K.km(h_re), K.km(h_ex), K.km(h_co)

    def whole(F, r):
        ks = np.nonzero((r >= MIN_RISK) & np.isfinite(F))[0]
        return (F[ks.max()], int(ks.max())) if ks.size else (np.nan, 0)
    naive_re = e_re.sum() / r_re[1] if r_re[1] else np.nan          # exit spells that came back
    naive_ex = e_ex.sum() / r_ex[1] if r_ex[1] else np.nan          # entry spells that left
    naive_co = e_co.sum() / r_co[1] if r_co[1] else np.nan
    return dict(h_re=h_re, h_ex=h_ex, h_co=h_co, r_re=r_re, r_ex=r_ex, r_co=r_co,
                F_re=F_re, F_ex=F_ex, F_co=F_co, W_re=whole(F_re, r_re), W_ex=whole(F_ex, r_ex),
                W_co=whole(F_co, r_co), naive_re=naive_re, naive_ex=naive_ex, naive_co=naive_co,
                pool_re=K.pooled_bins(e_re, r_re), pool_ex=K.pooled_bins(e_ex, r_ex),
                type1=o['type1'][0] / o['type1'][1], age1=e_re[1] / r_re[1])


def depth_ratios(o):
    R, E = o['risk_rda'], o['ev_rda']                  # [rep, depth7, age]
    res = dict(user=[], fine=[])
    for groups, key in ((K.DEPTH_USER, 'user'), ([[i] for i in range(R.shape[1])], 'fine')):
        for g in groups:
            r0, e0 = R[0, g].sum(), E[0, g].sum(); r1, e1 = R[1, g].sum(), E[1, g].sum()
            h0 = e0 / r0 if r0 else np.nan; h1 = e1 / r1 if r1 else np.nan
            res[key].append((h1 / h0 if h0 > 0 else np.nan, h1, h0, e1, r1, e0, r0))
    # pooled (crude) and Mantel-Haenszel rate ratio over depth7 x age strata
    r0, e0, r1, e1 = R[0], E[0], R[1], E[1]
    res['crude'] = (e1.sum() / r1.sum()) / (e0.sum() / r0.sum())
    T = r0 + r1
    with np.errstate(invalid='ignore', divide='ignore'):
        num = np.nansum(np.where(T > 0, e1 * r0 / T, 0)); den = np.nansum(np.where(T > 0, e0 * r1 / T, 0))
    res['MH'] = num / den if den > 0 else np.nan
    # MH over depth only (age pooled) and over age only
    for nm, ax in (('MH_depth', 1), ('MH_age', 0)):   # strata = depth (sum over age) | age
        a0, b0, a1, b1 = r0.sum(ax), e0.sum(ax), r1.sum(ax), e1.sum(ax)
        TT = a0 + a1
        with np.errstate(invalid='ignore', divide='ignore'):
            res[nm] = np.nansum(np.where(TT > 0, b1 * a0 / TT, 0)) / np.nansum(np.where(TT > 0, b0 * a1 / TT, 0))
    return res


def rho_check(o, Nmax_check=None, excl=2, min_count=200):
    """1/|a| vs histogram on the dense check line, bins away from caustics."""
    C = o['cycles']; rels = []; perN = {}
    for N in range(1, C + 1):
        ra = o['chk_rho_a'][N]; rh = o['chk_rho_hist'][N]; cb = o['chk_caus_bin'][N]
        cbd = cb.copy()
        for s in range(1, excl + 1):
            cbd |= np.roll(cb, s) | np.roll(cb, -s)
        ok = np.isfinite(ra) & ~cbd & (o['chk_count'][N] >= min_count)
        if ok.sum() >= 10:
            rel = np.abs(ra[ok] / rh[ok] - 1); rels.append(rel)
            perN[N] = (ok.sum(), np.median(rel), np.mean(rel < 0.02))
    allr = np.concatenate(rels) if rels else np.array([np.nan])
    return dict(n=allr.size, median=np.nanmedian(allr), frac2=np.mean(allr < 0.02), perN=perN)


# ============================================================ main
def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument('--out', default=os.path.join(here, 'out'))
    ap.add_argument('--report', default=os.path.join(os.path.dirname(here), 'reports'))
    a = ap.parse_args()
    figdir = os.path.join(a.report, 'fig'); os.makedirs(figdir, exist_ok=True)

    R = {m: [] for m in MODES}
    for p in sorted(glob.glob(os.path.join(a.out, 'run_*_s*.pkl.gz')),
                    key=lambda s: (s, int(re.search(r'_s(\d+)', s).group(1)))):
        m, s = re.search(r'run_(\w+?)_s(\d+)', os.path.basename(p)).groups()
        R[m].append(load(p))
    for m in MODES:
        R[m].sort(key=lambda o: o['seed'])
    iid = load(os.path.join(a.out, 'iid.pkl.gz'))
    mchk_path = os.path.join(a.out, 'matlab_checks.json')
    mchk = json.load(open(mchk_path)) if os.path.exists(mchk_path) else {}
    C = R['annealed'][0]['cycles']; N = np.arange(C + 1)
    nseed = {m: len(R[m]) for m in MODES}
    print('seeds:', nseed)

    T = {}           # all tables, as lists of rows
    figs = []        # (path, caption)

    def savefig(fig, name, cap):
        p = os.path.join(figdir, name + '.png'); fig.savefig(p); plt.close(fig); figs.append((p, cap))

    # ---------------------------------------------------------------- 4C (up front)
    rows = []
    fit4c = {}
    for m in MODES:
        F = [fit_4c(o) for o in R[m]]
        al = np.array([f[0] for f in F]); v = np.array([f[1] for f in F]); pr = np.array([f[2] for f in F])
        me = np.array([f[3] for f in F])
        fit4c[m] = (al, v, pr, me)
        rows.append([m, pm(*ms(al)), pm(*ms(v)), pm(*ms(pr))] +
                    ['%s (%.2f)' % (pm(*ms(me[:, i])), np.nanmean(me[:, i] / pr)) for i in range(len(S_FIT))])
    T['4C'] = (['mode', 'alpha', 'v', 'pred a²/2v'] + ['s=%g: rate (ratio)' % s for s in S_FIT], rows)
    rows = []
    for m in MODES:
        al, v, pr, me = (x[0] for x in fit4c[m])
        rows.append([m, '%.4f' % pr, '%.4f' % me[S_FIT.index(2)],
                     '%.4f–%.4f' % (np.nanmin(me), np.nanmax(me))])
    T['4C_seed0'] = (['mode (seed 0)', 'predicted', 'measured s=2', 'measured s=1.25..10'], rows)

    # ---------------------------------------------------------------- 1A caustics
    rows = []
    cf = {}
    for m in MODES:
        X = np.array([[caustic_fit(o, ri) for ri in range(3)] for o in R[m]])   # seeds x rows x 6
        cf[m] = X
        for ri, lab in enumerate(['row0-rows/8', 'row0', 'row0+rows/8']):
            rows.append([m, lab, pm(*ms(X[:, ri, 0])), pm(*ms(X[:, ri, 5])), pm(*ms(X[:, ri, 3])),
                         '%.2f' % np.nanmean(X[:, ri, 0] / X[:, ri, 3]),
                         pm(*ms(X[:, ri, 1])), pm(*ms(X[:, ri, 2])), '10–%d' % np.nanmin(X[:, ri, 4])])
    T['caustics'] = (['mode', 'launch row', 'caustic growth (fit)', 'line-length growth (exact)', 'alpha+beta/4', 'ratio caustic/pred',
                      'alpha', 'beta', 'fit N'], rows)

    # 1A(iv) correlations
    rows = []
    for m in MODES:
        for Nn in (5, 10, 20, 40, 60):
            cl = np.array([np.nanmean(o['line_corr_low'][:, Nn, 0]) for o in R[m]])
            cr = np.array([np.nanmean(o['line_corr_negr'][:, Nn, 0]) for o in R[m]])
            sl = np.array([np.nanmean(o['line_corr_low'][:, Nn, 1]) for o in R[m]])
            sr = np.array([np.nanmean(o['line_corr_negr'][:, Nn, 1]) for o in R[m]])
            rows.append([m, Nn, pm(*ms(cl), '{:+.3f}'), pm(*ms(sl), '{:+.3f}'),
                         pm(*ms(cr), '{:+.3f}'), pm(*ms(sr), '{:+.3f}')])
    T['corr'] = (['mode', 'N', 'r(rho, frac S≤2) Pearson', 'Spearman', 'r(rho, mean −r) Pearson', 'Spearman'], rows)

    # ---------------------------------------------------------------- Part 2
    P2 = {m: [part2(o) for o in R[m]] for m in MODES}
    rows = []
    for m in MODES:
        for bi, (lo, hi) in enumerate(K.AGE_BINS):
            lab = str(lo) if lo == hi else ('%d–%d' % (lo, hi) if hi < 10 ** 5 else '>%d' % (lo - 1))
            hr = np.array([p['pool_re'][bi][0] for p in P2[m]]); nr = np.array([p['pool_re'][bi][1] for p in P2[m]])
            he = np.array([p['pool_ex'][bi][0] for p in P2[m]]); ne = np.array([p['pool_ex'][bi][1] for p in P2[m]])
            rows.append([m, lab, pm(*ms(hr)), '%.0f' % nr.mean(), pm(*ms(he)), '%.0f' % ne.mean()])
    T['pooled'] = (['mode', 'age', 'h_re', 'at risk (mean)', 'h_exit (entrants)', 'at risk (mean)'], rows)

    rows = []
    KS = [1, 2, 3, 5, 10, 20, 40]
    for m in MODES:
        for d, Fk, Wk, nk in (('re-entry', 'F_re', 'W_re', 'naive_re'), ('exit (entrants)', 'F_ex', 'W_ex', 'naive_ex'),
                              ('exit (orig. cohort)', 'F_co', 'W_co', 'naive_co')):
            vals = []
            for k in KS:
                x = np.array([p[Fk][k] if k < p[Fk].size else np.nan for p in P2[m]])
                vals.append(pm(*ms(x), '{:.3f}') if np.isfinite(x).any() else '–')
            W = np.array([p[Wk][0] for p in P2[m]]); Wn = np.array([p[Wk][1] for p in P2[m]])
            nv = np.array([p[nk] for p in P2[m]])
            rows.append([m, d] + vals + ['%s (k=%d)' % (pm(*ms(W), '{:.3f}'), int(np.median(Wn))),
                                         pm(*ms(nv), '{:.3f}')])
    T['KM'] = (['mode', 'direction'] + ['F(%d)' % k for k in KS] + ['whole run', 'naive'], rows)

    rows = []
    D = {m: [depth_ratios(o) for o in R[m]] for m in MODES}
    for m in MODES:
        cr = np.array([d['crude'] for d in D[m]]); mh = np.array([d['MH'] for d in D[m]])
        mhd = np.array([d['MH_depth'] for d in D[m]]); mha = np.array([d['MH_age'] for d in D[m]])
        rows.append([m, 'all depths (crude)', pm(*ms(cr), '{:.2f}'), '', ''])
        for i, lab in enumerate(K.DEPTH_USER_LABELS):
            x = np.array([d['user'][i][0] for d in D[m]])
            e1 = sum(d['user'][i][3] for d in D[m]); r1 = sum(d['user'][i][4] for d in D[m])
            e0 = sum(d['user'][i][5] for d in D[m]); r0 = sum(d['user'][i][6] for d in D[m])
            rows.append([m, 'S_{N-1} ' + lab, pm(*ms(x), '{:.2f}'),
                         '%d/%d' % (e1 / nseed[m], r1 / nseed[m]), '%d/%d' % (e0 / nseed[m], r0 / nseed[m])])
        rows.append([m, 'MH, depth-adjusted', pm(*ms(mhd), '{:.2f}'), '', ''])
        rows.append([m, 'MH, age-adjusted', pm(*ms(mha), '{:.2f}'), '', ''])
        rows.append([m, 'MH, depth × age', pm(*ms(mh), '{:.2f}'), '', ''])
    T['depth'] = (['mode', 'stratum', 'h(returned before) / h(first)', 'events/at risk returned (per seed)',
                   'events/at risk first (per seed)'], rows)

    # ---------------------------------------------------------------- checks
    checks = []
    for m in MODES:
        ok_re = all(o['ev_re'].sum() == o['n_ent'].sum() for o in R[m])
        ok_ex = all(o['ev_ex'].sum() == o['n_ext'].sum() for o in R[m])
        ok_t1 = all(abs(p['age1'] - p['type1']) < 1e-12 for p in P2[m])
        checks.append(('Σ age re-entry events = Σ re-entries (haz_reent·count) [%s, %d seeds]' % (m, nseed[m]), ok_re,
                       '%d = %d (seed 0)' % (R[m][0]['ev_re'].sum(), R[m][0]['n_ent'].sum())))
        checks.append(('Σ age exit events = Σ exits (haz_gross·count) [%s]' % m, ok_ex,
                       '%d = %d (seed 0)' % (R[m][0]['ev_ex'].sum(), R[m][0]['n_ext'].sum())))
        checks.append(('pooled age-1 re-entry hazard = Type 1 (consecutive masks) [%s]' % m, ok_t1,
                       '%.6f = %.6f (seed 0)' % (P2[m][0]['age1'], P2[m][0]['type1'])))
        rc = rho_check(R[m][0])
        checks.append(('rho: 1/|a| vs histogram within 2%% away from caustics [%s, seed 0]' % m,
                       rc['frac2'] >= 0.95, '%d bins, median |rel| %.4f, %.1f%% within 2%%' % (
                           rc['n'], rc['median'], 100 * rc['frac2'])))
        v = sum(o['cond_r'][Nn]['violations_S_lt_er'] for o in R[m] for Nn in o['cond_r'])
        checks.append(('identity S ≥ e^r on every line parcel [%s]' % m, v == 0, '%d violations' % v))
    for k, v in mchk.items():
        checks.append((k, bool(v[0]), v[1]))
    T['checks'] = checks

    # ---------------------------------------------------------------- 4A tables
    rows = []
    for m in MODES:
        for Nn in (10, 20, 40, 60):
            P = np.array([o['P_le'][Nn] for o in R[m]])
            mm, ss = ms(P)
            rows.append([m, Nn] + ['%.2e' % x if x < 0.01 else '%.4f' % x for x in mm])
    T['Ple'] = (['mode', 'N'] + ['s=%g' % s for s in K.S_LIST], rows)
    rows = []
    for m in MODES:
        for Nn in (10, 30, 60):
            M = np.array([o['mom'][Nn] for o in R[m]]); Q = np.array([o['quant'][Nn] for o in R[m]])
            rows.append([m, Nn] + [pm(*ms(M[:, i]), '{:.3f}') for i in range(4)] +
                        ['%.2f' % x for x in Q.mean(0)])
    T['mom'] = (['mode', 'N', 'mean ln S', 'var ln S', 'skew', 'ex. kurt'] +
                ['q%d' % int(round(100 * q)) for q in K.QUANTS], rows)

    # ---------------------------------------------------------------- 4D
    n_iid = iid['n_samples']
    lse = {Nn: iid['Lq'][Nn] * Nn + np.log(n_iid) for Nn in iid['Lq']}
    Lq_diff = (lse[60] - lse[40]) / 20.0
    qs = iid['qs']
    Pi = iid['P_le']
    col2 = int(np.nonzero(K.S_LIST == 2)[0][0])
    sel = N >= 20
    al_i_pre = slope(N[sel], iid['meanlnS'][sel])
    dec_iid = [-slope(N[sel], np.log(Pi[sel, int(np.nonzero(K.S_LIST == s)[0][0])])) for s in S_FIT]
    dec_2d = np.array([fit_4c(o)[3] for o in R['annealed']])
    lam_e = iid['lam_edges']; lam_c = 0.5 * (lam_e[:-1] + lam_e[1:])
    dl = lam_e[1] - lam_e[0]
    Irate = {}
    for Nn in (20, 40, 60):
        h = iid['lam_hist'][Nn] / (n_iid * dl)
        Irate[Nn] = np.where(h > 0, -np.log(np.where(h > 0, h, 1)) / Nn, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        h40 = iid['lam_hist'][40] / (n_iid * dl); h60 = iid['lam_hist'][60] / (n_iid * dl)
        I_diff = np.where((h40 > 0) & (h60 > 0), -(np.log(h60) - np.log(h40)) / 20, np.nan)
    # I(0+) from the lambda-CDF, difference N 40 -> 60. Near S = 1, P(S <= s) ~ (ln s)^2
    # (sigma_max = 1 is a codimension-2 set), so at FIXED lambda P carries an extra
    # (N lambda)^2 and the finite-N slope is low by 2 ln(60/40)/20; corrected below.
    F40 = np.cumsum(iid['lam_hist'][40]) / n_iid; F60 = np.cumsum(iid['lam_hist'][60]) / n_iid
    I0_cdf = -(np.log(F60[0]) - np.log(F40[0])) / 20
    I0_cdf_corr = I0_cdf + 2 * np.log(60 / 40) / 20
    lse20 = iid['Lq'][20] * 20 + np.log(n_iid)
    Lq_diff2040 = (lse[40] - lse20) / 20
    with np.errstate(divide='ignore', invalid='ignore'):
        I_cdf = np.where(lam_c < al_i_pre, -(np.log(F60) - np.log(F40)) / 20,
                         -(np.log1p(-np.minimum(F60, 1 - 1e-16) ) - np.log1p(-np.minimum(F40, 1 - 1e-16))) / 20)
    I0_diff = dec_iid[S_FIT.index(2)]
    al_i = slope(N[sel], iid['meanlnS'][sel]); v_i = slope(N[sel], iid['varlnS'][sel])
    rows = [['i.i.d. model (1e7)', '%.4f' % al_i, '%.4f' % v_i, '%.4f' % (al_i ** 2 / (2 * v_i)),
             '%.4f' % dec_iid[S_FIT.index(2)], '%.4f–%.4f' % (min(dec_iid), max(dec_iid)),
             '%.4f → %.4f' % (I0_cdf, I0_cdf_corr), '%.4f / %.4f' % (-Lq_diff2040.min(), -Lq_diff.min())],
            ['2-D annealed (%d seeds)' % nseed['annealed'], pm(*ms(fit4c['annealed'][0])), pm(*ms(fit4c['annealed'][1])),
             pm(*ms(fit4c['annealed'][2])), pm(*ms(dec_2d[:, S_FIT.index(2)])),
             '%.4f–%.4f' % (np.nanmin(dec_2d.mean(0)), np.nanmax(dec_2d.mean(0))), '', '']]
    T['iid'] = (['model', 'alpha', 'v', 'Gaussian a²/2v', 'decay of P(S≤2) = I(0)', 'decay s=1.25..10',
                 'λ-CDF I(0+), N 40→60 → corrected for (ln S)²', '−min L(q): N 20→40 / 40→60'], rows)
    rows = []
    P2d = np.array([o['P_le'] for o in R['annealed']])
    for Nn in (5, 10, 20, 40, 60):
        r_ = [Nn]
        for s in (1.5, 2, 5, 10):
            c = int(np.nonzero(K.S_LIST == s)[0][0])
            mm, ss = ms(P2d[:, Nn, c])
            r_ += ['%.3e' % Pi[Nn, c], '%.3e ± %.1e' % (mm, ss), '%+.1f' % ((Pi[Nn, c] - mm) / max(ss / np.sqrt(len(P2d)), 1e-30))]
        rows.append(r_)
    T['iid_P'] = (['N'] + sum([['iid s=%g' % s, '2-D s=%g' % s, 'z'] for s in (1.5, 2, 5, 10)], []), rows)
    rows = []
    for i, q in enumerate(iid['TO_q']):
        j = np.argmin(np.abs(qs - q))
        if abs(qs[j] - q) < 1e-9:
            rows.append(['%.2f' % q, '%.4f' % iid['TO_L'][i], '%.4f' % Lq_diff[j], '%.4f' % iid['Lq'][60][j]])
    T['Lq'] = (['q', 'L(q) transfer operator (vector, exact)', 'L(q) MC, N 40→60', 'L(q) MC, (1/60) ln E'], rows)

    # ---------------------------------------------------------------- 4E
    rows = []
    for m in MODES:
        for fi, fct in enumerate([1.0] + list(K.METRIC_F)):
            rr = []
            for o in R[m]:
                P = o['P_le'] if fct == 1.0 else o['metric_P_le'][K.METRIC_F.tolist().index(fct)]
                mu = o['mom'][:, 0] if fct == 1.0 else o['metric_mean_lnS'][K.METRIC_F.tolist().index(fct)]
                vals = [slope(N[sel], mu[sel])]
                for s in (1.5, 2, 3):
                    c = int(np.nonzero(K.S_LIST == s)[0][0]); Pc = P[:, c]
                    ok = sel & (Pc * o['n_win'] >= 10)
                    vals.append(-slope(N[ok], np.log(Pc[ok])) if ok.sum() >= 5 else np.nan)
                    vals.append(Pc[C])
                rr.append(vals)
            rr = np.array(rr)
            rows.append([m, '%g' % fct, pm(*ms(rr[:, 0]))] +
                        sum([[pm(*ms(rr[:, 1 + 2 * i])), '%.2e' % np.nanmean(rr[:, 2 + 2 * i])] for i in range(3)], []))
    T['4E'] = (['mode', 'y scale f', 'alpha', 'rate S_THR=1.5', 'P(N=%d)' % C, 'rate S_THR=2', 'P(N=%d)' % C,
                'rate S_THR=3', 'P(N=%d)' % C], rows)

    # ---------------------------------------------------------------- 4F
    rows = []
    for m in MODES:
        for Nn in (10, 30, 60):
            g = lambda k: np.array([o['cond_r'][Nn][k] for o in R[m]])
            rows.append([m, Nn, pm(*ms(g('P_lowS')), '{:.4f}'), pm(*ms(g('P_lowS_given_lowr')), '{:.3f}'),
                         pm(*ms(g('P_lowr_given_lowS')), '{:.3f}'), pm(*ms(g('r_q10')), '{:.2f}'),
                         pm(*ms(g('P_lowS_given_r10')), '{:.3f}'), pm(*ms(g('P_r10_given_lowS')), '{:.3f}')])
    T['4F'] = (['mode', 'N', 'P(S≤2)', 'P(S≤2 | r≤ln2)', 'P(r≤ln2 | S≤2)', 'r 10% quantile',
                'P(S≤2 | r≤q10)', 'P(r≤q10 | S≤2)'], rows)

    # ================================================================ FIGURES
    o0 = {m: R[m][0] for m in MODES}
    # (1A-i) branched-flow images
    for m in MODES:
        fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), sharey=True)
        for ri, ax in enumerate(axs):
            img = o0[m]['line_rho_hist'][ri].copy(); img[img <= 0] = np.nan
            im = ax.imshow(img, origin='lower', aspect='auto', extent=[0, 1920, -0.5, C + 0.5],
                           norm=LogNorm(vmin=0.05, vmax=50), cmap='magma', interpolation='nearest')
            ax.set_title('launch row %s' % ['row0 − rows/8', 'row0', 'row0 + rows/8'][ri]); ax.set_xlabel('x (px)')
            ax.grid(False)
        axs[0].set_ylabel('cycle N')
        fig.colorbar(im, ax=axs, label='ρ(x,N)  (rays per px per launched px)', shrink=0.9)
        fig.suptitle('Branched flow — projected ray density, %s (seed 0, %d parcels/px, %d px bins)' % (m, 8, 4), y=1.02)
        savefig(fig, 'p1a_branched_%s' % m, 'Branched-flow image ρ(x,N), %s mode.' % m)
    # (1A-ii) manifold snapshots
    for m in MODES:
        sn = o0[m]['snaps']; Ns = [5, 10, 20, 40, 60]
        fig, axs = plt.subplots(2, 5, figsize=(13, 4.6), gridspec_kw=dict(height_ratios=[2.2, 1]), sharex=True)
        for j, Nn in enumerate(Ns):
            x, y, S, aa = sn[(1, Nn)]
            ax = axs[0, j]
            sc = ax.scatter(x, y, c=np.log10(S), s=0.15, cmap='viridis', vmin=0, vmax=8, rasterized=True)
            ax.set_title('N = %d' % Nn); ax.grid(False)
            caus = np.nonzero(aa * np.roll(aa, -1) < 0)[0]
            ax2 = axs[1, j]
            rh = o0[m]['line_rho_hist'][1, Nn]
            ax2.semilogy(np.arange(rh.size) * 4 + 2, np.where(rh > 0, rh, np.nan), color=INK, lw=0.6)
            ax2.plot(x[caus], np.full(caus.size, 0.03), '|', color='#e34948', ms=5, alpha=0.5)
            ax2.set_ylim(0.02, 200); ax2.set_xlabel('x (px)')
        axs[0, 0].set_ylabel('y (px, unwrapped)'); axs[1, 0].set_ylabel('ρ(x)  | caustics')
        fig.colorbar(sc, ax=axs[0, :], label='log10 S', shrink=0.9)
        fig.suptitle('Launched line from row0 in the phase plane, %s (seed 0); red ticks = caustics (sign changes of a)' % m, y=1.0)
        savefig(fig, 'p1a_manifold_%s' % m, 'Manifold snapshots coloured by S with ρ(x) and caustics, %s.' % m)
    # (1A-iii) caustic counts
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, m in zip(axs, MODES):
        E = np.array([o['caustics_exact'][1] for o in R[m]])
        band(ax, N, E, COL[m], label='exact (piecewise-affine line)', logy=True)
        for di, D_ in enumerate(o0[m]['caustic_dens']):
            ax.plot(N, o0[m]['caustics_dense'][di], ls=':', lw=1.0, color=INK2, alpha=0.4 + 0.2 * di,
                    label='grid %d/px (lower bound)' % D_ if di in (0, len(o0[m]['caustic_dens']) - 1) else None)
        rate, al, be, pred, Nmax, lrate = (np.nanmean(cf[m][:, 1, i]) for i in range(6))
        Em = np.nanmean(E, 0); k0 = 10
        if np.isfinite(Em[k0]):
            xx = np.arange(10, C + 1)
            ax.plot(xx, Em[k0] * np.exp(pred * (xx - k0)), color=INK, lw=1, ls='--',
                    label='α+β/4 = %.3f (pred.)' % pred)
            ax.plot(xx, Em[k0] * np.exp(rate * (xx - k0)), color=COL[m], lw=1, ls='-.',
                    label='fit = %.3f' % rate)
        ax.set_yscale('log'); ax.set_title(m); ax.set_xlabel('N'); ax.set_ylim(1, 1e8)
        ax.legend(loc='upper left')
    axs[0].set_ylabel('caustics on the row0 line')
    savefig(fig, 'p1a_caustics', 'Caustic count vs N (exact, mean ± 1 sd over seeds) against Kaplan Eq. (6) rate α+β/4.')
    # (1A-iv) correlations
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2), sharey=True)
    for m in MODES:
        Y1 = np.array([np.nanmean(o['line_corr_low'][:, :, 0], 0) for o in R[m]])
        Y2 = np.array([np.nanmean(o['line_corr_negr'][:, :, 0], 0) for o in R[m]])
        band(axs[0], N, Y1, COL[m], m); band(axs[1], N, Y2, COL[m], m)
    axs[0].set_title('corr over x of ρ with local fraction S ≤ 2'); axs[1].set_title('corr over x of ρ with local mean of −r')
    for ax in axs:
        ax.axhline(0, color=INK2, lw=0.6); ax.set_xlabel('N')
    axs[0].set_ylabel('Pearson r (mean of 3 rows)'); axs[0].legend()
    savefig(fig, 'p1a_corr', 'Correlation across x of ρ with low-stretch fraction and with −r.')
    # (1B) band projections
    for m in MODES:
        fig = plt.figure(figsize=(11, 5.2))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.6, 1])
        BX = np.array([o['band_x'] for o in R[m]])
        for b in range(3):
            ax = fig.add_subplot(gs[0, b]); img = BX[0, :, b, :].astype(float); img[img <= 0] = np.nan
            im = ax.imshow(img, origin='lower', aspect='auto', extent=[0, 1920, -0.5, C + 0.5],
                           norm=LogNorm(vmin=1, vmax=np.nanmax(BX[0]) + 1), cmap='magma', interpolation='nearest')
            ax.set_title(['S ≤ 2', '2 < S ≤ 5', 'S > 5'][b]); ax.grid(False); ax.set_xlabel('x (px)')
            if b == 0:
                ax.set_ylabel('cycle N')
        fig.colorbar(im, ax=fig.axes, label='tracers per 4-px bin', shrink=0.8)
        for j, Nn in enumerate((10, 30, 60)):
            ax = fig.add_subplot(gs[1, j])
            for b in range(3):
                band(ax, np.arange(BX.shape[-1]) * 4 + 2, BX[:, Nn, b, :] + 0.5, BANDC[b],
                     ['S ≤ 2', '2–5', '> 5'][b], lw=0.8, logy=True)
            ax.set_yscale('log'); ax.set_title('marginal at N = %d' % Nn); ax.set_xlabel('x (px)')
            if j == 0:
                ax.set_ylabel('tracers per bin (+0.5)'); ax.legend(ncol=3)
        fig.suptitle('Band projections of the tracer lattice onto x, %s (image seed 0; marginals mean ± 1 sd)' % m)
        savefig(fig, 'p1b_bands_%s' % m, 'Band projections, %s.' % m)
    # (2) hazards vs age per mode (at-risk counts on a second axis as requested)
    for m in MODES:
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.3))
        for ax, (hk, rk, lab) in zip(axs, (('h_re', 'r_re', 're-entry h_re(k)'), ('h_ex', 'r_ex', 'exit h_exit(k), entrants'),
                                          ('h_co', 'r_co', 'exit, original cohort'))):
            H = np.array([p[hk] for p in P2[m]]); Rk = np.array([p[rk] for p in P2[m]])
            k = np.arange(H.shape[1])
            band(ax, k[1:], H[:, 1:], COL[m], logy=True); ax.plot(k[1:], ms(H[:, 1:])[0], 'o', ms=2.5, color=COL[m])
            ax.set_yscale('log'); ax.set_xlabel('age k (cycles since last crossing)'); ax.set_title(lab)
            axr = ax.twinx(); axr.bar(k[1:], Rk[:, 1:].mean(0), color=INK2, alpha=0.15, width=0.8)
            axr.set_yscale('log'); axr.grid(False); axr.spines['right'].set_visible(True)
            axr.set_ylabel('at risk (mean per seed)', color=INK2)
        axs[0].set_ylabel('hazard')
        fig.suptitle('Age-structured hazards, %s (points with < %d at risk hidden)' % (m, MIN_RISK))
        savefig(fig, 'p2_hazard_%s' % m, 'Hazards vs age with at-risk counts, %s.' % m)
    # KM
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), sharey=True)
    for ax, m in zip(axs, MODES):
        for Fk, Wk, ls, lab in (('F_re', 'W_re', '-', 're-entry'), ('F_ex', 'W_ex', '--', 'exit (entrants)'),
                                ('F_co', 'W_co', ':', 'exit (orig. cohort)')):
            Fm = np.array([p[Fk] for p in P2[m]]); k = np.arange(Fm.shape[1])
            band(ax, k, Fm, COL[m], lab, ls=ls)
            mm = ms(Fm)[0]
            for kk in (1, 2):
                ax.plot(kk, mm[kk], 'o', color=COL[m], ms=4)
            W = np.mean([p[Wk][0] for p in P2[m]]); Wn = int(np.median([p[Wk][1] for p in P2[m]]))
            ax.plot(Wn, W, 's', color=COL[m], ms=5, mfc='white')
        ax.set_xscale('log'); ax.set_title(m + '  (dots: k=1,2; square: whole run)'); ax.set_xlabel('age k')
    axs[0].set_ylabel('Kaplan–Meier F(k)'); axs[0].legend()
    savefig(fig, 'p2_km', 'Kaplan–Meier cumulative probability of re-entry and exit.')
    # pooled bins bar chart
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.3), sharey=False)
    labs = ['1', '2', '3–5', '6–10', '11–20', '21–40', '>40']
    w = 0.26
    for ax, key, tt in ((axs[0], 'pool_re', 're-entry hazard'), (axs[1], 'pool_ex', 'exit hazard (entrants)')):
        for j, m in enumerate(MODES):
            Y = np.array([[b[0] for b in p[key]] for p in P2[m]]); mm, ss = ms(Y)
            ax.bar(np.arange(7) + (j - 1) * w, mm, w * 0.92, yerr=ss, color=COL[m], label=m,
                   error_kw=dict(lw=0.8, ecolor=INK2))
        ax.set_xticks(range(7)); ax.set_xticklabels(labs); ax.set_xlabel('age bin (cycles)'); ax.set_title(tt)
        ax.set_yscale('log')
    axs[0].legend(); axs[0].set_ylabel('pooled hazard')
    savefig(fig, 'p2_pooled_bars', 'Pooled age-bin hazards, both directions, three modes.')
    # cohort heat maps
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.5))
    for ax, m in zip(axs, MODES):
        rr = np.sum([o['risk_re'] for o in R[m]], 0); ee = np.sum([o['ev_re'] for o in R[m]], 0)
        H = np.where(rr >= MIN_RISK, ee / np.maximum(rr, 1), np.nan)
        im = ax.imshow(H, origin='lower', aspect='auto', norm=LogNorm(vmin=1e-3, vmax=0.5), cmap='viridis',
                       interpolation='nearest')
        ax.set_xlabel('age k'); ax.set_title(m); ax.grid(False)
    axs[0].set_ylabel('exit cycle (cohort)')
    fig.colorbar(im, ax=axs, label='re-entry hazard (seeds pooled)')
    savefig(fig, 'p2_cohort', 'Cohort heat map: re-entry hazard by exit cycle and age.')
    # depth-controlled ratios
    fig, ax = plt.subplots(figsize=(7.5, 3.3))
    for j, m in enumerate(MODES):
        Y = np.array([[d['fine'][i][0] for i in range(len(K.DEPTH_LABELS))] for d in D[m]]); mm, ss = ms(Y)
        ax.errorbar(np.arange(len(mm)) + (j - 1) * 0.15, mm, yerr=ss, fmt='o', color=COL[m], ms=4, label=m, lw=1)
        crude = np.mean([d['crude'] for d in D[m]]); mh = np.mean([d['MH'] for d in D[m]])
        ax.axhline(crude, color=COL[m], ls=':', lw=1); ax.axhline(mh, color=COL[m], ls='--', lw=1)
    ax.axhline(1, color=INK, lw=0.8)
    ax.set_xticks(range(len(K.DEPTH_LABELS))); ax.set_xticklabels(K.DEPTH_LABELS)
    ax.set_xlabel('depth S_{N−1}'); ax.set_ylabel('h(returned before) / h(first time)'); ax.set_yscale('log')
    ax.set_title('Returned-before effect at fixed depth (dotted: crude; dashed: Mantel–Haenszel depth×age)')
    ax.legend()
    savefig(fig, 'p2_depth', 'Depth-controlled re-entry ratio.')
    # 3-mode overlay
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.3))
    for m in MODES:
        H = np.array([p['h_re'] for p in P2[m]]); k = np.arange(H.shape[1])
        band(axs[0], k[1:], H[:, 1:], COL[m], m, logy=True)
        H = np.array([p['h_ex'] for p in P2[m]])
        band(axs[1], k[1:], H[:, 1:], COL[m], m, logy=True)
    for ax, tt in zip(axs, ('re-entry hazard', 'exit hazard (entrants)')):
        ax.set_yscale('log'); ax.set_xlabel('age k'); ax.set_title(tt)
    axs[0].legend(); axs[0].set_ylabel('hazard')
    savefig(fig, 'p2_overlay', 'Three modes overlaid: hazards vs age.')

    # (4A) distribution of S
    e = K.LOG10_EDGES; c10 = 0.5 * (e[:-1] + e[1:]); dl10 = e[1] - e[0]
    for m in MODES:
        Hs = np.array([o['hist_log10S'] / o['n_win'][:, None] / dl10 for o in R[m]])   # density per decade
        fig, axs = plt.subplots(1, 4, figsize=(14, 3.4))
        ax = axs[0]
        for j, Nn in enumerate((1, 3, 5, 10, 20, 30, 40, 60)):
            mm = Hs[:, Nn].mean(0)
            ax.fill_between(c10, j * 1.2, j * 1.2 + mm / max(mm.max(), 1e-12), color=COL[m], alpha=0.25, lw=0)
            ax.plot(c10, j * 1.2 + mm / max(mm.max(), 1e-12), color=COL[m], lw=0.9)
            ax.text(-0.3, j * 1.2 + 0.2, 'N=%d' % Nn, fontsize=7, ha='right', color=INK2)
        ax.set_xlim(-1.2, 10); ax.set_yticks([]); ax.set_xlabel('log10 S'); ax.set_title('ridge (normalised)')
        ax = axs[1]
        img = Hs.mean(0); img[img <= 0] = np.nan
        im = ax.imshow(img, origin='lower', aspect='auto', extent=[0, 13, -0.5, C + 0.5], norm=LogNorm(1e-5, 10),
                       cmap='magma', interpolation='nearest'); ax.grid(False)
        ax.set_xlabel('log10 S'); ax.set_ylabel('N'); ax.set_title('density per decade'); fig.colorbar(im, ax=ax)
        ax = axs[2]
        Q = np.array([o['quant'] for o in R[m]]).mean(0) / np.log(10)
        ax.fill_between(N, Q[:, 0], Q[:, 6], color=COL[m], alpha=0.12, lw=0, label='1–99%')
        ax.fill_between(N, Q[:, 1], Q[:, 5], color=COL[m], alpha=0.2, lw=0, label='10–90%')
        ax.fill_between(N, Q[:, 2], Q[:, 4], color=COL[m], alpha=0.35, lw=0, label='25–75%')
        ax.plot(N, Q[:, 3], color=COL[m], label='median'); ax.axhline(np.log10(2), color=INK, ls='--', lw=0.8)
        ax.set_xlabel('N'); ax.set_ylabel('log10 S'); ax.set_title('quantile fan'); ax.legend(loc='upper left')
        ax = axs[3]
        B = np.array([o['bandshare'] for o in R[m]]).mean(0)
        labs6 = ['≤2', '2–5', '5–10', '10–100', '100–1e4', '>1e4']
        cols6 = ['#2a78d6', '#1baf7a', '#eda100', '#eb6834', '#e34948', '#6250d6']
        ax.stackplot(N, B.T, colors=cols6, labels=labs6, edgecolor='white', linewidth=0.3)
        ax.set_xlim(0, C); ax.set_ylim(0, 1); ax.set_xlabel('N'); ax.set_title('band shares (S bands)')
        ax.legend(loc='upper right', ncol=2, fontsize=6.5)
        fig.suptitle('Distribution of S on the Eulerian window, %s (%d seeds)' % (m, nseed[m]))
        savefig(fig, 'p4a_%s' % m, 'Distribution of S, %s.' % m)
    # (4B) collapse
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), sharey=True)
    zz = np.linspace(-5, 5, 400)
    for ax, m in zip(axs, MODES):
        o = o0[m]
        for j, Nn in enumerate((10, 20, 30, 40, 60)):
            mu, v = o['mom'][Nn, 0], o['mom'][Nn, 1]
            z = (c10 * np.log(10) - mu) / np.sqrt(v)
            dens = o['hist_log10S'][Nn] / o['n_win'][Nn] / (dl10 * np.log(10)) * np.sqrt(v)
            ax.semilogy(z, np.where(dens > 0, dens, np.nan), lw=1, color=COL[m], alpha=0.35 + 0.13 * j, label='N=%d' % Nn)
        ax.semilogy(zz, np.exp(-zz ** 2 / 2) / np.sqrt(2 * np.pi), color=INK, ls='--', lw=1, label='N(0,1)')
        ax.set_ylim(1e-6, 1); ax.set_xlim(-5, 6); ax.set_xlabel('(ln S − mean)/sd'); ax.set_title(m)
    axs[0].set_ylabel('density'); axs[0].legend()
    savefig(fig, 'p4b_collapse', 'Collapse of (ln S − mean)/sd against a standard normal (seed 0).')
    # (4C) tail test
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4), sharey=True)
    for ax, m in zip(axs, MODES):
        al, v, pr, me = fit4c[m]
        for s, cc in zip((1.25, 2, 5, 10), ('#2a78d6', '#1baf7a', '#eb6834', '#6250d6')):
            col = int(np.nonzero(K.S_LIST == s)[0][0])
            Y = np.array([-np.log(np.where(o['P_le'][:, col] > 0, o['P_le'][:, col], np.nan)) for o in R[m]])
            band(ax, N, Y, cc, 's = %g' % s)
        mm = np.nanmean(np.array([-np.log(o['P_le'][:, 2]) for o in R[m]]), 0)
        xx = np.arange(20, C + 1)
        ax.plot(xx, mm[20] + np.mean(pr) * (xx - 20), color=INK, ls='--', lw=1, label='slope α²/2v = %.3f' % np.mean(pr))
        ax.set_title('%s: measured %.3f vs pred %.3f (s=2)' % (m, np.nanmean(me[:, 2]), np.mean(pr)))
        ax.set_xlabel('N')
    axs[0].set_ylabel('−ln P(S ≤ s)'); axs[0].legend()
    savefig(fig, 'p4c_tail', 'Kaplan log-normal tail test: −ln P(S≤s) vs N with the predicted slope.')
    # (4D) iid model
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
    ax = axs[0]
    for s, cc in zip((1.5, 2, 5, 10), ('#2a78d6', '#1baf7a', '#eb6834', '#6250d6')):
        col = int(np.nonzero(K.S_LIST == s)[0][0])
        band(ax, N, P2d[:, :, col], cc, '2-D s=%g' % s, logy=True)
        ax.plot(N, Pi[:, col], color=cc, ls='--', lw=1)
    ax.set_yscale('log'); ax.set_xlabel('N'); ax.set_ylabel('P(S ≤ s)'); ax.set_title('2-D annealed (bands) vs i.i.d. (dashed)')
    ax.legend(fontsize=6.5)
    ax = axs[1]
    for Nn, cc in zip((20, 40, 60), ('#a3c4ef', '#5b98e0', '#2a78d6')):
        ax.plot(lam_c, Irate[Nn] - np.nanmin(Irate[Nn]), color=cc, lw=1, label='N=%d, −(1/N) ln p' % Nn)
    ax.plot(lam_e[1:], I_cdf, color=INK, lw=1.2, label='CDF difference N 40→60')
    lnp = (lam_c - al_i) ** 2 / (2 * v_i)
    ax.plot(lam_c, lnp, color='#eb6834', ls='--', lw=1, label='log-normal (λ−α)²/2v')
    ax.axhline(np.nanmean(dec_2d[:, 2]), color='#1baf7a', ls=':', lw=1, label='2-D decay s=2')
    ax.set_xlim(0, 0.5); ax.set_ylim(0, 0.25); ax.set_xlabel('λ = ln S / N'); ax.set_ylabel('I(λ)')
    ax.set_title('rate function'); ax.legend(fontsize=6.5)
    ax = axs[2]
    ax.plot(qs, Lq_diff, 'o-', ms=3, color='#2a78d6', label='MC, N 40→60')
    ax.plot(qs, iid['Lq'][60], '-', color='#a3c4ef', lw=1, label='MC, (1/60) ln E')
    ax.plot(iid['TO_q'], iid['TO_L'], 's', ms=4, mfc='white', color=INK, label='transfer operator (vector)')
    ax.axhline(-I0_diff, color='#eb6834', ls='--', lw=1, label='−I(0) = −decay of P(S≤2)')
    ax.set_xlabel('q'); ax.set_ylabel('L(q)'); ax.set_title('generalised Lyapunov exponents'); ax.legend(fontsize=6.5)
    savefig(fig, 'p4d_iid', 'i.i.d. random-matrix model: P(S≤s,N), rate function, L(q).')
    # (4E) metric
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), sharey=True)
    for ax, m in zip(axs, MODES):
        band(ax, N, np.array([o['P_le'][:, 2] for o in R[m]]), INK, 'f = 1', logy=True)
        for fi, (fct, cc) in enumerate(zip(K.METRIC_F, ('#2a78d6', '#eb6834', '#6250d6'))):
            band(ax, N, np.array([o['metric_P_le'][fi, :, 2] for o in R[m]]), cc, 'f = %g' % fct, logy=True)
        ax.set_yscale('log'); ax.set_xlabel('N'); ax.set_title(m)
    axs[0].set_ylabel('P(S_f ≤ 2), S_f = σmax(D M D⁻¹)'); axs[0].legend()
    savefig(fig, 'p4e_metric', 'Metric rescaling y → f·y: threshold fractions move, decay rates do not.')
    # (4F) joint
    fig, axs = plt.subplots(3, 3, figsize=(10, 8.4))
    RB = np.arange(-6, 30.0001, 0.25); SB = np.arange(0, 30.0001, 0.25)
    for i, m in enumerate(MODES):
        for j, Nn in enumerate((10, 30, 60)):
            ax = axs[i, j]; H = sum(o['joint_r_lnS'][Nn] for o in R[m]).T.astype(float); H[H <= 0] = np.nan
            ax.imshow(H, origin='lower', aspect='auto', extent=[RB[0], RB[-1], SB[0], SB[-1]], norm=LogNorm(),
                      cmap='viridis', interpolation='nearest')
            ax.plot([0, 30], [0, 30], color='#e34948', lw=0.8, ls='--')
            ax.axhline(np.log(2), color='white', lw=0.6); ax.axvline(np.log(2), color='white', lw=0.6)
            ax.grid(False); ax.set_title('%s, N=%d' % (m, Nn))
            hi = max(4, 1.3 * Nn * 0.3); ax.set_xlim(-3, hi); ax.set_ylim(0, hi)
            if i == 2:
                ax.set_xlabel('r = ln |tangent|')
            if j == 0:
                ax.set_ylabel('ln S')
    fig.suptitle('Joint distribution of r and ln S on the launched lines (dashed: ln S = r; S ≥ e^r always)')
    savefig(fig, 'p4f_joint', 'Joint distribution of r vs ln S.')

    # ================================================================ write tables
    md = ['# Kick-drift LSVD — multi-seed analysis tables', '',
          'Seeds: ' + ', '.join('%s %d' % (m, nseed[m]) for m in MODES) + '; cycles = %d. Values are mean ± sd over seeds.' % C, '']
    def mdt(key, title):
        hdr, rows = T[key]
        md.extend(['## ' + title, '', '| ' + ' | '.join(hdr) + ' |', '|' + '---|' * len(hdr)])
        md.extend('| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows)
        md.append('')
    md.extend(['## Checks', ''] + ['- **%s** — %s: %s' % ('PASS' if ok else 'FAIL', n, d) for n, ok, d in checks] + [''])
    TITLES = [('4C', '4C. Kaplan log-normal tail test (fit N = 20..%d)' % C), ('4C_seed0', '4C reproduction, seed 0'),
              ('caustics', '1A. Caustic growth vs Kaplan Eq. (6)'), ('corr', '1A(iv). Correlations across x'),
              ('pooled', '2A. Pooled age-bin hazards'), ('KM', '2B. Kaplan–Meier F(k)'),
              ('depth', '2C. Returned-before effect at fixed depth'), ('Ple', '4A. P(S ≤ s)'),
              ('mom', '4A. Moments and quantiles of ln S'), ('iid', '4D. i.i.d. random-matrix model'),
              ('iid_P', '4D. P(S ≤ s, N): i.i.d. vs 2-D annealed (z = difference in standard errors)'),
              ('Lq', '4D. Generalised Lyapunov exponents'), ('4E', '4E. Robustness: threshold and metric'),
              ('4F', '4F. r vs ln S on the lines')]
    for k, t in TITLES:
        mdt(k, t)
    os.makedirs(a.report, exist_ok=True)
    open(os.path.join(a.report, 'kd_tables.md'), 'w').write('\n'.join(md))
    json.dump({'T': {k: v for k, v in T.items() if k != 'checks'}, 'checks': checks,
               'figs': figs, 'titles': TITLES, 'nseed': nseed, 'cycles': C},
              open(os.path.join(a.report, 'kd_tables.json'), 'w'), default=float, indent=0)
    print('\n'.join(md[:60]))
    return T, figs, checks


if __name__ == '__main__':
    main()
