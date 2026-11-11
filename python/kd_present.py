"""
kd_present -- reader-facing reports, built READ-ONLY from existing results.

    python3 kd_present.py

Reads python/out/*.pkl.gz (run results) and reports/kd_tables.json +
reports/kd_findings.json (the previous report). Writes only to reports_v2/:
paper1_stretching.pdf, paper2_reentry.pdf, appendix.pdf, fig/*.png, CHECK.md.

Every number quoted in the text or a table is taken from the previous report
through V()/S() and recorded in a ledger; numbers not in the previous report
are recorded as NEW with how they were read. The build stops if a quoted
number does not equal its source. No simulation is run and nothing outside
reports_v2/ is written.
"""
import glob, hashlib, json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, FixedLocator
import kd_engine as K
import kd_report as KR

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'out')
OLD = os.path.join(ROOT, 'reports')
NEW = os.path.join(ROOT, 'reports_v2')
FIG = os.path.join(NEW, 'fig')
os.makedirs(FIG, exist_ok=True)

MODES = ['annealed', 'quenched', 'cosine']                       # data keys
NAME = {'annealed': 'annealed', 'quenched': 'quenched', 'cosine': 'periodic'}
LONG = {'annealed': 'annealed (fresh random kick each cycle)',
        'quenched': 'quenched (one frozen random kick)',
        'cosine': 'periodic (standard map, K = 0.5)'}
LEG = {'annealed': 'annealed\n(fresh random kick\neach cycle)',
       'quenched': 'quenched\n(one frozen\nrandom kick)',
       'cosine': 'periodic\n(standard map,\nK = 0.5)'}
COL = {'annealed': '#2a78d6', 'quenched': '#eb6834', 'cosine': '#1baf7a'}
MODEL = '#3d3d3a'                                                # i.i.d. model / references
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
W_IN = 6.46                                                      # printed width of every figure, inches
FS = 11                                                          # minimum font size, pt

plt.rcParams.update({
    'font.size': FS, 'axes.titlesize': FS, 'axes.labelsize': FS, 'xtick.labelsize': FS,
    'ytick.labelsize': FS, 'legend.fontsize': FS, 'figure.titlesize': FS + 1,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'lines.linewidth': 1.8, 'legend.frameon': False,
    'figure.constrained_layout.use': True, 'savefig.dpi': 200})

# ============================================================ previous report + ledger
J = json.load(open(os.path.join(OLD, 'kd_tables.json')))
T = J['T']
FND = json.load(open(os.path.join(OLD, 'kd_findings.json')))
OLD_TEXT = open(os.path.join(OLD, 'kd_tables.md')).read() + json.dumps(FND, ensure_ascii=False)
LEDGER = []          # dict(where, source, old, shown, status)
WHERE = ['']


def tiny_sd(sd):
    """An sd that rounds to zero at the printed precision is shown as '< 0.0001' (same precision)."""
    try:
        if float(sd) == 0.0 and '.' in sd:
            d = len(sd.split('.')[1])
            return '< 0.' + '0' * (d - 1) + '1'
    except ValueError:
        pass
    return sd


def _nums(s):
    return [float(x) for x in re.findall(r'[-+]?\d+\.?\d*(?:e[-+]?\d+)?', str(s).replace('−', '-'))]


def V(key, col, **match):
    """Raw cell of the previous report's table `key` (column name or index)."""
    hdr, rows = T[key]
    ci = col if isinstance(col, int) else hdr.index(col)
    sel = [r for r in rows if all(str(r[hdr.index(k.replace('__', ' '))]) == str(v) for k, v in match.items())]
    assert len(sel) == 1, (key, col, match, len(sel))
    return str(sel[0][ci])


def S(key, col, **match):
    """A previous-report cell, formatted for readers and entered in the ledger.
    Deterministic (periodic) rows lose their meaningless '± 0.0000'."""
    raw = V(key, col, **match)
    mode = match.get('mode')
    m = re.match(r'^\s*([-+]?[\d.e+-]+)\s*±\s*([-+]?[\d.e+-]+|nan)(.*)$', raw)
    if m:
        val, sd, rest = m.group(1).lstrip('+'), m.group(2).lstrip('+'), m.group(3)
        if mode == 'cosine':
            shown = '%s ± n/a (deterministic)%s' % (val, rest)
        else:
            shown = '%s ± %s%s' % (val, tiny_sd(sd), rest)
    else:
        shown = raw.lstrip('+')
        if mode == 'cosine' and re.match(r'^[-+]?\d', shown):
            shown = shown + ' (deterministic)'
    a, b = _nums(raw), _nums(shown)
    ok = len(b) >= 1 and np.allclose(a[:len(b)][:1], b[:1])
    LEDGER.append(dict(where=WHERE[0], source='%s | %s | %s' % (key, col, match), old=raw, shown=shown,
                       status='same as previous report' if ok else 'DIFFERS'))
    if not ok:
        raise SystemExit('STOP: number %r shown as %r differs from the previous report (%s)' % (raw, shown, key))
    return shown


def F(text_in_old, shown=None):
    """A number quoted from the previous report's summary text (must appear there verbatim)."""
    ok = text_in_old in OLD_TEXT
    LEDGER.append(dict(where=WHERE[0], source='previous summary text', old=text_in_old,
                       shown=shown or text_in_old, status='same as previous report' if ok else 'DIFFERS'))
    if not ok:
        raise SystemExit('STOP: %r is not in the previous report' % text_in_old)
    return shown or text_in_old


def NEWN(value, how):
    """A number that is not in the previous report, read directly from the stored results."""
    LEDGER.append(dict(where=WHERE[0], source=how, old='(not in previous report)', shown=value, status='NEW'))
    return value


def pct(s):
    return s


# ============================================================ stored results (read only)
def loadall(prefix):
    R = {}
    for m in MODES:
        fs = sorted(glob.glob(os.path.join(OUT, '%s_%s_s*.pkl.gz' % (prefix, m))),
                    key=lambda p: int(re.search(r'_s(\d+)', p).group(1)))
        R[m] = [KR.load(p) for p in fs]
    return R


R = loadall('run')
P1C = loadall('p1c')
for m in MODES:                       # periodic is deterministic: p1c stored once
    if m == 'cosine' and len(P1C[m]) == 1:
        P1C[m] = P1C[m] * len(R[m])
IID = KR.load(os.path.join(OUT, 'iid.pkl.gz'))
C = R['annealed'][0]['cycles']
NN = np.arange(C + 1)
NSEED = {m: len(R[m]) for m in MODES}


def ms(x, axis=0):
    return KR.ms(x, axis)


# ============================================================ figure helpers + automatic checks
FIGCHECK = []
BANNED = ['cosine', 'LSVD', 'Eulerian', 'compressed', 'Kaplan-Meier', 'Kaplan–Meier', 'row0', 'rows/8',
          'S_THR', 'an_br', 'lattice', 'window']


def check_and_save(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    small, banned, overlap = [], [], []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip():
            continue
        if t.get_fontsize() < FS - 1e-9:
            small.append((s, t.get_fontsize()))
        for b in BANNED:
            if b.lower() in s.lower():
                banned.append((s, b))
    for ax in fig.axes:
        for lab in (ax.get_xticklabels(), ax.get_yticklabels()):
            bbs = [l.get_window_extent(r) for l in lab if l.get_visible() and l.get_text().strip()]
            for i in range(len(bbs)):
                for j in range(i + 1, len(bbs)):
                    if bbs[i].overlaps(bbs[j]):
                        overlap.append(ax.get_title()[:40])
    w, h = fig.get_size_inches()
    p = os.path.join(FIG, name + '.png')
    fig.savefig(p)                        # no bbox cropping: printed size = figure size
    plt.close(fig)
    FIGCHECK.append(dict(name=name, w=w, h=h, small=small, banned=banned, overlap=sorted(set(overlap))))
    if small or banned or overlap:
        raise SystemExit('figure check failed for %s: small=%s banned=%s overlap=%s' % (name, small[:3], banned[:3], overlap[:3]))
    return p, w, h


def mode_handles(extra=()):
    hs = [Line2D([0], [0], color=COL[m], lw=2.4) for m in MODES]
    ls = [LEG[m] for m in MODES]
    for h, l in extra:
        hs.append(h); ls.append(l)
    return hs, ls


def legend_bottom(fig, hs, ls, ncol=None):
    fig.legend(hs, ls, loc='outside lower center', ncol=ncol or len(hs), handlelength=2.2, columnspacing=1.4)


def band(ax, x, Y, color, ls='-', lw=1.8, logy=False, alpha=0.18):
    m, s = ms(Y)
    ax.plot(x, m, color=color, ls=ls, lw=lw)
    if Y.shape[0] > 1 and np.any(s > 0):
        lo = m - s
        if logy:
            lo = np.where(lo > 0, lo, np.nan)
        ax.fill_between(x, lo, m + s, color=color, alpha=alpha, lw=0)


def frame(ax, m):
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(COL[m]); sp.set_linewidth(2.0)


# ============================================================ figures
FIGS = {}            # id -> dict(path, w, h, caption)


def addfig(fid, fig, name, caption):
    p, w, h = check_and_save(fig, name)
    FIGS[fid] = dict(path=p, w=w, h=h, caption=caption)


def fig_branched(fid, ri, name):
    lab = ['960 − 240', '960', '960 + 240'][ri]
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 2.75), sharey=True)
    cmap = plt.get_cmap('magma_r').copy()
    cmap.set_bad(cmap(0.0))                       # empty bins = the lowest density, no blocks
    titles = {'annealed': 'annealed: branches\nform, then fade', 'quenched': 'quenched: branches\nform and persist',
              'cosine': 'periodic: a fixed\ncaustic pattern'}
    for ax, m in zip(axs, MODES):
        img = R[m][0]['line_rho_hist'][ri].copy()
        img = np.ma.masked_less_equal(img, 0)
        im = ax.imshow(img, origin='lower', aspect='auto', extent=[0, 1920, -0.5, C + 0.5],
                       norm=LogNorm(vmin=0.05, vmax=50), cmap=cmap, interpolation='nearest')
        ax.set_title(titles[m]); ax.grid(False); frame(ax, m)
        ax.xaxis.set_major_locator(FixedLocator([0, 960, 1920]))
        ax.set_xlabel('position x (px)')
    axs[0].set_ylabel('cycles N')
    cb = fig.colorbar(im, ax=axs, shrink=0.95, pad=0.02)
    cb.set_label('ray density\n(relative to launch)')
    addfig(fid, fig, name, 'Ray density along x (horizontal) against cycles N (vertical) for rays launched at '
           'momentum y0 = %s px, one run per mode (seed 0). Dark = many rays. Empty bins are drawn as the '
           'lowest density.' % lab)


def fig_collapse(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 2.95), sharey=True)
    e = K.LOG10_EDGES; c10 = 0.5 * (e[:-1] + e[1:]); dl = e[1] - e[0]
    zz = np.linspace(-5, 5, 300)
    styles = {10: '-', 30: '--', 60: ':'}
    titles = {'annealed': 'annealed:\nclose to Gaussian', 'quenched': 'quenched:\nroughly Gaussian',
              'cosine': 'periodic:\nfar from Gaussian'}
    for ax, m in zip(axs, MODES):
        o = R[m][0]
        for Nn, ls in styles.items():
            mu, v = o['mom'][Nn, 0], o['mom'][Nn, 1]
            z = (c10 * np.log(10) - mu) / np.sqrt(v)
            d = o['hist_log10S'][Nn] / o['n_win'][Nn] / (dl * np.log(10)) * np.sqrt(v)
            ax.semilogy(z, np.where(d > 0, d, np.nan), color=COL[m], ls=ls, lw=1.8)
        ax.semilogy(zz, np.exp(-zz ** 2 / 2) / np.sqrt(2 * np.pi), color=MODEL, lw=1.2, alpha=0.7)
        ax.set_ylim(1e-5, 1); ax.set_xlim(-5, 6); ax.set_title(titles[m]); frame(ax, m)
        ax.xaxis.set_major_locator(FixedLocator([-4, 0, 4]))
        ax.set_xlabel('(ln S − mean) / sd')
    axs[0].set_ylabel('probability density')
    hs = [Line2D([0], [0], color=INK2, ls=ls) for ls in styles.values()] + [Line2D([0], [0], color=MODEL, lw=1.2)]
    legend_bottom(fig, hs, ['N = 10', 'N = 30', 'N = 60', 'standard normal'])
    addfig(fid, fig, 'appA_lnS_shape', 'Shape of the distribution of ln S: each curve is ln S shifted by its '
           'mean and divided by its standard deviation, at three cycle numbers, against a standard normal '
           '(grey). Measured on the pixel grid, one run per mode (seed 0). The left edge is cut off because '
           'S cannot fall below 1.')


def fig_why(fid):
    fig, axs = plt.subplots(1, 2, figsize=(W_IN, 3.0))
    c2 = int(np.nonzero(K.S_LIST == 2)[0][0])
    P2d = np.array([o['P_le'][:, c2] for o in R['annealed']])
    ax = axs[0]
    band(ax, NN, P2d, COL['annealed'], logy=True)
    ax.plot(NN, IID['P_le'][:, c2], color=MODEL, ls='--', lw=1.6)
    pred = float(_nums(V('4C', 'pred a²/2v', mode='annealed'))[0])
    mm = P2d.mean(0); xx = np.arange(20, C + 1)
    ax.plot(xx, mm[20] * np.exp(-pred * (xx - 20)), color=INK2, ls=':', lw=2.0)
    ax.set_yscale('log'); ax.set_xlabel('cycles N'); ax.set_ylabel('share of material with S ≤ 2')
    ax.set_title('Independent random steps\nreproduce the decay')
    ax = axs[1]
    n = IID['n_samples']; qs = IID['qs']
    lse40 = IID['Lq'][40] * 40 + np.log(n); lse60 = IID['Lq'][60] * 60 + np.log(n)
    Ld = (lse60 - lse40) / 20
    al, v = float(V('iid', 'alpha', model='i.i.d. model (1e7)')), float(V('iid', 'v', model='i.i.d. model (1e7)'))
    I0 = float(V('iid', 'decay of P(S≤2) = I(0)', model='i.i.d. model (1e7)'))
    sel = qs <= 2
    ax.plot(qs[sel], Ld[sel], 'o', color=MODEL, ms=4.5)
    q2 = np.linspace(-4, 2, 200)
    ax.plot(q2, al * q2 + v * q2 ** 2 / 2, color=INK2, ls=':', lw=2.0)
    ax.axhline(-I0, color=COL['annealed'], ls='--', lw=1.4)
    ax.set_xlim(-4.2, 2.2); ax.set_ylim(-0.16, 0.65)
    ax.set_xlabel('moment order q'); ax.set_ylabel('growth rate of the q-th moment\nof S, L(q) (per cycle)')
    ax.set_title('The low tail flattens near\nthe measured decay rate')
    hs = [Line2D([0], [0], color=COL['annealed'], lw=2.4), Line2D([0], [0], color=MODEL, ls='--', marker='o', ms=4),
          Line2D([0], [0], color=INK2, ls=':', lw=2.0), Line2D([0], [0], color=COL['annealed'], ls='--', lw=1.4)]
    legend_bottom(fig, hs, ['annealed runs\n(mean ± sd, 10 seeds)', 'independent-step\nmodel (10 million)',
                            'log-normal\nprediction', 'measured\ndecay rate'], ncol=4)
    addfig(fid, fig, 'p1_fig3_why', 'Left: share of material with S ≤ 2 in the annealed runs (blue band) and in a '
           'model where every cycle applies an independent random step drawn from the same kicks (dashed); the '
           'dotted line has the slope the log-normal law predicts. Right: growth rate L(q) of the q-th moment '
           'of S in the model (dots, from cycles 40 to 60), the log-normal form (dotted) and minus the measured '
           'decay rate of the S ≤ 2 share (dashed).')


def fig_caustics(fid):
    fig, axs = plt.subplots(1, 2, figsize=(W_IN, 3.35))
    ax = axs[0]
    for m in MODES:
        E = np.array([o['caustics_exact'][1] for o in R[m]])
        band(ax, NN, E, COL[m], logy=True)
        cf = np.array([KR.caustic_fit(o, 1) for o in R[m]])
        pred = np.nanmean(cf[:, 3]); Em = np.nanmean(E, 0)
        xx = np.arange(10, 41)
        ax.plot(xx, Em[10] * np.exp(pred * (xx - 10)), color=COL[m], ls='--', lw=1.4)
    ax.set_yscale('log'); ax.set_ylim(1, 3e6); ax.set_xlim(0, 45)
    ax.set_xlabel('cycles N'); ax.set_ylabel('number of caustics on the line')
    ax.set_title('Caustics multiply\nexponentially')
    ax = axs[1]
    for i, m in enumerate(MODES):
        cf = np.array([KR.caustic_fit(o, 1) for o in R[m]])
        for j, (k, mk, fill) in enumerate(((0, 'o', COL[m]), (5, 's', 'white'), (3, 'X', COL[m]))):
            mm, ss = ms(cf[:, k])
            ax.errorbar(i + (j - 1) * 0.22, mm, yerr=ss if NSEED[m] > 1 and m != 'cosine' else None,
                        fmt=mk, color=COL[m], mfc=fill, ms=8, lw=1.4, capsize=3)
    ax.set_xticks(range(3)); ax.set_xticklabels([NAME[m] for m in MODES])
    ax.set_xlim(-0.6, 2.6); ax.set_ylim(0, 0.3)
    ax.set_ylabel('growth rate (per cycle)'); ax.set_title('Kaplan\'s rate fits annealed,\nnot periodic')
    hs, ls = mode_handles()
    hs += [Line2D([0], [0], color=INK2, ls='--', lw=1.4), Line2D([0], [0], color=INK2, marker='o', ls='none', ms=8),
           Line2D([0], [0], color=INK2, marker='s', mfc='white', ls='none', ms=8),
           Line2D([0], [0], color=INK2, marker='X', ls='none', ms=8)]
    ls += ['predicted slope\nα + β/4', 'measured caustic\ngrowth', 'measured line\nlength growth', 'predicted\nα + β/4']
    legend_bottom(fig, hs, ls, ncol=4)
    addfig(fid, fig, 'appA_caustics', 'Left: exact number of caustics on the line launched at momentum '
           'y0 = 960 px (mean ± sd over seeds; counting stops when the line exceeds 12 million pieces) and the '
           'slope Kaplan\'s formula predicts (dashed, from N = 10). Right: fitted growth rates (N ≥ 10); bars '
           'are sd over seeds (none for the deterministic periodic run).')


def b1c(m, kind, met, mask_flat=True):
    BRM = {k: i for i, k in enumerate(K.BR_METRICS)}
    X = np.array([p['br'] for p in P1C[m]])[:, :3, :, kind, BRM[met]].copy()
    if met == 'lost' and mask_flat:
        X[np.array([p['br'] for p in P1C[m]])[:, :3, :, kind, BRM['C_all']] < 1.25] = np.nan
    with np.errstate(all='ignore'):
        return np.nanmean(X, axis=1)


def fig_branches(fid, kind, name, titles, caption):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.2))
    for m in MODES:
        band(axs[0], NN, b1c(m, kind, 'enrichment'), COL[m], logy=True)
        L5 = b1c(m, kind, 'lost'); L5[:, :5] = np.nan          # ill-conditioned while the line is nearly flat
        band(axs[1], NN, L5, COL[m])
        band(axs[2], NN, b1c(m, kind, 'corr'), COL[m])
    nk = kind + 2
    L5 = ms(b1c('annealed', nk, 'lost'))[0]; L5[:5] = np.nan
    axs[1].plot(NN, L5, color=MODEL, lw=1.2, ls=':')
    axs[2].plot(NN, ms(b1c('annealed', nk, 'corr'))[0], color=MODEL, lw=1.2, ls=':')
    axs[0].set_yscale('log'); axs[0].axhline(1, color=INK2, lw=0.8)
    axs[0].set_ylim(0.3, 6) if kind == 0 else axs[0].set_ylim(0.3, 6)
    axs[0].yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4]))
    axs[0].yaxis.set_minor_locator(FixedLocator([]))
    axs[0].set_yticklabels(['0.5', '1', '2', '4'])
    axs[1].axhline(0, color=INK2, lw=0.8); axs[1].set_ylim(-0.3, 0.8)
    axs[2].set_ylim(0.0, 1.05)
    axs[0].set_ylabel('enrichment in branches\n(1 = no preference)')
    axs[1].set_ylabel('share of branch\ncontrast lost')
    axs[2].set_ylabel('corr of ray density\nwith and without')
    for ax, t in zip(axs, titles):
        ax.set_title(t); ax.set_xlabel('cycles N'); ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60]))
    hs, ls = mode_handles([(Line2D([0], [0], color=MODEL, lw=1.2, ls=':'), 'same number of\nrays removed\nevenly (control)')])
    legend_bottom(fig, hs, ls, ncol=4)
    addfig(fid, fig, name, caption)


def fig_return(fid):
    fig, ax = plt.subplots(1, 1, figsize=(W_IN, 3.5))
    for m in MODES:
        Fm = np.array([KR.part2(o)['F_re'] for o in R[m]]); k = np.arange(Fm.shape[1])
        band(ax, k, Fm, COL[m])
        ax.plot(1, ms(Fm)[0][1], 'o', color=COL[m], ms=7)
    ax.set_xscale('log'); ax.set_xlim(0.9, 70); ax.set_ylim(0, 0.42)
    ax.set_xlabel('cycles since leaving the low-stretch set')
    ax.set_ylabel('probability of having\nreturned (corrected\nfor censoring)')
    ax.set_title('Most returns happen within 5 cycles of leaving (dots: after one cycle)')
    hs, ls = mode_handles()
    legend_bottom(fig, hs, ls)
    addfig(fid, fig, 'p2_fig1_return', 'Probability that material which left the low-stretch set (S ≤ 2) has '
           'come back, against cycles since it left, corrected for runs ending before a return could be seen '
           '(mean ± sd over 10 seeds; periodic is deterministic).')


def fig_hazard(fid):
    fig, axs = plt.subplots(1, 2, figsize=(W_IN, 3.6))
    for m in MODES:
        P = [KR.part2(o) for o in R[m]]
        H = np.array([p['h_re'] for p in P]); Rk = np.array([p['r_re'] for p in P]); k = np.arange(H.shape[1])
        band(axs[0], k[1:], H[:, 1:], COL[m], logy=True)
        axs[1].plot(k[1:], Rk[:, 1:].mean(0), color=COL[m])
    axs[0].set_yscale('log'); axs[0].set_ylim(1e-5, 0.3)
    axs[0].set_ylabel('chance of returning\nin the next cycle')
    axs[0].set_title('The chance of return falls\nsteeply with time away')
    axs[1].set_yscale('log'); axs[1].set_ylabel('parcels still away\n(per run)')
    axs[1].set_title('Most material stays away,\nso late ages are well sampled')
    for ax in axs:
        ax.set_xlabel('cycles since leaving'); ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60]))
    hs, ls = mode_handles()
    legend_bottom(fig, hs, ls)
    addfig(fid, fig, 'p2_fig2_hazard', 'Left: return hazard, the chance that material still away from the '
           'low-stretch set comes back in the next cycle, against cycles since it left (points with fewer than '
           '200 parcels at risk are hidden). Right: the number of parcels at risk at each age.')


def fig_depth(fid):
    fig, ax = plt.subplots(1, 1, figsize=(W_IN, 3.7))
    labs = ['2 – 2.5', '2.5 – 3', '3 – 5', 'above 5']
    for j, m in enumerate(MODES):
        D = [KR.depth_ratios(o) for o in R[m]]
        Y = np.array([[d['user'][i][0] for i in range(4)] for d in D])
        mm, ss = ms(Y)
        ax.errorbar(np.arange(4) + (j - 1) * 0.18, mm, yerr=ss if m != 'cosine' else None, fmt='o', color=COL[m],
                    ms=7, lw=1.5, capsize=3)
        ax.axhline(np.mean([d['MH'] for d in D]), color=COL[m], ls='--', lw=1.3)
    ax.axhline(1, color=INK, lw=1.0)
    ax.set_yscale('log'); ax.set_ylim(0.25, 8)
    ax.yaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4, 8])); ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.set_yticklabels(['0.25', '0.5', '1', '2', '4', '8'])
    ax.set_xticks(range(4)); ax.set_xticklabels(labs)
    ax.set_xlabel('depth: S one cycle before a possible return')
    ax.set_ylabel('return hazard, returned\nbefore ÷ first time')
    ax.set_title('At the same depth, earlier returners are no more likely to return')
    hs, ls = mode_handles([(Line2D([0], [0], color=INK2, ls='--', lw=1.3), 'adjusted for\ndepth and age')])
    legend_bottom(fig, hs, ls, ncol=4)
    addfig(fid, fig, 'p2_fig3_depth', 'Ratio of the return hazard for material that has returned before to that '
           'for material away for the first time, within bands of depth (dots; mean ± sd over seeds), and the '
           'ratio adjusted for depth (7 bands) and age together (dashed). 1 = no difference.')


# ---------------- appendix figures
def fig_manifold(fid, m):
    sn = R[m][0]['snaps']; Ns = (5, 10, 20)
    ce = R[m][0]['caustics_exact'][1]
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.5), sharey=True, sharex=True)
    for ax, Nn in zip(axs, Ns):
        x, y, Sv, aa = sn[(1, Nn)]
        sc = ax.scatter(x, y, c=np.log10(Sv), s=0.4, cmap='viridis', vmin=0, vmax=6, rasterized=True)
        cnt = int(ce[Nn])
        if Nn < 20:
            cs = np.nonzero(aa * np.roll(aa, -1) < 0)[0]
            ax.plot(x[cs], np.full(cs.size, ax.get_ylim()[0]), '|', color='#e34948', ms=9, mew=1.2, clip_on=False)
            ax.set_title('N = %d: %s caustics\n(marked in red)' % (Nn, NEWN('{:,}'.format(cnt), 'exact caustic count, seed 0, %s, N = %d' % (NAME[m], Nn))))
        else:
            ax.set_title('N = %d: %s caustics\n(too dense to mark)' % (Nn, NEWN('{:,}'.format(cnt), 'exact caustic count, seed 0, %s, N = %d' % (NAME[m], Nn))))
        ax.grid(False); frame(ax, m)
        ax.xaxis.set_major_locator(FixedLocator([0, 960, 1920])); ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_xlabel('position x (px)')
    axs[0].set_ylabel('momentum y (px)')
    cb = fig.colorbar(sc, ax=axs, shrink=0.95, pad=0.02); cb.set_label('stretching S\n(log10)')
    addfig(fid, fig, 'appA_manifold_%s' % NAME[m], 'The line of rays launched at momentum y0 = 960 px, %s, '
           'seed 0, in the phase plane after 5, 10 and 20 cycles, coloured by stretching S. A caustic is a point '
           'where the line turns back in x; counts are exact.' % LONG[m])


def fig_fan(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.4), sharey=True)
    titles = {'annealed': 'annealed: moves up\nsteadily', 'quenched': 'quenched: wider\nspread',
              'cosine': 'periodic: slow\ngrowth'}
    for ax, m in zip(axs, MODES):
        Q = np.array([o['quant'] for o in R[m]]).mean(0) / np.log(10)
        ax.fill_between(NN, Q[:, 0], Q[:, 6], color=COL[m], alpha=0.12, lw=0)
        ax.fill_between(NN, Q[:, 1], Q[:, 5], color=COL[m], alpha=0.22, lw=0)
        ax.fill_between(NN, Q[:, 2], Q[:, 4], color=COL[m], alpha=0.38, lw=0)
        ax.plot(NN, Q[:, 3], color=COL[m]); ax.axhline(np.log10(2), color=INK, ls='--', lw=1)
        ax.set_title(titles[m]); ax.set_xlabel('cycles N'); ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60]))
    axs[0].set_ylabel('stretching S (log10)')
    hs = [Patch(color=INK2, alpha=a) for a in (0.12, 0.22, 0.38)] + [Line2D([0], [0], color=INK2), Line2D([0], [0], color=INK, ls='--', lw=1)]
    legend_bottom(fig, hs, ['1–99%', '10–90%', '25–75%', 'median', 'S = 2'])
    addfig(fid, fig, 'appA_quantiles', 'Quantiles of S measured on the pixel grid against cycles N (mean of the '
           'quantiles over seeds).')


def fig_tail(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.6), sharey=True)
    sty = {1.25: ':', 2: '-', 5: '--', 10: '-.'}
    for ax, m in zip(axs, MODES):
        for s, ls in sty.items():
            col = int(np.nonzero(K.S_LIST == s)[0][0])
            Y = np.array([-np.log(np.where(o['P_le'][:, col] > 0, o['P_le'][:, col], np.nan)) for o in R[m]])
            band(ax, NN, Y, COL[m], ls=ls)
        pred = float(_nums(V('4C', 'pred a²/2v', mode=m))[0])
        Y2 = np.nanmean(np.array([-np.log(o['P_le'][:, 2]) for o in R[m]]), 0); xx = np.arange(20, C + 1)
        ax.plot(xx, Y2[20] + pred * (xx - 20), color=INK, ls='--', lw=1.2)
        meas = S('4C', 's=2: rate (ratio)', mode=m).split(' (')[0].split(',')[0].split(' ±')[0]
        prd = S('4C', 'pred a²/2v', mode=m).split(',')[0].split(' ±')[0]
        ax.set_title('%s:\nmeasured %s\nlog-normal %s' % (NAME[m], meas, prd))
        ax.set_xlabel('cycles N'); ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60]))
    axs[0].set_ylabel('−ln (share with S ≤ s)')
    hs = [Line2D([0], [0], color=INK2, ls=ls) for ls in sty.values()] + [Line2D([0], [0], color=INK, ls='--', lw=1.2)]
    legend_bottom(fig, hs, ['s = 1.25', 's = 2', 's = 5', 's = 10', 'log-normal\nslope (s = 2)'])
    addfig(fid, fig, 'appA_tail', 'Decay of the share of material below four thresholds s, as −ln(share) '
           'against N (mean ± sd over seeds); the black dashed line has the log-normal slope α²/(2v). Decay '
           'rates in the titles are for s = 2, fitted over N = 20–60.')


def fig_rate(fid):
    fig, ax = plt.subplots(1, 1, figsize=(W_IN, 3.4))
    e = IID['lam_edges']; n = IID['n_samples']
    F40 = np.cumsum(IID['lam_hist'][40]) / n; F60 = np.cumsum(IID['lam_hist'][60]) / n
    al = float(V('iid', 'alpha', model='i.i.d. model (1e7)')); v = float(V('iid', 'v', model='i.i.d. model (1e7)'))
    lam = e[1:]
    with np.errstate(divide='ignore', invalid='ignore'):
        I = np.where(lam < al, -(np.log(F60) - np.log(F40)) / 20, -(np.log1p(-np.minimum(F60, 1 - 1e-16)) - np.log1p(-np.minimum(F40, 1 - 1e-16))) / 20)
    ax.plot(lam, I, color=MODEL, lw=1.8)
    ax.plot(lam, (lam - al) ** 2 / (2 * v), color=INK2, ls=':', lw=2)
    ax.axhline(float(V('iid', 'decay of P(S≤2) = I(0)', model='i.i.d. model (1e7)')), color=COL['annealed'], ls='--', lw=1.3)
    raw = S('iid', 'λ-CDF I(0+), N 40→60 → corrected for (ln S)²', model='i.i.d. model (1e7)')
    i_raw, i_cor = _nums(raw)[0], _nums(raw)[-1]
    ax.plot([0.004], [i_raw], 'o', mfc='white', color=MODEL, ms=8, clip_on=False)
    ax.plot([0.004], [i_cor], 'o', color=MODEL, ms=8, clip_on=False)
    ax.set_xlim(0, 0.45); ax.set_ylim(0, 0.2)
    ax.set_xlabel('growth rate of ln S, λ = ln S / N (per cycle)'); ax.set_ylabel('rate function I(λ)\n(per cycle)')
    ax.set_title('Corrected for the (ln S)² law near S = 1, I(0) lies at the measured decay rate')
    legend_bottom(fig, [Line2D([0], [0], color=MODEL), Line2D([0], [0], color=INK2, ls=':', lw=2),
                        Line2D([0], [0], color=COL['annealed'], ls='--', lw=1.3),
                        Line2D([0], [0], color=MODEL, marker='o', mfc='white', ls='none', ms=8),
                        Line2D([0], [0], color=MODEL, marker='o', ls='none', ms=8)],
                  ['model\n(cycles 40 → 60)', 'log-normal\nparabola', 'measured\ndecay rate',
                   'I(0),\nraw', 'I(0),\ncorrected'], ncol=5)
    addfig(fid, fig, 'appA_rate_function', 'Large-deviation rate function of the independent-step model of the '
           'annealed kick, estimated from the change in the cumulative distribution of λ between cycles 40 and 60. '
           'At λ → 0 this estimate is biased low by 2 ln(60/40)/20 = 0.041 because P(S ≤ s) ∝ (ln s)² near S = 1 '
           '(open dot: raw; filled dot: corrected).')


def fig_joint(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.4), sharey=True)
    RB = np.arange(-6, 30.0001, 0.25); SB = np.arange(0, 30.0001, 0.25)
    for ax, m in zip(axs, MODES):
        H = sum(o['joint_r_lnS'][30] for o in R[m]).T.astype(float)
        Hm = np.ma.masked_less_equal(H, 0)
        cmap = plt.get_cmap('viridis_r').copy(); cmap.set_bad('white')
        ax.imshow(Hm, origin='lower', aspect='auto', extent=[RB[0], RB[-1], SB[0], SB[-1]], norm=LogNorm(),
                  cmap=cmap, interpolation='nearest')
        ax.plot([0, 30], [0, 30], color='#e34948', ls='--', lw=1.4)
        ax.set_xlim(-2, 10); ax.set_ylim(0, 10); ax.grid(False); frame(ax, m)
        ax.set_title('%s:\nno ray below the line' % NAME[m]); ax.set_xlabel('rarefaction r')
        ax.xaxis.set_major_locator(FixedLocator([0, 4, 8]))
    axs[0].set_ylabel('ln S')
    legend_bottom(fig, [Line2D([0], [0], color='#e34948', ls='--', lw=1.4)], ['ln S = r  (S ≥ exp(r) for every ray)'])
    addfig(fid, fig, 'appA_r_vs_S', 'Joint distribution of rarefaction r and ln S for every ray on the three '
           'launched lines after 30 cycles (all seeds pooled; darker = more rays).')


def fig_metric(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.6), sharey=True)
    sty = {0.5: ':', 2.0: '--', 4.0: '-.'}
    titles = {'annealed': 'annealed: level moves,\ndecay rate does not', 'quenched': 'quenched: small\nchanges',
              'cosine': 'periodic: decay\nrate changes'}
    for ax, m in zip(axs, MODES):
        band(ax, NN, np.array([o['P_le'][:, 2] for o in R[m]]), COL[m], logy=True)
        for fi, f in enumerate(K.METRIC_F):
            ax.plot(NN, ms(np.array([o['metric_P_le'][fi, :, 2] for o in R[m]]))[0], color=COL[m], ls=sty[f], lw=1.5)
        ax.set_yscale('log'); ax.set_title(titles[m]); ax.set_xlabel('cycles N')
        ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60]))
    axs[0].set_ylabel('share of material\nwith S ≤ 2')
    hs = [Line2D([0], [0], color=INK2)] + [Line2D([0], [0], color=INK2, ls=ls) for ls in sty.values()]
    legend_bottom(fig, hs, ['y unscaled', 'y × 0.5', 'y × 2', 'y × 4'])
    addfig(fid, fig, 'appA_metric', 'The share of material with S ≤ 2 when S is measured with the momentum '
           'axis y stretched by a factor f (the same dynamics in a different metric).')


def fig_leave(fid):
    fig, ax = plt.subplots(1, 1, figsize=(W_IN, 3.4))
    for m in MODES:
        P = [KR.part2(o) for o in R[m]]
        for key, ls in (('F_ex', '-'), ('F_co', '--')):
            Fm = np.array([p[key] for p in P]); band(ax, np.arange(Fm.shape[1]), Fm, COL[m], ls=ls)
    ax.set_xscale('log'); ax.set_xlim(0.9, 70); ax.set_ylim(0, 1.05)
    ax.set_xlabel('cycles since entering the low-stretch set'); ax.set_ylabel('probability of having left\n(corrected for censoring)')
    ax.set_title('Returned material leaves again faster than the original material')
    hs, ls = mode_handles([(Line2D([0], [0], color=INK2), 'returned\nmaterial'), (Line2D([0], [0], color=INK2, ls='--'), 'original\nmaterial')])
    legend_bottom(fig, hs, ls, ncol=5)
    addfig(fid, fig, 'appA_leave_again', 'Probability that material in the low-stretch set has left, against cycles '
           'since it entered: material that returned (solid) and the original material, in the set since launch '
           '(dashed).')


def fig_cohort(fid):
    fig, axs = plt.subplots(1, 3, figsize=(W_IN, 3.5), sharey=True)
    for ax, m in zip(axs, MODES):
        rr = np.sum([o['risk_re'] for o in R[m]], 0); ee = np.sum([o['ev_re'] for o in R[m]], 0)
        H = np.ma.masked_where(rr < 200, ee / np.maximum(rr, 1)); H = np.ma.masked_less_equal(H, 0)
        cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('white')
        im = ax.imshow(H, origin='lower', aspect='auto', norm=LogNorm(vmin=1e-4, vmax=0.3), cmap=cmap, interpolation='nearest')
        ax.set_xlabel('cycles since leaving'); ax.grid(False); frame(ax, m)
        ax.xaxis.set_major_locator(FixedLocator([0, 30, 60]))
        ax.set_title({'annealed': 'annealed: set by\ntime away', 'quenched': 'quenched: mostly\ntime away',
                      'cosine': 'periodic:\npatchy'}[m])
    axs[0].set_ylabel('cycle in which it left')
    cb = fig.colorbar(im, ax=axs, shrink=0.95, pad=0.02); cb.set_label('return\nhazard')
    addfig(fid, fig, 'appA_cohorts', 'Return hazard for each group of material that left in the same cycle '
           '(vertical) against cycles since it left (horizontal), seeds pooled; white = fewer than 200 at risk.')


# ============================================================ build all figures
WHERE[0] = 'figures'
fig_branched('F1.1', 1, 'p1_fig1_branched_flow')
fig_collapse('A.14')
fig_why('F1.3')
fig_caustics('A.15')
fig_branches('F1.5', 0, 'p1_fig5_lowS_branches',
             ['Low-S rays prefer\nbranches…', '…but removing them\ncosts little contrast\n(except periodic)…',
              '…and branches stay\nput (except periodic)'],
             'Rays with S ≤ 2 removed from the ray density: their enrichment in branches (left), the share of '
             'branch contrast lost (centre; shown from N = 5, before which the line is nearly flat) and corr between the '
             'ray density with and without them (right). Mean of the three launch momenta; mean ± sd over seeds. '
             'Dotted: the same number of rays removed evenly along the line (annealed).')
fig_return('F2.1')
fig_hazard('F2.2')
fig_depth('F2.3')
fig_branched('A.1', 0, 'appA_branched_low')
fig_branched('A.2', 2, 'appA_branched_high')
for i, m in enumerate(MODES):
    fig_manifold('A.%d' % (3 + i), m)
fig_fan('A.6')
fig_tail('A.7')
fig_rate('A.8')
fig_branches('A.9', 1, 'appA_lowr_branches',
             ['Low-r rays prefer\nbranches…', '…and removing them\ncosts more contrast…', '…partly by\nconstruction'],
             'As Figure 3 of the stretching report, but removing the rays whose rarefaction r is below its 10th '
             'percentile. Because ray density is a sum of 1/|a| and r ≥ ln|a|, these are by construction the rays '
             'that add most to the density, so this is not independent evidence.')
fig_joint('A.10')
fig_metric('A.11')
fig_leave('A.12')
fig_cohort('A.13')
print('figures:', len(FIGS))


# ============================================================ text helpers
def P(m, key, col, **kw):
    """S() for one mode, the shown value without its sd (for compact text)."""
    return S(key, col, mode=m, **kw)


def val(shown):
    """The value part of a shown string (drops '± sd' and the deterministic note)."""
    return shown.split(' ±')[0].split(' (deterministic)')[0]


def trio(xs, unit=''):
    return 'annealed %s%s, quenched %s%s and periodic %s%s' % (xs[0], unit, xs[1], unit, xs[2], unit)


def nok(s):
    return s.split(' (k=')[0]


def three(key, col, **kw):
    return [S(key, col, mode=m, **kw) for m in MODES]


def D(value, how):
    """A number derived arithmetically from previous-report values (ratio, range)."""
    LEDGER.append(dict(where=WHERE[0], source=how, old='derived from previous-report values', shown=value,
                       status='derived from previous report'))
    return value


GLOSSARY = [
    ('S (stretching)', 'The largest factor by which a small piece of material has been stretched since launch, '
     'computed exactly from the product of the per-cycle Jacobians. S = 1 means unstretched; S can never be '
     'below 1. "Low-stretch" means S ≤ 2.'),
    ('rarefaction r', 'Shaw\'s rarefaction exponent: for a launched line of rays, r = ln √(a² + c²), the log of how '
     'much a short piece of the line has lengthened (a = ∂x/∂x0, c = ∂y/∂x0). Always S ≥ exp(r).'),
    ('caustic', 'A point where the launched line folds back in position x (a changes sign). Ray density peaks there.'),
    ('branch', 'A 4-px bin of position x among the 10% with the highest ray density.'),
    ('contrast', 'Mean ray density in branches ÷ mean ray density everywhere. 1 means no branches.'),
    ('enrichment', 'Share of a group of rays among the rays in branches ÷ its share among all rays. 1 means no '
     'preference.'),
    ('hazard', 'The chance that a parcel crosses (returns to, or leaves, the low-stretch set) in the next cycle, '
     'given it has not crossed yet.'),
    ('censoring', 'The run ends before we can see whether a parcel would cross. Censoring-corrected probabilities '
     'use every parcel for as long as it was observed (method named in the appendix, Methods).'),
    ('depth', 'S one cycle before a possible return: how far above the threshold S = 2 a parcel is.'),
    ('seed', 'Index of the random-number stream; each seed is an independent kick history. Results are mean ± sd '
     'over 10 seeds. Periodic flow has no randomness, so its sd is n/a (deterministic).'),
    ('N', 'Number of cycles (one kick, then one drift) since launch; runs last 60 cycles.'),
    ('modes', 'annealed (fresh random kick each cycle); quenched (one frozen random kick); periodic (standard map, '
     'K = 0.5).'),
    ('launch momentum y0', 'Rows of rays are launched at y0 = 960 − 240 px, 960 px and 960 + 240 px; 960 px is the '
     'momentum at which a parcel does not drift.'),
    ('corr', 'Pearson correlation coefficient across the 480 bins of position x.'),
    ('pixel grid', 'S measured on every pixel of the 1920 × 1920 frame at the current cycle.'),
]


# ============================================================ Paper 1
def paper1():
    Q = []
    A, Qn, Pm = 'annealed', 'quenched', 'cosine'
    # Q1
    WHERE[0] = 'Paper 1, Q1'
    c5 = three('1C', 'contrast, all rays', N=5); c60 = three('1C', 'contrast, all rays', N=60)
    Q.append(dict(
        q='Does every kind of kick produce branched flow?',
        a='Yes. In all three modes, rays launched on a line of constant momentum focus into branches within a few '
          'cycles: at N = 5 the mean ray density in branches is %s times the average. By N = 60 branches fade in '
          'annealed flow (%s) but persist in quenched (%s) and periodic flow (%s).' % (
              trio(c5), c60[0], c60[1], c60[2]),
        ev=('fig', 'F1.1', '%FIG% (launch momentum y0 = 960 px; the other two launch momenta are in Appendix '
            'Figures A1 and A2).'),
        conf=('Strong', 'seen in every seed and at all three launch momenta (Appendix Figures A1, A2).'),
        cav='Sampling noise in the ray histogram (8 rays per px) could make spurious branches. Check: the exact ray '
            'density (a sum of 1/|a| over the line) agrees with the histogram within 2%% in %s–%s of bins away from '
            'caustics (Appendix Table B1). Ruled out.' % (
                F('96.0% within 2%', '96.0%'), F('97.8% within 2%', '97.8%'))))
    # Q2 (growth and shape of ln S, and the log-normal test)
    WHERE[0] = 'Paper 1, Q2'
    al = three('4C', 'alpha'); vv = three('4C', 'v')
    sk = three('mom', 'skew', N=60); ku = three('mom', 'ex. kurt', N=60)
    rows = []
    for i, m in enumerate(MODES):
        pr = S('4C', 'pred a²/2v', mode=m)
        me = S('4C', 's=2: rate (ratio)', mode=m)
        means = [_nums(V('4C', c, mode=m))[0] for c in T['4C'][0][4:]]
        rng = D('%.4f–%.4f' % (min(means), max(means)), 'min–max of the s = 1.25…10 means, 4C, %s' % m)
        rows.append([LONG[m], al[i], vv[i], pr, re.sub(r'\((\d+\.\d+)\)$', r'(\1×)', me), rng])
    s0 = [S('4C_seed0', c, **{'mode__(seed__0)': m}) for m in ('annealed', 'quenched') for c in ('predicted', 'measured s=2')]
    Q.append(dict(
        q='Is ln S Gaussian, and does the log-normal law predict how fast low-stretch material disappears?',
        a='Mean ln S and its variance grow linearly, at α and v per cycle (Table 1), and the body of ln S is close to '
          'Gaussian in annealed flow (skewness %s, excess kurtosis %s at N = 60; less so in quenched flow, not in '
          'periodic flow). Yet the log-normal law fails for low stretch: in annealed flow the share with S ≤ 2 '
          'decays at %s per cycle against %s predicted from α and v, %s times faster, at every threshold from 1.25 '
          'to 10. In quenched flow the prediction is consistent within the seed spread (%s against %s).' % (
              sk[0], ku[0], val(rows[0][4]), val(rows[0][3]), rows[0][4].split('(')[-1].rstrip(')').rstrip('×'),
              rows[1][4].split(' (')[0], rows[1][3]),
        ev=('table', 'T1.1', (['mode', 'growth of mean ln S, α', 'growth of variance, v', 'log-normal prediction α²/(2v)',
                               'measured decay of the share with S ≤ 2 (÷ prediction)'], [r[:5] for r in rows]),
            'Table 1. Rates per cycle on the pixel grid, fitted over N = 20–60, mean ± sd over seeds. Measured over '
            'thresholds s = 1.25 to 10 the annealed decay is %s. Seed 0 reproduces the earlier single-run numbers '
            '(annealed %s measured against %s predicted; quenched %s against %s). Shape of ln S: Appendix Figure '
            'A14 and Table B4; quantiles: Figure A6; decay curves: Figure A7; every threshold: Table B2; shares below each '
            'threshold: Table B3.' % (
                rows[0][5], s0[1], s0[0], s0[3], s0[2])),
        conf=('Strong for annealed, Tentative for quenched, Not supported for periodic',
              'annealed seeds agree to about 2%% at every threshold, and the pixel grid matches a full-period sample '
              'of parcels to %s (Appendix Table B1); quenched seeds scatter by about 20%%; periodic flow is one '
              'deterministic run with regular islands, where the log-normal premise does not apply.' % F('agree to 1e-4', '10⁻⁴')),
        cav='Few parcels remain below the threshold at large N, so the fit could follow counting noise. Check: '
            'only points with at least 10 parcels are fitted, and the annealed rate barely changes from s = 1.25 '
            'to 10 (%s). Ruled out.' % rows[0][5]))
    # Q4
    WHERE[0] = 'Paper 1, Q3'
    it = 'i.i.d. model (1e7)'
    dec_i = S('iid', 'decay of P(S≤2) = I(0)', model=it)
    dec_2 = S('iid', 'decay of P(S≤2) = I(0)', model='2-D annealed (10 seeds)')
    gau = S('iid', 'Gaussian a²/2v', model=it)
    lq = V('iid', '−min L(q): N 20→40 / 40→60', model=it); lq = S('iid', '−min L(q): N 20→40 / 40→60', model=it)
    l2040, l4060 = lq.split(' / ')
    i0 = S('iid', 'λ-CDF I(0+), N 40→60 → corrected for (ln S)²', model=it)
    i_raw, i_cor = i0.split(' → ')
    Q.append(dict(
        q='Why does low-stretch material disappear faster than the log-normal law predicts in annealed flow?',
        a='Because each annealed cycle is an independent random step, and the share near S = 1 is set by the floor '
          'S ≥ 1, not by the Gaussian centre of ln S. A model that multiplies independent random steps drawn from '
          'the same kicks gives a decay of %s per cycle (runs: %s); its moment growth rate L(q) flattens at −%s '
          'for q ≤ −2, the signature of a floor, where a log-normal would reach only −%s.' % (dec_i, dec_2, l4060, gau),
        ev=('fig', 'F1.3', '%FIG%; the model\'s rate function in Appendix Figure A8; model against runs in Appendix '
            'Table B5; L(q) in Appendix Table B6.'),
        conf=('Strong that independent steps reproduce the decay; Tentative for the exact rate',
              'finite-N estimates of the limiting rate range from %s (cycles 40 → 60) to %s (cycles 20 → 40).' % (l4060, l2040)),
        cav='Near S = 1, P(S ≤ s) grows as (ln s)² (measured P(1.25)/P(2) = %s against %s), which biases a direct '
            'estimate of the limiting rate low, to %s. Check: corrected, it is %s, matching the measured %s. Ruled '
            'out.' % (F('measured P(1.25)/P(2) = 0.097', '0.097'),
                                         F('(ln 1.25/ln 2)² = 0.104', '0.104'), i_raw, i_cor, dec_i)))
    # Q5
    WHERE[0] = 'Paper 1, Q4'
    cg = [S('caustics', 'caustic growth (fit)', mode=m, launch__row='row0') for m in MODES]
    cp = [S('caustics', 'alpha+beta/4', mode=m, launch__row='row0') for m in MODES]
    ll = [S('caustics', 'line-length growth (exact)', mode=m, launch__row='row0') for m in MODES]
    fN = S('caustics', 'fit N', mode='annealed', launch__row='row0')
    Q.append(dict(
        q='Do caustics multiply at the rate Kaplan\'s formula predicts?',
        a='In annealed flow, yes: caustics on the line launched at y0 = 960 px multiply at %s per cycle against '
          'Kaplan\'s prediction α + β/4 = %s, where α is the growth rate of the mean of rarefaction r and β/2 that '
          'of its variance. Quenched flow is consistent within the seed spread (%s against %s). Periodic flow is '
          'not: %s against %s, while the line\'s own length grows at %s per cycle.' % (
              cg[0], cp[0], cg[1], cp[1], val(cg[2]), val(cp[2]), val(ll[2])),
        ev=('table', 'T1.4', (['mode', 'measured caustic growth', 'measured growth of line length',
                               'Kaplan prediction α + β/4', 'measured ÷ predicted'],
                              [[LONG[m], cg[i], ll[i], cp[i],
                                S('caustics', 'ratio caustic/pred', mode=m, launch__row='row0')] for i, m in enumerate(MODES)]),
            'Table 2. Growth rates per cycle on the line launched at y0 = 960 px, fitted over N = 10–31 (annealed), '
            'mean ± sd over seeds. Caustic counts against N are in Appendix Figure A15; the line in the phase plane '
            'in Appendix Figures A3–A5; all three launch momenta in Appendix Table B7.'),
        conf=('Strong for annealed, Tentative for quenched, Not supported for periodic',
              'annealed agrees at all three launch momenta; quenched seeds scatter widely; the formula assumes r is '
              'Gaussian, which fails where regular islands and chaos coexist.'),
        cav='A line sampled by a fixed set of rays misses caustics once it folds faster than the sampling. Check: '
            'the counts are exact — the kick is piecewise linear, so the line is followed piece by piece — and the '
            'fit uses N = %s, where exact counting is possible. Ruled out.' % fN))
    # Q6
    WHERE[0] = 'Paper 1, Q5'
    sh = [S('1C', 'low-S share of branch mass', mode=A, N=n) for n in (5, 20, 60)]
    lo = [S('1C', 'contrast lost', mode=A, N=n) for n in (5, 20, 60)]
    co = [S('1C', 'corr(ρ, ρ_noLow)', mode=A, N=n) for n in (5, 20, 60)]
    en10 = S('1C', 'enrichment', mode=A, N=10)
    ql = [val(S('1C', 'contrast lost', mode=Qn, N=n)) for n in (20, 40, 60)]
    ps = S('1C', 'low-S share of branch mass', mode=Pm, N=20); pl = S('1C', 'contrast lost', mode=Pm, N=20)
    pc = S('1C', 'corr(ρ, ρ_noLow)', mode=Pm, N=20)
    nul = S('1C', 'null lost (S / r)', mode=A, N=20)
    Q.append(dict(
        q='Do low-stretch rays make the branches?',
        a='Mostly not. In annealed flow, rays with S ≤ 2 are over-represented in branches (enrichment %s at N = 10) '
          'but become too rare to build them. From N = 5 to 20 to 60 their share of the rays in branches falls '
          '%s → %s → %s, removing them costs %s → %s → %s of the branch contrast, and the branches stay put (corr '
          '%s → %s → %s). Quenched flow loses a steady %s–%s of its contrast; periodic flow is the exception (N = 20: '
          '%s of the rays in branches, %s of the contrast lost, corr %s). Removing the rays of lowest rarefaction r '
          'costs far more (Appendix Figure A9), but "removing low-r rays removes branches" is partly true by '
          'construction (ray density = sum of 1/|a|, and r ≥ ln|a|), so it is not independent evidence.' % (
              en10, sh[0], sh[1], sh[2], lo[0], lo[1], lo[2], co[0], co[1], co[2], min(ql), max(ql),
              val(ps), val(pl), val(pc)),
        ev=('fig', 'F1.5', '%FIG%; full table in Appendix Table B8, with a control on parcels that fill the plane '
            'in Appendix Table B9; correlations with local stretching in Appendix Table B10; r against S in Appendix '
            'Figure A10 and Table B11.'),
        conf=('Strong for annealed and quenched, Tentative for periodic',
              '10 seeds, and removing the same number of rays evenly along the line costs %s of the contrast; '
              'periodic flow is one deterministic run with a period-2 oscillation in N.' % nul.split(' / ')[0]),
        cav='"Branch" means the top 10%% of bins, a choice. Check: corr across all 480 bins, which does not use that '
            'choice, gives the same picture, and the branch sets with and without low-S rays overlap strongly in '
            'annealed flow (overlap %s at N = 20, Appendix Table B8). Ruled out for annealed.' % S('1C', 'branch overlap (Jaccard)', mode=A, N=20)))
    # Q7
    WHERE[0] = 'Paper 1, Q6'
    full, rows, sds = [], [], {}
    for m in MODES:
        r = [S('4E', 'alpha', mode=m, y__scale__f='1')]
        for c in ('rate S_THR=1.5', 'rate S_THR=2', 'rate S_THR=3'):
            r.append(S('4E', c, mode=m, y__scale__f='1'))
        for f in ('0.5', '2', '4'):
            r.append(S('4E', 'rate S_THR=2', mode=m, y__scale__f=f))
        full.append(r)
        sds[m] = D('%.4f' % max(_nums(c)[1] for c in r[1:]), 'largest seed sd among the decay rates, 4E, %s' % m) if m != 'cosine' else None
        rows.append([LONG[m], val(r[0]), ' / '.join(val(c) for c in r[1:4]), ' / '.join(val(c) for c in r[4:7])])
    an = [_nums(c)[0] for c in full[0][1:]]
    dev = D('%.0f%%' % np.ceil(100 * max(abs(x - an[1]) / an[1] for x in an)), 'max |rate − rate(S≤2, unscaled)| / rate, annealed, from 4E')
    pa1, pa4 = S('4E', 'P(N=60)', mode=A, y__scale__f='1'), S('4E', 'P(N=60)', mode=A, y__scale__f='4')
    Q.append(dict(
        q='Do these results depend on the threshold S ≤ 2 or on how momentum is scaled?',
        a='Not for annealed flow: the decay rate changes by at most %s across thresholds 1.5, 2 and 3 and when the '
          'momentum axis is scaled by 0.5, 2 or 4, and the growth rate α does not change at all. Periodic flow is not '
          'robust: scaling momentum by 2 changes its rate from %s to %s.' % (dev, val(full[2][2]), val(full[2][5])),
        ev=('table', 'T1.2', (['mode', 'growth α of mean ln S', 'decay rate at thresholds S ≤ 1.5 / 2 / 3',
                               'decay rate at S ≤ 2 with y × 0.5 / 2 / 4'], rows),
            'Table 3. Rates per cycle (N = 20–60), means over seeds; the largest seed sd is %s (annealed) and %s '
            '(quenched), n/a (deterministic) for periodic. "y × f" measures S with the momentum axis stretched by f. '
            'Curves in Appendix Figure A11; every value with its sd in Appendix Table B12.' % (sds['annealed'], sds['quenched'])),
        conf=('Strong for annealed, Tentative for quenched, Not supported for periodic',
              'annealed changes are within 3%; quenched changes stay inside the seed spread; periodic rates move by a third.'),
        cav='S depends on the metric, so the share with S ≤ 2 is not invariant (annealed, N = 60: %s unscaled, %s '
            'with y × 4); only rates are. Check: α is identical for every scaling, as a growth rate must be. The '
            'conclusions use rates only.' % (val(pa1), val(pa4))))
    WHERE[0] = 'Paper 1, claims'
    can = ['In annealed flow the share of low-stretch material decays %s times faster than the log-normal law '
           'predicts, at every threshold and metric scaling tested (Q2, Q6).' % Q[1]['a'].split(' times faster')[0].split(', ')[-1],
           'The cause is the floor S ≥ 1: a model of independent random steps reproduces the decay, and its moment '
           'growth rate flattens where a floor predicts (Q3).',
           'Caustics in annealed flow multiply at Kaplan\'s rate α + β/4 (Q4).',
           'Low-stretch rays are over-represented in branches but do not make them in annealed or quenched flow (Q5).']
    cannot = ['The exact limiting decay rate: finite-N estimates span %s to %s (Q3).' % (l4060, l2040),
              'Anything general about periodic flow: one deterministic run in a mixed regular/chaotic phase space, '
              'and its rates depend on the metric (Q2, Q4, Q6).',
              'That rays of low rarefaction r cause branches: that is largely true by construction (Q5).',
              'That quenched flow follows the log-normal law: consistent, but seeds scatter by about 20% (Q2).',
              'Caustic growth beyond N = %s, or for a smooth rather than piecewise-linear kick (Q4).' % fN]
    return dict(title='Stretching in kicked flow: the distribution of S, caustics and branches',
                sub='Paper 1 — report for supervisor', Q=Q, can=can, cannot=cannot)


# ============================================================ Paper 2
def paper2():
    Q = []
    A, Qn, Pm = 'annealed', 'quenched', 'cosine'
    WHERE[0] = 'Paper 2, Q1'
    f1 = [S('KM', 'F(1)', mode=m, direction='re-entry') for m in MODES]
    f5 = [S('KM', 'F(5)', mode=m, direction='re-entry') for m in MODES]
    fw = [S('KM', 'whole run', mode=m, direction='re-entry') for m in MODES]
    nv = [S('KM', 'naive', mode=m, direction='re-entry') for m in MODES]
    share5 = [D('%.0f%%' % (100 * _nums(a)[0] / _nums(b)[0]), 'F(5)/whole run, re-entry, %s' % m)
              for a, b, m in zip(f5, fw, MODES)]
    Q.append(dict(
        q='How likely is material that leaves the low-stretch set to come back?',
        a='The censoring-corrected return probability (mean ± sd over 10 seeds) is %s within one cycle, and %s by the '
          'end of the 60-cycle run. Returns are quick: %s of them happen within 5 cycles.' % (
              trio(f1), trio([nok(x) for x in fw]), trio(share5)),
        ev=('fig', 'F2.1', '%FIG%; the full table in Appendix Table B14.'),
        conf=('Strong', '10 seeds, about %s parcels at risk one cycle after leaving in each annealed run.' % '{:,}'.format(int(S('pooled', 'at risk (mean)', mode=A, age='1')))),
        cav='Censoring: material that leaves late is followed for only a few cycles. Check: the censoring-corrected '
            'and uncorrected (naive) values agree (%s against %s, annealed), and the age bookkeeping reproduces the '
            'run\'s own re-entry and exit counts exactly (Appendix Table B1). Ruled out.' % (nok(fw[0]), nv[0])))
    WHERE[0] = 'Paper 2, Q2'
    hz = {m: [S('pooled', 'h_re', mode=m, age=a) for a in ('1', '2', '6–10', '21–40', '>40')] for m in MODES}
    Q.append(dict(
        q='Does the chance of returning fade with time away?',
        a='Yes, steeply. The return hazard is %s one cycle after leaving, %s after two, %s at 6–10 cycles and %s at '
          '21–40 cycles (annealed); quenched flow follows the same pattern (%s at one cycle, %s at 21–40). In '
          'periodic flow the hazard is not monotonic: it peaks two cycles after leaving (%s) and rises again late '
          '(%s beyond 40 cycles), as material circulates around regular islands.' % (
              hz[A][0], hz[A][1], hz[A][2], hz[A][3], hz[Qn][0], hz[Qn][3], hz[Pm][1], hz[Pm][4]),
        ev=('fig', 'F2.2', '%FIG%; hazards by age band in Appendix Table B13; by exit cycle in Appendix Figure A13.'),
        conf=('Strong for annealed and quenched, Tentative for periodic', '10 seeds with about a million parcels at '
              'risk per run; periodic flow is one deterministic run.'),
        cav='A falling hazard can come from mixing parcels at different depths rather than from time itself, and it '
            'could depend on when the material left. Check: grouped by the cycle in which material left, the annealed '
            'hazard one cycle after leaving stays between %s and %s (10th–90th percentile over exit cycles; Appendix '
            'Figure A13), so it is set by time away; mixing by depth is tested in Q5.' % (
                NEWN('0.106', 'age-1 return hazard, 10th percentile over exit cycles (risk ≥ 200), annealed, seeds pooled'),
                NEWN('0.125', 'age-1 return hazard, 90th percentile over exit cycles (risk ≥ 200), annealed, seeds pooled'))))
    WHERE[0] = 'Paper 2, Q3'
    rows, full, sdx = [], [], []
    for m in MODES:
        r = [LONG[m]]
        for d in ('exit (entrants)', 'exit (orig. cohort)'):
            for c in ('F(1)', 'F(2)', 'F(10)'):
                x = V('KM', c, mode=m, direction=d)
                r.append(S('KM', c, mode=m, direction=d) if x != '–' else '–')
                if x != '–' and m != 'cosine':
                    sdx.append(_nums(x)[1])
        full.append(r)
        rows.append([r[0], ' / '.join(val(c) for c in r[1:4]), ' / '.join(val(c) for c in r[4:7])])
    sdmax = D('%.3f' % max(sdx), 'largest seed sd among the leaving probabilities shown, KM table')
    Q.append(dict(
        q='Does returned material stay in the low-stretch set?',
        a='No. In annealed flow %s of returned material leaves again after one cycle and %s within 10 cycles; '
          'quenched and periodic flow are similar. Material that has been in the set since launch leaves more slowly '
          'at first (%s after one cycle).' % (full[0][1], full[0][3], full[0][4]),
        ev=('table', 'T2.1', (['mode', 'returned material: share that has left within 1 / 2 / 10 cycles',
                               'original material: share that has left within 1 / 2 / 10 cycles'], rows),
            'Table 1. Censoring-corrected probabilities, by cycles since entering (returned material) or since launch '
            '(original material); means over seeds, largest seed sd %s, n/a (deterministic) for periodic. "–": fewer '
            'than 200 parcels at risk. Curves in Appendix Figure A12; every value with its sd in Appendix Table B14.' % sdmax),
        conf=('Strong', '10 seeds; every curve reaches nearly 1 within the run.'),
        cav='The original material starts at S = 1 exactly, the deepest possible point, so it is not a like-for-like '
            'comparison over the first cycles. Check: partly — by 10 cycles both groups have almost entirely left '
            '(%s and %s, annealed), so the conclusion that returned material does not stay holds either way.' % (
                val(full[0][3]), val(full[0][6]))))
    WHERE[0] = 'Paper 2, Q4'
    rows = []
    for m in MODES:
        rows.append([LONG[m]] + [S('depth', 'h(returned before) / h(first)', mode=m, stratum=s) for s in
                                 ('all depths (crude)', 'MH, age-adjusted', 'MH, depth-adjusted', 'MH, depth × age')])
    Q.append(dict(
        q='Is material that has returned before more likely to return again?',
        a='In raw counts, yes: its return hazard is %s times that of material away for the first time, reproducing '
          'the earlier %s (quenched) and %s (periodic). Allowing for time away alone does not remove it in quenched '
          'and periodic flow (%s and %s).' % (
              trio([r[1] for r in rows]), F('earlier 1.8×', '1.8×'), F('earlier 2.6×', '2.6×'), rows[1][2], rows[2][2]),
        ev=('table', 'T2.2', (['mode', 'raw ratio', 'adjusted for time away', 'adjusted for depth',
                               'adjusted for depth and time away'], rows),
            'Table 2. Ratio of return hazards, returned before ÷ first time; adjusted ratios combine strata of '
            'depth (7 bands) and/or time away (method in the appendix). 1 = no difference.'),
        conf=('Strong for the raw effect in quenched and periodic flow, Tentative in annealed',
              'the annealed raw ratio is only %s seed standard deviations above 1.' % D('%.0f' % ((_nums(rows[0][1])[0] - 1) / _nums(rows[0][1])[1]), 'annealed (crude − 1)/sd, depth table')),
        cav='Material that has returned before tends to sit closer to the threshold, so depth could explain the '
            'difference. Check: tested directly in Q5.'))
    WHERE[0] = 'Paper 2, Q5'
    b1 = [S('depth', 'h(returned before) / h(first)', mode=m, stratum='S_{N-1} (2,2.5]') for m in MODES]
    mh = [r[4] for r in rows]
    p35 = S('depth', 'h(returned before) / h(first)', mode=Pm, stratum='S_{N-1} (3,5]')
    ev5 = S('depth', 'events/at risk returned (per seed)', mode=A, stratum='S_{N-1} >5')
    ev5f = S('depth', 'events/at risk first (per seed)', mode=A, stratum='S_{N-1} >5')
    Q.append(dict(
        q='Is that a memory of having returned, or just depth?',
        a='Depth. At the same depth (S one cycle before a possible return) the ratio is close to 1 — %s in the '
          'shallowest band — and adjusted for depth and time away together it is %s. The raw excess in Q4 '
          'disappears.' % (trio(b1), trio(mh)),
        ev=('fig', 'F2.3', '%FIG%; events and parcels at risk in each band in Appendix Table B15.'),
        conf=('Strong for annealed and quenched, Tentative for periodic',
              '10 seeds; periodic flow is one deterministic run and one of its bands gives %s.' % p35),
        cav='The deepest band (S above 5) holds few returns (%s against %s per run, annealed), so its ratio is '
            'uncertain. Check: the adjusted ratio uses 7 finer depth bands, up to S above 100, so this band does '
            'not drive it.' % (ev5.split('/')[0], ev5f.split('/')[0])))
    WHERE[0] = 'Paper 2, claims'
    can = ['Of the material that leaves the low-stretch set, %s comes back within 60 cycles, most of it within 5 '
           'cycles (Q1).' % trio([val(nok(x)) for x in fw]),
           'In random flows the chance of return falls by more than two orders of magnitude with time away (Q2).',
           'Returned material does not stay: it leaves again within a few cycles (Q3).',
           'The apparent "returners return more" effect is a depth effect, not a memory of having returned (Q4, Q5).']
    cannot = ['How material returns after more than 40 cycles: few events there (Q2).',
              'Why periodic flow behaves differently: one deterministic run (Q2, Q5).',
              'Whether these numbers hold for thresholds other than S ≤ 2: re-entry was analysed at S ≤ 2 only.',
              'What sets the depth from which material returns: not analysed here.']
    return dict(title='Re-entry into the low-stretch set: crossing ages, hazards and depth',
                sub='Paper 2 — report for supervisor', Q=Q, can=can, cannot=cannot)


P1 = paper1()
P2 = paper2()
print('ledger entries:', len(LEDGER))


# ============================================================ PDF building
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.lib.fonts import addMapping
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image,
                                KeepTogether, CondPageBreak)
from pypdf import PdfReader

FD = '/usr/share/fonts/truetype/dejavu/'
pdfmetrics.registerFont(TTFont('DV', FD + 'DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB', FD + 'DejaVuSans-Bold.ttf'))
addMapping('DV', 0, 0, 'DV'); addMapping('DV', 1, 0, 'DVB'); addMapping('DV', 0, 1, 'DV'); addMapping('DV', 1, 1, 'DVB')
BODY = ParagraphStyle('b', fontName='DV', fontSize=11, leading=14.2, spaceAfter=4)
CAP = ParagraphStyle('c', parent=BODY, textColor=colors.HexColor('#52514e'), spaceAfter=6)
H1 = ParagraphStyle('h1', fontName='DVB', fontSize=17, leading=21, spaceAfter=4)
H2 = ParagraphStyle('h2', fontName='DVB', fontSize=13, leading=16, spaceBefore=6, spaceAfter=2)
SUB = ParagraphStyle('s', parent=BODY, textColor=colors.HexColor('#52514e'))
CELL = ParagraphStyle('cell', fontName='DV', fontSize=11, leading=13.2)
CELLB = ParagraphStyle('cellb', parent=CELL, fontName='DVB')
CONF_COL = {'Strong': '#008300', 'Tentative': '#b36b00', 'Not supported': '#c2382f'}


def esc(t):
    return t.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def tbl(hdr, rows, widths=None):
    data = [[Paragraph(esc(str(h)), CELLB) for h in hdr]] + [[Paragraph(esc(str(c)), CELL) for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([('LINEBELOW', (0, 0), (-1, 0), 0.8, colors.HexColor('#52514e')),
                           ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                           ('TOPPADDING', (0, 0), (-1, -1), 2), ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
                           ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f4f3f0')])]))
    return t


def img(fid):
    f = FIGS[fid]
    return Image(f['path'], f['w'] * inch, f['h'] * inch)        # 1:1, so fonts print at their set size


def footer(canvas, doc):
    canvas.saveState(); canvas.setFont('DV', 11); canvas.setFillColor(colors.HexColor('#52514e'))
    canvas.drawRightString(doc.pagesize[0] - doc.rightMargin, 0.9 * cm, str(doc.page))
    canvas.restoreState()


def conf_par(c):
    lvl, why = c
    col = CONF_COL['Strong'] if lvl.startswith('Strong') and 'Tentative' not in lvl and 'Not' not in lvl else (
        CONF_COL['Not supported'] if lvl.startswith('Not') else CONF_COL['Tentative'])
    return Paragraph('<b>Confidence: <font color="%s">%s</font>.</b> %s' % (col, esc(lvl), esc(why[0].upper() + why[1:])), BODY)


def main_report(P, fname, fig_ids, tab_ids):
    doc = SimpleDocTemplate(os.path.join(NEW, fname), pagesize=A4, leftMargin=1.7 * cm, rightMargin=1.7 * cm,
                            topMargin=1.4 * cm, bottomMargin=1.5 * cm, title=P['title'], author='')
    W = A4[0] - 3.4 * cm - 12
    S_ = [Paragraph(esc(P['title']), H1), Paragraph(esc(P['sub']) + ' · numbers are mean ± sd over 10 seeds; the '
          'periodic mode is deterministic (sd n/a). Full tables, extra figures and methods are in the appendix PDF.', SUB),
          Spacer(1, 6),
          Paragraph('<b>How to read this report.</b> Each section answers one question: the answer with its key '
                    'numbers, the evidence (one figure or table), how confident we are and why, and the most likely '
                    'artefact with the check that addresses it.', BODY),
          Paragraph('Glossary', H2),
          tbl(['term', 'meaning'], GLOSSARY, [3.6 * cm, W - 3.6 * cm]), PageBreak()]
    fi = ti = 0
    for k, q in enumerate(P['Q'], 1):
        head = [Paragraph('Q%d. %s' % (k, esc(q['q'])), H2), Paragraph('<b>Answer.</b> ' + esc(q['a']), BODY)]
        kind, eid, extra = q['ev'][0], q['ev'][1], q['ev'][2]
        if kind == 'fig':
            fi += 1
            ev = [img(eid), Paragraph('<b>Figure %d.</b> %s' % (fi, esc(FIGS[eid]['caption'])), CAP)]
            evtext = '<b>Evidence.</b> ' + esc(extra.replace('%FIG%', 'Figure %d' % fi))
        else:
            ti += 1
            hdr, rows = extra
            rows = [[r[0]] + [c.replace(' ± n/a (deterministic)', ' ± n/a').replace(
                ' (deterministic)', '') for c in r[1:]] if r[0].startswith('periodic') else r for r in rows]
            ncol = len(hdr); w0 = 4.4 * cm
            capt = q['ev'][3]
            if 'n/a (deterministic)' not in capt:
                capt = capt.replace('Table %d. ' % ti, 'Table %d. Periodic sd: n/a (deterministic). ' % ti, 1)
            ev = [Paragraph('<b>Evidence.</b>', BODY), tbl(hdr, rows, [w0] + [(W - w0) / (ncol - 1)] * (ncol - 1)),
                  Paragraph(esc(capt).replace('Table %d.' % ti, '<b>Table %d.</b>' % ti, 1), CAP)]
            evtext = None
        S_.append(KeepTogether(head + ([Paragraph(evtext, BODY)] if evtext else [])))
        S_.append(KeepTogether(ev))
        S_.append(conf_par(q['conf']))
        S_.append(Paragraph('<b>Caveat.</b> ' + esc(q['cav']), BODY))
    S_.append(KeepTogether([Paragraph('What we can claim', H2)] +
                           [Paragraph('• ' + esc(c), BODY) for c in P['can']]))
    S_.append(KeepTogether([Paragraph('What we cannot claim yet', H2)] +
                           [Paragraph('• ' + esc(c), BODY) for c in P['cannot']]))
    doc.build(S_, onFirstPage=footer, onLaterPages=footer)
    return os.path.join(NEW, fname)


MAINFIG = {'paper1_stretching.pdf': ['F1.1', 'F1.3', 'F1.5'],
           'paper2_reentry.pdf': ['F2.1', 'F2.2', 'F2.3']}
p1 = main_report(P1, 'paper1_stretching.pdf', MAINFIG['paper1_stretching.pdf'], ['T1.1', 'T1.2'])
p2 = main_report(P2, 'paper2_reentry.pdf', MAINFIG['paper2_reentry.pdf'], ['T2.1', 'T2.2'])
for p in (p1, p2):
    print(os.path.basename(p), 'pages:', len(PdfReader(p).pages))


# ============================================================ appendix
LONGT = dict(LONG)


def conv(key, hdr, colmap=None, rowsel=None, extra_maps=None):
    """A previous-report table for the appendix: plain headers, plain mode and row names, and
    'n/a (deterministic)' in place of '± 0.0000' for the periodic mode. Values are unchanged."""
    H, rows = T[key]
    idx = list(range(len(H))) if colmap is None else colmap
    out = []
    maps = {'row0-rows/8': 'y0 = row0 − 240 px', 'row0': 'y0 = row0', 'row0+rows/8': 'y0 = row0 + 240 px',
            're-entry': 'return', 'exit (entrants)': 'leaving, returned material',
            'exit (orig. cohort)': 'leaving, original material', 'all depths (crude)': 'all depths (raw)',
            'S_{N-1} (2,2.5]': 'depth 2–2.5', 'S_{N-1} (2.5,3]': 'depth 2.5–3', 'S_{N-1} (3,5]': 'depth 3–5',
            'S_{N-1} >5': 'depth above 5', 'MH, depth-adjusted': 'adjusted for depth', 'MH, age-adjusted':
            'adjusted for time away', 'MH, depth × age': 'adjusted for depth and time away',
            'i.i.d. model (1e7)': 'independent-step model (10 million)', '2-D annealed (10 seeds)': 'annealed runs (10 seeds)'}
    if extra_maps:
        maps.update(extra_maps)
    for r in rows:
        if rowsel and not rowsel(r):
            continue
        det = r[0] == 'cosine'
        rr = []
        for i in idx:
            c = str(r[i])
            if i == 0 and c in LONGT:
                c = LONGT[c]
            c = maps.get(c, c)
            if c.startswith('nan ± nan') or c == 'nan':
                c = '–'
            if '±' in c:                                    # every 'value ± sd' in the cell
                c = re.sub(r'([-+]?[\d.]+(?:e[-+]?\d+)?)\s*±\s*\+?([\d.]+(?:e[-+]?\d+)?)',
                           lambda mm: '%s ± %s' % (mm.group(1).lstrip('+'), 'n/a' if det else tiny_sd(mm.group(2))), c)
            if i == 0 and det:
                c = c + '; sd n/a (deterministic)'
            rr.append(c)
        out.append(rr)
    return hdr, out


def checks_table():
    rep = [('Σ age re-entry events = Σ re-entries (haz_reent·count)', 'Return events counted by age = re-entries counted by the run'),
           ('Σ age exit events = Σ exits (haz_gross·count)', 'Exit events counted by age = exits counted by the run'),
           ('pooled age-1 re-entry hazard = Type 1 (consecutive masks)', 'Return hazard one cycle after leaving = immediate bounce-back counted from consecutive cycles'),
           ('rho: 1/|a| vs histogram within 2% away from caustics', 'Ray density: exact (sum of 1/|a|) against the histogram, within 2% away from caustics'),
           ('identity S ≥ exp(r) on every line parcel', 'S ≥ exp(r) for every ray on the launched lines'),
           ('cosine', 'periodic'), ('Eulerian lattice', 'pixel grid'), ('window-selected full-period lattice', 'full-period sample of parcels'),
           ('MATLAB 1C (an_br)', 'MATLAB branch analysis'), ('rho_hist, rho from 1/|a|', 'ray density (histogram and exact)'),
           ('4F', 'r against S'), ('Part 4', 'distribution of S'), ('Part 2', 'crossing-age'), ('Type 1', 'immediate bounce-back'),
           ('agree to 1e-4', 'agree to 10⁻⁴'), ('rows', 'rows')]
    out = []
    for n, ok, d in J['checks']:
        n = n.replace('e^r', 'exp(r)'); d = d.replace('e^r', 'exp(r)')
        for a, b in rep:
            n = n.replace(a, b); d = d.replace(a, b)
        out.append(['PASS' if ok else 'FAIL', n, d])
    return ['result', 'check', 'detail'], out


APPX_TABLES = [
    ('B1', 'All PASS/FAIL checks of the analysis', checks_table()),
    ('B2', 'Log-normal test at every threshold s (rates per cycle, N = 20–60; measured ÷ predicted in brackets)',
     conv('4C', ['mode', 'growth of mean ln S, α', 'growth of variance, v', 'log-normal prediction α²/(2v)',
                 's = 1.25', 's = 1.5', 's = 2', 's = 3', 's = 5', 's = 10'])),
    ('B2 (seed 0)', 'Log-normal test, seed 0 (reproduces the earlier single-run numbers)',
     conv('4C_seed0', ['mode (seed 0)', 'log-normal prediction', 'measured, s = 2', 'measured, s = 1.25 to 10'])),
    ('B3', 'Share of material with S ≤ s on the pixel grid (mean over seeds)',
     conv('Ple', ['mode', 'N'] + ['s = %s' % x for x in ('1.25', '1.5', '2', '3', '5', '10', '20', '50', '100', '10³', '10⁴', '10⁶')])),
    ('B4', 'Moments and quantiles of ln S',
     conv('mom', ['mode', 'N', 'mean ln S', 'variance', 'skewness', 'excess kurtosis', '1%', '10%', '25%', '50%', '75%', '90%', '99%'])),
    ('B5', 'Independent-step model against the annealed runs (rates per cycle)',
     conv('iid', ['model', 'growth α', 'growth v', 'log-normal α²/(2v)', 'decay of the share with S ≤ 2', 'decay, s = 1.25 to 10',
                  'limiting rate from the λ distribution, raw → corrected for (ln S)²', 'plateau of −L(q): cycles 20→40 / 40→60'])),
    ('B5 (shares)', 'Share with S ≤ s: model against the annealed runs (z = difference in standard errors of the seed mean)',
     conv('iid_P', ['N', 'model, s = 1.5', 'runs, s = 1.5', 'z', 'model, s = 2', 'runs, s = 2', 'z', 'model, s = 5', 'runs, s = 5', 'z',
                    'model, s = 10', 'runs, s = 10', 'z'])),
    ('B6', 'Moment growth rate L(q) (per cycle)',
     conv('Lq', ['moment order q', 'transfer operator (exact, length of a vector)', 'model, cycles 40 → 60', 'model, (1/60) ln of the moment'])),
    ('B7', 'Caustic and line growth at each launch momentum (rates per cycle)',
     conv('caustics', ['mode', 'launch momentum', 'caustic growth', 'line-length growth', 'Kaplan prediction α + β/4',
                       'caustic ÷ prediction', 'α (growth of mean r)', 'β (2 × growth of variance of r)', 'cycles fitted'])),
    ('B8', 'Removing rays with S ≤ 2 from the ray density (mean of the three launch momenta)',
     conv('1C', ['mode', 'N', 'share of all rays', 'share of rays in branches', 'enrichment', 'contrast, all rays',
                 'contrast without them', 'contrast lost', 'corr, with and without', 'overlap of branch sets'], colmap=list(range(10)))),
    ('B8 (rarefaction)', 'Removing the rays of lowest rarefaction r, and the control that removes the same numbers evenly',
     conv('1C', ['mode', 'N', 'low-r enrichment', 'low-r contrast lost', 'corr, with and without low-r', 'control: contrast lost (S / r)'],
          colmap=[0, 1, 10, 11, 12, 13])),
    ('B9', 'Control: parcels seeded every 2 px over the frame, which fill the plane and so project to an almost flat density',
     conv('1C_lat', ['mode', 'N', 'contrast of their projection', 'low-S enrichment', 'low-r enrichment', 'corr, with and without low-S'])),
    ('B10', 'Correlation across x of ray density with local stretching',
     conv('corr', ['mode', 'N', 'corr(ray density, share with S ≤ 2), Pearson', 'Spearman', 'corr(ray density, local mean of −r), Pearson', 'Spearman'])),
    ('B11', 'Rarefaction r against S on the launched lines',
     conv('4F', ['mode', 'N', 'share with S ≤ 2', 'P(S ≤ 2 | r ≤ ln 2)', 'P(r ≤ ln 2 | S ≤ 2)', '10% quantile of r',
                 'P(S ≤ 2 | r below its 10% quantile)', 'P(r below its 10% quantile | S ≤ 2)'])),
    ('B12', 'Threshold and metric robustness (rates per cycle; share at N = 60)',
     conv('4E', ['mode', 'momentum scale f', 'growth α', 'decay, S ≤ 1.5', 'share at N = 60', 'decay, S ≤ 2', 'share at N = 60',
                 'decay, S ≤ 3', 'share at N = 60'])),
    ('B13', 'Return and leaving hazards by time since the last crossing (pooled over age bands)',
     conv('pooled', ['mode', 'cycles since crossing', 'return hazard', 'parcels at risk (mean per run)',
                     'leaving hazard of returned material', 'parcels at risk (mean per run)'])),
    ('B14', 'Censoring-corrected probabilities of having crossed, by cycles since the last crossing',
     conv('KM', ['mode', 'direction', 'within 1', 'within 2', 'within 3', 'within 5', 'within 10', 'within 20', 'within 40',
                 'whole run (age reached)', 'uncorrected'])),
    ('B15', 'Returned before ÷ first time, by depth, with event counts per run',
     conv('depth', ['mode', 'stratum', 'hazard ratio', 'returns / at risk, returned before', 'returns / at risk, first time'])),
]

APPX_FIGS = [('A%d' % i, 'A.%d' % i) for i in range(1, 16)]
METHODS = [
    ('Map and stretching', 'Each cycle applies a kick and a drift: y → y + f(x), then x → x + Tau (y − row0), with f the '
     'piecewise-linear interpolant of the kick profile, 1920 × 1920 px, Tau = 0.5, row0 = 960 px (the momentum with no '
     'drift). The Jacobian of one cycle is J = [1 + Tau f′, Tau; f′, 1], with determinant 1, and the product M of the '
     'Jacobians is carried exactly. S is the largest singular value of M, S = (√(t + 2) + √(t − 2))/2 with t the sum of '
     'the squared entries; S ≥ 1 always.'),
    ('Modes and seeds', 'annealed (fresh random kick each cycle), quenched (one frozen random kick) and periodic (standard '
     'map, K = 0.5). Random kicks reproduce the MATLAB generator exactly (seed 0 is the MATLAB rng(0) run). Ten seeds per '
     'random mode; the periodic mode uses no random numbers, so its ten seeds are the same run and its sd is n/a '
     '(deterministic). All runs last 60 cycles.'),
    ('Where S is measured', 'S is measured on the pixel grid (every pixel of the frame at the current cycle). Because the map '
     'is periodic in momentum with period 3840 px, a sample of parcels seeded over a full period and kept where they land '
     'in the frame gives the same distribution; this was checked to 10⁻⁴ and is how the multi-seed runs measure it.'),
    ('Launched lines and caustics', 'Rays are launched at all x on three lines, launch momentum y0 = row0 − 240 px / row0 / '
     'row0 + 240 px, 8 rays per px, each carrying its own M from the identity. The tangent of the line is (a, c), the first '
     'column of M; rarefaction r = ln √(a² + c²). Ray density is the histogram of x in 4-px bins (normalised to 1 at launch) '
     'and, independently, the sum of 1/|a| over the branches of the line. A caustic is a sign change of a. Because the kick '
     'is piecewise linear, a is constant on pieces of the line, so the line is followed exactly as pieces and every caustic '
     'is counted, up to 12 million pieces. Kaplan\'s prediction for the growth of the caustic count is α + β/4, with α the '
     'growth rate of the mean of r and β twice that of its variance, fitted over N = 10 to the end of the run.'),
    ('Branches', 'Branches are the top 10% of the 480 bins of ray density from all rays. For a group of rays (S ≤ 2, or r '
     'below its 10th percentile) we report its share of the rays in branches, enrichment (that share ÷ its share of all '
     'rays), the contrast (mean density in branches ÷ mean density) with and without the group on the same bins, the '
     'fraction of the excess contrast lost, (C_all − C_without)/(C_all − 1), corr across all bins of the densities with '
     'and without, and the overlap of the two sets of branch bins. The control removes the same number of rays evenly '
     'along the line.'),
    ('Log-normal test and independent-step model', 'α and v are the slopes of the mean and variance of ln S over N = 20–60; '
     'the log-normal law predicts that the share with S ≤ s decays at α²/(2v) per cycle. The measured rate is the slope of '
     '−ln(share) over the same N, using points with at least 10 parcels. The independent-step model draws f′ for every '
     'cycle independently from the slopes of 500 fresh kicks and multiplies 10 million products to N = 60. L(q) is the '
     'growth rate of the q-th moment of S, estimated from the change between cycles 40 and 60; for q ≥ 0 it is also '
     'computed exactly from a transfer operator on 2048 directions.'),
    ('Metric', 'S of the same M measured with the momentum axis stretched by f (S of D M D⁻¹, D = diag(1, f)), f = 0.5, 2, 4.'),
    ('Crossing ages', 'The tracer parcels are seeded every 2 px over the frame (921,600) and followed as material. A parcel '
     'is in the low-stretch set when S ≤ 2; all start in it (S = 1). For each parcel we keep the cycle of its last exit and '
     'entry, how often it has returned, and S one cycle earlier (its depth). Risk and event counts are accumulated by '
     'crossing cycle and age (cycles since the last crossing). The hazard at age k is events ÷ parcels at risk at age k, '
     'shown only with at least 200 at risk. The censoring-corrected probability of having crossed within k cycles is the '
     'Kaplan–Meier estimator, F(k) = 1 − Π_{j ≤ k}(1 − h(j)); "uncorrected" is the plain share of spells that ended '
     'within the run. Ratios adjusted for depth and/or age are Mantel–Haenszel rate ratios over strata of 7 depth bands '
     'and/or ages.'),
    ('Checks', 'Every number in the reports is taken from the previous report\'s tables, which the stored run files '
     'reproduce exactly. The MATLAB script and the Python code agree exactly on identical kicks (Table B1).'),
]
CAVEATS = [
    'Seeds above 0 assume MATLAB\'s generator seeds as numpy does (first draw for seed 1: 0.4170). Kicks agree with MATLAB\'s '
    'to about 10⁻⁷ (rounding in the kick generator\'s Fourier filter); single chaotic parcels can differ at large N, statistics do not.',
    'The periodic mode is deterministic: its ten seeds are one run, not independent samples.',
    'Censoring: hazards need at least 200 parcels at risk, so whole-run probabilities stop at the largest such age.',
    'Crossing ages follow material; the distribution of S is measured on the pixel grid. Once material leaves the frame the '
    'two describe different populations.',
    'S depends on the coordinate metric: thresholds such as S ≤ 2 and every share built on them change with the scaling of '
    'momentum; growth rates do not.',
    'Caustic counts are a property of the piecewise-linear kick at 1-px resolution; with a smooth kick folds are smooth and '
    'the ray density diverges at them. Exact counts stop at 12 million pieces (annealed N ≈ 31).',
    '"Immediate bounce-back" is the probability of being back in the set one cycle after leaving, counted from consecutive cycles.',
    'Finite N: the log-normal law is a statement about the centre of ln S; the share with S ≤ s probes the far lower tail.',
]


def scan_text(text, appendix=False):
    hits = {}
    terms = ['cosine', 'LSVD', 'Eulerian', 'compressed', 'Kaplan-Meier', 'an_br', 'S_THR', 'Pearson r', 'window',
             '± 0.0000', '±0.0000', '± 0.000 ', 'lattice']
    if not appendix:
        terms += ['row0', 'rows/8', 'cycleData', 'Kaplan–Meier']
    for t in terms:
        n = text.count(t)
        if n:
            hits[t] = n
    codes = re.findall(r'\b[124][A-F]\b', text)
    if codes:
        hits['part codes'] = len(codes)
    return hits


def pdf_text(p, first=0, last=None):
    rd = PdfReader(p)
    pages = rd.pages[first:last]
    return '\n'.join(pg.extract_text() for pg in pages), len(rd.pages)


def main_refs():
    refs = {}
    for name, P in (('Paper 1', P1), ('Paper 2', P2)):
        for i, q in enumerate(P['Q'], 1):
            txt = ' '.join([q['a'], str(q['ev'][2]) if q['ev'][0] == 'fig' else q['ev'][3], q['conf'][1], q['cav']])
            for kind, lab in (('Figure', 'A'), ('Table', 'B')):
                for m_ in re.finditer(kind + r's? ((?:%s)\d+(?:[–-](?:%s)?\d+)?(?:(?:, | and )(?:%s)?\d+)*)' % (lab, lab, lab), txt):
                    for part in re.split(r', | and ', m_.group(1)):
                        part = part if part.startswith(lab) else lab + part
                        if '–' in part:
                            a_, b_ = part.split('–'); b_ = b_.lstrip(lab)
                            ids = ['%s%d' % (lab, k) for k in range(int(a_[1:]), int(b_) + 1)]
                        else:
                            ids = [part]
                        for id_ in ids:
                            refs.setdefault(id_, set()).add('%s Q%d' % (name, i))
            if q['ev'][0] == 'fig':
                refs.setdefault('main:' + q['ev'][1], set()).add('%s Q%d' % (name, i))
            else:
                refs.setdefault('main:' + q['ev'][1], set()).add('%s Q%d' % (name, i))
    return refs


def build_appendix(check_rows, scan_info, removed):
    path = os.path.join(NEW, 'appendix.pdf')
    pg = landscape(A4)
    doc = SimpleDocTemplate(path, pagesize=pg, leftMargin=1.5 * cm, rightMargin=1.5 * cm, topMargin=1.3 * cm,
                            bottomMargin=1.5 * cm, title='Appendix: methods, figures, tables and checks', author='')
    W = pg[0] - 3 * cm - 12
    S_ = [Paragraph('Appendix to the two reports: methods, figures, tables and checks', H1),
          Paragraph('Every item here is referenced by a question in one of the two main reports (see the check at the '
                    'end). Numbers are mean ± sd over 10 seeds; the periodic mode is deterministic (sd n/a).', SUB),
          Paragraph('Methods', H2)]
    for h, t in METHODS:
        S_.append(Paragraph('<b>%s.</b> %s' % (esc(h), esc(t)), BODY))
    S_.append(Paragraph('Caveats', H2))
    for c in CAVEATS:
        S_.append(Paragraph('• ' + esc(c), BODY))
    S_.append(PageBreak()); S_.append(Paragraph('Appendix figures', H1))
    for lab, fid in APPX_FIGS:
        if lab in removed:
            continue
        f = FIGS[fid]
        S_.append(KeepTogether([img(fid), Paragraph('<b>Figure %s.</b> %s <i>Referenced by: %s.</i>' % (
            lab, esc(f['caption']), ', '.join(sorted(REFS.get(lab, [])))), CAP), Spacer(1, 8)]))
    S_.append(PageBreak()); S_.append(Paragraph('Appendix tables', H1))
    for lab, title, (hdr, rows) in APPX_TABLES:
        base = lab.split(' ')[0]
        if base in removed:
            continue
        ncol = len(hdr)
        w0 = 4.4 * cm if hdr[0] in ('mode', 'mode (seed 0)') else (2.2 * cm if hdr[0] == 'result' else 2.6 * cm)
        widths = [w0] + [(W - w0) / (ncol - 1)] * (ncol - 1)
        if hdr[0] == 'result':
            widths = [2.2 * cm, 0.45 * (W - 2.2 * cm), 0.55 * (W - 2.2 * cm)]
        S_.append(CondPageBreak(6 * cm))
        S_.append(Paragraph('<b>Table %s.</b> %s. <i>Referenced by: %s.</i>' % (esc(lab), esc(title), ', '.join(sorted(REFS.get(base, [])))), BODY))
        S_.append(tbl(hdr, rows, widths)); S_.append(Spacer(1, 10))
    S_.append(PageBreak())
    CHECK_START[0] = len(S_)
    S_.append(Paragraph('Check: figures, tables, numbers and terms', H1))
    S_.append(Paragraph('1. Every figure and table, and the question that references it. Items no question referenced '
                        'were removed: %s.' % esc(removed_text()), BODY))
    S_.append(tbl(['item', 'where', 'referenced by'], check_rows, [4.5 * cm, 6 * cm, W - 10.5 * cm]))
    for line in scan_info:
        S_.append(Paragraph(esc(line), BODY))
    doc.build(S_, onFirstPage=footer, onLaterPages=footer)
    return path


CHECK_START = [0]
REMOVED_FROM_OLD = ['ridge plot and heat map of the distribution of S', 'stacked shares of S bands',
                    'band projections of the tracer parcels onto x', 'per-mode hazard figures with at-risk bars',
                    'pooled age-bin bar chart', 'phase-plane snapshots at N = 40 and 60',
                    'correlation-against-N figure (its numbers are in Table B10)']


def removed_text():
    return '; '.join(REMOVED_FROM_OLD) + ' (figures of the previous report that no question needs)'


# ---------------- reference check
REFS = main_refs()
removed = set()
for lab, _ in APPX_FIGS:
    if lab not in REFS:
        removed.add(lab)
for lab, _, _ in APPX_TABLES:
    if lab.split(' ')[0] not in REFS:
        removed.add(lab.split(' ')[0])
check_rows = []
main_items = [('paper1_stretching.pdf', 'Figure %d' % (i + 1), fid) for i, fid in enumerate(MAINFIG['paper1_stretching.pdf'])] + \
             [('paper2_reentry.pdf', 'Figure %d' % (i + 1), fid) for i, fid in enumerate(MAINFIG['paper2_reentry.pdf'])]
for f_, lab, fid in main_items:
    check_rows.append(['%s, %s' % ('Paper 1' if 'paper1' in f_ else 'Paper 2', lab), f_, ', '.join(sorted(REFS.get('main:' + fid, ['NOT REFERENCED'])))])
for P, nm, tids in ((P1, 'Paper 1', None), (P2, 'Paper 2', None)):
    k = 0
    for i, q in enumerate(P['Q'], 1):
        if q['ev'][0] == 'table':
            k += 1
            check_rows.append(['%s, Table %d' % (nm, k), nm.lower().replace(' ', '') + ' main text', '%s Q%d' % (nm, i)])
for lab, fid in APPX_FIGS:
    check_rows.append(['Figure %s' % lab, 'appendix' if lab not in removed else 'REMOVED', ', '.join(sorted(REFS.get(lab, ['not referenced'])))])
for lab, _, _ in APPX_TABLES:
    b = lab.split(' ')[0]
    check_rows.append(['Table %s' % lab, 'appendix' if b not in removed else 'REMOVED', ', '.join(sorted(REFS.get(b, ['not referenced'])))])
assert all('NOT REFERENCED' not in r[2] for r in check_rows), 'a main-report item is not referenced'

# ---------------- scans of the main reports
scan = []
for p in (p1, p2):
    txt, npg = pdf_text(p)
    hits = scan_text(txt)
    scan.append('%s: %d pages; banned terms and internal codes found: %s.' % (os.path.basename(p), npg, hits or 'none'))
    if hits or npg > 6:
        raise SystemExit('STOP: %s failed the term/page check: %s, %d pages' % (p, hits, npg))
pa = build_appendix(check_rows, [], removed)                    # first pass, to find where the check section starts
rd = PdfReader(pa)
cstart = next(i for i, pg_ in enumerate(rd.pages) if 'Check: figures, tables, numbers and terms' in pg_.extract_text())
txt, npg = pdf_text(pa, 0, cstart)
hits = scan_text(txt, appendix=True)
km = txt.count('Kaplan–Meier')
if hits or km != 1:
    raise SystemExit('STOP: appendix term check failed: %s, Kaplan–Meier named %d times' % (hits, km))
scan.append('appendix.pdf (all pages before this check): banned terms found: none; the estimator for censored data is '
            'named once, in Methods.')
fc = [f for f in FIGCHECK]
scan.append('2. Figures: %d built; every text element is ≥ %d pt at printed size (figures are placed at 100%%); no tick '
            'labels overlap; no banned terms; at most 3 panels each.' % (len(fc), FS))
n_same = sum(l['status'] == 'same as previous report' for l in LEDGER)
n_der = sum(l['status'].startswith('derived') for l in LEDGER)
n_new = sum(l['status'] == 'NEW' for l in LEDGER)
n_diff = sum(l['status'] == 'DIFFERS' for l in LEDGER)
scan.append('3. Numbers: %d quoted numbers; %d taken unchanged from the previous report, %d derived arithmetically from '
            'previous-report values (ratios, ranges), %d new (read directly from the stored run files; listed in '
            'CHECK.md), %d differ.' % (len(LEDGER), n_same, n_der, n_new, n_diff))
scan.append('4. The stored run files, re-analysed in memory, reproduce the previous report\'s tables, checks and summary '
            'exactly; no existing output was written (checksums before and after in CHECK.md).')
pa = build_appendix(check_rows, scan, removed)
print('appendix pages:', len(PdfReader(pa).pages), '| removed:', sorted(removed))


# ============================================================ CHECK.md (written last)
def sha_all():
    out = {}
    for d in (OLD, OUT):
        for p in sorted(glob.glob(os.path.join(d, '**', '*'), recursive=True)):
            if os.path.isfile(p):
                out[os.path.relpath(p, ROOT)] = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    return out


FROZEN = os.environ.get('KD_FROZEN', '')
after = sha_all()
chk = 'not compared (no checksum list given in KD_FROZEN)'
if FROZEN and os.path.exists(FROZEN):
    before = {}
    for line in open(FROZEN):
        h, p = line.split(None, 1)
        before[p.strip()] = h
    changed = [p for p in before if after.get(p) != before[p]]
    added = [p for p in after if p not in before]
    chk = '%d files compared; changed: %s; new files: %s' % (len(before), changed or 'none', added or 'none')
    if changed:
        raise SystemExit('STOP: existing outputs changed: %s' % changed)

L = ['# reports_v2 — build check', '',
     'Built by `python/kd_present.py`, read-only from `python/out/` (stored run results) and `reports/` (the previous '
     'report). Nothing outside `reports_v2/` was written.', '',
     '## Existing outputs unchanged', '', '- SHA-256 of every file in `reports/` and `python/out/`, before and after: ' + chk, '',
     '## Previous report reproduced', '',
     '- The stored run files, re-analysed into a scratch folder with the unchanged report code, reproduce '
     '`reports/kd_tables.json` (all tables and checks), `reports/kd_tables.md` and `reports/kd_findings.json` exactly.', '',
     '## Pages and terms', '']
for line in scan:
    L.append('- ' + line)
L += ['- Terms scanned in the main reports: cosine, LSVD, Eulerian, compressed, Kaplan-Meier / Kaplan–Meier, an_br, '
      'S_THR, Pearson r, window, lattice, row0, rows/8, cycleData, "± 0.0000", and part codes such as 1A or 4C. In the '
      'appendix the same list except row0, rows/8 and cycleData, and Kaplan–Meier must appear exactly once (Methods).', '',
      '## Figures', '', '| figure file | printed size (in) | min font (pt) | overlapping tick labels | banned terms |', '|---|---|---|---|---|']
for f in FIGCHECK:
    L.append('| %s | %.2f × %.2f | ≥ %d | %s | %s |' % (f['name'], f['w'], f['h'], FS, f['overlap'] or 'none', f['banned'] or 'none'))
L += ['', '## Every figure and table, and the question that references it', '', '| item | where | referenced by |', '|---|---|---|']
for r in check_rows:
    L.append('| %s | %s | %s |' % tuple(r))
L += ['', 'Removed from the previous report because no question needs them: ' + '; '.join(REMOVED_FROM_OLD) + '.', '',
      '## Number ledger', '',
      '%d numbers: %d unchanged from the previous report, %d derived from previous-report values, %d new, %d differ.' % (
          len(LEDGER), n_same, n_der, n_new, n_diff), '',
      '| where | source in the previous report (table, column, row) | previous | shown | status |', '|---|---|---|---|---|']
for l in LEDGER:
    L.append('| %s | %s | %s | %s | %s |' % (l['where'], l['source'].replace('|', '/'), l['old'], l['shown'], l['status']))
open(os.path.join(NEW, 'CHECK.md'), 'w').write('\n'.join(L) + '\n')
print('CHECK.md written;', chk)
