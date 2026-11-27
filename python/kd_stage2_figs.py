"""kd_stage2_figs -- figures for the stage-2 questions (at most 4 panels each)."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm, Normalize
import kd_engine as E
import kd_stage2 as S2
import kd_stage2_analysis as A

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 11.5, 'figure.titlesize': 14})
MODE_COL = {'annealed': '#1f77b4', 'quenched': '#ff7f0e', 'cosine': '#2ca02c'}
MODE_NAME = {'annealed': 'annealed (new random kick every cycle)', 'quenched': 'quenched (one random kick, repeated)',
             'cosine': 'periodic (sine kick)'}
FATE_COL = ['#2ca02c', '#ff7f0e', '#7f7f7f']
FATE_LAB = ['never above S = 2', 'left S <= 2 and came back', 'left S <= 2 for good']
_t = np.linspace(0, 1, 256)[:, None]
CMAP_S = ListedColormap(np.hstack([0.0 * (1 - _t) ** 2 + 0.129 * 2 * _t * (1 - _t) + 0.993 * _t ** 2,
                                   0.25 * (1 - _t) ** 2 + 0.566 * 2 * _t * (1 - _t) + 0.906 * _t ** 2,
                                   0.35 * (1 - _t) ** 2 + 0.692 * 2 * _t * (1 - _t) + 0.144 * _t ** 2]))


def _save(fig, folder, name, dpi=110):
    os.makedirs(folder, exist_ok=True)
    p = os.path.join(folder, name)
    fig.savefig(p, dpi=dpi)
    plt.close(fig)
    return p


def _seedband(ax, x, Y, col, label=None, lw=1.6):
    """mean over seeds as a line, +- sd as a band (single seed: just the line)"""
    Y = np.atleast_2d(Y)
    m = np.nanmean(Y, 0)
    ax.plot(x, m, color=col, lw=lw, label=label)
    if Y.shape[0] > 1:
        s = np.nanstd(Y, 0, ddof=1)
        ax.fill_between(x, m - s, m + s, color=col, alpha=0.2, lw=0)


# ------------------------------------------------------------------ the starting picture (core.py)
def start_image(R=S2.R):
    u = R / 41.0
    row_of = lambda yd: (41 - yd) * u
    img = np.ones((R, R, 3), np.float32)
    rgb = {'red': (1, 0, 0), 'green': (0, 1, 0), 'blue': (0, 0, 1)}
    for k, (a, b) in S2.BAND_EDGES.items():
        img[int(round(a)):int(round(b))] = rgb[k]
    yy, xx = np.mgrid[0:R, 0:R].astype(float)
    for i in range(42):
        for j in range(42):
            cx = i * u; cy = row_of(j)
            x0, x1 = int(max(cx - 6, 0)), int(min(cx + 7, R)); y0, y1 = int(max(cy - 6, 0)), int(min(cy + 7, R))
            if x1 <= x0 or y1 <= y0:
                continue
            m = (xx[y0:y1, x0:x1] - cx) ** 2 + (yy[y0:y1, x0:x1] - cy) ** 2 <= 4.4 ** 2
            img[y0:y1, x0:x1][m] = 0
    return img


_IMG0 = None


def picture_and_M(kicks, N, Tau, row0, t_last=1.0, window=None):
    """The exact picture at time (N-1) + t_last: every pixel traced back to its
    starting point (nearest pixel of the starting picture, as core.py), and the
    monodromy M of that material. window = (x0, x1, y0, y1) restricts to a crop."""
    global _IMG0
    if _IMG0 is None:
        _IMG0 = start_image()
    R = S2.R
    x0, x1, y0, y1 = window if window else (0, R, 0, R)
    ys, xs = np.mgrid[y0:y1, x0:x1].astype(float)
    x = xs.ravel(); y = ys.ravel()
    # backwards: the partial drift of cycle N, its kick, then whole cycles
    x = x - t_last * Tau * (y - row0)
    dp = kicks[N - 1]; xf = np.floor(x); i0 = np.mod(xf.astype(np.int64), R)
    y = y - (dp[i0] + (dp[(i0 + 1) % R] - dp[i0]) * (x - xf))
    for k in range(N - 2, -1, -1):
        dp = kicks[k]
        x = x - Tau * (y - row0); xf = np.floor(x); i0 = np.mod(xf.astype(np.int64), R)
        y = y - (dp[i0] + (dp[(i0 + 1) % R] - dp[i0]) * (x - xf))
    pic = _IMG0[np.mod(np.rint(y).astype(int), R), np.mod(np.rint(x).astype(int), R)].reshape(ys.shape + (3,))
    st = E.new_state(x, y)
    for k in range(N - 1):
        E.step(st, kicks[k], Tau, row0)
    # kick of cycle N, then the partial drift
    dp = kicks[N - 1]; xx = st['x']; xf = np.floor(xx); i0 = np.mod(xf.astype(np.int64), R)
    fp = dp[(i0 + 1) % R] - dp[i0]
    st['y'] = st['y'] + dp[i0] + fp * (xx - xf)
    st['c'] = st['c'] + fp * st['a']; st['d'] = st['d'] + fp * st['b']
    st['x'] = xx + t_last * Tau * (st['y'] - row0)
    st['a'] = st['a'] + t_last * Tau * st['c']; st['b'] = st['b'] + t_last * Tau * st['d']
    S = S2.sigma(st['a'], st['b'], st['c'], st['d']).reshape(ys.shape)
    th, _ = S2.orient(st['a'], st['b'], st['c'], st['d'])
    return pic, S, th.reshape(ys.shape)


def overlay_S(pic, S, alpha=0.28, smax=15):
    tt = np.clip(np.log(np.maximum(S, 1)) / np.log(smax), 0, 1)
    return (1 - alpha) * pic + alpha * CMAP_S(tt)[..., :3]


def overlay_theta(pic, S, th, amax=0.75, smax=15):
    col = plt.get_cmap('twilight')((th + 90.0) / 180.0)[..., :3]
    al = amax * np.clip(np.log(np.maximum(S, 1)) / np.log(smax), 0, 1)[..., None]
    return (1 - al) * pic + al * col


# ================================================================== Part 1A
def fig1a(res, mode, folder):
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    t, p, w = res['scatter_true'], res['scatter_prod'], res['scatter_wrong']
    ax = axs[0]
    ax.loglog(p, t, ',', color='k', alpha=0.4)
    lim = [1, max(p.max(), t.max()) * 1.2]
    ax.plot(lim, lim, 'r-', lw=1.5, label='y = x: the upper bound S15 <= S10 x S(11-15)')
    ax.set_xlabel('S10 x S(11-15)   (stretch to cycle 10) x (stretch of cycles 11-15 alone)')
    ax.set_ylabel('true S15 (one monodromy through all 15 cycles)')
    ax.set_title('Keeping M: the product only bounds the truth from above')
    ax.legend(loc='upper left', fontsize=9)
    ax = axs[1]
    ax.loglog(w, t, ',', color='k', alpha=0.4)
    ax.axhline(2, color='b', lw=1); ax.axvline(2, color='b', lw=1)
    ax.plot(lim, lim, 'r-', lw=1)
    ax.set_xlabel('S(11-15): M reset to the identity at cycle 10 (the WRONG way)')
    ax.set_ylabel('true S15')
    ax.set_title('Forgetting M at cycle 10: %d of %d parcels (%.1f%%)\nchange low (S <= 2) / not-low class'
                 % (res['n_changed'], res['n_parcels'], 100 * res['n_changed'] / res['n_parcels']))
    fig.suptitle('Does the past matter? %s: true stretch vs restarting the count at cycle 10' % MODE_NAME[mode])
    fig.tight_layout()
    return _save(fig, folder, '1A_restart_vs_true_stretch.png')


# ================================================================== Part 1B
def fig1b(b_list, mode, folder):
    """b_list: part1b results per seed"""
    fig, axs = plt.subplots(2, 2, figsize=(14, 9.5))
    axs = axs.ravel()
    Ns = b_list[0]['N']
    mk = ['o', 's', '^']
    for gi in range(3):
        ax = axs[gi]
        for si in range(3):
            if mode == 'annealed' and si == 2:
                pass
            v = np.array([b['rr'][:, gi, si] for b in b_list])
            lo = np.array([b['lo'][:, gi, si] for b in b_list]); hi = np.array([b['hi'][:, gi, si] for b in b_list])
            xs = Ns + (si - 1) * 0.35
            m = np.nanmean(v, 0)
            if v.shape[0] > 1:
                e = 2.776 * np.nanstd(v, 0, ddof=1) / np.sqrt(v.shape[0])
                ax.errorbar(xs, m, yerr=e, fmt=mk[si], capsize=3, label='controlling for ' + A.STRATA[si])
            else:
                ax.errorbar(xs, m, yerr=[m - lo[0], hi[0] - m], fmt=mk[si], capsize=3,
                            label='controlling for ' + A.STRATA[si])
        ax.axhline(1, color='k', lw=1); ax.axhspan(0.95, 1.05, color='g', alpha=0.12, lw=0)
        ax.set_xticks(Ns); ax.set_xlabel('cycle N (current state)')
        ax.set_ylabel('risk ratio of S <= 2 at N+1\n(group 1 / group 0, Mantel-Haenszel)')
        ax.set_title(A.HIST_GROUPS[gi])
        if gi == 0:
            ax.legend(fontsize=8.5)
    ax = axs[3]
    for gi in range(3):
        v = np.array([b['dmed'][:, gi, 1] for b in b_list])
        ax.plot(Ns + (gi - 1) * 0.3, np.nanmean(v, 0), mk[gi] + '-', label=A.HIST_GROUPS[gi])
    ax.axhline(0, color='k', lw=1)
    ax.set_xticks(Ns); ax.set_xlabel('cycle N')
    ax.set_ylabel('difference in median ln(S_N+1 / S_N)\n(group 1 - group 0, within S + orientation bins)')
    ax.set_title('Next-step stretch growth, same comparison'); ax.legend(fontsize=8)
    fig.suptitle('Does the past matter beyond the current M? %s\n(1 = no memory; green band = 1 +/- 0.05; bars: %s)'
                 % (MODE_NAME[mode], '95% interval over seeds' if len(b_list) > 1 else '95% interval within the run'))
    fig.tight_layout()
    return _save(fig, folder, '1B_memory_test.png')


# ================================================================== Part 1C
def fig1c(c_list, folder):
    if not c_list or c_list[0]['t0'].size == 0:
        return None
    fig, ax = plt.subplots(figsize=(9, 6.5))
    cols = plt.get_cmap('viridis')(np.linspace(0, 0.85, c_list[0]['t0'].size))
    for i, t0 in enumerate(c_list[0]['t0']):
        P = np.array([c['P'][i] for c in c_list])
        ax.loglog(A.S_GRID, P.mean(0), '-', color=cols[i], lw=2, label='window of kicks %d-%d' % (t0 + 1, t0 + 10))
        for p in P:
            ax.loglog(A.S_GRID, p, 'o', ms=3, color=cols[i], alpha=0.6)
    ax.set_xlabel('s'); ax.set_ylabel('P(S <= s) after 10 kicks, every pixel')
    ax.set_title('Does the start cycle matter? Annealed: 10-kick windows starting at different cycles\n'
                 '(dots = individual seeds)')
    ax.grid(alpha=0.3, which='both'); ax.legend()
    fig.tight_layout()
    return _save(fig, folder, '1C_start_cycle.png')


# ================================================================== Part 2A
def fig2a_rate(a_by_mode, folder):
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    for m, lst in a_by_mode.items():
        R_ = np.array([a['rate'] for a in lst])
        _seedband(ax, A.S_GRID, R_, MODE_COL[m], MODE_NAME[m])
    ax.set_xscale('log'); ax.axvline(2, color='k', ls=':', lw=1)
    ax.set_xlabel('threshold s'); ax.set_ylabel('decay rate of P(S <= s) per cycle\n(slope of -ln P over cycles 10 to the end)')
    ax.set_title('Low stretch without the number 2: how fast the share below s disappears, for every s')
    ax.grid(alpha=0.3, which='both'); ax.legend(fontsize=9)
    fig.tight_layout()
    return _save(fig, folder, '2A_decay_rate_vs_threshold.png')


def fig2a_heat(a_list, mode, folder):
    P = np.mean([a['P'] for a in a_list], 0)
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    Pp = np.where(P > 0, P, np.nan)
    im = ax.pcolormesh(A.S_GRID, np.arange(P.shape[0]), Pp, norm=LogNorm(1e-5, 1), cmap='viridis', shading='nearest')
    ax.set_xscale('log'); ax.set_xlabel('threshold s'); ax.set_ylabel('cycle N')
    cb = fig.colorbar(im, ax=ax); cb.set_label('P(S <= s)')
    cs = ax.contour(A.S_GRID, np.arange(P.shape[0]), np.log10(np.maximum(P, 1e-12)), [-4, -3, -2, -1], colors='w', linewidths=0.8)
    ax.clabel(cs, fmt=lambda v: '10^%d' % v, fontsize=8)
    ax.axvline(2, color='w', ls=':', lw=1)
    ax.set_title('Share of the picture with S <= s, for every s and every cycle: %s' % MODE_NAME[mode])
    fig.tight_layout()
    return _save(fig, folder, '2A_share_below_s_heat_map.png')


# ================================================================== Part 2B
def fig2b(b_by_mode, folder):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.8))
    ls = {0.01: ':', 0.05: '-', 0.10: '--'}
    for m, lst in b_by_mode.items():
        sq = np.array([b['sq'] for b in lst])
        Ns = np.arange(sq.shape[1])
        for qi, q in enumerate(A.Q_LIST):
            _seedband(axs[0], Ns, np.log(sq[:, :, qi]), MODE_COL[m], '%s, q = %g%%' % (m, 100 * q) if qi == 1 or q == 0.01 else None)
            axs[0].lines[-1].set_linestyle(ls[q])
        pers = np.array([b['pers'] for b in lst])
        for ki, k in enumerate((1, 3, 5)):
            _seedband(axs[1], Ns, pers[:, 1, :, ki], MODE_COL[m], '%s, %d cycle(s) later' % (m, k) if True else None)
            axs[1].lines[-1].set_linestyle(['-', '--', ':'][ki])
        en = np.array([b['branch_enrich'] for b in lst])
        for qi, q in enumerate(A.Q_LIST):
            _seedband(axs[2], np.arange(1, en.shape[1] + 1), en[:, :, qi], MODE_COL[m], '%s, q = %g%%' % (m, 100 * q))
            axs[2].lines[-1].set_linestyle(ls[q])
    axs[0].set_xlabel('cycle N'); axs[0].set_ylabel('ln s_q(N)')
    axs[0].set_title('Relative threshold: the q-quantile of S\n(dotted 1%, solid 5%, dashed 10%)')
    axs[0].legend(fontsize=8)
    axs[1].axhline(0.05, color='k', ls=':', lw=1)
    axs[1].set_xlabel('cycle N'); axs[1].set_ylabel('fraction still in the least-stretched 5%')
    axs[1].set_title('Persistence of the least-stretched 5%\n(dotted black line = chance)'); axs[1].legend(fontsize=7.5)
    axs[1].set_ylim(0, 1)
    axs[2].axhline(1, color='k', ls=':', lw=1)
    axs[2].set_xlabel('cycle N'); axs[2].set_ylabel('enrichment = share / q')
    axs[2].set_title('Is branch material (amount >= 2) among the\nleast-stretched q% of the picture? (1 = by chance)')
    axs[2].legend(fontsize=7.5)
    fig.suptitle('Low stretch without the number 2: a threshold that moves with the run')
    fig.tight_layout()
    return _save(fig, folder, '2B_relative_threshold.png')


# ================================================================== Part 3B
def fig3b(b3, mode, folder):
    if b3['N'].size == 0:
        return None
    n = b3['N'].size
    fig, axs = plt.subplots(1, n, figsize=(4.6 * n, 5.6), squeeze=False)
    e = A.LOGS_BINS_3B
    for i in range(n):
        ax = axs[0, i]
        im = ax.pcolormesh(np.arange(-90, 91, 15), e, b3['heat'][i], cmap='magma', vmin=-1, vmax=2.5, shading='flat')
        ax.set_xlabel('orientation now (deg; 0 = upright, + = leaning right)')
        if i == 0:
            ax.set_ylabel('log10 S')
        r2 = b3['r2'][i]
        ax.set_title('cycle %d\nexplained: S %.0f%%, orientation %.0f%%,\nS + orientation %.0f%%, + launch direction %.0f%%'
                     % (b3['N'][i], *(100 * np.array(r2))), fontsize=9.5)
        ax.set_ylim(0, 4.5)
    cb = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.85); cb.set_label('median log10 brightness (1/|a|)')
    fig.suptitle('Brightness of band material, 1/|a|: what sets it? %s' % MODE_NAME[mode])
    return _save(fig, folder, '3B_brightness_vs_stretch_and_orientation.png')


# ================================================================== Part 3C
def fig3c(c_by_mode, folder):
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.8))
    ax = axs[0]
    cen = np.arange(6) * 15 + 7.5
    ax.plot(cen, A.REENTRY_REF, 'ks--', ms=8, mfc='w', label='random-matrix model (reference)')
    for m, lst in c_by_mode.items():
        i8 = [list(c['N']).index(8) for c in lst if 8 in c['N']]
        if not i8:
            continue
        P = np.array([c['p_unsigned'][list(c['N']).index(8)] for c in lst])
        n = np.array([c['n_unsigned'][list(c['N']).index(8)] for c in lst]).sum(0)
        m_ = P.mean(0)
        se = np.sqrt(m_ * (1 - m_) / np.maximum(n, 1))
        ax.errorbar(cen + (1 if m == 'quenched' else 0), m_, yerr=2 * se, fmt='o-', color=MODE_COL[m], capsize=3,
                    label='%s, 2-D tracers' % m)
    ax.set_xticks(np.arange(0, 91, 15)); ax.set_xlabel('orientation now, degrees from upright (either side)')
    ax.set_ylabel('P(S <= 2 at N+1 | 2 < S <= 2.5 at N)'); ax.set_title('Cycle 8: re-entry vs orientation')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    ax = axs[1]
    cs = np.arange(-90, 90, 15) + 7.5
    for m, lst in c_by_mode.items():
        c = lst[0]
        for j, N in enumerate(c['N']):
            P = np.mean([cc['p_signed'][j] for cc in lst], 0)
            ax.plot(cs, P, ['-', '--', ':', '-.'][j % 4], color=MODE_COL[m], label='%s, cycle %d' % (m, N))
    ax.set_xticks(np.arange(-90, 91, 30)); ax.set_xlabel('orientation now (deg; + = leaning right)')
    ax.set_ylabel('P(S <= 2 at N+1)'); ax.set_title('Signed orientation, all cycles tested')
    ax.legend(fontsize=7.5, ncol=2); ax.grid(alpha=0.3)
    fig.suptitle('Re-entry into S <= 2 depends on which way the stretched material points')
    fig.tight_layout()
    return _save(fig, folder, '3C_reentry_vs_orientation.png')


# ================================================================== Part 3D
def fig3d(r, folder):
    Ns = [N for N in (3, 8, 13) if N <= r.N]
    if not Ns:
        return []
    Tau, row0 = r.prm.Tau, r.prm.row0
    imgs = [picture_and_M(r.kicks, N, Tau, row0) for N in Ns]
    out = []
    for kind in ('stretch', 'orientation'):
        fig, axs = plt.subplots(1, len(Ns), figsize=(6.2 * len(Ns), 6.6), squeeze=False)
        for i, (N, (pic, S, th)) in enumerate(zip(Ns, imgs)):
            ax = axs[0, i]
            im = overlay_S(pic, S) if kind == 'stretch' else overlay_theta(pic, S, th)
            ax.imshow(np.clip(im, 0, 1), interpolation='nearest')
            ax.contour(S, [2], colors='w', linewidths=0.6)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title('cycle %d' % N)
        if kind == 'orientation':
            cax = fig.add_axes([0.92, 0.25, 0.012, 0.5])
            sm = plt.cm.ScalarMappable(norm=Normalize(-90, 90), cmap='twilight')
            cb = fig.colorbar(sm, cax=cax); cb.set_label('orientation of the stretch now (deg)\n0 = upright, +-90 = flat')
            cb.set_ticks([-90, -45, 0, 45, 90])
            fig.suptitle('%s: overlay colour = which way the material is stretched, strength = log S; white line S = 2'
                         % MODE_NAME[r.mode])
            fig.subplots_adjust(left=0.02, right=0.9, top=0.88, bottom=0.03, wspace=0.04)
        else:
            fig.suptitle('%s: the picture with the stretch overlay (S, log scale); white line S = 2' % MODE_NAME[r.mode])
            fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.03, wspace=0.04)
        out.append(_save(fig, folder, '3D_picture_%s.png' % kind, dpi=90))
    return out


# ================================================================== Part 4A
def fig4a(a4_by_mode, folder):
    fig, axs = plt.subplots(1, 3, figsize=(19, 6.2))
    for panel, m in ((0, 'annealed'), (2, 'quenched')):
        if m not in a4_by_mode:
            continue
        d = a4_by_mode[m]
        ax = axs[panel]
        Ns = np.arange(d['lnS_pick'].shape[0])
        ax.fill_between(Ns, d['q'][0], d['q'][2], color='0.85', label='10-90% of all parcels')
        ax.plot(Ns, d['q'][1], 'k-', lw=2, label='median of all parcels')
        for f in (2, 1, 0):
            sel = d['fate'] == f
            ax.plot(Ns, d['lnS_pick'][:, sel], color=FATE_COL[f], lw=0.5, alpha=0.6)
            ax.plot([], [], color=FATE_COL[f], label='%s (%d)' % (FATE_LAB[f], sel.sum()))
        ax.axhline(np.log(2), color='b', ls=':', lw=1)
        ax.set_xlabel('kick number N'); ax.set_ylabel('ln S'); ax.set_title('%s: 300 random pieces' % MODE_NAME[m])
        ax.legend(fontsize=8, loc='upper left')
    if 'annealed' in a4_by_mode:
        d = a4_by_mode['annealed']
        ax = axs[1]
        dl = np.diff(d['lnS_pick'], axis=0)
        Ns = np.arange(1, dl.shape[0] + 1)
        ax.plot(Ns, dl, color='0.5', lw=0.4, alpha=0.5)
        ax.plot([], [], color='0.5', label='ln(S_N / S_N-1) of each piece')
        lj = d['lnJ']
        ax.fill_between(Ns, -np.quantile(lj, 0.9, 1), np.quantile(lj, 0.9, 1), color='r', alpha=0.12,
                        label='+-ln sigma_max(J_N): the most one kick can stretch or un-stretch (90th pct)')
        ax.plot(Ns, np.median(lj, 1), 'r-', lw=1.5, label='ln sigma_max(J_N), median')
        ax.plot(Ns, np.median(dl, 1), 'k-', lw=2, label='median change')
        ax.axhline(0, color='k', lw=0.8)
        ax.set_xlabel('kick number N'); ax.set_ylabel('change of ln S at this kick')
        ax.set_title('annealed: change per kick (below 0 = that kick un-stretched it)')
        ax.legend(fontsize=7.5, loc='lower right')
    fig.suptitle('How the stretch of individual pieces of material grows, and sometimes shrinks, kick by kick')
    fig.tight_layout()
    return _save(fig, folder, '4A_stretch_kick_by_kick.png')


# ================================================================== Part 4B
def fig4b(b4_by_mode, sstar, folder):
    paths = []
    c = 10 ** (0.5 * (S2.LOGS_EDGES[1:] + S2.LOGS_EDGES[:-1]))
    # (i) where cusps are born
    fig, axs = plt.subplots(1, len(b4_by_mode), figsize=(7.5 * len(b4_by_mode), 5.6), squeeze=False)
    for i, (m, lst) in enumerate(b4_by_mode.items()):
        ax = axs[0, i]
        hc = np.sum([b['hist_cusps'] for b in lst], 0); hm = np.sum([b['mat_at_cusps'] for b in lst], 0)
        ax.step(c, hc / hc.sum() / 0.05, where='mid', color=MODE_COL[m], lw=2, label='S of the material at each cusp')
        ax.step(c, hm / hm.sum() / 0.05, where='mid', color='k', lw=1.2,
                label='S of all centre-line material at the same moments')
        ax.axvline(2, color='b', ls=':', lw=1, label='S = 2')
        ax.set_xscale('log'); ax.set_xlim(1, 1e4); ax.set_xlabel('S'); ax.set_ylabel('density per decade')
        ax.set_title('%s: %d cusps' % (m, int(hc.sum()))); ax.legend(fontsize=8.5)
    fig.suptitle('Are cusps born in low- or high-stretch material? (band centre lines, exact)')
    fig.tight_layout(); paths.append(_save(fig, folder, '4B_i_stretch_at_cusps.png'))
    # (ii) cusp rate per unit material
    fig, ax = plt.subplots(figsize=(9.5, 6.3))
    for m, lst in b4_by_mode.items():
        e = lst[0]['rate_edges']; cc = 10 ** (0.5 * (e[1:] + e[:-1]))
        R_ = np.array([b['rate'] for b in lst])
        _seedband(ax, cc, R_, MODE_COL[m], m)
        s5 = np.nanmedian([b['s5'] for b in lst])
        ax.axvline(s5, color=MODE_COL[m], ls='--', lw=1, label='%s: 5%% quantile of S (median over cycles 10+)' % m)
    ax.axvline(2, color='b', ls=':', lw=1.2, label='S = 2')
    ax.axvline(sstar, color='m', ls='-.', lw=1.2, label='S* = xi / w = %.2f (hypothesis)' % sstar)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1, 1e4)
    ax.set_xlabel('S of the material'); ax.set_ylabel('cusps per launch px of line per cycle')
    ax.set_title('How often material of a given stretch gives birth to a cusp'); ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which='both')
    fig.tight_layout(); paths.append(_save(fig, folder, '4B_ii_cusp_rate_vs_stretch.png'))
    # (iii) branches vs S at their cusp
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    for j, m in enumerate([mm for mm in ('annealed', 'quenched') if mm in b4_by_mode]):
        br = np.vstack([b['branches'] for b in b4_by_mode[m]])
        ok = np.isfinite(br[:, 7])
        for i, (col, lab) in enumerate(((2, 'lifetime (cycles)'), (5, 'mean amount over the branch'))):
            ax = axs[i, j]
            yv = br[ok, 2] - br[ok, 1] if col == 2 else br[ok, 5]
            ax.semilogx(br[ok, 7], yv, 'o', color=MODE_COL[m], ms=5, alpha=0.7)
            if ok.sum() >= 6:
                ed = np.quantile(np.log10(br[ok, 7]), np.linspace(0, 1, 6))
                bi = np.clip(np.digitize(np.log10(br[ok, 7]), ed[1:-1]), 0, 4)
                med = [np.median(yv[bi == k]) if np.any(bi == k) else np.nan for k in range(5)]
                ax.semilogx(10 ** (0.5 * (ed[1:] + ed[:-1])), med, 'k-s', lw=2, label='binned medians')
            ax.axvline(2, color='b', ls=':'); ax.axvline(sstar, color='m', ls='-.')
            ax.set_xlabel('S of the material at the branch\'s cusp'); ax.set_ylabel(lab)
            ax.set_title('%s: %d of %d branches matched to a cusp' % (m, ok.sum(), br.shape[0]))
            ax.legend(fontsize=8)
    fig.suptitle('Do branches born in low-stretch material live longer or shine brighter?')
    fig.tight_layout(); paths.append(_save(fig, folder, '4B_iii_branches_vs_stretch_at_cusp.png'))
    # (iv) S profile around cusps
    fig, axs = plt.subplots(1, len(b4_by_mode), figsize=(7.5 * len(b4_by_mode), 5.6), squeeze=False)
    for i, (m, lst) in enumerate(b4_by_mode.items()):
        ax = axs[0, i]
        pr = np.concatenate([b['profiles'] for b in lst if b['profiles'].size])
        for li, lag in enumerate(S2.PROFILE_LAGS):
            v = pr[:, li, :]
            ok = np.all(np.isfinite(v), 1)
            if ok.sum() < 10:
                continue
            ax.semilogy(S2.PROFILE_U, np.median(v[ok], 0), lw=2, label='%d cycle(s) after the cusp (%d cusps)' % (lag, ok.sum()))
        ax.axhline(2, color='b', ls=':'); ax.set_xlabel('distance along the line from the cusp (launch px)')
        ax.set_ylabel('median S'); ax.set_title(m); ax.legend(fontsize=8.5)
    fig.suptitle('Stretch along the line around a cusp, at the cusp and later')
    fig.tight_layout(); paths.append(_save(fig, folder, '4B_iv_stretch_profile_around_cusps.png'))
    return paths


def fig4b_example(r, b4, folder, tc_target=2.44, xc_target=765.0):
    """Worked example: the cusp nearest (2.44, x 765) on the red band."""
    cu = b4['cusps']
    red = cu[cu[:, 7] == 0]
    if red.size == 0:
        return None
    i = np.argmin(((red[:, 0] - tc_target) / 0.1) ** 2 + (A._circ_dx(red[:, 1], xc_target) / 30) ** 2)
    tc, xc, yc = red[i, 0], red[i, 1], red[i, 2]
    N = int(np.ceil(tc)); t = tc - (N - 1)
    Tau, row0 = r.prm.Tau, r.prm.row0
    win = (int(xc) - 180, int(xc) + 180, int(yc) - 130, int(yc) + 130)
    pic, S, th = picture_and_M(r.kicks, N, Tau, row0, t_last=t, window=win)
    # the red centre line at that moment, coloured by S
    P = S2.line_new(S2.band_centre('red'))
    for k in range(N - 1):
        P, _ = S2.line_step(P, r.kicks[k], Tau, row0)
    _, D = S2.line_step(P, r.kicks[N - 1], Tau, row0)
    a, b, c, d = S2.drift_M(D, t)
    xl = D['xl'] + t * Tau * (D['yk'] - row0)
    xr = xl + a * D['L']; yl = D['yk']; yr = yl + c * D['L']
    Sl = S2.sigma(a, b, c, d)
    sel = (np.mod(xl - win[0], S2.R) < 360) & (yl > win[2]) & (yl < win[3])
    fig = plt.figure(figsize=(12, 12.5))
    g = fig.add_gridspec(2, 1, height_ratios=[1, 0.8], hspace=0.25)
    ax = fig.add_subplot(g[0])
    ax.imshow(np.clip(overlay_S(pic, S), 0, 1), extent=[win[0], win[1], win[3], win[2]], interpolation='nearest')
    ax.contour(np.arange(win[0], win[1]), np.arange(win[2], win[3]), S, [2], colors='w', linewidths=0.7)
    from matplotlib.collections import LineCollection
    x0 = win[0] + np.mod(xl[sel] - win[0], S2.R); x1 = x0 + (xr[sel] - xl[sel])
    segs = np.stack([np.stack([x0, yl[sel]], -1), np.stack([x1, yr[sel]], -1)], 1)
    lc = LineCollection(segs, cmap='plasma', norm=LogNorm(1, 50), linewidths=3)
    lc.set_array(Sl[sel]); ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, fraction=0.04); cb.set_label('S along the red band\'s centre line')
    ax.plot(xc, yc, '*', ms=26, mfc='cyan', mec='k', mew=2)
    ax.set_xlim(win[0], win[1]); ax.set_ylim(win[3], win[2])
    ax.set_xlabel('x (px)'); ax.set_ylabel('y (px)')
    ax.set_title('Cycle %.3f, red band, %s: the cusp (star) where the centre line first stands exactly upright\n'
                 'S there = %.2f' % (tc, r.mode, red[i, 4]))
    ax = fig.add_subplot(g[1])
    fine = np.array(r.rec['fine'])[:, 0, :]
    tt = np.arange(fine.shape[0]) / S2.NSUB
    xcs = (np.arange(fine.shape[1]) + 0.5) * S2.BIN
    im = ax.pcolormesh(xcs, tt, np.maximum(fine, 0.05), norm=LogNorm(0.05, 50), cmap='magma', shading='nearest')
    bl, lab = S2.branches(fine)
    ax.contour(xcs, tt, (lab > 0).astype(float), [0.5], colors='cyan', linewidths=1)
    ax.plot(xc, tc, '*', ms=22, mfc='cyan', mec='k', mew=2)
    ax.set_xlim(xc - 350, xc + 350); ax.set_ylim(max(0, tc - 1.5), min(tt[-1], tc + 5))
    ax.set_xlabel('x (px)'); ax.set_ylabel('time (cycles)')
    cb = fig.colorbar(im, ax=ax, fraction=0.04); cb.set_label('amount of the red band per 4-px strip')
    ax.set_title('The red band\'s fine-time map near the cusp: branches outlined in cyan')
    return _save(fig, folder, '4B_worked_example_cusp.png', dpi=90)
