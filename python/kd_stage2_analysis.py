"""
kd_stage2_analysis -- the stage-2 questions, computed from finished runs.

Every function takes a finished kd_stage2.Stage2Run (or the pieces it needs)
and returns plain numpy arrays / scalars, so results stack over seeds.
"""
import numpy as np
import kd_engine as E
import kd_stage2 as S2

S_GRID = np.geomspace(1.1, 1e4, 40)                 # Part 2A thresholds
Q_LIST = (0.01, 0.05, 0.10)                          # Part 2B quantiles
MEM_N = (5, 10, 15, 20)                              # Part 1B cycles
MEM_S_EDGES = 2.0 ** (np.arange(1, 8) / 3.0)         # 8 bins of S: <1.26, ..., 2 is an edge, >5.04
THETA_EDGES = np.arange(-90.0, 90.01, 15.0)          # 12 signed bins of 15 deg
REENTRY_REF = np.array([0.04, 0.25, 0.31, 0.30, 0.21, 0.12])   # random-matrix model, annealed, N = 8
S_THR = 2.0


def theta_bin(th):
    """signed 15-deg bins over (-90, 90]: 0 .. 11"""
    return np.clip(np.ceil((np.asarray(th) + 90.0) / 15.0).astype(int) - 1, 0, 11)


def S_hist(r, shifted=False):
    return np.array(r.rec['S_T2' if shifted else 'S_T'])


# ====================================================================== Part 1A
def part1a(mode, seed, n_first=10, n_more=5):
    """Checkpoint after n_first cycles, resume n_more, against one straight run.
    Then the WRONG way: reset M to I at cycle n_first and run n_more."""
    tot = n_first + n_more
    run = S2.Stage2Run(mode, seed, tot).run(n_first)
    state = run.save_state()
    run.run(n_more)
    res = S2.Stage2Run.resume(state, n_more)
    diffs = {}

    def md(a, b):
        return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float)))) if np.size(a) else 0.0
    for key in ('x', 'y', 'a', 'b', 'c', 'd'):
        diffs['tracers ' + key] = md(run.T[key], res.T[key])
        diffs['tracers shifted ' + key] = md(run.T2[key], res.T2[key])
        diffs['band material ' + key] = max(md(p[key], q[key]) for p, q in zip(run.B, res.B))
        diffs['centre-line pieces ' + key] = max(md(p[key], q[key]) for p, q in zip(run.Lp, res.Lp))
    diffs['tracers S'] = md(run.rec['S_T'][-1], res.rec['S_T'][-1])
    diffs['band amounts'] = md(run.rec['amount'][-1], res.rec['amount'][-1])
    diffs['fine-time maps'] = md(np.array(run.rec['fine']), np.array(res.rec['fine']))
    diffs['kicks'] = md(np.array(run.kicks), np.array(res.kicks))
    cu_a = np.vstack([np.vstack(c) for c in run.rec['cusps']])
    cu_b = np.vstack([np.vstack(c) for c in res.rec['cusps']])
    diffs['cusps'] = md(cu_a, cu_b) if cu_a.shape == cu_b.shape else np.inf
    # the wrong way: forget M at cycle n_first
    T10 = state['T']
    wr = E.new_state(T10['x'], T10['y'])
    for k in range(n_first, tot):
        E.step(wr, run.kicks[k], run.prm.Tau, run.prm.row0)
    S_true = S2.sigma(run.T['a'], run.T['b'], run.T['c'], run.T['d'])
    S10 = S2.sigma(T10['a'], T10['b'], T10['c'], T10['d'])
    S_wrong = S2.sigma(wr['a'], wr['b'], wr['c'], wr['d'])
    prod = S10 * S_wrong
    viol = int(np.sum(S_true > prod * (1 + 1e-12)))
    low_t = S_true <= S_THR; low_w = S_wrong <= S_THR
    sub = np.random.default_rng(7).choice(S_true.size, 20000, replace=False)
    return dict(diffs=diffs, n_parcels=int(S_true.size), bound_violations=viol,
                n_low_true=int(low_t.sum()), n_low_wrong=int(low_w.sum()),
                n_changed=int(np.sum(low_t != low_w)), n_wrong_says_low=int(np.sum(low_w & ~low_t)),
                n_wrong_says_high=int(np.sum(~low_w & low_t)),
                scatter_true=S_true[sub], scatter_prod=prod[sub], scatter_wrong=S_wrong[sub],
                frac_at_bound=float(np.mean(S_true >= 0.99 * prod)))


# ====================================================================== Part 1B
def _mh_rr(sid, g, y):
    """Mantel-Haenszel risk ratio (group g / not g) of y, over strata sid, with
    the Greenland-Robins 95% interval."""
    m = sid.max() + 1
    n1 = np.bincount(sid[g], minlength=m).astype(float)
    n0 = np.bincount(sid[~g], minlength=m).astype(float)
    a = np.bincount(sid[g & y], minlength=m).astype(float)
    c = np.bincount(sid[~g & y], minlength=m).astype(float)
    n = n1 + n0
    ok = (n1 > 0) & (n0 > 0)
    n1, n0, a, c, n = n1[ok], n0[ok], a[ok], c[ok], n[ok]
    Rn = np.sum(a * n0 / n); Sd = np.sum(c * n1 / n)
    if Rn <= 0 or Sd <= 0:
        return np.nan, np.nan, np.nan, int(np.sum(n1)), int(np.sum(n0))
    rr = Rn / Sd
    P = np.sum((n1 * n0 * (a + c) - a * c * n) / n ** 2)
    se = np.sqrt(max(P, 0) / (Rn * Sd))
    return rr, rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se), int(np.sum(n1)), int(np.sum(n0))


def _grouped_median(sid, v):
    o = np.lexsort((v, sid))
    s, vv = sid[o], v[o]
    m = sid.max() + 1
    cnt = np.bincount(s, minlength=m)
    start = np.r_[0, np.cumsum(cnt)[:-1]]
    med = np.full(m, np.nan)
    ok = cnt > 0
    lo = start[ok] + (cnt[ok] - 1) // 2
    hi = start[ok] + cnt[ok] // 2
    med[ok] = 0.5 * (vv[lo] + vv[hi])
    return med, cnt


def _strat_median_diff(sid, g, v, nmin=5):
    """pooled difference of median v (group 1 - group 0) within strata, weights n1 n0 / n"""
    m = sid.max() + 1
    m1, c1 = _grouped_median(sid[g], v[g]) if g.any() else (np.full(m, np.nan), np.zeros(m))
    m0, c0 = _grouped_median(sid[~g], v[~g]) if (~g).any() else (np.full(m, np.nan), np.zeros(m))
    m1 = np.r_[m1, np.full(m - m1.size, np.nan)]; c1 = np.r_[c1, np.zeros(m - c1.size)]
    m0 = np.r_[m0, np.full(m - m0.size, np.nan)]; c0 = np.r_[c0, np.zeros(m - c0.size)]
    ok = (c1 >= nmin) & (c0 >= nmin)
    if not ok.any():
        return np.nan
    w = c1[ok] * c0[ok] / (c1[ok] + c0[ok])
    return float(np.sum(w * (m1[ok] - m0[ok])) / np.sum(w))


HIST_GROUPS = ('re-entered before vs never', 'S two cycles ago above vs below the bin median',
               '1 vs >= 3 cycles since the last crossing')
STRATA = ('S', 'S + orientation', 'S + orientation + position')


def part1b(r):
    S = S_hist(r)
    C = S.shape[0] - 1
    low = S <= S_THR
    out = dict(N=[], rr=[], lo=[], hi=[], n1=[], n0=[], dmed=[])
    for N in MEM_N:
        if N + 1 > C or N not in r.rec['xy_T']:
            continue
        Sn, Sx = S[N], S[N + 1]
        yb = Sx <= S_THR
        yc = np.log(Sx / Sn)
        sb = np.digitize(Sn, MEM_S_EDGES)                       # 0..7
        tb = theta_bin(r.rec['th_T'][N])                        # 0..11
        x, y = r.rec['xy_T'][N]
        pb = (np.clip((np.mod(y, S2.R) // 240).astype(int), 0, 7) * 8 +
              np.clip((x // 240).astype(int), 0, 7))            # 0..63
        sids = [sb, sb * 12 + tb, (sb * 12 + tb) * 64 + pb]
        # history groups
        ent = (~low[:-1]) & low[1:]                             # entered the low set at k = 1..C
        g_reent = ent[:N].any(0)
        med = np.array([np.median(S[N - 2][sb == k]) if np.any(sb == k) else np.nan for k in range(8)])
        g_prev = S[N - 2] > med[sb]
        cross = low[1:N + 1] != low[:N]                         # crossing at k = 1..N
        anyc = cross.any(0)
        lastk = np.where(anyc, N - np.argmax(cross[::-1], axis=0), 0)
        age = N - lastk + 1
        groups = [(g_reent, np.ones_like(g_reent)), (g_prev, np.ones_like(g_prev)),
                  (age == 1, anyc & ((age == 1) | (age >= 3)))]
        rr = np.full((3, 3), np.nan); lo = rr.copy(); hi = rr.copy(); dm = rr.copy()
        n1 = np.zeros((3, 3)); n0 = np.zeros((3, 3))
        for gi, (g, use) in enumerate(groups):
            for si, sid in enumerate(sids):
                rr[gi, si], lo[gi, si], hi[gi, si], n1[gi, si], n0[gi, si] = _mh_rr(sid[use], g[use], yb[use])
                dm[gi, si] = _strat_median_diff(sid[use], g[use], yc[use])
        for k, v in zip(('N', 'rr', 'lo', 'hi', 'n1', 'n0', 'dmed'), (N, rr, lo, hi, n1, n0, dm)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


# ====================================================================== Part 1C
def window_windows(C):
    return [t0 for t0 in (0, 10, 20) if t0 + 10 <= C] if C >= 30 else [t0 for t0 in (0, 5) if t0 + 10 <= C]


def part1c(r, step=1, block=1_000_000):
    """Annealed: P(S <= s) of a 10-kick window starting at t0, every pixel."""
    if r.mode != 'annealed':
        return dict(t0=np.zeros(0), P=np.zeros((0, S_GRID.size)), note='identical by construction (same kick every cycle)')
    C = len(r.kicks)
    t0s = window_windows(C)
    g = np.arange(0, S2.R, step, dtype=float)
    X, Y = np.meshgrid(g, g)
    X = X.ravel(); Y = Y.ravel()
    P = np.zeros((len(t0s), S_GRID.size))
    for i, t0 in enumerate(t0s):
        cnt = np.zeros(S_GRID.size)
        for s0 in range(0, X.size, block):
            st = E.new_state(X[s0:s0 + block], Y[s0:s0 + block])
            for k in range(t0, t0 + 10):
                E.step(st, r.kicks[k], r.prm.Tau, r.prm.row0)
            Sv = S2.sigma(st['a'], st['b'], st['c'], st['d'])
            cnt += (Sv[:, None] <= S_GRID[None, :]).sum(0)
        P[i] = cnt / X.size
    return dict(t0=np.array(t0s), P=P, note='')


# ====================================================================== Part 1D
def headline(S, C):
    """A few headline numbers of one tracer population."""
    P2 = (S <= S_THR).mean(1)
    out = {}
    for N in (5, 10, 15, 20, 30):
        if N <= C:
            out['share with S <= 2 at cycle %d' % N] = P2[N]
    Ns = np.arange(10, C + 1)
    if Ns.size >= 3 and np.all(P2[Ns] > 0):
        out['decay rate of P(S <= 2), cycles 10-%d' % C] = -np.polyfit(Ns, np.log(P2[Ns]), 1)[0]
    if C >= 10:
        out['5% quantile of S at cycle 10'] = np.quantile(S[10], 0.05)
        q = np.quantile(S[10], 0.05); qn = np.quantile(S[11], 0.05) if C >= 11 else np.nan
        if C >= 11:
            inn = S[10] <= q
            out['least-stretched 5% at cycle 10 still in it at 11'] = np.mean(S[11][inn] <= qn)
    return out


def part1d(r):
    S = S_hist(r); S2_ = S_hist(r, shifted=True)
    C = S.shape[0] - 1
    return headline(S, C), headline(S2_, C)


# ====================================================================== Part 2A / 2B
def part2a(r):
    S = S_hist(r)
    C = S.shape[0] - 1
    P = np.stack([(S[N][:, None] <= S_GRID[None, :]).mean(0) for N in range(C + 1)])
    n = S.shape[1]
    Ns = np.arange(10, C + 1)
    rate = np.full(S_GRID.size, np.nan)
    for j in range(S_GRID.size):
        p = P[Ns, j]
        ok = p * n >= 20
        if ok.sum() >= 3:
            rate[j] = -np.polyfit(Ns[ok], np.log(p[ok]), 1)[0]
    return dict(P=P, rate=rate)


def part2b(r):
    S = S_hist(r)
    C = S.shape[0] - 1
    sq = np.stack([np.quantile(S[N], Q_LIST) for N in range(C + 1)])       # [N, q]
    pers = np.full((len(Q_LIST), C + 1, 3), np.nan)                        # [q, N, k in (1,3,5)]
    for qi in range(len(Q_LIST)):
        for N in range(1, C + 1):
            inn = S[N] <= sq[N, qi]
            for ki, k in enumerate((1, 3, 5)):
                if N + k <= C:
                    pers[qi, N, ki] = np.mean(S[N + k][inn] <= sq[N + k, qi])
    bl = np.array(r.rec['branch_link'])                                    # rows N = 1..C
    share = bl[:, 4:7] / np.maximum(bl[:, 3:4], 1)
    enrich = share / np.array(Q_LIST)[None, :]
    return dict(sq=sq, pers=pers, branch_share=share, branch_enrich=enrich, branch_n=bl[:, 3])


# ====================================================================== Part 3A
def part3a(r, n=100_000):
    rng = np.random.default_rng(11)
    pools = [r.T] + list(r.B)
    A = np.concatenate([p['a'] for p in pools]); B = np.concatenate([p['b'] for p in pools])
    Cc = np.concatenate([p['c'] for p in pools]); D = np.concatenate([p['d'] for p in pools])
    i = rng.choice(A.size, n, replace=False)
    a, b, c, d = A[i], B[i], Cc[i], D[i]
    M = np.stack([np.stack([a, b], -1), np.stack([c, d], -1)], -2)
    U, s, Vt = np.linalg.svd(M)
    al, be = S2.axis_angles(a, b, c, d)
    al_s = np.arctan2(U[:, 1, 0], U[:, 0, 0]); be_s = np.arctan2(Vt[:, 0, 1], Vt[:, 0, 0])
    err_u = np.abs(np.mod(al - al_s + np.pi / 2, np.pi) - np.pi / 2)
    err_v = np.abs(np.mod(be - be_s + np.pi / 2, np.pi) - np.pi / 2)
    S = S2.sigma(a, b, c, d)
    u1 = np.stack([np.cos(al), np.sin(al)], -1); v1 = np.stack([np.cos(be), np.sin(be)], -1)
    Mv1 = np.einsum('nij,nj->ni', M, v1)
    sg = np.sign(np.einsum('ni,ni->n', u1, Mv1)); sg[sg == 0] = 1
    u1 = u1 * sg[:, None]
    u2 = np.stack([-u1[:, 1], u1[:, 0]], -1); v2 = np.stack([-np.sin(be), np.cos(be)], -1)
    Mv2 = np.einsum('nij,nj->ni', M, v2)
    sg2 = np.sign(np.einsum('ni,ni->n', u2, Mv2)); sg2[sg2 == 0] = 1
    u2 = u2 * sg2[:, None]
    a_rec = S * u1[:, 0] * v1[:, 0] + (1 / S) * u2[:, 0] * v2[:, 0]
    rel = np.abs(a_rec - a) / S
    S_svd_err = float(np.max(np.abs(S - s[:, 0]) / s[:, 0]))
    return dict(max_angle_err_u=float(err_u.max()), max_angle_err_v=float(err_v.max()),
                max_identity_rel=float(rel.max()), max_S_rel=S_svd_err, min_S=float(S.min()), n=n)


# ====================================================================== Part 3B
LOGS_BINS_3B = np.linspace(0, 6, 25)


def _cv_r2(keys, y, rng):
    """2-fold cross-validated share of the variance of y explained by the mean of
    y in each bin (keys = integer bin ids)."""
    half = rng.random(y.size) < 0.5
    sse = 0.0
    m = keys.max() + 1
    for tr in (half, ~half):
        te = ~tr
        s = np.bincount(keys[tr], weights=y[tr], minlength=m)
        n = np.bincount(keys[tr], minlength=m)
        mu = np.where(n > 0, s / np.maximum(n, 1), y[tr].mean())
        sse += np.sum((y[te] - mu[keys[te]]) ** 2)
    return 1 - sse / np.sum((y - y.mean()) ** 2)


def part3b(r):
    rng = np.random.default_rng(3)
    out = dict(N=[], heat=[], r2=[], counts=[])
    for N in (3, 5, 10, 15):
        if N not in r.rec['band_snap']:
            continue
        sn = r.rec['band_snap'][N]
        a = np.concatenate([s['a'] for s in sn]).astype(float)
        S = np.concatenate([s['S'] for s in sn]).astype(float)
        th = np.concatenate([s['th'] for s in sn]).astype(float)
        ph = np.concatenate([s['ph'] for s in sn]).astype(float)
        ok = a != 0
        y = np.log10(1 / np.abs(a[ok])); lS = np.log10(S[ok]); th = th[ok]; ph = ph[ok]
        sb = np.clip(np.digitize(lS, LOGS_BINS_3B) - 1, 0, LOGS_BINS_3B.size - 2)
        tb = theta_bin(th)
        hid = sb * 12 + tb
        med, cnt = _grouped_median(hid, y)
        nbS = LOGS_BINS_3B.size - 1
        med = np.r_[med, np.full(nbS * 12 - med.size, np.nan)].reshape(nbS, 12)
        cnt = np.r_[cnt, np.zeros(nbS * 12 - cnt.size)].reshape(nbS, 12)
        med[cnt < 20] = np.nan
        qS = np.quantile(lS, np.linspace(0, 1, 25)[1:-1])
        kS = np.digitize(lS, qS)                                  # 24 quantile bins of S
        kT = np.clip(((th + 90) // 7.5).astype(int), 0, 23)       # 24 bins of 7.5 deg
        kP = np.clip(((ph + 90) // 7.5).astype(int), 0, 23)
        r2 = [_cv_r2(kS, y, rng), _cv_r2(kT, y, rng), _cv_r2(kS * 24 + kT, y, rng),
              _cv_r2((kS * 24 + kT) * 24 + kP, y, rng)]
        out['N'].append(N); out['heat'].append(med); out['r2'].append(r2); out['counts'].append(cnt)
    return {k: np.array(v) for k, v in out.items()}


# ====================================================================== Part 3C
def part3c(r):
    S = S_hist(r)
    C = S.shape[0] - 1
    out = dict(N=[], p_signed=[], n_signed=[], p_unsigned=[], n_unsigned=[])
    for N in (5, 8, 10, 15):
        if N + 1 > C:
            continue
        sel = (S[N] > 2.0) & (S[N] <= 2.5)
        y = S[N + 1][sel] <= S_THR
        th = r.rec['th_T'][N][sel].astype(float)
        tb = theta_bin(th)
        n_s = np.bincount(tb, minlength=12); k_s = np.bincount(tb, weights=y, minlength=12)
        ub = np.clip((np.abs(th) // 15).astype(int), 0, 5)
        n_u = np.bincount(ub, minlength=6); k_u = np.bincount(ub, weights=y, minlength=6)
        out['N'].append(N)
        out['p_signed'].append(k_s / np.maximum(n_s, 1)); out['n_signed'].append(n_s)
        out['p_unsigned'].append(k_u / np.maximum(n_u, 1)); out['n_unsigned'].append(n_u)
    return {k: np.array(v) for k, v in out.items()}


# ====================================================================== Part 4A
def part4a(r):
    S = S_hist(r)
    lnS = np.log(S)
    pk = r.pick
    lp = lnS[:, pk]
    high = S[:, pk] > S_THR
    ever = high.any(0)
    first = np.argmax(high, 0)
    back = np.array([np.any(~high[f:, i]) if ever[i] else False for i, f in enumerate(first)])
    fate = np.where(~ever, 0, np.where(back, 1, 2))        # never above / left and came back / left for good
    return dict(lnS_pick=lp, fate=fate, lnJ=np.array(r.rec['lnJ_pick']),
                q=np.quantile(lnS, [0.1, 0.5, 0.9], axis=1))


# ====================================================================== Part 4B
def _circ_dx(a, b, cols=S2.R):
    return np.abs(np.mod(a - b + cols / 2, cols) - cols / 2)


def part4b(r, s5):
    """Cusps of the three exact centre lines, matched to branches of the same band."""
    C = r.N
    fine = np.array(r.rec['fine'])                          # [time, band, strip]
    tt = np.arange(fine.shape[0]) / S2.NSUB
    cus, brs, prof = [], [], []
    mat = np.zeros(S2.LOGS_EDGES.size - 1)
    mat_cyc = np.zeros((C, S2.LOGS_EDGES.size - 1))
    for k in range(3):
        cu = np.vstack(r.rec['cusps'][k]) if r.rec['cusps'][k] else np.zeros((0, 8))
        bl, _ = S2.branches(fine[:, k, :])
        mh = np.array(r.rec['mat_hist'][k])
        mat += mh.sum(0); mat_cyc += mh
        # match each cusp to the nearest branch start within 0.5 cycle and 20 px
        mb = np.full(cu.shape[0], -1)
        if cu.shape[0] and bl.shape[0]:
            dt = np.abs(cu[:, 0:1] - bl[None, :, 1]); dx = _circ_dx(cu[:, 1:2], bl[None, :, 3])
            dd = (dt / 0.5) ** 2 + (dx / 20.0) ** 2
            dd[(dt > 0.5) | (dx > 20)] = np.inf
            j = np.argmin(dd, 1)
            mb = np.where(np.isfinite(dd[np.arange(cu.shape[0]), j]), j, -1)
        # each branch: its nearest matched cusp
        bS = np.full(bl.shape[0], np.nan); bth = bS.copy(); bt = bS.copy()
        for jb in range(bl.shape[0]):
            ii = np.nonzero(mb == jb)[0]
            if ii.size:
                dd2 = ((cu[ii, 0] - bl[jb, 1]) / 0.5) ** 2 + (_circ_dx(cu[ii, 1], bl[jb, 3]) / 20.0) ** 2
                i0 = ii[np.argmin(dd2)]
                bS[jb], bth[jb], bt[jb] = cu[i0, 4], cu[i0, 5], cu[i0, 0]
        cus.append(np.column_stack([cu[:, :7], np.full(cu.shape[0], k), mb]))
        brs.append(np.column_stack([bl, np.full(bl.shape[0], k), bS, bth, bt]) if bl.shape[0] else np.zeros((0, 10)))
        for ci, pr in r.rec['prof'][k].items():
            prof.append(pr)
    cus = np.vstack(cus) if cus else np.zeros((0, 9))
    brs = np.vstack(brs) if brs else np.zeros((0, 10))
    prof = np.array(prof) if prof else np.zeros((0, len(S2.PROFILE_LAGS), S2.PROFILE_U.size))
    # (i) cusps vs all centre-line material at the same moments
    ncyc = np.bincount(np.minimum(np.floor(cus[:, 0]).astype(int), C - 1), minlength=C)[:C] if cus.size else np.zeros(C)
    w = mat_cyc / np.maximum(mat_cyc.sum(1, keepdims=True), 1e-300)
    mat_at_cusps = (w * ncyc[:, None]).sum(0)
    lg_c = np.log10(cus[:, 4]) if cus.size else np.zeros(0)
    hc = np.histogram(lg_c, S2.LOGS_EDGES)[0].astype(float)
    # (ii) cusp rate per launch px per cycle, against S (coarser bins: 0.1 in log10 S)
    edges = np.arange(0, 6.01, 0.1)
    fold = lambda h: np.add.reduceat(h[:(edges.size - 1) * 2], np.arange(0, (edges.size - 1) * 2, 2))
    rate = fold(hc) / np.maximum(fold(mat), 1e-300)
    rate[fold(mat) < 50] = np.nan                          # less than 50 launch-px cycles of material
    return dict(cusps=cus, branches=brs, profiles=prof, mat=mat, mat_at_cusps=mat_at_cusps,
                hist_cusps=hc, rate_edges=edges, rate=rate, n_cusps=int(cus.shape[0]),
                n_matched=int(np.sum(cus[:, 8] >= 0)) if cus.size else 0,
                n_branches=int(brs.shape[0]), n_branches_with_cusp=int(np.sum(np.isfinite(brs[:, 7]))) if brs.size else 0,
                s5=s5)


# ====================================================================== Part 2C
def part2c(r, b4, sq):
    """The S* hypothesis: S* = xi / w with xi = lc*cols (kick correlation length)
    and w = band thickness."""
    xi = r.prm.lc * r.prm.cols
    w = S2.band_rows('red').size
    Sst = xi / w
    # (a) cusp rate just below vs just above S*
    e = b4['rate_edges']; c = 10 ** (0.5 * (e[1:] + e[:-1]))
    below = (c >= Sst / 1.5) & (c < Sst); above = (c >= Sst) & (c < Sst * 1.5)
    ra = np.nanmean(b4['rate'][above]) / np.nanmean(b4['rate'][below])
    # (b) median branch lifetime above minus below a cut on S at the branch's cusp
    br = b4['branches']
    life = br[:, 2] - br[:, 1]
    Sc = br[:, 7]
    ok = np.isfinite(Sc)
    tc = br[:, 9]
    s5_at = np.where(ok, np.interp(np.nan_to_num(tc), np.arange(sq.shape[0]), sq[:, 1]), np.nan)
    cuts = {'S* = %.2f' % Sst: np.full(Sc.shape, Sst), 'S = 2': np.full(Sc.shape, 2.0),
            '5% quantile of S at that moment': s5_at}
    sep = {}
    for nm, cv in cuts.items():
        lo_ = ok & (Sc <= cv); hi_ = ok & (Sc > cv)
        sep[nm] = (np.median(life[hi_]) - np.median(life[lo_]) if lo_.sum() >= 3 and hi_.sum() >= 3 else np.nan,
                   int(lo_.sum()), int(hi_.sum()))
    return dict(xi=xi, w=w, S_star=Sst, rate_ratio_above_below=ra, sep=sep)
