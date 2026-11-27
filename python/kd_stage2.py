"""
kd_stage2 -- stage-2 engine: one run that can be checkpointed and resumed.

Same map, kicks and random stream as the MATLAB script and kd_engine.py:
  one step   y' = y + f(x);  x' = x + Tau*(y' - row0)   (f = piecewise-linear dp)
  M          [a b; c d], rows (x, y);  J = [1+Tau*f', Tau; f', 1]
  S          sigma_max(M)

A run carries
  * the tracer grid (TRACER_STEP px over the 1920 x 1920 frame) and the same
    grid shifted by half a spacing, each parcel with its own M;
  * the material of the red, green and blue bands (1 px in y, 0.5 px in x)
    with M, giving amount(x) per 4-px strip every cycle and at NSUB moments
    inside every drift (the fine-time maps the branches are found on);
  * each band's centre line, tracked EXACTLY as a list of affine pieces (the
    kick is piecewise linear, so the line is piecewise affine): caustics and
    cusps are exact, at every cycle, with no sampling.
The kicks are drawn cycle by cycle from the MATLAB stream (kick, then the dot
jitter the MATLAB script draws), so save_state() / resume() can stop and
continue a run anywhere and reproduce it bit for bit.
"""
import copy
import numpy as np
import kd_engine as E

R = 1920
NSUB = 20                    # in-between moments per cycle, fine-time maps
BIN = 4                      # px per strip of amount(x)
MERGE_PX = 0.5               # cusp bookkeeping tolerance, px of CURRENT arc length
BAND_NAMES = ('red', 'green', 'blue')
_u41 = R / 41.0
BAND_EDGES = {'red': ((41 - 23.5) * _u41, (41 - 22.5) * _u41),
              'green': ((41 - 20.5) * _u41, (41 - 19.5) * _u41),
              'blue': ((41 - 17.5) * _u41, (41 - 16.5) * _u41)}      # as core.py
LOGS_EDGES = np.arange(0.0, 12.0 + 1e-9, 0.05)                         # log10 S bins (material vs S)
PROFILE_U = np.arange(-40, 41, dtype=float)                           # launch px around a cusp
PROFILE_LAGS = (0, 1, 2, 4)                                           # cycles after the cusp
PROFILE_MAX = 2500                                                    # cusps per line that get profiles


def band_rows(name):
    a, b = BAND_EDGES[name]
    return np.arange(np.ceil(a), np.floor(b) + 1.0)                    # 0-based y, as core.py


def band_centre(name):
    a, b = BAND_EDGES[name]
    return 0.5 * (a + b)


# ------------------------------------------------------------------ geometry of M
def sigma(a, b, c, d):
    return E.sigma_max(a, b, c, d)


def axis_angles(a, b, c, d):
    """alpha = angle of u1 (long axis now) from +x, beta = angle of v1 (launch
    direction stretched most) from +x; both in (-pi/2, pi/2], (x, y) data coords."""
    alpha = 0.5 * np.arctan2(2 * (a * c + b * d), (a * a + b * b) - (c * c + d * d))
    beta = 0.5 * np.arctan2(2 * (a * b + c * d), (a * a + c * c) - (b * b + d * d))
    return alpha, beta


def theta_now_deg(alpha):
    """Orientation of the long axis in the picture: 0 = upright (along the row
    axis y), positive = leaning right (top end to the right; y points down the
    picture), in (-90, 90]."""
    ad = np.degrees(alpha)
    return np.where(ad > 0, ad - 90.0, ad + 90.0)


def phi_launch_deg(beta):
    """Angle of v1 from the horizontal as seen in the picture (positive = rising
    to the right, i.e. counter-clockwise on screen), in (-90, 90]."""
    p = -np.degrees(beta)
    return np.where(p <= -90.0, p + 180.0, p)


def orient(a, b, c, d):
    al, be = axis_angles(a, b, c, d)
    return theta_now_deg(al), phi_launch_deg(be)


# ------------------------------------------------------------------ kicks
class KickStream:
    """MATLAB's kicks, cycle by cycle: kick (annealed) then the dot jitter draw."""

    def __init__(self, prm):
        self.prm = prm
        self.rs = E.matlab_stream(prm.seed)
        self.frozen = prm.K * E.randpot1d_four(prm.cf_x, self.rs) if prm.mode == 'quenched' else None
        self.cos = None
        if prm.mode == 'cosine':
            x = np.arange(prm.cols)
            self.cos = (prm.K_STD * prm.cols / (2 * np.pi * prm.Tau)) * np.sin(2 * np.pi * x / prm.cols)

    def next(self):
        p = self.prm
        if p.mode == 'annealed':
            dp = p.K * E.randpot1d_four(p.cf_x, self.rs)
        elif p.mode == 'quenched':
            dp = self.frozen
        else:
            dp = self.cos
        self.rs.random_sample(2 * p.np_dots)
        return dp


def fp_at(x, dp):
    n = dp.size
    i0 = np.mod(np.floor(x).astype(np.int64), n)
    return dp[(i0 + 1) % n] - dp[i0]


# ------------------------------------------------------------------ exact centre line
def line_new(y_row, cols=R):
    n = cols
    u = np.arange(n, dtype=float)
    return dict(u=u, L=np.ones(n), x=u.copy(), y=np.full(n, float(y_row)), a=np.ones(n),
                b=np.zeros(n), c=np.zeros(n), d=np.ones(n))


def line_step(P, dp, Tau, row0, cols=R):
    """One cycle of a piecewise-affine line (kd_engine.caustics_exact, with b and
    d carried too). Returns the new pieces and the drift record D: for every
    piece, M at the start of the drift (a0, b0, ck, dk -- the kick changes c
    and d only) and at its end (a1, b1); during the drift M(t) = [a0 + t(a1-a0),
    b0 + t(b1-b0); ck, dk] exactly."""
    n = cols
    x, y, a, b, c, d, u, L = P['x'], P['y'], P['a'], P['b'], P['c'], P['d'], P['u'], P['L']
    xR = x + a * L
    lo = np.minimum(x, xR); hi = np.maximum(x, xR)
    k0 = np.floor(lo).astype(np.int64) + 1
    k1 = np.ceil(hi).astype(np.int64) - 1
    ncut = np.maximum(k1 - k0 + 1, 0)
    tot = int(ncut.sum()) + x.size
    reps = ncut + 1
    pid = np.repeat(np.arange(x.size), reps)
    j = np.arange(tot) - np.repeat(np.cumsum(reps) - reps, reps)
    ap = a[pid]
    asc = ap >= 0
    kint = np.where(asc, k0[pid] + j - 1, k1[pid] - j + 1)
    xl = np.where(j == 0, x[pid], kint.astype(float))
    s_l = np.where(j == 0, 0.0, (xl - x[pid]) / np.where(ap != 0, ap, 1))
    last = j == reps[pid] - 1
    s_r = np.empty(tot); s_r[:-1] = s_l[1:]; s_r[last] = L[pid[last]]
    u_n = u[pid] + s_l; L_n = s_r - s_l
    y_n = y[pid] + c[pid] * s_l
    a_n = ap; b_n = b[pid]; c_n = c[pid]; d_n = d[pid]
    xm = xl + 0.5 * a_n * L_n
    i0 = np.mod(np.floor(xm).astype(np.int64), n)
    fp = dp[np.mod(i0 + 1, n)] - dp[i0]
    cell = np.floor(xm)
    y_k = y_n + dp[i0] + fp * (xl - cell)
    g = 1 + Tau * fp
    a1 = g * a_n + Tau * c_n; b1 = g * b_n + Tau * d_n
    ck = fp * a_n + c_n; dk = fp * b_n + d_n
    x1 = xl + Tau * (y_k - row0)
    keep = L_n > 0
    D = dict(u=u_n[keep], L=L_n[keep], xl=xl[keep], yk=y_k[keep], a0=a_n[keep], b0=b_n[keep],
             ck=ck[keep], dk=dk[keep], a1=a1[keep], b1=b1[keep])
    Pn = dict(u=D['u'], L=D['L'], x=x1[keep], y=D['yk'], a=D['a1'], b=D['b1'], c=D['ck'], d=D['dk'])
    return Pn, D


def drift_M(D, t, idx=None):
    s = slice(None) if idx is None else idx
    a = D['a0'][s] + t * (D['a1'][s] - D['a0'][s])
    b = D['b0'][s] + t * (D['b1'][s] - D['b0'][s])
    return a, b, D['ck'][s], D['dk'][s]


def _circ_groups(gap_next, merge):
    """Members of a closed loop with gap_next[j] = distance from member j to
    member j+1 (last -> first wraps). Groups break where the gap >= merge.
    Returns the member order (starting at a group boundary) and group starts."""
    m = gap_next.size
    brk = gap_next >= merge
    if not brk.any():
        order = np.arange(m)
        return order, np.array([0]), np.array([m])
    j0 = np.nonzero(brk)[0][-1]
    order = np.r_[np.arange(j0 + 1, m), np.arange(0, j0 + 1)]
    bo = brk[order]
    starts = np.r_[0, np.nonzero(bo[:-1])[0] + 1]
    ends = np.r_[starts[1:], m]
    return order, starts, ends


def line_cusps(D, N, Tau, row0, cols=R, merge=MERGE_PX):
    """Cusps born during cycle N's drift: runs of pieces whose a changes sign
    in this drift, lying away from every caustic present at the start of the
    drift (so the stretch ends with two caustics more), each born at the first
    t at which a piece of the run reaches a = 0. Distances (merging runs split
    by flicker, nearness to an existing caustic) are measured in px of CURRENT
    arc length. Returns rows [time, x, y, launch_u, S, theta_now, phi_launch,
    piece index]."""
    a0, a1, ck = D['a0'], D['a1'], D['ck']
    n = a0.size
    flip = (a0 > 0) != (a1 > 0)
    if not flip.any():
        return np.zeros((0, 8))
    L = D['L']
    l0 = L * np.hypot(a0, ck)                       # piece lengths at the start of the drift
    l1 = L * np.hypot(a1, ck)                       # and at its end
    pos0 = np.r_[0.0, np.cumsum(l0)]; tot0 = pos0[-1]
    pos1 = np.r_[0.0, np.cumsum(l1)]; tot1 = pos1[-1]
    s0 = a0 > 0
    kb = np.nonzero(s0 != np.roll(s0, -1))[0]       # start caustics, between kb and kb+1
    cpos = np.sort(np.mod(pos0[kb + 1], tot0)) if kb.size else np.zeros(0)
    fi = np.nonzero(flip)[0]
    gap = np.r_[pos1[fi[1:]] - pos1[fi[:-1] + 1], tot1 - pos1[fi[-1] + 1] + pos1[fi[0]]]
    order, st, en = _circ_groups(gap, merge)
    fo = fi[order]
    # run extent along the line at the start of the drift
    first = fo[st]; lastp = fo[en - 1]
    lo_ = pos0[first] - merge; hi_ = pos0[lastp + 1] + merge
    if cpos.size:
        span = np.minimum(hi_ - lo_, tot0)
        lo_m = np.mod(lo_, tot0)
        hi_m = lo_m + span
        c1 = np.searchsorted(cpos, np.minimum(hi_m, tot0), 'right') - np.searchsorted(cpos, lo_m, 'left')
        c2 = np.where(hi_m > tot0, np.searchsorted(cpos, hi_m - tot0, 'right'), 0)
        near = (c1 + c2) > 0
    else:
        near = np.zeros(st.size, bool)
    tstar = a0[fo] / (a0[fo] - a1[fo])
    tstar = np.where((tstar > 0) & (tstar <= 1), tstar, np.inf)
    gmin = np.minimum.reduceat(tstar, st)
    gid = np.repeat(np.arange(st.size), en - st)
    hit = tstar == gmin[gid]
    posn = np.where(hit, np.arange(tstar.size), tstar.size)
    first_hit = np.minimum.reduceat(posn, st)
    ok = (~near) & np.isfinite(gmin)
    if not ok.any():
        return np.zeros((0, 8))
    p = fo[first_hit[ok]]
    t = gmin[ok]
    a_, b_, c_, d_ = drift_M(D, t, p)
    S_ = sigma(a_, b_, c_, d_)
    th, ph = orient(a_, b_, c_, d_)
    x = np.mod(D['xl'][p] + t * Tau * (D['yk'][p] - row0), cols)   # a(t) = 0: the piece is vertical
    y = D['yk'][p] + c_ * D['L'][p] / 2
    u = np.mod(D['u'][p] + D['L'][p] / 2, cols)
    return np.column_stack([N - 1 + t, x, y, u, S_, th, ph, p])


def line_S_at(D, t, uq, cols=R):
    """S of the line at launch positions uq, at time t of this cycle's drift."""
    uq = np.mod(uq, cols)
    idx = np.searchsorted(D['u'], uq, side='right') - 1
    idx = np.clip(idx, 0, D['u'].size - 1)
    return sigma(*drift_M(D, t, idx))


# ------------------------------------------------------------------ band amounts, branches
def strip_hist(x, cols=R, binw=BIN, w=None):
    nb = cols // binw
    j = np.minimum((np.mod(x, cols) // binw).astype(np.int64), nb - 1)
    return np.bincount(j, weights=w, minlength=nb)[:nb].astype(float)


def branches(A, thr=2.0, minlife=1.0, nsub=NSUB, binw=BIN, cols=R):
    """Connected regions (8-connectivity, x periodic) of a fine-time map
    A[time, strip] with A >= thr, lasting at least minlife cycles. Rows:
    [id, start time, end time, start x, peak amount, mean amount]; and the
    label map."""
    from scipy import ndimage as ndi
    T, nb = A.shape
    lab, n = ndi.label(A >= thr, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros((0, 6)), np.zeros_like(lab)
    par = np.arange(n + 1)

    def root(i):
        while par[i] != i:
            par[i] = par[par[i]]; i = par[i]
        return i
    rr = np.nonzero(lab[:, -1])[0]
    for t in rr:
        for t2 in (t - 1, t, t + 1):
            if 0 <= t2 < T and lab[t2, 0]:
                r1, r2 = root(lab[t, -1]), root(lab[t2, 0])
                if r1 != r2:
                    par[max(r1, r2)] = min(r1, r2)
    roots = np.array([root(i) for i in range(n + 1)])
    lab = roots[lab]
    xc = (np.arange(nb) + 0.5) * binw
    objs = ndi.find_objects(lab)
    rows, keepids = [], []
    tt = np.arange(T) / nsub
    for k in np.unique(lab[lab > 0]):
        ii, jj = np.nonzero(lab == k)
        t0, t1 = tt[ii.min()], tt[ii.max()]
        if t1 - t0 < minlife:
            continue
        sel = ii == ii.min()
        ang = 2 * np.pi * xc[jj[sel]] / cols
        xs = np.mod(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) * cols / (2 * np.pi), cols)
        v = A[ii, jj]
        rows.append([0, t0, t1, xs, v.max(), v.mean()]); keepids.append(k)
    if not rows:
        return np.zeros((0, 6)), np.zeros_like(lab)
    rows = np.array(rows)
    o = np.argsort(rows[:, 1], kind='stable')
    rows = rows[o]; rows[:, 0] = np.arange(1, len(o) + 1)
    out = np.zeros_like(lab)
    for new, k in enumerate(np.array(keepids)[o], 1):
        out[lab == k] = new
    return rows, out


# ------------------------------------------------------------------ the run
class Stage2Run:
    """One (mode, seed) run. step() advances one cycle; save_state() returns a
    deep copy of everything; Stage2Run.resume(state, k) continues it."""

    def __init__(self, mode, seed, cycles, tracer_step=4, record=True, lines=True,
                 snap_N=(3, 5, 10, 15, 20), pick_seed=12345, n_pick=300):
        prm = E.Params(cycles=cycles, mode=mode, seed=seed)
        self.prm = prm
        self.mode, self.seed, self.cycles = mode, seed, cycles
        self.N = 0
        self.ks = KickStream(prm)
        self.kicks = []
        st = float(tracer_step)
        g = np.arange(0, R, st)
        X, Y = np.meshgrid(g, g)                       # row-major: y slow, x fast
        self.T = E.new_state(X.ravel(), Y.ravel())
        self.T2 = E.new_state(X.ravel() + st / 2, Y.ravel() + st / 2)   # shifted by half a spacing
        self.tracer_step = tracer_step
        self.B = []
        for nm in BAND_NAMES:
            Xb, Yb = np.meshgrid(np.arange(0, R, 0.5), band_rows(nm))
            self.B.append(E.new_state(Xb.ravel(), Yb.ravel()))
        self.n0 = [strip_hist(b['x']) for b in self.B]
        self.lines = lines
        self.Lp = [line_new(band_centre(nm)) for nm in BAND_NAMES] if lines else None
        self.record = record
        self.snap_N = tuple(snap_N)
        self.pick = np.sort(np.random.default_rng(pick_seed).choice(self.T['x'].size, n_pick, replace=False))
        self.rec = self._new_record() if record else None

    # ---------------------------------------------------------- recording
    def _new_record(self):
        nb = R // BIN
        r = dict(S_T=[sigma(**_m(self.T))], th_T=[np.zeros(self.T['x'].size, np.float32)],
                 S_T2=[sigma(**_m(self.T2))], lnJ_pick=[], xy_T={}, band_snap={},
                 amount=[np.ones((3, nb))], amount_low=[np.ones((3, nb))], fine=[np.ones((3, nb))],
                 branch_link=[], cusps=[[] for _ in BAND_NAMES], prof=[{} for _ in BAND_NAMES],
                 prof_req=[{} for _ in BAND_NAMES], n_cusps=[0, 0, 0], mat_hist=[[] for _ in BAND_NAMES],
                 caus_count=[[0] for _ in BAND_NAMES], pieces=[[R] for _ in BAND_NAMES])
        r['th_T'][0][:] = 0.0
        return r

    # ---------------------------------------------------------- one cycle
    def step(self):
        prm = self.prm
        Tau, row0 = prm.Tau, prm.row0
        dp = self.ks.next()
        self.kicks.append(dp)
        N = self.N + 1
        rec = self.rec
        if rec is not None:
            fpp = fp_at(self.T['x'][self.pick], dp)
            rec['lnJ_pick'].append(np.log(sigma(1 + Tau * fpp, np.full_like(fpp, Tau), fpp, np.ones_like(fpp))))
        E.step(self.T, dp, Tau, row0)
        E.step(self.T2, dp, Tau, row0)
        fine = np.zeros((NSUB, 3, R // BIN))
        amt = np.zeros((3, R // BIN)); amt_lo = np.zeros((3, R // BIN))
        for k, b in enumerate(self.B):
            xprev = b['x'].copy()
            E.step(b, dp, Tau, row0)
            for q in range(1, NSUB + 1):
                fine[q - 1, k] = strip_hist(xprev + (q / NSUB) * Tau * (b['y'] - row0)) / self.n0[k]
            amt[k] = strip_hist(b['x']) / self.n0[k]
            Sb = sigma(**_m(b))
            amt_lo[k] = strip_hist(b['x'], w=(Sb <= 2.0)) / self.n0[k]
        self.N = N
        if rec is not None:
            ST = sigma(**_m(self.T))
            th, _ = orient(**_m(self.T))
            rec['S_T'].append(ST); rec['th_T'].append(th.astype(np.float32))
            rec['S_T2'].append(sigma(**_m(self.T2)))
            rec['amount'].append(amt); rec['amount_low'].append(amt_lo)
            rec['fine'].extend(list(fine))
            if N in self.snap_N:
                rec['xy_T'][N] = (np.mod(self.T['x'], R).astype(np.float32), self.T['y'].astype(np.float32))
                snaps = []
                for k, b in enumerate(self.B):
                    th_b, ph_b = orient(**_m(b))
                    snaps.append(dict(a=b['a'].astype(np.float32), S=sigma(**_m(b)).astype(np.float32),
                                      th=th_b.astype(np.float32), ph=ph_b.astype(np.float32)))
                rec['band_snap'][N] = snaps
            # Part 2B branch link: band material over strips with amount >= 2, vs the
            # least-stretched q% of the whole frame (tracers)
            qs = np.quantile(ST, [0.01, 0.05, 0.10])
            tot = np.zeros(4)
            for k, b in enumerate(self.B):
                on = amt[k][np.minimum((np.mod(b['x'], R) // BIN).astype(np.int64), R // BIN - 1)] >= 2.0
                Sb = sigma(**_m(b))[on]
                tot += [on.sum(), (Sb <= qs[0]).sum(), (Sb <= qs[1]).sum(), (Sb <= qs[2]).sum()]
            rec['branch_link'].append(np.r_[qs, tot])
        if self.lines:
            for k in range(3):
                Pn, D = line_step(self.Lp[k], dp, Tau, row0)
                self.Lp[k] = Pn
                if rec is not None:
                    cu = line_cusps(D, N, Tau, row0)
                    rec['cusps'][k].append(cu)
                    # material of the line by S (launch px), start and end of the drift
                    lg0 = np.log10(sigma(*drift_M(D, 0.0)))
                    lg1 = np.log10(sigma(*drift_M(D, 1.0)))
                    h = 0.5 * (np.histogram(lg0, LOGS_EDGES, weights=D['L'])[0] +
                               np.histogram(lg1, LOGS_EDGES, weights=D['L'])[0])
                    rec['mat_hist'][k].append(h)
                    a1 = D['a1']
                    rec['caus_count'][k].append(int(np.count_nonzero((a1 > 0) != (np.roll(a1, -1) > 0))))
                    rec['pieces'][k].append(int(a1.size))
                    # S profiles around cusps: queue each kept cusp's requests by the
                    # cycle they fall in, then serve every request due in this cycle
                    n_before = rec['n_cusps'][k]
                    for i, row in enumerate(cu):
                        ci = n_before + i
                        if not _profile_kept(ci):
                            continue
                        rec['prof'][k][ci] = np.full((len(PROFILE_LAGS), PROFILE_U.size), np.nan, np.float32)
                        for li, lag in enumerate(PROFILE_LAGS):
                            tt = row[0] + lag
                            cyc = int(np.ceil(tt)) if tt > np.floor(tt) else int(tt)
                            rec['prof_req'][k].setdefault(cyc, []).append((ci, li, tt - (cyc - 1), row[3]))
                    rec['n_cusps'][k] += len(cu)
                    for ci, li, tf, uc in rec['prof_req'][k].pop(N, []):
                        rec['prof'][k][ci][li] = line_S_at(D, tf, uc + PROFILE_U).astype(np.float32)
        return self

    # ---------------------------------------------------------- checkpoint
    def save_state(self):
        """Every tracked parcel's x, y, a, b, c, d (tracers, shifted tracers, band
        material, centre-line pieces), the kick history, the random-number
        state and the records so far -- a deep copy."""
        st = dict(N=self.N, mode=self.mode, seed=self.seed, cycles=self.cycles,
                  rng=self.ks.rs.get_state(), frozen=self.ks.frozen, kicks=[k.copy() for k in self.kicks],
                  T=copy.deepcopy(self.T), T2=copy.deepcopy(self.T2), B=copy.deepcopy(self.B),
                  Lp=copy.deepcopy(self.Lp), rec=copy.deepcopy(self.rec), pick=self.pick.copy(),
                  tracer_step=self.tracer_step, lines=self.lines, snap_N=self.snap_N)
        return st

    @classmethod
    def resume(cls, state, k):
        """Continue a saved run for k more cycles."""
        obj = cls.__new__(cls)
        obj.prm = E.Params(cycles=state['cycles'], mode=state['mode'], seed=state['seed'])
        obj.mode, obj.seed, obj.cycles = state['mode'], state['seed'], state['cycles']
        obj.N = state['N']
        obj.ks = KickStream.__new__(KickStream)
        obj.ks.prm = obj.prm
        obj.ks.rs = np.random.RandomState()
        obj.ks.rs.set_state(state['rng'])
        obj.ks.frozen = state['frozen']
        obj.ks.cos = None
        if obj.mode == 'cosine':
            x = np.arange(obj.prm.cols)
            obj.ks.cos = (obj.prm.K_STD * obj.prm.cols / (2 * np.pi * obj.prm.Tau)) * np.sin(2 * np.pi * x / obj.prm.cols)
        obj.kicks = [kk.copy() for kk in state['kicks']]
        obj.T = copy.deepcopy(state['T']); obj.T2 = copy.deepcopy(state['T2'])
        obj.B = copy.deepcopy(state['B']); obj.Lp = copy.deepcopy(state['Lp'])
        obj.lines = state['lines']; obj.snap_N = state['snap_N']
        obj.tracer_step = state['tracer_step']
        obj.n0 = [strip_hist(b0) for b0 in _band_launch_x()]
        obj.rec = copy.deepcopy(state['rec']); obj.record = obj.rec is not None
        obj.pick = state['pick'].copy()
        for _ in range(k):
            obj.step()
        return obj

    def run(self, k):
        for _ in range(k):
            self.step()
        return self


def _band_launch_x():
    out = []
    for nm in BAND_NAMES:
        Xb, _ = np.meshgrid(np.arange(0, R, 0.5), band_rows(nm))
        out.append(Xb.ravel())
    return out


def _profile_kept(i):
    """Which cusps of a line get an S profile: every one up to PROFILE_MAX, then
    an evenly thinned sample (deterministic, no random numbers)."""
    if i < PROFILE_MAX:
        return True
    return (i * 2654435761) % 1000 < 200


def _m(s):
    return dict(a=s['a'], b=s['b'], c=s['c'], d=s['d'])
