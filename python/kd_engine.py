"""
kd_engine -- kick-drift dynamics and analyses, matching the MATLAB script.

Conventions (identical to "Kick-drift LSVD"):
  one step     Yc = y + f(x);  Xc = x + Tau*(Yc - row0)
  f            piecewise-linear interpolant of dp_x (period cols in x)
  M            [a b; c d], rows (x, y); per-step J = [1+Tau*fp, Tau; fp, 1]
  S            sigma_max(M) = (sqrt(t+2) + sqrt(t-2))/2, t = a^2+b^2+c^2+d^2
  S_THR        2

Randomness reproduces MATLAB exactly: rng(seed,'twister') is
np.random.RandomState(5489) for seed 0 and RandomState(seed) otherwise (MT19937
init_genrand). The stream is consumed in the script's order -- per cycle the
kick (959 uniforms, annealed only) then the dot-lattice jitter (2 x 25600) --
so seed 0 gives the MATLAB run's kicks bit-for-bit.

Geometry. The map is exactly periodic in y with period P = cols/Tau (a shift
of P moves the next drift by exactly cols). One forward sweep over a lattice
seeded on a FULL period therefore yields
  * the MATLAB tracer lattice      : the subset seeded in the window [0, rows)
  * the Eulerian window distribution: parcels whose y mod P lands in the
    window at step N -- an exact, uniform sample of what Sfield_snapshot
    measures, in O(N) instead of O(N^2)
  * the full phase-space distribution: all parcels.
"""
import numpy as np


# ============================================================ parameters
class Params:
    def __init__(self, cycles=60, mode='annealed', seed=0, K_STD=0.50,
                 S_THR=2.0, lag_step=2, line_rows=None, line_dens=8,
                 line_bin=4, band_bin=4, caustic_dens=(64, 512), exact_max_pieces=1.2e7):
        self.rows = self.cols = 1920
        self.lc, self.K0, self.Tau0, self.Lx, self.Lp = 0.05, 0.025, 0.5, 1.0, 1.0
        self.K = self.K0 / self.Lp * self.rows
        self.Tau = (self.Tau0 * self.Lp * self.cols) / (self.Lx * self.rows)
        self.row0 = self.rows / 2.0
        self.P = self.cols / self.Tau                  # exact momentum period in y
        self.cycles, self.mode, self.seed = cycles, mode, seed
        self.K_STD, self.S_THR = K_STD, S_THR
        self.lag_step = lag_step
        self.line_rows = (line_rows if line_rows is not None else
                          [self.row0 - self.rows / 8, self.row0, self.row0 + self.rows / 8])
        self.line_dens, self.line_bin, self.band_bin = line_dens, line_bin, band_bin
        self.caustic_dens = caustic_dens
        self.exact_max_pieces = exact_max_pieces
        self.np_dots = len(range(6, self.rows, 12)) * len(range(6, self.cols, 12))
        xg = np.linspace(-0.5, 0.5, self.cols)
        c_x = np.exp(-(xg ** 2) / self.lc ** 2)
        self.cf_x = np.sqrt(np.abs(np.fft.fft(np.fft.ifftshift(c_x))))


S_LIST = np.array([1.25, 1.5, 2, 3, 5, 10, 20, 50, 100, 1e3, 1e4, 1e6])
LOG10_EDGES = np.arange(0, 13.0 + 1e-9, 0.05)          # 260 bins + overflow
QUANTS = np.array([0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
BAND_EDGES = np.array([2, 5, 10, 100, 1e4])             # S<=2, 2-5, 5-10, 10-100, 100-1e4, >1e4
METRIC_F = np.array([0.5, 2.0, 4.0])
DEPTH_EDGES = np.array([2.0, 2.5, 3.0, 5.0, 10.0, 20.0, 100.0, np.inf])
DEPTH_LABELS = ['(2,2.5]', '(2.5,3]', '(3,5]', '(5,10]', '(10,20]', '(20,100]', '>100']
DEPTH_USER = [[0], [1], [2], [3, 4, 5, 6]]              # the four requested bins: ..., >5
DEPTH_USER_LABELS = ['(2,2.5]', '(2.5,3]', '(3,5]', '>5']
AGE_BINS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 40), (41, 10 ** 6)]


# ============================================================ randomness
def matlab_stream(seed):
    """rng(seed,'twister') in MATLAB."""
    return np.random.RandomState(5489 if seed == 0 else seed)


def randpot1d_four(cf, rs):
    """MATLAB randpot1d_four, drawing from stream rs."""
    n = cf.size
    p = (2 * rs.random_sample(n // 2 - 1) - 1) * np.pi
    phase = np.zeros(n)
    phase[1:n // 2] = p
    phase[n // 2 + 1:] = -p[::-1]
    cff = cf * np.exp(1j * phase)
    cff[0] = 0
    u = np.real(np.sqrt(n) * np.fft.ifft(cff))
    return u / max(u.std(ddof=1), 1e-12)


def make_kicks(prm):
    """The kicks a MATLAB run with randomSeed = prm.seed applies, per cycle."""
    rs = matlab_stream(prm.seed)
    if prm.mode == 'cosine':
        x = np.arange(prm.cols)
        dc = (prm.K_STD * prm.cols / (2 * np.pi * prm.Tau)) * np.sin(2 * np.pi * x / prm.cols)
        return [dc] * prm.cycles
    frozen = prm.K * randpot1d_four(prm.cf_x, rs) if prm.mode == 'quenched' else None
    out = []
    for _ in range(prm.cycles):
        out.append(prm.K * randpot1d_four(prm.cf_x, rs) if prm.mode == 'annealed' else frozen)
        rs.random_sample(2 * prm.np_dots)                # Smap_dots_snapshot jitter
    return out


# ============================================================ dynamics
def sigma_max(a, b, c, d):
    t = a * a + b * b + c * c + d * d
    return (np.sqrt(t + 2) + np.sqrt(np.maximum(t - 2, 0))) / 2


def step(st, dp, Tau, row0):
    """In-place step of a state dict with x, y, a, b, c, d (monodromy first)."""
    x = st['x']
    n = dp.size
    xf = np.floor(x)
    i0 = np.mod(xf.astype(np.int64), n)
    d0 = dp[i0]
    fp = dp[np.mod(i0 + 1, n)] - d0
    g = 1 + Tau * fp
    a, b, c, d = st['a'], st['b'], st['c'], st['d']
    a2 = g * a + Tau * c
    b2 = g * b + Tau * d
    st['c'] = fp * a + c
    st['d'] = fp * b + d
    st['a'] = a2
    st['b'] = b2
    st['y'] = st['y'] + d0 + fp * (x - xf)
    st['x'] = x + Tau * (st['y'] - row0)


def new_state(x, y):
    x = np.asarray(x, float).ravel().copy()
    y = np.asarray(y, float).ravel().copy()
    return dict(x=x, y=y, a=np.ones_like(x), b=np.zeros_like(x),
                c=np.zeros_like(x), d=np.ones_like(x))


# ============================================================ helpers
def pearson(u, v):
    u = u - u.mean(); v = v - v.mean()
    den = np.sqrt((u * u).sum() * (v * v).sum())
    return float((u * v).sum() / den) if den > 0 else np.nan


def spearman(u, v):
    """Spearman's rho with tied values given their average rank (many bins of
    the low-stretch fraction are exactly 0)."""
    from scipy.stats import rankdata
    return pearson(rankdata(u), rankdata(v))


def bincount2(i, j, w, shape):
    flat = np.bincount(i * shape[1] + j, weights=w, minlength=shape[0] * shape[1])
    return flat.reshape(shape)


def moments(lnS):
    m = lnS.mean()
    z = lnS - m
    v = (z * z).mean()
    sk = (z ** 3).mean() / v ** 1.5 if v > 0 else np.nan
    ku = (z ** 4).mean() / v ** 2 - 3 if v > 0 else np.nan
    return m, v, sk, ku


def rho_from_a(xs, a, bin_w, cols, max_span_bins=2):
    """Ray density at bin centres from the exact Jacobian, rho = sum 1/|a|.

    xs: unwrapped x_N along the line in launch order (cyclic). Every segment
    between neighbours is followed, and each bin centre it crosses gets 1/|a|
    linearly interpolated at the crossing. A segment longer than
    max_span_bins bins is not resolved (the line folds between the samples):
    the bins it touches are flagged unresolved and returned as NaN.
    """
    nb = int(round(cols / bin_w))
    x0, x1 = xs, np.roll(xs, -1)
    a0, a1 = a, np.roll(a, -1)
    x1 = x1 + np.round((x0 - x1) / cols) * cols          # nearest image on the ring
    lo, hi = np.minimum(x0, x1), np.maximum(x0, x1)
    c0 = bin_w / 2.0
    m_lo = np.floor((lo - c0) / bin_w).astype(np.int64) + 1
    m_hi = np.floor((hi - c0) / bin_w).astype(np.int64)
    ncross = np.maximum(m_hi - m_lo + 1, 0)
    rho = np.zeros(nb)
    long_seg = (hi - lo) > max_span_bins * bin_w
    bad = _mark_spans(m_lo[long_seg] - 1, m_hi[long_seg] + 1, nb)   # unresolved: poison span
    ok = (~long_seg) & (ncross > 0)
    sel = np.nonzero(ok)[0]
    reps = ncross[sel]
    seg = np.repeat(sel, reps)
    offs = np.arange(reps.sum()) - np.repeat(np.cumsum(reps) - reps, reps)
    m = m_lo[seg] + offs
    xc = c0 + m * bin_w
    dx = x1[seg] - x0[seg]
    tt = np.where(dx != 0, (xc - x0[seg]) / np.where(dx != 0, dx, 1), 0.5)
    ac = a0[seg] + tt * (a1[seg] - a0[seg])
    np.add.at(rho, np.mod(m, nb), 1.0 / np.maximum(np.abs(ac), 1e-300))
    # bins whose segment contains a caustic (a changes sign inside it)
    cs = (a0 * a1 < 0) & ~long_seg
    caus = _mark_spans(m_lo[cs] - 1, m_hi[cs] + 1, nb)
    rho[bad] = np.nan
    return rho, caus


def _mark_spans(lo, hi, nb):
    """Boolean mask of ring bins covered by any integer interval [lo, hi]."""
    mark = np.zeros(nb, bool)
    if lo.size == 0:
        return mark
    full = (hi - lo + 1) >= nb
    if full.any():
        mark[:] = True
        return mark
    lo_m = np.mod(lo, nb); ln = hi - lo + 1
    diff = np.zeros(2 * nb + 1)
    np.add.at(diff, lo_m, 1); np.add.at(diff, lo_m + ln, -1)
    cov = np.cumsum(diff)[:2 * nb]
    return (cov[:nb] + cov[nb:]) > 0


def r_of(st):
    return 0.5 * np.log(st['a'] ** 2 + st['c'] ** 2)


# ============================================================ one run
def run(prm, keep_snapshots=False, snapshot_N=(5, 10, 20, 40, 60), verbose=False, kicks=None):
    """All per-run analyses of Parts 1, 2 and 4. Returns a dict of arrays.

    kicks: optional list of dp profiles to use instead of make_kicks(prm), e.g.
    a MATLAB run's dp_x_hist, so both implementations see identical kicks."""
    C, cols, rows, Tau, row0, thr = prm.cycles, prm.cols, prm.rows, prm.Tau, prm.row0, prm.S_THR
    kicks = make_kicks(prm) if kicks is None else [np.asarray(k, float).ravel() for k in kicks]
    out = dict(mode=prm.mode, seed=prm.seed, cycles=C)

    # -------- full-period lattice; the window-seeded subset is the MATLAB tracer lattice
    ys = np.arange(row0 - prm.P / 2, row0 + prm.P / 2, prm.lag_step)
    xs = np.arange(0, cols, prm.lag_step, dtype=float)
    YY, XX = np.meshgrid(ys, xs, indexing='ij')
    st = new_state(XX, YY)
    inwin0 = (YY.ravel() >= 0) & (YY.ravel() < rows)              # MATLAB tracer subset
    trk = np.nonzero(inwin0)[0]
    nT = trk.size
    ybase = row0 - prm.P / 2

    # -------- launched lines (Part 1A)
    L = dict()
    x_line = np.arange(cols * prm.line_dens) / prm.line_dens
    lines = [new_state(x_line, np.full_like(x_line, yr)) for yr in prm.line_rows]
    nrow = len(lines)
    cdens = list(prm.caustic_dens)
    cl_lines = [new_state(np.arange(cols * D) / D, np.full(cols * D, prm.row0)) for D in cdens]
    nbL = int(round(cols / prm.line_bin))

    # -------- storage
    S_EDGES = LOG10_EDGES
    nbin = S_EDGES.size - 1
    out['hist_log10S'] = np.zeros((C + 1, nbin)); out['overflow'] = np.zeros(C + 1)
    out['P_le'] = np.full((C + 1, S_LIST.size), np.nan)
    out['P_le_full'] = np.full((C + 1, S_LIST.size), np.nan)     # whole period
    out['n_win'] = np.zeros(C + 1)
    out['mom'] = np.full((C + 1, 4), np.nan)
    out['quant'] = np.full((C + 1, QUANTS.size), np.nan)
    out['bandshare'] = np.full((C + 1, BAND_EDGES.size + 1), np.nan)
    out['metric_P_le'] = np.full((METRIC_F.size, C + 1, S_LIST.size), np.nan)
    out['metric_hist'] = np.zeros((METRIC_F.size, C + 1, nbin))
    out['metric_mean_lnS'] = np.full((METRIC_F.size, C + 1), np.nan)

    nbB = int(round(cols / prm.band_bin))
    out['band_x'] = np.zeros((C + 1, 3, nbB))

    out['line_rho_hist'] = np.zeros((nrow, C + 1, nbL))
    out['line_rho_a'] = np.full((nrow, C + 1, nbL), np.nan)
    out['line_caus_bin'] = np.zeros((nrow, C + 1, nbL), bool)
    out['line_caustics'] = np.zeros((nrow, C + 1))
    out['line_mean_r'] = np.zeros((nrow, C + 1)); out['line_var_r'] = np.zeros((nrow, C + 1))
    out['line_mean_lnS'] = np.zeros((nrow, C + 1))
    out['line_corr_low'] = np.full((nrow, C + 1, 2), np.nan)       # pearson, spearman
    out['line_corr_negr'] = np.full((nrow, C + 1, 2), np.nan)
    out['caustic_dens'] = np.array(cdens)
    out['caustics_dense'] = np.zeros((len(cdens), C + 1))
    out['dense_mean_r'] = np.zeros((len(cdens), C + 1)); out['dense_var_r'] = np.zeros((len(cdens), C + 1))
    chk_i = int(np.argmin(np.abs(np.array(cdens) - 512)))      # rho check on the ~512/px line
    out['chk_dens'] = cdens[chk_i]
    out['chk_rho_hist'] = np.zeros((C + 1, nbL)); out['chk_count'] = np.zeros((C + 1, nbL))
    out['chk_rho_a'] = np.full((C + 1, nbL), np.nan); out['chk_caus_bin'] = np.zeros((C + 1, nbL), bool)
    RB = np.arange(-6, 30.0001, 0.25); SB = np.arange(0, 30.0001, 0.25)
    out['joint_r_lnS'] = {}
    out['cond_r'] = {}
    snaps = {}

    # -------- Part 2 state (window-seeded subset)
    low_prev = np.ones(nT, bool)                 # S_0 = 1: everyone starts in the set
    S_prev = np.ones(nT)
    exit_cyc = np.full(nT, -1, np.int32)
    entry_cyc = np.zeros(nT, np.int32)           # 0 = original cohort
    n_ret = np.zeros(nT, np.int32)
    risk_re = np.zeros((C + 1, C + 1)); ev_re = np.zeros((C + 1, C + 1))     # [exit cycle, age]
    risk_ex = np.zeros((C + 1, C + 1)); ev_ex = np.zeros((C + 1, C + 1))     # [entry cycle, age]
    nd = DEPTH_EDGES.size - 1
    risk_rda = np.zeros((2, nd, C + 1)); ev_rda = np.zeros((2, nd, C + 1))   # [returned before?, depth S_{N-1}, age]
    n_ent = np.zeros(C + 1); n_ext = np.zeros(C + 1)
    type1_num = 0.0; type1_den = 0.0
    ext_prev = np.zeros(nT, bool)
    haz_f = np.zeros(C + 1); haz_f[0] = 1.0

    def part4_row(N, S, a, b, c, d, win):
        Sw = S[win]
        lnS = np.log(Sw)
        l10 = lnS / np.log(10)
        out['n_win'][N] = Sw.size
        h, _ = np.histogram(l10, bins=S_EDGES)
        out['hist_log10S'][N] = h
        out['overflow'][N] = np.count_nonzero(l10 >= S_EDGES[-1])
        out['P_le'][N] = (Sw[:, None] <= S_LIST[None, :]).mean(0)
        out['P_le_full'][N] = (S[:, None] <= S_LIST[None, :]).mean(0)
        out['mom'][N] = moments(lnS)
        out['quant'][N] = np.quantile(lnS, QUANTS)
        idx = np.searchsorted(BAND_EDGES, Sw, side='left')        # S<=2 -> 0, etc.
        out['bandshare'][N] = np.bincount(idx, minlength=BAND_EDGES.size + 1) / Sw.size
        aw, bw, cw, dw = a[win], b[win], c[win], d[win]
        for i, fct in enumerate(METRIC_F):                          # D M D^-1, D = diag(1,f)
            Sf = sigma_max(aw, bw / fct, cw * fct, dw)
            out['metric_P_le'][i, N] = (Sf[:, None] <= S_LIST[None, :]).mean(0)
            out['metric_hist'][i, N], _ = np.histogram(np.log10(Sf), bins=S_EDGES)
            out['metric_mean_lnS'][i, N] = np.log(Sf).mean()

    # N = 0 rows
    S0 = sigma_max(st['a'], st['b'], st['c'], st['d'])
    part4_row(0, S0, st['a'], st['b'], st['c'], st['d'], inwin0)

    for N in range(0, C + 1):
        if N > 0:
            dp = kicks[N - 1]
            step(st, dp, Tau, row0)
            for ln in lines:
                step(ln, dp, Tau, row0)
            for ln in cl_lines:
                step(ln, dp, Tau, row0)
            S = sigma_max(st['a'], st['b'], st['c'], st['d'])
            yfd = ybase + np.mod(st['y'] - ybase, prm.P)
            win = (yfd >= 0) & (yfd < rows)
            part4_row(N, S, st['a'], st['b'], st['c'], st['d'], win)

            # ---------------- Part 2 + 1B on the tracer subset
            St = S[trk]
            low = St <= thr
            haz_f[N] = low.mean()
            ent = low & ~low_prev
            ext = low_prev & ~low
            n_ent[N] = ent.sum(); n_ext[N] = ext.sum()
            # Type 1: left at N-1, back at N -- counted straight from the masks
            type1_num += np.count_nonzero(ext_prev & low)
            type1_den += np.count_nonzero(ext_prev)
            # re-entry risk: out at N-1
            o = ~low_prev
            if o.any():
                e = exit_cyc[o]; k = N - e; ev = low[o].astype(float)
                risk_re += bincount2(e, k, None, risk_re.shape)
                ev_re += bincount2(e, k, ev, ev_re.shape)
                db = np.clip(np.searchsorted(DEPTH_EDGES, S_prev[o], side='left') - 1, 0, nd - 1)
                rep = (n_ret[o] >= 1).astype(np.int64)
                flat = (rep * nd + db) * (C + 1) + k
                risk_rda += np.bincount(flat, minlength=risk_rda.size).reshape(risk_rda.shape)
                ev_rda += np.bincount(flat, weights=ev, minlength=ev_rda.size).reshape(ev_rda.shape)
            # exit risk: low at N-1
            if low_prev.any():
                g = entry_cyc[low_prev]; k = N - g; ev = (~low[low_prev]).astype(float)
                risk_ex += bincount2(g, k, None, risk_ex.shape)
                ev_ex += bincount2(g, k, ev, ev_ex.shape)
            exit_cyc[ext] = N
            entry_cyc[ent] = N
            n_ret[ent] += 1
            ext_prev = ext
            low_prev = low
            S_prev = St
        else:
            St = np.ones(nT)
        # 1B: band histograms of CURRENT x of the tracer parcels
        xm = np.mod(st['x'][trk], cols)
        band = np.where(St <= 2, 0, np.where(St <= 5, 1, 2))
        xb = np.minimum((xm / prm.band_bin).astype(np.int64), nbB - 1)
        out['band_x'][N] = np.bincount(band * nbB + xb, minlength=3 * nbB).reshape(3, nbB)

        # ---------------- Part 1A: lines
        for ri, ln in enumerate(lines):
            a, c = ln['a'], ln['c']
            Sl = sigma_max(ln['a'], ln['b'], ln['c'], ln['d'])
            r = r_of(ln)
            xm = np.mod(ln['x'], cols)
            xb = np.minimum((xm / prm.line_bin).astype(np.int64), nbL - 1)
            cnt = np.bincount(xb, minlength=nbL)
            out['line_rho_hist'][ri, N] = cnt / (prm.line_dens * prm.line_bin)
            rho_a, cb = rho_from_a(ln['x'], a, prm.line_bin, cols)
            out['line_rho_a'][ri, N] = rho_a
            out['line_caus_bin'][ri, N] = cb
            out['line_caustics'][ri, N] = np.count_nonzero(a * np.roll(a, -1) < 0)
            out['line_mean_r'][ri, N] = r.mean(); out['line_var_r'][ri, N] = r.var()
            out['line_mean_lnS'][ri, N] = np.log(Sl).mean()
            ok = cnt >= 5
            if N > 0 and ok.sum() > 10:
                flow = np.bincount(xb, weights=(Sl <= thr).astype(float), minlength=nbL)[ok] / cnt[ok]
                negr = np.bincount(xb, weights=-r, minlength=nbL)[ok] / cnt[ok]
                rh = out['line_rho_hist'][ri, N][ok]
                if flow.std() > 0:
                    out['line_corr_low'][ri, N] = (pearson(rh, flow), spearman(rh, flow))
                out['line_corr_negr'][ri, N] = (pearson(rh, negr), spearman(rh, negr))
            if keep_snapshots and N in snapshot_N:
                snaps[(ri, N)] = (xm.astype(np.float32), ln['y'].astype(np.float32),
                                  Sl.astype(np.float32), a.astype(np.float32))
        for di, ln in enumerate(cl_lines):
            a = ln['a']
            out['caustics_dense'][di, N] = np.count_nonzero(a * np.roll(a, -1) < 0)
            r = r_of(ln)
            out['dense_mean_r'][di, N] = r.mean(); out['dense_var_r'][di, N] = r.var()
            if di == chk_i:       # rho check line: exact 1/|a| vs histogram, same bins
                D = cdens[di]
                xb = np.minimum((np.mod(ln['x'], cols) / prm.line_bin).astype(np.int64), nbL - 1)
                cnt = np.bincount(xb, minlength=nbL)
                out['chk_rho_hist'][N] = cnt / (D * prm.line_bin)
                out['chk_count'][N] = cnt
                out['chk_rho_a'][N], out['chk_caus_bin'][N] = rho_from_a(ln['x'], a, prm.line_bin, cols)
        # ---------------- 4F: joint (r, ln S) on the lines
        if N in (10, 30, 60):
            r_all = np.concatenate([r_of(ln) for ln in lines])
            lS = np.concatenate([np.log(sigma_max(ln['a'], ln['b'], ln['c'], ln['d'])) for ln in lines])
            H, _, _ = np.histogram2d(r_all, lS, bins=[RB, SB])
            out['joint_r_lnS'][N] = H
            lowS = lS <= np.log(thr)
            lowr = r_all <= np.log(thr)
            q10 = np.quantile(r_all, 0.10)
            lowr10 = r_all <= q10
            out['cond_r'][N] = dict(
                P_lowS_given_lowr=lowS[lowr].mean() if lowr.any() else np.nan,
                P_lowr_given_lowS=lowr[lowS].mean() if lowS.any() else np.nan,
                P_lowS_given_r10=lowS[lowr10].mean(),
                P_r10_given_lowS=lowr10[lowS].mean() if lowS.any() else np.nan,
                r_q10=q10, P_lowS=lowS.mean(), violations_S_lt_er=int(np.count_nonzero(lS < r_all - 1e-9)))
        if verbose and N % 10 == 0:
            print(f'  [{prm.mode} s{prm.seed}] N={N}', flush=True)

    # exact caustic counts on every launched row (piecewise-affine line)
    out['caustics_exact'] = np.full((nrow, C + 1), np.nan)
    out['exact_pieces'] = np.zeros((nrow, C + 1))
    for ri, yr in enumerate(prm.line_rows):
        out['caustics_exact'][ri], out['exact_pieces'][ri] = caustics_exact(
            prm, kicks, yr, C, max_pieces=prm.exact_max_pieces)
    out.update(risk_re=risk_re, ev_re=ev_re, risk_ex=risk_ex, ev_ex=ev_ex,
               risk_rda=risk_rda, ev_rda=ev_rda,
               n_ent=n_ent, n_ext=n_ext, type1=(type1_num, type1_den), haz_f=haz_f,
               n_tracers=nT, kicks_first=np.asarray(kicks[0]))
    if keep_snapshots:
        out['snaps'] = snaps
    return out


# ============================================================ hazards / KM
def hazard(ev, risk, min_risk=0):
    h = np.where(risk > 0, ev / np.maximum(risk, 1), np.nan)
    if min_risk:
        h = np.where(risk >= min_risk, h, np.nan)
    return h


def km(h_age):
    """F(k) = 1 - prod_{j<=k} (1 - h(j)), ages from 1. NaN hazards stop the curve."""
    hh = np.asarray(h_age, float)
    F = np.full(hh.size, np.nan)
    surv = 1.0
    for k in range(1, hh.size):
        if not np.isfinite(hh[k]):
            break
        surv *= (1 - hh[k])
        F[k] = 1 - surv
    return F


def pooled_bins(ev_age, risk_age):
    res = []
    for lo, hi in AGE_BINS:
        sl = slice(lo, min(hi, ev_age.size - 1) + 1)
        r = risk_age[sl].sum(); e = ev_age[sl].sum()
        res.append((e / r if r > 0 else np.nan, r))
    return res


# ============================================================ i.i.d. model (4D)
def fp_pool(prm, n_kicks=500, seed=12345):
    rs = np.random.RandomState(seed)
    pool = []
    for _ in range(n_kicks):
        dp = prm.K * randpot1d_four(prm.cf_x, rs)
        pool.append(np.roll(dp, -1) - dp)
    return np.concatenate(pool)


def iid_model(prm, n_samples=10 ** 7, Nmax=60, chunk=10 ** 6, seed=777, pool=None,
              Nkeep=(20, 40, 60), qs=np.linspace(-4, 4, 33)):
    """Exact i.i.d. 2x2 product for the annealed kick, sampled per parcel."""
    if pool is None:
        pool = fp_pool(prm)
    rs = np.random.RandomState(seed)
    Tau = prm.Tau
    P_le = np.zeros((Nmax + 1, S_LIST.size))
    mom = np.zeros((Nmax + 1, 2))        # sum lnS, sum lnS^2
    lam_hist = {N: np.zeros(400) for N in Nkeep}
    lam_edges = np.linspace(0, 1.0, 401)
    # log E[S^q]: running logsumexp per q, per N in Nkeep
    lse = {N: np.full(qs.size, -np.inf) for N in Nkeep}
    top_share = {N: np.zeros(qs.size) for N in Nkeep}
    nch = int(np.ceil(n_samples / chunk))
    P_chunks = np.zeros((nch, Nmax + 1, S_LIST.size))
    done = 0
    for ci in range(nch):
        m = min(chunk, n_samples - done)
        a = np.ones(m); b = np.zeros(m); c = np.zeros(m); d = np.ones(m)
        P_chunks[ci, 0] = 1.0 * m
        for N in range(1, Nmax + 1):
            fp = pool[rs.randint(0, pool.size, m)]
            g = 1 + Tau * fp
            a2 = g * a + Tau * c; b2 = g * b + Tau * d
            c = fp * a + c; d = fp * b + d
            a = a2; b = b2
            S = sigma_max(a, b, c, d)
            P_chunks[ci, N] = (S[:, None] <= S_LIST[None, :]).sum(0)
            lnS = np.log(S)
            mom[N, 0] += lnS.sum(); mom[N, 1] += (lnS ** 2).sum()
            if N in Nkeep:
                lam_hist[N] += np.histogram(lnS / N, bins=lam_edges)[0]
                for qi, q in enumerate(qs):
                    z = q * lnS
                    zm = z.max()
                    s = np.log(np.exp(z - zm).sum()) + zm
                    lse[N][qi] = np.logaddexp(lse[N][qi], s)
        done += m
    P_le = P_chunks.sum(0) / n_samples
    P_chunk_frac = P_chunks / chunk
    meanlnS = mom[:, 0] / n_samples
    varlnS = mom[:, 1] / n_samples - meanlnS ** 2
    Lq = {N: (lse[N] - np.log(n_samples)) / N for N in Nkeep}
    return dict(P_le=P_le, P_le_sd_chunk=P_chunk_frac.std(0, ddof=1), n_chunks=nch,
                meanlnS=meanlnS, varlnS=varlnS, lam_hist=lam_hist, lam_edges=lam_edges,
                Lq=Lq, qs=qs, n_samples=n_samples)


def transfer_operator_L(prm, pool, qs, n_theta=2048, n_fp=301):
    """Generalised Lyapunov exponent of the VECTOR norm, L(q) = ln rho(T_q).

    T_q f(theta) = E_fp[ |J u_theta|^q f(theta') ], theta in [0, pi). Exact for
    sigma_max when q >= 0 (all norms share the exponential rate); for q < 0
    sigma_max >= 1 makes its L(q) saturate, which this operator does not see.
    """
    Tau = prm.Tau
    w, e = np.histogram(pool, bins=n_fp)
    fpv = 0.5 * (e[:-1] + e[1:]); w = w / w.sum()
    keep = w > 0; fpv, w = fpv[keep], w[keep]
    th = (np.arange(n_theta) + 0.5) * np.pi / n_theta
    ux, uy = np.cos(th), np.sin(th)
    rows_i, cols_j, vals_n, ws = [], [], [], []
    for f, wf in zip(fpv, w):
        vx = (1 + Tau * f) * ux + Tau * uy
        vy = f * ux + uy
        nrm = np.hypot(vx, vy)
        tp = np.mod(np.arctan2(vy, vx), np.pi)
        pos = tp / (np.pi / n_theta) - 0.5
        j0 = np.floor(pos).astype(int); fr = pos - j0
        rows_i.append(np.arange(n_theta)); cols_j.append(np.mod(j0, n_theta))
        vals_n.append((nrm, 1 - fr, wf)); ws.append(None)
        rows_i.append(np.arange(n_theta)); cols_j.append(np.mod(j0 + 1, n_theta))
        vals_n.append((nrm, fr, wf)); ws.append(None)
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import eigs
    I = np.concatenate(rows_i); Jc = np.concatenate(cols_j)
    NRM = np.concatenate([v[0] for v in vals_n])
    LIN = np.concatenate([v[1] for v in vals_n])
    WF = np.concatenate([np.full(v[0].size, v[2]) for v in vals_n])
    Ls = []
    for q in qs:
        T = csr_matrix((WF * LIN * NRM ** q, (I, Jc)), shape=(n_theta, n_theta))
        val = eigs(T, k=1, which='LM', return_eigenvectors=False)
        Ls.append(np.log(np.abs(val[0])))
    return np.array(Ls)


# ============================================================ exact caustic count (1A)
def caustics_exact(prm, kicks, y_row, Nmax, max_pieces=3e7):
    """Exact caustic count on a launched line, N = 0..Nmax.

    The kick f is the piecewise-linear interpolant of dp, so one step is affine
    inside every unit cell [i, i+1) of x and the N-step map is piecewise affine
    along the line. The tangent component a = dx_N/dx0 is therefore constant on
    each affine piece, and the line is tracked exactly as a list of pieces
    (left end in launch and phase space, a, c). Every step a piece is cut where
    its image crosses an integer x; caustics are the sign changes of a between
    consecutive pieces (each is a turning point of x_N(x0); in this model the
    fold is a corner, so rho jumps there instead of diverging as d^-1/2).
    Returns counts (NaN once the piece budget is exceeded) and piece numbers.
    """
    cols, Tau, row0 = prm.cols, prm.Tau, prm.row0
    n = cols
    # N = 0: one piece per cell, a = 1, c = 0
    u = np.arange(n, dtype=float); L = np.ones(n)            # launch left end, launch length
    x = u.copy(); y = np.full(n, float(y_row)); a = np.ones(n); c = np.zeros(n)
    counts = np.full(Nmax + 1, np.nan); npc = np.zeros(Nmax + 1)
    counts[0] = 0; npc[0] = n
    for N in range(1, Nmax + 1):
        dp = kicks[N - 1]
        # cut every piece where its image crosses an integer x
        xR = x + a * L
        lo = np.minimum(x, xR); hi = np.maximum(x, xR)
        k0 = np.floor(lo).astype(np.int64) + 1                # first integer strictly inside
        k1 = np.ceil(hi).astype(np.int64) - 1                 # last integer strictly inside
        ncut = np.maximum(k1 - k0 + 1, 0)
        tot = int(ncut.sum()) + x.size
        if tot > max_pieces:
            return counts, npc
        reps = ncut + 1
        pid = np.repeat(np.arange(x.size), reps)
        j = np.arange(tot) - np.repeat(np.cumsum(reps) - reps, reps)   # sub-index in piece
        # breakpoints in image x, ordered along the launch direction
        asc = a[pid] >= 0
        kint = np.where(asc, k0[pid] + j - 1, k1[pid] - j + 1)        # integer at the sub-piece's left end
        xl = np.where(j == 0, x[pid], kint.astype(float))
        s_l = np.where(j == 0, 0.0, (xl - x[pid]) / np.where(a[pid] != 0, a[pid], 1))
        # right end of each sub-piece in launch distance
        last = j == reps[pid] - 1
        s_r = np.empty(tot); s_r[:-1] = s_l[1:]; s_r[last] = L[pid[last]]
        u_n = u[pid] + s_l; L_n = s_r - s_l
        y_n = y[pid] + c[pid] * s_l
        a_n = a[pid]; c_n = c[pid]
        # cell of each sub-piece from its midpoint
        xm = xl + 0.5 * a_n * L_n
        i0 = np.mod(np.floor(xm).astype(np.int64), n)
        fp = dp[np.mod(i0 + 1, n)] - dp[i0]
        cell = np.floor(xm)
        fx = dp[i0] + fp * (xl - cell)
        y_n = y_n + fx
        x_n = xl + Tau * (y_n - row0)
        a2 = (1 + Tau * fp) * a_n + Tau * c_n
        c2 = fp * a_n + c_n
        # drop zero-length slivers from exact-integer endpoints
        keep = L_n > 0
        u, L, x, y, a, c = u_n[keep], L_n[keep], x_n[keep], y_n[keep], a2[keep], c2[keep]
        counts[N] = np.count_nonzero(a * np.roll(a, -1) < 0)
        npc[N] = x.size
    return counts, npc


def integrate_line(x0, y_row, kicks, N, prm):
    """Direct N-step integration of launch points (x0, y_row)."""
    st = new_state(x0, np.full(np.size(x0), float(y_row)))
    for k in range(N):
        step(st, kicks[k], prm.Tau, prm.row0)
    return st


# ============================================================ Part 1C: do low-stretch rays make the branches?
BR_KINDS = ['lowS', 'lowR', 'nullS', 'nullR']
BR_METRICS = ['share_branch', 'share_all', 'enrichment', 'C_all', 'C_no', 'lost', 'C_resel',
              'corr', 'jaccard', 'n_removed']
BR_POPS = ['line row0-rows/8', 'line row0', 'line row0+rows/8', 'tracer lattice']


def even_pick(n, k):
    """k of n indices spread evenly (deterministic, no random numbers): the null
    removal. Rounds half up, as floor(x + 0.5), so MATLAB gives the same set."""
    if k <= 0:
        return np.zeros(0, np.int64)
    return np.floor(np.linspace(0, n - 1, k) + 0.5).astype(np.int64)


def branch_stats(xm, S, r, nb, binw, thr, frac=0.10):
    """Part 1C for one population of rays at one N.

    rho = rays per x bin; branches B = the top `frac` of bins of rho (all rays).
    For each removal set -- S <= thr, r <= its 10th percentile, and the same
    numbers of rays removed evenly along the population (nulls) -- returns
    [share of branch ray mass in the set, its overall share, enrichment,
     contrast of rho on B, contrast of rho_no on the SAME B, fraction of the
     excess contrast lost (C_all - C_no)/(C_all - 1), contrast of rho_no on
     its own top bins, Pearson corr(rho, rho_no) across x, Jaccard overlap of
     the two branch sets, number removed].
    """
    n = xm.size
    xb = np.minimum((xm / binw).astype(np.int64), nb - 1)
    cnt = np.bincount(xb, minlength=nb).astype(float)
    nbr = max(int(np.floor(frac * nb + 0.5)), 1)
    B = np.argsort(-cnt, kind='stable')[:nbr]
    C_all = cnt[B].mean() / cnt.mean()
    q10 = np.quantile(r, 0.10)
    lowS = S <= thr
    lowR = r <= q10
    masks = [lowS, lowR]
    for m in (lowS, lowR):
        nm = np.zeros(n, bool); nm[even_pick(n, int(m.sum()))] = True
        masks.append(nm)
    res = np.full((len(masks), len(BR_METRICS)), np.nan)
    setB = set(B.tolist())
    for k, m in enumerate(masks):
        nrem = int(m.sum())
        crem = np.bincount(xb[m], minlength=nb).astype(float)
        cno = cnt - crem
        share_b = crem[B].sum() / cnt[B].sum()
        share_a = nrem / n
        enr = share_b / share_a if share_a > 0 else np.nan
        C_no = cno[B].mean() / cno.mean() if cno.sum() > 0 else np.nan
        lost = (C_all - C_no) / (C_all - 1) if C_all > 1 else np.nan
        B2 = np.argsort(-cno, kind='stable')[:nbr]
        C_resel = cno[B2].mean() / cno.mean() if cno.sum() > 0 else np.nan
        corr = pearson(cnt, cno)
        s2 = set(B2.tolist())
        jac = len(setB & s2) / len(setB | s2)
        res[k] = [share_b, share_a, enr, C_all, C_no, lost, C_resel, corr, jac, nrem]
    return res


def run_part1c(prm, kicks=None):
    """Part 1C alone: the launched lines and the window-seeded tracer lattice
    (the MATLAB tracer lattice), stepped with the run's kicks. Same numbers as
    the 'br' array of run()."""
    C, cols, rows, Tau, row0 = prm.cycles, prm.cols, prm.rows, prm.Tau, prm.row0
    kicks = make_kicks(prm) if kicks is None else [np.asarray(k, float).ravel() for k in kicks]
    x_line = np.arange(cols * prm.line_dens) / prm.line_dens
    pops = [new_state(x_line, np.full_like(x_line, yr)) for yr in prm.line_rows]
    ys = np.arange(0, rows, prm.lag_step, dtype=float); xs = np.arange(0, cols, prm.lag_step, dtype=float)
    XX, YY = np.meshgrid(xs, ys)                 # MATLAB's tracer order: column-major, y fastest
    pops.append(new_state(XX.ravel(order='F'), YY.ravel(order='F')))
    nbL = int(round(cols / prm.line_bin))
    br = np.full((len(pops), C + 1, len(BR_KINDS), len(BR_METRICS)), np.nan)
    for N in range(C + 1):
        if N > 0:
            for p in pops:
                step(p, kicks[N - 1], Tau, row0)
        for pi, p in enumerate(pops):
            S = sigma_max(p['a'], p['b'], p['c'], p['d'])
            br[pi, N] = branch_stats(np.mod(p['x'], cols), S, r_of(p), nbL, prm.line_bin, prm.S_THR)
    return dict(mode=prm.mode, seed=prm.seed, cycles=C, br=br)
