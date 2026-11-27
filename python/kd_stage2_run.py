"""
Stage 2: does the past matter, low stretch without the number 2, orientation,
two new plots. All analyses in Python on the MATLAB map and kicks.

    python3 kd_stage2_run.py --cycles 15 --seeds 0 --modes annealed quenched     # first run
    python3 kd_stage2_run.py --cycles 30 --seeds 0 1 2 3 4                       # full run

Figures -> ../figures_stage2/<mode>/ (cross-mode figures in ../figures_stage2/all_modes/),
data -> ../stage2_<mode>.npz, numbers for the report -> ../figures_stage2/summary.json.
"""
import argparse, json, os, time
import numpy as np
import kd_stage2 as S2
import kd_stage2_analysis as A
import kd_stage2_figs as F

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
FIGDIR = os.path.join(ROOT, 'figures_stage2')


class Timer:
    def __init__(self):
        self.t = {}

    def __call__(self, name, fn, *a, **k):
        t0 = time.time()
        out = fn(*a, **k)
        dt = time.time() - t0
        self.t[name] = self.t.get(name, 0.0) + dt
        print('    %-46s %7.1f s' % (name, dt), flush=True)
        return out


def flatten(prefix, obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            flatten('%s/%s' % (prefix, k), v, out)
    elif isinstance(obj, (list, tuple)) and obj and isinstance(obj[0], dict):
        for i, v in enumerate(obj):
            flatten('%s/%d' % (prefix, i), v, out)
    else:
        try:
            out[prefix] = np.asarray(obj)
        except Exception:
            out[prefix] = np.asarray(str(obj))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cycles', type=int, default=30)
    ap.add_argument('--seeds', type=int, nargs='+', default=[0, 1, 2, 3, 4])
    ap.add_argument('--modes', nargs='+', default=['annealed', 'quenched', 'cosine'])
    args = ap.parse_args()
    C = args.cycles
    T = Timer()
    res = {m: {} for m in args.modes}
    runs0 = {}
    for m in args.modes:
        for sd in args.seeds:
            print('== %s, seed %d, %d cycles' % (m, sd, C), flush=True)
            r = T('run (tracers, bands, exact centre lines)', lambda: S2.Stage2Run(m, sd, C).run(C))
            out = {}
            out['1A'] = T('1A checkpoint 10 + resume 5 vs 15 straight', A.part1a, m, sd)
            out['1B'] = T('1B memory test', A.part1b, r)
            out['1C'] = T('1C start cycle (annealed, every pixel)', A.part1c, r)
            out['1D'] = T('1D headline numbers, shifted grid', A.part1d, r)
            out['2A'] = T('2A P(S <= s) for 40 thresholds', A.part2a, r)
            out['2B'] = T('2B quantiles, persistence, branch link', A.part2b, r)
            out['3A'] = T('3A orientation checks against svd', A.part3a, r)
            out['3B'] = T('3B brightness vs S and orientation', A.part3b, r)
            out['3C'] = T('3C re-entry vs orientation', A.part3c, r)
            out['4A'] = T('4A stretch kick by kick', A.part4a, r)
            s5 = float(np.median(out['2B']['sq'][10:, 1])) if C >= 10 else float(np.median(out['2B']['sq'][1:, 1]))
            out['4B'] = T('4B cusps on the exact centre lines', A.part4b, r, s5)
            out['2C'] = T('2C the S* hypothesis', A.part2c, r, out['4B'], out['2B']['sq'])
            out['P1'] = out['2A']['P'][1]
            res[m][sd] = out
            if sd == args.seeds[0]:
                runs0[m] = r
            else:
                del r

    # ------------------------------------------------------------ figures
    print('== figures', flush=True)
    t0 = time.time()
    allf = os.path.join(FIGDIR, 'all_modes')
    sd0 = args.seeds[0]
    sstar = res[args.modes[0]][sd0]['2C']['S_star']
    for m in args.modes:
        fo = os.path.join(FIGDIR, m)
        L = [res[m][s] for s in args.seeds]
        F.fig1a(res[m][sd0]['1A'], m, fo)
        F.fig1b([x['1B'] for x in L], m, fo)
        if m == 'annealed':
            F.fig1c([x['1C'] for x in L], fo)
        F.fig2a_heat([x['2A'] for x in L], m, fo)
        F.fig3b(res[m][sd0]['3B'], m, fo)
        T('3D pictures (%s)' % m, F.fig3d, runs0[m], fo)
    F.fig2a_rate({m: [res[m][s]['2A'] for s in args.seeds] for m in args.modes}, allf)
    F.fig2b({m: [res[m][s]['2B'] for s in args.seeds] for m in args.modes}, allf)
    F.fig3c({m: [res[m][s]['3C'] for s in args.seeds] for m in args.modes}, allf)
    F.fig4a({m: res[m][sd0]['4A'] for m in args.modes}, allf)
    F.fig4b({m: [res[m][s]['4B'] for s in args.seeds] for m in args.modes}, sstar, allf)
    if 'annealed' in args.modes:
        F.fig4b_example(runs0['annealed'], res['annealed'][sd0]['4B'], allf)
    print('    %-46s %7.1f s' % ('figures', time.time() - t0))

    # ------------------------------------------------------------ checks
    print('\n================ PASS / FAIL ================')
    checks = []

    def chk(name, ok, detail):
        checks.append((name, bool(ok), detail))
        print('  %s  %s  (%s)' % ('PASS' if ok else 'FAIL', name, detail))
    for m in args.modes:
        for s in args.seeds:
            d = res[m][s]['1A']['diffs']
            worst = max(d.values())
            chk('Python resume bit-identical, %s seed %d (10 + 5 == 15)' % (m, s), worst == 0,
                'max |diff| = %g over positions, M, S, band amounts, fine maps, cusps' % worst)
    for m in args.modes:
        for s in args.seeds:
            a3 = res[m][s]['3A']
            chk('angles vs svd and the identity for a, %s seed %d' % (m, s),
                a3['max_angle_err_u'] < 1e-8 and a3['max_angle_err_v'] < 1e-8 and a3['max_identity_rel'] < 1e-10,
                'max angle error u1 %.1e, v1 %.1e rad; identity %.1e relative to S; %d parcels'
                % (a3['max_angle_err_u'], a3['max_angle_err_v'], a3['max_identity_rel'], a3['n']))
    if 'annealed' in args.modes:
        bl = [res['annealed'][s]['1B'] for s in args.seeds]
        rr = np.nanmean([b['rr'][:, :, 1] for b in bl], 0)          # [N, group], S + orientation
        Ns = bl[0]['N']
        worst = np.nanmax(np.abs(rr - 1))
        det = '; '.join('cycle %d: %s' % (N, ', '.join('%.3f' % v for v in rr[i])) for i, N in enumerate(Ns))
        chk('annealed memory ratio within 1 +- 0.05 once orientation is included', worst <= 0.05,
            'risk ratios for the 3 history groups, controlling for S + orientation: ' + det)
        c3 = [res['annealed'][s]['3C'] for s in args.seeds]
        if 8 in c3[0]['N']:
            j = list(c3[0]['N']).index(8)
            P = np.mean([c['p_unsigned'][j] for c in c3], 0)
            n = np.sum([c['n_unsigned'][j] for c in c3], 0)
            se = np.sqrt(P * (1 - P) / np.maximum(n, 1))
            tol = np.maximum(0.05, 2 * se)
            chk('re-entry vs orientation near the random-matrix reference (annealed, cycle 8)',
                np.all(np.abs(P - A.REENTRY_REF) <= tol),
                '2-D: %s vs reference %s (tolerance max(0.05, 2 SE) per bin)'
                % (', '.join('%.2f' % v for v in P), ', '.join('%.2f' % v for v in A.REENTRY_REF)))
    if 'annealed' in args.modes and 'quenched' in args.modes:
        for s in args.seeds:
            d = np.max(np.abs(res['annealed'][s]['P1'] - res['quenched'][s]['P1']))
            chk('P(S <= s) at cycle 1 identical in annealed and quenched, seed %d' % s, d == 0, 'max |diff| = %g' % d)
    for m in args.modes:
        worst = 0.0
        for s in args.seeds:
            pass
        r = runs0[m]
        am = np.array(r.rec['amount']); fm = np.array(r.rec['fine'])
        worst = max(np.max(np.abs(am.mean(-1) - 1)), np.max(np.abs(fm.mean(-1) - 1)))
        chk('band amounts average to 1 at every cycle and moment, %s' % m, worst < 1e-9, 'max |mean - 1| = %.1e' % worst)
    nfail = sum(1 for c in checks if not c[1])
    print('%d checks, %d failed.' % (len(checks), nfail))

    # ------------------------------------------------------------ save
    os.makedirs(FIGDIR, exist_ok=True)
    for m in args.modes:
        flat = {}
        for s in args.seeds:
            o = dict(res[m][s])
            o['1A'] = {k: v for k, v in o['1A'].items() if not k.startswith('scatter')}
            o['1D'] = {'main': o['1D'][0], 'shifted': o['1D'][1]}
            flatten('seed%d' % s, o, flat)
        np.savez_compressed(os.path.join(ROOT, 'stage2_%s.npz' % m), **flat)
    summ = dict(cycles=C, seeds=args.seeds, modes=args.modes, checks=[list(c) for c in checks],
                time=T.t)
    for m in args.modes:
        L = [res[m][s] for s in args.seeds]
        summ[m] = dict(
            part1a={k: (v if not isinstance(v, np.ndarray) else None) for k, v in L[0]['1A'].items() if not k.startswith('scatter')},
            part1b=dict(N=L[0]['1B']['N'].tolist(), rr=np.nanmean([x['1B']['rr'] for x in L], 0).tolist(),
                        lo=L[0]['1B']['lo'].tolist(), hi=L[0]['1B']['hi'].tolist(),
                        rr_sd=(np.nanstd([x['1B']['rr'] for x in L], 0, ddof=1).tolist() if len(L) > 1 else None),
                        dmed=np.nanmean([x['1B']['dmed'] for x in L], 0).tolist(),
                        n1=L[0]['1B']['n1'].tolist(), n0=L[0]['1B']['n0'].tolist()),
            part1c=dict(t0=L[0]['1C']['t0'].tolist(), P=np.mean([x['1C']['P'] for x in L], 0).tolist() if L[0]['1C']['P'].size else [],
                        P_at_2=[[float(np.interp(np.log(2), np.log(A.S_GRID), p)) for p in x['1C']['P']] for x in L]),
            part1d=dict(main=[{k: float(v) for k, v in x['1D'][0].items()} for x in L],
                        shifted=[{k: float(v) for k, v in x['1D'][1].items()} for x in L]),
            part2a=dict(rate=np.nanmean([x['2A']['rate'] for x in L], 0).tolist(),
                        rate_sd=(np.nanstd([x['2A']['rate'] for x in L], 0, ddof=1).tolist() if len(L) > 1 else None),
                        P2=[float(np.mean([x['2A']['P'][N, np.argmin(abs(A.S_GRID - 2))] for x in L])) for N in range(C + 1)]),
            part2b=dict(sq=np.mean([x['2B']['sq'] for x in L], 0).tolist(),
                        pers=np.nanmean([x['2B']['pers'] for x in L], 0).tolist(),
                        enrich=np.nanmean([x['2B']['branch_enrich'] for x in L], 0).tolist(),
                        branch_n=L[0]['2B']['branch_n'].tolist()),
            part2c={k: (v if not isinstance(v, dict) else {kk: list(vv) for kk, vv in v.items()})
                    for k, v in L[0]['2C'].items()},
            part3a=L[0]['3A'],
            part3b=dict(N=L[0]['3B']['N'].tolist(), r2=np.mean([x['3B']['r2'] for x in L], 0).tolist()),
            part3c=dict(N=L[0]['3C']['N'].tolist(), p_unsigned=np.mean([x['3C']['p_unsigned'] for x in L], 0).tolist(),
                        n_unsigned=np.sum([x['3C']['n_unsigned'] for x in L], 0).tolist(),
                        p_signed=np.mean([x['3C']['p_signed'] for x in L], 0).tolist()),
            part4a=dict(fates=np.bincount(L[0]['4A']['fate'], minlength=3).tolist()),
            part4b={k: L[0]['4B'][k] for k in ('n_cusps', 'n_matched', 'n_branches', 'n_branches_with_cusp', 's5')},
        )
        b4 = L[0]['4B']
        lg = np.log10(b4['cusps'][:, 4]) if b4['cusps'].size else np.zeros(0)
        summ[m]['part4b'].update(median_S_cusp=float(10 ** np.median(lg)) if lg.size else None,
                                 frac_cusps_S_le_2=float(np.mean(lg <= np.log10(2))) if lg.size else None)
        mat = b4['mat_at_cusps']; cc = 10 ** (0.5 * (S2.LOGS_EDGES[1:] + S2.LOGS_EDGES[:-1]))
        cdf = np.cumsum(mat) / max(mat.sum(), 1e-300)
        summ[m]['part4b']['median_S_material'] = float(cc[np.searchsorted(cdf, 0.5)])
        summ[m]['part4b']['frac_material_S_le_2'] = float(cdf[np.searchsorted(cc, 2.0) - 1])
    with open(os.path.join(FIGDIR, 'summary.json'), 'w') as f:
        json.dump(summ, f, indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    print('\nRun time per part (s):')
    for k, v in T.t.items():
        print('  %-50s %8.1f' % (k, v))


if __name__ == '__main__':
    main()
