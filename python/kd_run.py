"""
kd_run -- multi-seed driver for kd_engine.

    python3 kd_run.py                          # 3 modes x seeds 0..9 + i.i.d. model
    python3 kd_run.py --modes annealed --seeds 0 1 2 --workers 3

Each run is pickled to <out>/run_<mode>_s<seed>.pkl.gz, the i.i.d. model to
<out>/iid.pkl.gz. Cosine mode uses no random numbers in its dynamics, so all
its seeds are identical; it is integrated once and copied (sd = 0 by
construction, stated in the report).
"""
import argparse, gzip, os, pickle, time
from multiprocessing import Pool
import numpy as np
import kd_engine as K


def _save(obj, path):
    with gzip.open(path, 'wb', compresslevel=4) as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def load(path):
    with gzip.open(path, 'rb') as f:
        return pickle.load(f)


def one(args):
    mode, seed, cycles, out = args
    path = os.path.join(out, f'run_{mode}_s{seed}.pkl.gz')
    if os.path.exists(path):
        return path, 0.0
    t = time.time()
    first = seed == 0
    prm = K.Params(cycles=cycles, mode=mode, seed=seed,
                   caustic_dens=(64, 512, 4096) if first else (64, 512),
                   exact_max_pieces=1.2e7)
    res = K.run(prm, keep_snapshots=first)
    _save(res, path)
    return path, time.time() - t


def p1c_job(args):
    mode, seed, cycles, out = args
    path = os.path.join(out, f'p1c_{mode}_s{seed}.pkl.gz')
    if not os.path.exists(path):
        _save(K.run_part1c(K.Params(cycles=cycles, mode=mode, seed=seed)), path)
    return path, 0.0


def iid_job(out, n_samples=10 ** 7, cycles=60):
    path = os.path.join(out, 'iid.pkl.gz')
    if os.path.exists(path):
        return path, 0.0
    t = time.time()
    prm = K.Params(cycles=cycles)
    pool = K.fp_pool(prm, n_kicks=500)
    res = K.iid_model(prm, n_samples=n_samples, Nmax=cycles, pool=pool)
    res['pool_moments'] = dict(rms_a=float((prm.Tau * pool).std()),
                               kurt=float(((pool - pool.mean()) ** 4).mean() / pool.var() ** 2 - 3))
    qpos = np.linspace(0, 4, 17)
    res['TO_q'] = qpos
    res['TO_L'] = K.transfer_operator_L(prm, pool, qpos, n_theta=2048, n_fp=301)
    _save(res, path)
    return path, time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--modes', nargs='+', default=['annealed', 'quenched', 'cosine'])
    ap.add_argument('--seeds', nargs='+', type=int, default=list(range(10)))
    ap.add_argument('--cycles', type=int, default=60)
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(__file__), 'out'))
    ap.add_argument('--no-iid', action='store_true')
    ap.add_argument('--only-part1c', action='store_true', help='just the Part 1C branch statistics')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    for m in a.modes:
        for s in (a.seeds if m != 'cosine' else [0]):
            jobs.append((m, s, a.cycles, a.out))
    jobs.sort(key=lambda j: j[1] != 0)           # the heavier seed-0 runs first
    with Pool(a.workers, maxtasksperchild=1) as P:   # Part 1C (seconds per run)
        for p, _ in P.imap_unordered(p1c_job, jobs, chunksize=1):
            print(p, flush=True)
    if a.only_part1c:
        return
    if not a.no_iid:
        p, dt = iid_job(a.out, cycles=a.cycles)
        print(f'iid model -> {p} ({dt:.0f}s)', flush=True)
    # one task per worker process, so each run's memory is returned to the OS
    with Pool(a.workers, maxtasksperchild=1) as P:
        for p, dt in P.imap_unordered(one, jobs, chunksize=1):
            print(f'{p} ({dt:.0f}s)', flush=True)
    if 'cosine' in a.modes:                       # identical for every seed
        src = load(os.path.join(a.out, 'run_cosine_s0.pkl.gz'))
        for s in a.seeds:
            if s != 0:
                d = dict(src); d['seed'] = s; d.pop('snaps', None)
                _save(d, os.path.join(a.out, f'run_cosine_s{s}.pkl.gz'))


if __name__ == '__main__':
    main()
