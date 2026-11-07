"""Runs all experiments for the HALO paper. Evaluation seeds 1..10 (tuning used seed 101 only)."""
import json, sys, os, time
from multiprocessing import Pool
import sim

SEEDS = list(range(1, 11))
BASE = dict(n_users=500, t_warm=1200.0, t_meas=3600.0, hap_lambda=4.0)
OUT = "results"
os.makedirs(OUT, exist_ok=True)


def job(args):
    tag, seed, over, schemes = args
    P = dict(BASE, **over)
    t0 = time.time()
    r = sim.run(seed, P, schemes=schemes)
    return tag, seed, r, time.time() - t0


def build():
    J = []
    # E1 main comparison
    for sd in SEEDS:
        J.append(("E1", sd, {}, sim.SCHEMES))
    # E2 LA size sweep
    for la in [2, 3, 6]:
        for sd in SEEDS:
            J.append((f"E2_la{la}", sd, dict(la_rows=la, la_cols=la), ["GSM", "TSP", "VAHO+TSP", "HALO"]))
    # E3 call-rate sweep (MO = MT = half of total)
    for tot in [1.0, 6.0]:
        for sd in SEEDS:
            J.append((f"E3_calls{tot:g}", sd, dict(lam_mo=tot/2/3600, lam_mt=tot/2/3600), ["GSM", "TSP", "VAHO+TSP", "HALO"]))
    # E4 mobility mixes
    mixes = {"ped": [(0.40, 0.0, 0.0), (0.50, 1.2, 0.30), (0.08, 10.0, 0.08), (0.02, 25.0, 0.03)],
             "veh": [(0.20, 0.0, 0.0), (0.20, 1.2, 0.30), (0.40, 10.0, 0.08), (0.20, 25.0, 0.03)]}
    for k, mx in mixes.items():
        for sd in SEEDS:
            J.append((f"E4_{k}", sd, dict(classes=mx), ["GSM", "VAHO", "TSP", "VAHO+TSP", "HALO"]))
    # E5 HAP delay-penalty sensitivity (paging frontier)
    for lam in [0.0, 2.0, 8.0, 16.0]:
        for sd in SEEDS[:5]:
            J.append((f"E5_lam{lam:g}", sd, dict(hap_lambda=lam), ["HALO"]))
    return J


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    J = [j for j in build() if (only is None or j[0].startswith(only))]
    done = set()
    fn = os.path.join(OUT, "raw.jsonl")
    if os.path.exists(fn):
        for line in open(fn):
            d = json.loads(line); done.add((d["tag"], d["seed"]))
    J = [j for j in J if (j[0], j[1]) not in done]
    print(f"{len(J)} jobs to run", flush=True)
    with Pool(2) as pool, open(fn, "a") as f:
        for tag, seed, r, el in pool.imap_unordered(job, J):
            f.write(json.dumps(dict(tag=tag, seed=seed, res=r)) + "\n"); f.flush()
            print(f"{tag} seed={seed} {el:.0f}s", flush=True)
    print("ALL DONE", flush=True)
