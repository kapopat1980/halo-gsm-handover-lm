import json, collections, numpy as np, os
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = collections.defaultdict(dict)   # R[tag][seed] = res
for line in open("results/raw.jsonl"):
    d = json.loads(line); R[d["tag"]][d["seed"]] = d["res"]

KEYS = ["lau_per_uh", "lau_post", "ilr_per_uh", "page_cells_per_call", "page_miss_rate", "page_delay_extra",
        "lm_cost_per_uh", "ho_per_callmin", "pp_ratio", "drop_rate", "poorq_frac"]


def agg(tag, scheme, key):
    xs = np.array([R[tag][s][scheme][key] for s in sorted(R[tag]) if scheme in R[tag][s]], float)
    n = len(xs); m = xs.mean(); sd = xs.std(ddof=1) if n > 1 else 0
    h = stats.t.ppf(0.975, n - 1)*sd/np.sqrt(n) if n > 1 else 0
    return m, h, xs


def paired(tag, a, b, key):
    xa, xb = agg(tag, a, key)[2], agg(tag, b, key)[2]
    try:
        p = stats.wilcoxon(xa, xb).pvalue
    except ValueError:
        p = float("nan")
    return (xa.mean() - xb.mean())/xb.mean()*100 if xb.mean() else float("nan"), p


out = {}
for tag in sorted(R):
    sch = list(next(iter(R[tag].values())).keys())
    out[tag] = {s: {k: agg(tag, s, k)[:2] for k in KEYS} for s in sch}
    out[tag]["_n"] = len(R[tag])
json.dump(out, open("results/summary.json", "w"), indent=1, default=float)

# pairwise tests for E1
tests = {}
if "E1" in R:
    for b in ["GSM", "VAHO", "TSP", "VAHO+TSP"]:
        tests[b] = {k: paired("E1", "HALO", b, k) for k in KEYS}
json.dump(tests, open("results/tests_E1.json", "w"), indent=1, default=float)

# print E1 table
if "E1" in out:
    print("E1 (n=%d)" % out["E1"]["_n"])
    for s in out["E1"]:
        if s == "_n": continue
        print(f"{s:11s} " + " ".join(f"{k}={out['E1'][s][k][0]:.4g}±{out['E1'][s][k][1]:.2g}" for k in KEYS))
    for b, d in tests.items():
        print("HALO vs", b, {k: (round(v[0], 1), float('%.3g' % v[1])) for k, v in d.items()})

# --------------- figures ---------------
os.makedirs("figs", exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "axes.grid.axis": "y", "grid.alpha": 0.3, "grid.linewidth": 0.5})
COL = {"GSM": "#7f7f7f", "VAHO": "#9e7bb5", "TSP": "#4c8bbf", "VAHO+TSP": "#2a5d8a", "HALO": "#c0392b",
       "HALO-noILR": "#e59866", "HALO-noHAP": "#f0b27a", "HALO-noDAH": "#f5cba7"}

if "E1" in out:
    sch = ["GSM", "VAHO", "TSP", "VAHO+TSP", "HALO"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for a, (k, lab) in zip(ax, [("lm_cost_per_uh", "LM signalling cost\n(msg / user·h)"),
                                ("page_cells_per_call", "Cells paged per MT call"),
                                ("page_miss_rate", "First-step paging miss ratio")]):
        m = [out["E1"][s][k][0] for s in sch]; h = [out["E1"][s][k][1] for s in sch]
        a.bar(range(len(sch)), m, yerr=h, color=[COL[s] for s in sch], capsize=2, width=0.7)
        a.set_xticks(range(len(sch))); a.set_xticklabels(sch, rotation=40, ha="right"); a.set_title(lab, fontsize=8.5)
    fig.tight_layout(); fig.savefig("figs/fig_lm.png", dpi=300); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for a, (k, lab, sc) in zip(ax, [("ho_per_callmin", "Handovers per call-minute", 1),
                                    ("pp_ratio", "Ping-pong ratio (%)", 100),
                                    ("poorq_frac", "Call time with C/I < 9 dB (%)", 100)]):
        m = [out["E1"][s][k][0]*sc for s in sch]; h = [out["E1"][s][k][1]*sc for s in sch]
        a.bar(range(len(sch)), m, yerr=h, color=[COL[s] for s in sch], capsize=2, width=0.7)
        a.set_xticks(range(len(sch))); a.set_xticklabels(sch, rotation=40, ha="right"); a.set_title(lab, fontsize=8.5)
    fig.tight_layout(); fig.savefig("figs/fig_ho.png", dpi=300); plt.close(fig)

    # ablation: relative to HALO
    ab = ["HALO", "HALO-noILR", "HALO-noHAP", "HALO-noDAH", "GSM"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
    for a, (k, lab, sc) in zip(ax, [("lau_post", "Post-call LAUs (count / run)", 1),
                                    ("page_cells_per_call", "Cells paged per MT call", 1),
                                    ("pp_ratio", "Ping-pong ratio (%)", 100)]):
        m = [out["E1"][s][k][0]*sc for s in ab]; h = [out["E1"][s][k][1]*sc for s in ab]
        a.bar(range(len(ab)), m, yerr=h, color=[COL[s] for s in ab], capsize=2, width=0.7)
        a.set_xticks(range(len(ab))); a.set_xticklabels(ab, rotation=40, ha="right"); a.set_title(lab, fontsize=8.5)
    fig.tight_layout(); fig.savefig("figs/fig_ablation.png", dpi=300); plt.close(fig)

# LA-size sweep
las = [(2, "E2_la2"), (3, "E2_la3"), (4, "E1"), (6, "E2_la6")]
if all(t in out for _, t in las):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
    for s in ["GSM", "TSP", "VAHO+TSP", "HALO"]:
        xs = [l*l for l, _ in las]
        m = [out[t][s]["lm_cost_per_uh"][0] for _, t in las]; h = [out[t][s]["lm_cost_per_uh"][1] for _, t in las]
        ax[0].errorbar(xs, m, yerr=h, marker="o", ms=3.5, color=COL[s], label=s, capsize=2, lw=1.4)
        m = [out[t][s]["page_cells_per_call"][0] for _, t in las]
        ax[1].plot(xs, m, marker="o", ms=3.5, color=COL[s], label=s, lw=1.4)
    ax[0].set_xlabel("Cells per location area"); ax[0].set_ylabel("LM cost (msg / user·h)")
    ax[1].set_xlabel("Cells per location area"); ax[1].set_ylabel("Cells paged per MT call")
    ax[0].legend(frameon=False, fontsize=7.5)
    fig.tight_layout(); fig.savefig("figs/fig_la.png", dpi=300); plt.close(fig)

# paging frontier
lams = [(0, "E5_lam0"), (2, "E5_lam2"), (4, "E1"), (8, "E5_lam8"), (16, "E5_lam16")]
if all(t in out for _, t in lams):
    fig, a = plt.subplots(figsize=(3.6, 2.8))
    x = [out[t]["HALO"]["page_cells_per_call"][0] for _, t in lams]; y = [out[t]["HALO"]["page_miss_rate"][0]*100 for _, t in lams]
    a.plot(x, y, "-o", color=COL["HALO"], ms=4, label="HALO (HAP), λ = 0…16")
    for (l, _), xi, yi in zip(lams, x, y): a.annotate(f"λ={l}", (xi, yi), textcoords="offset points", xytext=(4, 3), fontsize=7)
    a.plot(out["E1"]["TSP"]["page_cells_per_call"][0], out["E1"]["TSP"]["page_miss_rate"][0]*100, "s", color=COL["TSP"], label="TSP (last cell first)")
    a.plot(out["E1"]["GSM"]["page_cells_per_call"][0], 0, "^", color=COL["GSM"], label="GSM (LA-wide)")
    a.set_xlabel("Cells paged per MT call"); a.set_ylabel("First-step miss ratio (%)"); a.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig("figs/fig_frontier.png", dpi=300); plt.close(fig)
print("figures written")
