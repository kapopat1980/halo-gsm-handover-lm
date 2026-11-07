import json
S = json.load(open("results/summary.json")); T = json.load(open("results/tests_E1.json"))
N = {}
def g(tag, s, k, i=0): return S[tag][s][k][i]
def f(x, d=2): return f"{x:.{d}f}"
def pct(a, b): return (a - b)/b*100

rows = []
for s in ["GSM", "VAHO", "TSP", "VAHO+TSP", "HALO", "HALO-noILR", "HALO-noHAP", "HALO-noDAH"]:
    lab = s.replace("HALO-no", "HALO−")
    rows.append([lab,
        f"{f(g('E1',s,'lau_per_uh'))} ± {f(g('E1',s,'lau_per_uh',1))}",
        f"{g('E1',s,'lau_post'):.0f} ± {g('E1',s,'lau_post',1):.0f}",
        f"{f(g('E1',s,'page_cells_per_call'))} ± {f(g('E1',s,'page_cells_per_call',1))}",
        f"{g('E1',s,'page_miss_rate')*100:.1f} ± {g('E1',s,'page_miss_rate',1)*100:.1f}",
        f"{g('E1',s,'lm_cost_per_uh'):.1f} ± {g('E1',s,'lm_cost_per_uh',1):.1f}",
        f"{f(g('E1',s,'ho_per_callmin'),3)} ± {f(g('E1',s,'ho_per_callmin',1),3)}",
        f"{g('E1',s,'pp_ratio')*100:.1f} ± {g('E1',s,'pp_ratio',1)*100:.1f}",
        f"{g('E1',s,'poorq_frac')*100:.2f} ± {g('E1',s,'poorq_frac',1)*100:.2f}",
        f"{g('E1',s,'drop_rate')*100:.2f} ± {g('E1',s,'drop_rate',1)*100:.2f}"])
N["main_rows"] = rows
N["main_note"] = ("Post-call LAU counts are totals per replication (500 users, 1 h). C_LM uses w_LU = 10, w_ILR = 2, w_P = 1. "
                  "Paired two-sided Wilcoxon signed-rank tests (ten pairs; minimum attainable p = 0.002): HALO differs from every baseline with p = 0.002 in LAUs per user-hour, post-call LAUs, cells paged, C_LM and ping-pong ratio; "
                  "its first-step miss ratio is lower than TSP and VAHO+TSP with p = 0.002 (GSM and VAHO page the whole LA in one step and have no misses); "
                  "handovers per call-minute are lower than GSM (p = 0.002) and VAHO (p = 0.049); differences in drop rate and in time with C/I < 9 dB are not significant (p > 0.05).")
for s in ["GSM", "VAHO", "TSP", "VAHO+TSP", "HALO", "HALO-noILR", "HALO-noHAP", "HALO-noDAH"]:
    k = s.replace("+", "p").replace("-", "_")
    for m in ["lau_per_uh", "lau_post", "page_cells_per_call", "page_miss_rate", "lm_cost_per_uh", "ho_per_callmin", "pp_ratio", "drop_rate", "poorq_frac", "ilr_per_uh", "page_delay_extra"]:
        N[f"{k}_{m}"] = S["E1"][s][m][0]

# convenient formatted values
E = lambda s, m: S["E1"][s][m][0]
N.update({
 "cost_gsm": f"{E('GSM','lm_cost_per_uh'):.1f}", "cost_tsp": f"{E('TSP','lm_cost_per_uh'):.1f}", "cost_vt": f"{E('VAHO+TSP','lm_cost_per_uh'):.1f}", "cost_halo": f"{E('HALO','lm_cost_per_uh'):.1f}",
 "cost_red_gsm": f"{-pct(E('HALO','lm_cost_per_uh'), E('GSM','lm_cost_per_uh')):.1f}",
 "cost_red_vt": f"{-pct(E('HALO','lm_cost_per_uh'), E('VAHO+TSP','lm_cost_per_uh')):.1f}",
 "post_gsm": f"{E('GSM','lau_post'):.0f}", "post_halo": f"{E('HALO','lau_post'):.0f}",
 "post_red": f"{-pct(E('HALO','lau_post'), E('GSM','lau_post')):.0f}",
 "lau_gsm": f"{E('GSM','lau_per_uh'):.2f}", "lau_halo": f"{E('HALO','lau_per_uh'):.2f}", "lau_red": f"{-pct(E('HALO','lau_per_uh'), E('GSM','lau_per_uh')):.1f}",
 "ilr_uh": f"{E('HALO','ilr_per_uh'):.2f}",
 "pg_gsm": f"{E('GSM','page_cells_per_call'):.0f}", "pg_tsp": f"{E('TSP','page_cells_per_call'):.2f}", "pg_halo": f"{E('HALO','page_cells_per_call'):.2f}",
 "pg_red_gsm": f"{-pct(E('HALO','page_cells_per_call'), E('GSM','page_cells_per_call')):.1f}",
 "pg_red_tsp": f"{-pct(E('HALO','page_cells_per_call'), E('TSP','page_cells_per_call')):.1f}",
 "miss_tsp": f"{E('TSP','page_miss_rate')*100:.1f}", "miss_halo": f"{E('HALO','page_miss_rate')*100:.1f}",
 "miss_red": f"{-pct(E('HALO','page_miss_rate'), E('TSP','page_miss_rate')):.0f}", "miss_red_vt": f"{-pct(E('HALO','page_miss_rate'), E('VAHO+TSP','page_miss_rate')):.0f}",
 "dly_tsp": f"{E('TSP','page_delay_extra'):.2f}", "dly_halo": f"{E('HALO','page_delay_extra'):.2f}",
 "ho_gsm": f"{E('GSM','ho_per_callmin'):.3f}", "ho_vaho": f"{E('VAHO','ho_per_callmin'):.3f}", "ho_halo": f"{E('HALO','ho_per_callmin'):.3f}",
 "ho_red_gsm": f"{-pct(E('HALO','ho_per_callmin'), E('GSM','ho_per_callmin')):.1f}",
 "pp_gsm": f"{E('GSM','pp_ratio')*100:.1f}", "pp_vaho": f"{E('VAHO','pp_ratio')*100:.1f}", "pp_halo": f"{E('HALO','pp_ratio')*100:.1f}",
 "pp_red_gsm": f"{-pct(E('HALO','pp_ratio'), E('GSM','pp_ratio')):.1f}", "pp_red_vaho": f"{-pct(E('HALO','pp_ratio'), E('VAHO','pp_ratio')):.1f}",
 "drop_gsm": f"{E('GSM','drop_rate')*100:.2f}", "drop_halo": f"{E('HALO','drop_rate')*100:.2f}",
 "pq_gsm": f"{E('GSM','poorq_frac')*100:.2f}", "pq_halo": f"{E('HALO','poorq_frac')*100:.2f}",
 "cost_noilr": f"{E('HALO-noILR','lm_cost_per_uh'):.1f}", "cost_nohap": f"{E('HALO-noHAP','lm_cost_per_uh'):.1f}",
 "pp_nodah": f"{E('HALO-noDAH','pp_ratio')*100:.1f}",
})
# idle share of LAUs
# sweep rows
sw = []
for la, tag in [(2, "E2_la2"), (3, "E2_la3"), (4, "E1"), (6, "E2_la6")]:
    c = {s: S[tag][s]["lm_cost_per_uh"][0] for s in ["GSM", "TSP", "VAHO+TSP", "HALO"]}
    pg = {s: S[tag][s]["page_cells_per_call"][0] for s in ["GSM", "TSP", "HALO"]}
    sw.append([str(la*la), f"{c['GSM']:.1f}", f"{c['TSP']:.1f}", f"{c['VAHO+TSP']:.1f}", f"{c['HALO']:.1f}",
               f"{pct(c['HALO'], c['GSM']):+.1f}%", f"{pct(c['HALO'], c['VAHO+TSP']):+.1f}%",
               f"{pg['GSM']:.0f} / {pg['TSP']:.2f} / {pg['HALO']:.2f}"])
    N[f"sw{la}"] = c; N[f"swpg{la}"] = pg
N["sweep_rows"] = sw
gmin = min(((la, N[f"sw{la}"]["GSM"]) for la in [2, 3, 4, 6]), key=lambda x: x[1])
hmin = min(((la, N[f"sw{la}"]["HALO"]) for la in [2, 3, 4, 6]), key=lambda x: x[1])
vmin = min(((la, N[f"sw{la}"]["VAHO+TSP"]) for la in [2, 3, 4, 6]), key=lambda x: x[1])
N.update({"gmin_la": f"{gmin[0]*gmin[0]}", "gmin": f"{gmin[1]:.1f}", "hmin_la": f"{hmin[0]*hmin[0]}", "hmin": f"{hmin[1]:.1f}",
          "vmin_la": f"{vmin[0]*vmin[0]}", "vmin": f"{vmin[1]:.1f}",
          "opt_red_gsm": f"{-pct(hmin[1], gmin[1]):.1f}", "opt_red_vt": f"{-pct(hmin[1], vmin[1]):.1f}",
          "la36_red_gsm": f"{-pct(N['sw6']['HALO'], N['sw6']['GSM']):.1f}", "la4_red_gsm": f"{-pct(N['sw2']['HALO'], N['sw2']['GSM']):.1f}",
          "la36_red_vt": f"{-pct(N['sw6']['HALO'], N['sw6']['VAHO+TSP']):.1f}", "la4_red_vt": f"{-pct(N['sw2']['HALO'], N['sw2']['VAHO+TSP']):.1f}"})
for k in ["sw2", "sw3", "sw4", "sw6", "swpg2", "swpg3", "swpg4", "swpg6"]: N.pop(k)

sc = []
for lab, tag in [("Low traffic (1 call/h)", "E3_calls1"), ("Default (1.5 calls/h)", "E1"), ("High traffic (6 calls/h)", "E3_calls6"),
                 ("Pedestrian-dominant mix", "E4_ped"), ("Vehicular-dominant mix", "E4_veh")]:
    a, b = S[tag]["VAHO+TSP"], S[tag]["HALO"]
    sc.append([lab, f"{a['lm_cost_per_uh'][0]:.1f}", f"{b['lm_cost_per_uh'][0]:.1f} ({pct(b['lm_cost_per_uh'][0], a['lm_cost_per_uh'][0]):+.1f}%)",
               f"{a['page_miss_rate'][0]*100:.1f}", f"{b['page_miss_rate'][0]*100:.1f}",
               f"{a['pp_ratio'][0]*100:.1f}", f"{b['pp_ratio'][0]*100:.1f}"])
    key = tag.replace("E3_calls", "c").replace("E4_", "")
    N[f"sc_{key}_cost"] = f"{pct(b['lm_cost_per_uh'][0], a['lm_cost_per_uh'][0]):+.1f}"
    N[f"sc_{key}_pp"] = f"{pct(b['pp_ratio'][0], a['pp_ratio'][0]):+.1f}"
    N[f"sc_{key}_costgsm"] = f"{pct(b['lm_cost_per_uh'][0], S[tag]['GSM']['lm_cost_per_uh'][0]):+.1f}"
N["scen_rows"] = sc
# frontier
for lam, tag in [(0, "E5_lam0"), (2, "E5_lam2"), (8, "E5_lam8"), (16, "E5_lam16")]:
    N[f"fr{lam}_pg"] = f"{S[tag]['HALO']['page_cells_per_call'][0]:.2f}"; N[f"fr{lam}_miss"] = f"{S[tag]['HALO']['page_miss_rate'][0]*100:.1f}"
# VAHO+TSP vs GSM in E4 veh: VAHO effect
N["vaho_veh_pp"] = f"{S['E4_veh']['VAHO']['pp_ratio'][0]*100:.1f}"; N["gsm_veh_pp"] = f"{S['E4_veh']['GSM']['pp_ratio'][0]*100:.1f}"; N["halo_veh_pp"] = f"{S['E4_veh']['HALO']['pp_ratio'][0]*100:.1f}"
N["drop_veh_gsm"] = f"{S['E4_veh']['GSM']['drop_rate'][0]*100:.2f}"; N["drop_veh_halo"] = f"{S['E4_veh']['HALO']['drop_rate'][0]*100:.2f}"; N["drop_veh_vaho"] = f"{S['E4_veh']['VAHO']['drop_rate'][0]*100:.2f}"
# cost-weight sensitivity computed from raw counts (E1): ratio w_LU/w_P in {5, 20}
import collections
R = collections.defaultdict(dict)
for line in open("results/raw.jsonl"):
    d = json.loads(line)
    if d["tag"] == "E1": R[d["seed"]] = d["res"]
def cost(r, wl, wi, wp): return (wl*r["lau"] + wi*r["ilr"] + wp*r["page_cells"])/r["user_hours"]
for wl in [5, 20]:
    for b in ["GSM", "VAHO+TSP"]:
        h = sum(cost(R[s]["HALO"], wl, wl/5, 1) for s in R)/len(R); x = sum(cost(R[s][b], wl, wl/5, 1) for s in R)/len(R)
        N[f"w{wl}_{b.replace('+','p')}"] = f"{-pct(h, x):.1f}"
json.dump(N, open("numbers.json", "w"), indent=1, default=float)
print({k: v for k, v in N.items() if not isinstance(v, list)})
