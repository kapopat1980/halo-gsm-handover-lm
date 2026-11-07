import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import sim
base=dict(n_users=400, t_warm=300, t_meas=2400)
def show(tag,r): print(f"{tag:40s} HO/min={r['ho_per_callmin']:.3f} pp={r['pp_ratio']:.3f} drop={r['drop_rate']:.4f} poorq={r['poorq_frac']:.4f} LAU/uh={r['lau_per_uh']:.2f} clsLAU={r['cls_lau']} post={r['lau_post']} ilr={r['ilr']} pg={r['page_cells_per_call']:.2f} miss={r['page_miss_rate']:.3f}",flush=True)
for hom in [4,6,8]:
    r=sim.run(101, dict(base, hom_std=hom), schemes=["GSM","TSP","HALO"])
    show(f"GSM hom={hom}", r["GSM"]); 
    if hom==4: show("TSP",r["TSP"]); show("HALO default",r["HALO"])
for hh,hl,g,tg in [(6,4,3,10),(6,4,6,15),(8,4,4,15),(7,3,4,15)]:
    r=sim.run(101, dict(base, hom_hi=hh, hom_lo=hl, pp_guard=g, t_guard=tg), schemes=["HALO"])["HALO"]
    show(f"HALO hi={hh} lo={hl} g={g} tg={tg}", r)
