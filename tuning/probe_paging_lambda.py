import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import sim,sys
base=dict(n_users=400, t_warm=300, t_meas=2400)
lam=float(sys.argv[1]) if len(sys.argv)>1 else 8
r=sim.run(101, dict(base, hap_lambda=lam), schemes=["GSM","VAHO","TSP","HALO"])
for n,x in r.items(): print(f"{n:6s} HO/min={x['ho_per_callmin']:.3f} pp={x['pp_ratio']:.3f} poorq={x['poorq_frac']:.4f} drop={x['drop_rate']:.4f} LAU/uh={x['lau_per_uh']:.2f} post={x['lau_post']} ilr={x['ilr']} pg={x['page_cells_per_call']:.2f} miss={x['page_miss_rate']:.3f} cost={x['lm_cost_per_uh']:.2f}",flush=True)
