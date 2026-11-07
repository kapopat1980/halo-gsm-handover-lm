import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import sim, sys, numpy as np
for hom in [2,4,6,8]:
    for n in [1,2,4]:
        r = sim.run(101, dict(n_users=400, t_warm=300, t_meas=2400, hom_std=hom, n_std=n), schemes=["GSM"])["GSM"]
        print(f"hom={hom} n={n} HO/min={r['ho_per_callmin']:.3f} pp={r['pp_ratio']:.3f} drop={r['drop_rate']:.4f} poorq={r['poorq_frac']:.3f} clsHO={r['cls_ho']} clsPP={r['cls_pp']} clsDrop={r['cls_drop']}/{r['cls_calls']} LAU/uh={r['lau_per_uh']:.2f} clsLAU={r['cls_lau']}", flush=True)
