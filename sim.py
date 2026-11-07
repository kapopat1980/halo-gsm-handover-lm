"""
HALO simulator: Handover-Anchored Location Optimisation for GSM networks.

Time-stepped system-level simulator of a GSM network (omni cells, hexagonal
layout on a torus, reuse-4 frequency plan, COST-231 Hata path loss,
Gudmundson-correlated log-normal shadowing). All schemes run in lockstep on
identical radio, mobility and traffic realisations (common random numbers),
so differences between schemes are attributable to the schemes alone.

Schemes
  GSM     : standard GSM: PBGT handover (fixed margin), post-call LAU, LA-wide paging
  VAHO    : velocity-aware handover margin using the TRUE user speed (oracle) + GSM LM
  TSP     : GSM handover + two-step paging (last known cell, then LA)
  VAHO+TSP: strongest combination of prior techniques
  HALO    : proposed (ILR + HAP + DAH)
  HALO-noILR, HALO-noHAP, HALO-noDAH : ablations
"""
import numpy as np, math, json, sys, time

SCHEMES = ["GSM", "VAHO", "TSP", "VAHO+TSP", "HALO", "HALO-noILR", "HALO-noHAP", "HALO-noDAH"]

DEFAULT = dict(
    R=12, C=12, ISD=1500.0, la_rows=4, la_cols=4,       # 144 cells, LA = 4x4 cells
    n_users=500, dt=0.5, t_warm=900.0, t_meas=3600.0,
    f_mhz=900.0, hb=30.0, hm=1.5, eirp=55.0, noise=-113.0, load=0.25,
    shadow_sigma=8.0, shadow_dcorr=50.0, shadow_site_corr=0.5,
    tau_lev=2.0, tau_idle=8.0, tau_q=2.0, t_resel_pen=15.0, resel_pen=5.0,
    # mobility classes: fraction, mean speed (m/s), heading diffusion (rad/sqrt(s))
    classes=[(0.40, 0.0, 0.0), (0.30, 1.2, 0.30), (0.22, 10.0, 0.08), (0.08, 25.0, 0.03)],
    lam_mo=1.5/3600, lam_mt=1.5/3600, hold=120.0,
    # GSM radio resource / handover parameters
    hom_std=4.0, n_std=2, q_ho=9.0, q_drop=6.0, t_rlt=8.0, t_pp=10.0,
    crh=4.0, hyst_intra=2.0, t3212=1800.0,
    # VAHO / DAH speed-to-margin mapping
    hom_hi=6.0, hom_lo=4.0, v1=3.0, v2=15.0, n_hi=3, n_lo=1,
    pp_guard=6.0, t_guard=15.0,
    # HAP
    hap_sigma0=0.55, hap_beta=0.6, hap_lambda=8.0, hap_vprior=1.0, hap_tmax=3600.0,
    mr_valid=1800.0,
    # signalling weights (messages)
    w_lau=10.0, w_ilr=2.0, w_page=1.0, t3113=5.0,
    msc_split=True,
)


def hata(d_m, f, hb, hm):
    a_hm = (1.1*math.log10(f) - 0.7)*hm - (1.56*math.log10(f) - 0.8)
    dk = np.maximum(d_m, 20.0)/1000.0
    return 69.55 + 26.16*math.log10(f) - 13.82*math.log10(hb) - a_hm + (44.9 - 6.55*math.log10(hb))*np.log10(dk)


class Net:
    def __init__(s, P):
        R, C, ISD = P["R"], P["C"], P["ISD"]
        assert R % 4 == 0 and C % 2 == 0
        s.W, s.H = C*ISD, R*ISD*math.sqrt(3)/2
        rr, cc = np.meshgrid(np.arange(R), np.arange(C), indexing="ij")
        rr, cc = rr.ravel(), cc.ravel()
        s.B = R*C
        s.bx = cc*ISD + (rr % 2)*ISD/2
        s.by = rr*ISD*math.sqrt(3)/2
        q = cc - (rr - (rr & 1))//2               # axial coordinates (odd-r)
        s.color = (q % 2) + 2*(rr % 2)              # valid reuse-4 colouring, torus-consistent
        lr, lc = P["la_rows"], P["la_cols"]
        s.la = (rr//lr)*(C//lc) + (cc//lc)
        s.n_la = (R//lr)*(C//lc)
        s.la_cells = [np.where(s.la == k)[0] for k in range(s.n_la)]
        # MSC areas: left / right half of the LAs
        s.msc = ((cc//lc) >= (C//lc)//2).astype(int) if P["msc_split"] else np.zeros(s.B, int)
        s.la_msc = np.array([s.msc[s.la_cells[k][0]] for k in range(s.n_la)])
        dx, dy = s.tor(s.bx[:, None] - s.bx[None, :], s.by[:, None] - s.by[None, :])
        s.cdist = np.hypot(dx, dy)
        s.ISD = ISD

    def tor(s, dx, dy):
        dx = dx - s.W*np.round(dx/s.W)
        dy = dy - s.H*np.round(dy/s.H)
        return dx, dy


def run(seed, P=None, schemes=SCHEMES, verbose=False):
    P = dict(DEFAULT, **(P or {}))
    net = Net(P)
    rng = np.random.default_rng(seed)
    N, B, dt = P["n_users"], net.B, P["dt"]
    T_end = P["t_warm"] + P["t_meas"]
    nsteps = int(round(T_end/dt))
    S = len(schemes)
    sidx = {n: i for i, n in enumerate(schemes)}

    # ---------------- users / mobility ----------------
    fr = np.array([c[0] for c in P["classes"]])
    cls = rng.choice(len(fr), size=N, p=fr/fr.sum())
    vmean = np.array([c[1] for c in P["classes"]])[cls]
    hdiff = np.array([c[2] for c in P["classes"]])[cls]
    px, py = rng.uniform(0, net.W, N), rng.uniform(0, net.H, N)
    th = rng.uniform(0, 2*np.pi, N)
    vel = np.maximum(vmean*(1 + 0.1*rng.standard_normal(N)), 0)

    # ---------------- traffic (shared) ----------------
    def poisson_times(lam):
        out = []
        for _ in range(N):
            ts, t = [], rng.exponential(1/lam)
            while t < T_end:
                ts.append(t); t += rng.exponential(1/lam)
            out.append(ts)
        return out
    mo, mt = poisson_times(P["lam_mo"]), poisson_times(P["lam_mt"])
    ev = []   # (step, user, type, hold)
    for u in range(N):
        for t in mo[u]: ev.append((int(t/dt), u, 0, rng.exponential(P["hold"])))
        for t in mt[u]: ev.append((int(t/dt), u, 1, rng.exponential(P["hold"])))
    ev.sort()
    ev_ptr = 0

    # ---------------- radio state (shared) ----------------
    def rss_now():
        dx, dy = net.tor(px[:, None] - net.bx[None, :], py[:, None] - net.by[None, :])
        d = np.hypot(dx, dy)
        return P["eirp"] - hata(d, P["f_mhz"], P["hb"], P["hm"]) + sh
    # shadowing = common (user-specific, all sites) + independent per link; inter-site correlation = shadow_site_corr
    cw, iw = math.sqrt(P["shadow_site_corr"]), math.sqrt(1 - P["shadow_site_corr"])
    shc = (P["shadow_sigma"]*rng.standard_normal((N, 1))).astype(np.float32)
    shi = (P["shadow_sigma"]*rng.standard_normal((N, B))).astype(np.float32)
    sh = cw*shc + iw*shi
    rss = rss_now()
    lev = rss.copy(); lev_i = rss.copy()
    a_lev, a_idle, a_q = math.exp(-dt/P["tau_lev"]), math.exp(-dt/P["tau_idle"]), math.exp(-dt/P["tau_q"])
    onehot = np.zeros((B, 4), np.float32); onehot[np.arange(B), net.color] = 1
    noise_mw = 10**(P["noise"]/10)

    # ---------------- per-scheme state ----------------
    camp0 = np.argmax(lev_i, 1)
    camp = np.tile(camp0, (S, 1))
    reg_la = net.la[camp].copy()
    serv = camp.copy()
    in_call = np.zeros((S, N), bool)
    call_end = np.zeros((S, N)); call_start = np.zeros((S, N))
    qavg = np.zeros((S, N)); rlt = np.zeros((S, N))
    hocnt = np.zeros((S, N), int)
    prev_cell = -np.ones((S, N), int); last_ho = -1e9*np.ones((S, N)); entry_t = np.zeros((S, N))
    t_lastmm = np.tile(-rng.uniform(0, P["t3212"], N), (S, 1))   # T3212 reference (desynchronised)
    # last known contact (network side)
    lc_cell = camp.copy(); lc_t = np.zeros((S, N))
    # mobility record (HALO): speed estimate, heading unit vector, validity time
    mr_v = -np.ones((S, N)); mr_ux = np.zeros((S, N)); mr_uy = np.zeros((S, N)); mr_hasdir = np.zeros((S, N), bool)
    mr_t = -1e9*np.ones((S, N))
    call_first_cell = np.zeros((S, N), int); call_v = -np.ones((S, N))
    cc_n = np.zeros((S, N)); cc_hit = np.zeros((S, N))   # contact-continuity index (HAP)

    def contact(s, u, t, cell, ref_cell, ref_t):
        if ref_t > 0 and t - ref_t >= 30.0:
            cc_n[s, u] += 1; cc_hit[s, u] += (cell == ref_cell)

    M = {n: dict(lau=0, lau_post=0, lau_idle=0, lau_per=0, ilr=0, page_cells=0, pages=0, page_miss=0,
                 ho=0, pp=0, calls=0, drops=0, call_time=0.0, poorq_time=0.0, busy=0,
                 inter_la_ho=0, cls_ho=[0]*4, cls_pp=[0]*4, cls_lau=[0]*4, cls_drop=[0]*4, cls_calls=[0]*4) for n in schemes}

    is_halo = np.array([n.startswith("HALO") for n in schemes])
    use_ilr = np.array([n in ("HALO", "HALO-noHAP", "HALO-noDAH") for n in schemes])
    use_hap = np.array([n in ("HALO", "HALO-noILR", "HALO-noDAH") for n in schemes])
    use_dah = np.array([n in ("HALO", "HALO-noILR", "HALO-noHAP") for n in schemes])
    use_vaho = np.array([n in ("VAHO", "VAHO+TSP") for n in schemes])
    use_tsp = np.array([n in ("TSP", "VAHO+TSP") for n in schemes])
    ISD = net.ISD

    def margin_from_speed(v):
        v = np.asarray(v, float)
        f = np.clip((v - P["v1"])/(P["v2"] - P["v1"]), 0, 1)
        hom = P["hom_hi"] - f*(P["hom_hi"] - P["hom_lo"])
        nreq = np.where(v <= P["v1"], P["n_hi"], np.where(v >= P["v2"], P["n_lo"], P["n_std"]))
        return hom, nreq

    hocond = np.zeros((S, N), int)
    last_resel = -1e9*np.ones((S, N))
    t_meas0 = P["t_warm"]

    def measuring(t): return t >= t_meas0

    def do_lau(s, u, t, kind):
        if is_halo[s] and kind == "idle" and lc_t[s, u] > 0:
            dtl = t - lc_t[s, u]
            dx, dy = net.tor(net.bx[camp[s, u]] - net.bx[lc_cell[s, u]], net.by[camp[s, u]] - net.by[lc_cell[s, u]])
            dist = math.hypot(dx, dy)
            if dtl >= 30.0:
                mr_v[s, u] = max(dist - 0.5*ISD, 0.0)/dtl; mr_t[s, u] = t
                if dist > 1: mr_ux[s, u], mr_uy[s, u], mr_hasdir[s, u] = dx/dist, dy/dist, True
                else: mr_hasdir[s, u] = False
        contact(s, u, t, camp[s, u], lc_cell[s, u], lc_t[s, u])
        reg_la[s, u] = net.la[camp[s, u]]
        lc_cell[s, u] = camp[s, u]; lc_t[s, u] = t; t_lastmm[s, u] = t
        if measuring(t):
            m = M[schemes[s]]; m["lau"] += 1; m["lau_" + kind] += 1; m["cls_lau"][cls[u]] += 1

    def page(s, u, t):
        """returns nothing; accounts paging cost for MT call to idle user u in scheme s"""
        cells = net.la_cells[reg_la[s, u]]
        target = camp[s, u]
        first = None
        if use_tsp[s]:
            if net.la[lc_cell[s, u]] == reg_la[s, u]:
                first = np.array([lc_cell[s, u]])
        elif use_hap[s]:
            e = t - lc_t[s, u]
            c0 = lc_cell[s, u]
            if e <= P["hap_tmax"] and net.la[c0] == reg_la[s, u]:
                fresh = mr_v[s, u] >= 0 and (t - mr_t[s, u]) <= P["mr_valid"]
                v = mr_v[s, u] if fresh else P["hap_vprior"]
                cx, cy = net.bx[c0], net.by[c0]
                if fresh and mr_hasdir[s, u]:
                    mx, my = cx + v*e*mr_ux[s, u], cy + v*e*mr_uy[s, u]
                    sig = P["hap_sigma0"]*ISD + P["hap_beta"]*v*e
                else:
                    mx, my = cx, cy
                    sig = P["hap_sigma0"]*ISD + v*e
                # location distribution over the LA's cells: "still in anchor" mass + Gaussian displacement
                pi_stay = (cc_hit[s, u] + 1.0)/(cc_n[s, u] + 2.0)   # learned contact-continuity index
                dx, dy = net.tor(net.bx[cells] - mx, net.by[cells] - my)
                g = np.exp(-0.5*(dx*dx + dy*dy)/(sig*sig)); g = g/g.sum()
                p = (1 - pi_stay)*g + pi_stay*(cells == c0)
                order = np.argsort(-p); cp = np.cumsum(p[order]); L = len(cells)
                ks = np.arange(1, L + 1)
                # expected cost of a two-step plan with step-1 = top-k cells (+ delay penalty per miss)
                ecost = ks + (1 - cp)*(L - ks) + P["hap_lambda"]*(1 - cp)
                k = int(ks[np.argmin(ecost)])
                if k < L:
                    first = cells[order[:k]]
        if measuring(t):
            m = M[schemes[s]]; m["pages"] += 1
            if first is None:
                m["page_cells"] += len(cells)
            elif target in first:
                m["page_cells"] += len(first)
            else:
                m["page_cells"] += len(cells); m["page_miss"] += 1   # step 1 + remaining cells
        contact(s, u, t, target, lc_cell[s, u], lc_t[s, u])
        lc_cell[s, u] = target; lc_t[s, u] = t

    def start_call(s, u, t, hold):
        in_call[s, u] = True
        call_start[s, u] = t; call_end[s, u] = t + hold
        serv[s, u] = camp[s, u]; entry_t[s, u] = t
        qavg[s, u] = 20.0; rlt[s, u] = 0; hocnt[s, u] = 0; hocond[s, u] = 0
        prev_cell[s, u] = -1; last_ho[s, u] = -1e9
        call_first_cell[s, u] = serv[s, u]
        # DAH starts from mobility record if fresh
        call_v[s, u] = mr_v[s, u] if (t - mr_t[s, u] <= P["mr_valid"]) else -1
        lc_cell[s, u] = serv[s, u]; lc_t[s, u] = t
        if measuring(t): M[schemes[s]]["calls"] += 1; M[schemes[s]]["cls_calls"][cls[u]] += 1

    def end_call(s, u, t, dropped):
        in_call[s, u] = False
        c = serv[s, u]
        camp[s, u] = c
        dur = t - call_start[s, u]
        # --- mobility record update (network side, from handover history) ---
        if is_halo[s]:
            if hocnt[s, u] == 0:
                vest = min(0.25*ISD/max(dur, 1.0), 2.0)
                if call_v[s, u] >= 0: vest = min(vest, call_v[s, u])
                mr_v[s, u] = vest; mr_hasdir[s, u] = False
            else:
                mr_v[s, u] = max(call_v[s, u], 0.0)
                f0 = call_first_cell[s, u]
                dx, dy = net.tor(net.bx[c] - net.bx[f0], net.by[c] - net.by[f0])
                nrm = math.hypot(dx, dy)
                if nrm > 1:
                    mr_ux[s, u], mr_uy[s, u], mr_hasdir[s, u] = dx/nrm, dy/nrm, True
                else:
                    mr_hasdir[s, u] = False
            mr_t[s, u] = t
        contact(s, u, t, c, call_first_cell[s, u], call_start[s, u])
        lc_cell[s, u] = c; lc_t[s, u] = t; t_lastmm[s, u] = t
        if net.la[c] != reg_la[s, u]:
            same_msc = net.la_msc[net.la[c]] == net.la_msc[reg_la[s, u]]
            if use_ilr[s] and (not dropped) and same_msc:
                reg_la[s, u] = net.la[c]
                if measuring(t): M[schemes[s]]["ilr"] += 1
            else:
                do_lau(s, u, t, "post")
        if measuring(t) and dropped: M[schemes[s]]["drops"] += 1; M[schemes[s]]["cls_drop"][cls[u]] += 1

    rngs = np.random.default_rng(seed + 10_000)    # radio randomness stream
    t0 = time.time()
    for k in range(1, nsteps + 1):
        t = k*dt
        # ---- mobility (Gauss-Markov heading, constant class speed) ----
        th += hdiff*math.sqrt(dt)*rng.standard_normal(N)
        px = (px + vel*np.cos(th)*dt) % net.W
        py = (py + vel*np.sin(th)*dt) % net.H
        # ---- shadowing (Gudmundson AR(1)) and RSS ----
        rho = np.exp(-vel*dt/P["shadow_dcorr"]).astype(np.float32)[:, None]
        inn = np.sqrt(1 - rho*rho)*P["shadow_sigma"]
        shc = rho*shc + inn*rngs.standard_normal((N, 1), dtype=np.float32)
        shi = rho*shi + inn*rngs.standard_normal((N, B), dtype=np.float32)
        sh = cw*shc + iw*shi
        rss = rss_now()
        lev = a_lev*lev + (1 - a_lev)*rss
        lev_i = a_idle*lev_i + (1 - a_idle)*rss
        lin = np.power(10.0, rss/10.0, dtype=np.float32)
        Icol = lin @ onehot
        top2 = np.argpartition(-lev, 2, axis=1)[:, :2]
        l0, l1 = lev[np.arange(N), top2[:, 0]], lev[np.arange(N), top2[:, 1]]
        best1 = np.where(l0 >= l1, top2[:, 0], top2[:, 1]); best2 = np.where(l0 >= l1, top2[:, 1], top2[:, 0])
        best_i = np.argmax(lev_i, 1)
        meas = measuring(t)
        ar = np.arange(N)
        for s in range(S):
            name = schemes[s]; m = M[name]
            ic = in_call[s]
            # ---------- dedicated mode ----------
            if ic.any():
                u = ar[ic]; sv = serv[s, u]
                Sl = lin[u, sv]
                I = P["load"]*(Icol[u, net.color[sv]] - Sl) + noise_mw
                cir = 10*np.log10(Sl/I)
                qavg[s, u] = a_q*qavg[s, u] + (1 - a_q)*cir
                if meas:
                    m["call_time"] += dt*len(u); m["poorq_time"] += dt*np.count_nonzero(qavg[s, u] < P["q_ho"])
                # radio link timeout
                bad = qavg[s, u] < P["q_drop"]
                rlt[s, u] = np.where(bad, rlt[s, u] + dt, 0.0)
                # handover candidate & margins
                cand = np.where(best1[u] == sv, best2[u], best1[u])
                diff = lev[u, cand] - lev[u, sv]
                if use_vaho[s]:
                    hom, nreq = margin_from_speed(vel[u])
                elif use_dah[s]:
                    cv = call_v[s, u]
                    hom, nreq = margin_from_speed(np.where(cv >= 0, cv, 8.0))
                    hom = np.where(cv >= 0, hom, P["hom_std"]); nreq = np.where(cv >= 0, nreq, P["n_std"])
                    guard = (cand == prev_cell[s, u]) & (t - last_ho[s, u] < P["t_guard"])
                    hom = hom + np.where(guard, P["pp_guard"], 0.0)
                else:
                    hom = np.full(len(u), P["hom_std"]); nreq = np.full(len(u), P["n_std"])
                cond = (diff > hom) | ((qavg[s, u] < P["q_ho"]) & (diff > 0.0))
                hocond[s, u] = np.where(cond, hocond[s, u] + 1, 0)
                fire = hocond[s, u] >= nreq
                for j in np.nonzero(fire)[0]:
                    uu = u[j]; old = serv[s, uu]; new = cand[j]
                    if meas:
                        m["ho"] += 1; m["cls_ho"][cls[uu]] += 1
                        if new == prev_cell[s, uu] and t - last_ho[s, uu] <= P["t_pp"]: m["pp"] += 1; m["cls_pp"][cls[uu]] += 1
                        if net.la[new] != net.la[old]: m["inter_la_ho"] += 1
                    # DAH speed estimate from dwell time (network-observable only)
                    if is_halo[s]:
                        # net-displacement estimator: robust to ping-pong (returning
                        # to an earlier cell yields ~zero displacement -> low speed)
                        el = max(t - call_start[s, uu], 10.0)
                        disp = net.cdist[call_first_cell[s, uu], new]
                        call_v[s, uu] = max(disp - 0.5*ISD, 0.0)/el if disp > 0 else 0.0
                    prev_cell[s, uu] = old; serv[s, uu] = new; last_ho[s, uu] = t; entry_t[s, uu] = t
                    hocnt[s, uu] += 1; hocond[s, uu] = 0; rlt[s, uu] = 0
                    lc_cell[s, uu] = new; lc_t[s, uu] = t
                drop = rlt[s, u] >= P["t_rlt"]
                for j in np.nonzero(drop)[0]:
                    end_call(s, u[j], t, True)
            # ---------- idle mode ----------
            idl = ~in_call[s]
            u = ar[idl]
            cp = camp[s, u]; bi = best_i[u]
            hy = np.where(net.la[bi] != net.la[cp], P["crh"], P["hyst_intra"])
            hy = hy + np.where(t - last_resel[s, u] < P["t_resel_pen"], P["resel_pen"], 0.0)  # GSM 15 s rule
            resel = (bi != cp) & (lev_i[u, bi] > lev_i[u, cp] + hy)
            camp[s, u[resel]] = bi[resel]; last_resel[s, u[resel]] = t
            chg = u[resel][net.la[bi[resel]] != reg_la[s, u[resel]]]
            for uu in chg: do_lau(s, uu, t, "idle")
            per = u[(t - t_lastmm[s, u]) >= P["t3212"]]
            for uu in per: do_lau(s, uu, t, "per")
            # planned call ends
            ce = ar[in_call[s] & (call_end[s] <= t)]
            for uu in ce: end_call(s, uu, t, False)
        # ---- traffic events ----
        while ev_ptr < len(ev) and ev[ev_ptr][0] <= k:
            _, u, typ, hold = ev[ev_ptr]; ev_ptr += 1
            for s in range(S):
                if in_call[s, u]:
                    if measuring(t): M[schemes[s]]["busy"] += 1
                    continue
                if typ == 1: page(s, u, t)
                else:
                    contact(s, u, t, camp[s, u], lc_cell[s, u], lc_t[s, u])
                    lc_cell[s, u] = camp[s, u]; lc_t[s, u] = t
                start_call(s, u, t, hold)
        if verbose and k % 1200 == 0:
            print(f"  seed {seed} t={t:.0f}s  {time.time()-t0:.1f}s", flush=True)

    uh = N*P["t_meas"]/3600.0
    out = {}
    for n in schemes:
        m = M[n]
        r = dict(m)
        r["user_hours"] = uh
        r["lau_per_uh"] = m["lau"]/uh
        r["ilr_per_uh"] = m["ilr"]/uh
        r["page_cells_per_call"] = m["page_cells"]/max(m["pages"], 1)
        r["page_miss_rate"] = m["page_miss"]/max(m["pages"], 1)
        r["page_delay_extra"] = r["page_miss_rate"]*P["t3113"]
        r["lm_cost_per_uh"] = (P["w_lau"]*m["lau"] + P["w_ilr"]*m["ilr"] + P["w_page"]*m["page_cells"])/uh
        cm = m["call_time"]/60.0
        r["ho_per_callmin"] = m["ho"]/max(cm, 1e-9)
        r["pp_ratio"] = m["pp"]/max(m["ho"], 1)
        r["drop_rate"] = m["drops"]/max(m["calls"], 1)
        r["poorq_frac"] = m["poorq_time"]/max(m["call_time"], 1e-9)
        out[n] = r
    return out


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    t0 = time.time()
    res = run(seed, dict(n_users=200, t_warm=300, t_meas=1200), verbose=True)
    for n, r in res.items():
        print(f"{n:11s} idle={r['lau_idle']} post={r['lau_post']} per={r['lau_per']} LAU/uh={r['lau_per_uh']:.2f} ILR/uh={r['ilr_per_uh']:.2f} pg/call={r['page_cells_per_call']:.2f} "
              f"miss={r['page_miss_rate']:.2f} cost/uh={r['lm_cost_per_uh']:.1f} HO/min={r['ho_per_callmin']:.3f} "
              f"pp={r['pp_ratio']:.3f} drop={r['drop_rate']:.3f} poorq={r['poorq_frac']:.3f} calls={r['calls']}")
    print("elapsed", time.time() - t0)
