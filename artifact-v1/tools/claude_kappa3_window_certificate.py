"""Chantier 1: rigorous window certificate for the cumulant floor
    kappa_3(theta) >= -1/theta^2
on a theta-window, using python-flint (arb/acb ball arithmetic + acb_calc
rigorous integration). Claude branch, 2026-07-10.

Objects:
  V(s) = e^s Phi(e^s),  Phi(u) = 2 sum_m (2 pi^2 m^4 e^{9u/2} - 3 pi m^2 e^{5u/2}) e^{-pi m^2 e^{2u}}
  I_k(theta; c) = int (s-c)^k e^{theta s} V(s) ds   (k = 0..4, centered at c)
  kappa_3 = a3 - 3 a2 a1 + 2 a1^3,  a_k = I_k/I_0   (exact for any center c)

Certification scheme on a segment [ta, tb]:
  1. kappa_3(ta) enclosed by rigorous integrals (truncation tails bounded by
     elementary closed forms; m-tail of Phi bounded by 4 pi^2 e^{4.5 Re u} * 1e-300).
  2. Lipschitz: |kappa_3(t) - kappa_3(ta)| <= 3 * MU4BOX * (tb-ta), where
     MU4BOX >= sup_box mu_4 via the monotone theta-split bound
       mu_4(t) <= [S+(s>=0 at tb) + S-(s<0 at ta)] / [S+(s>=0 at ta) + S-(s<0 at tb)]  (numerators with (s-c)^4).
  3. Accept segment iff lower(kappa_3(ta)) - 3*MU4BOX*(tb-ta) >= -1/tb^2.

All comparisons are certified ball comparisons (True only if provable).
"""
import sys, math, json, time
from flint import arb, acb, ctx

ctx.prec = 192

PI = arb.pi()
MTERMS = 15
# elementary m-tail bound constant: sum_{m>15} m^4 e^{-pi m^2} < 1e-300 (first term e^{-804})
MTAIL = arb("1e-295")

def phi_ball(u):
    s = acb(0)
    for m in range(1, MTERMS + 1):
        mm = m * m
        s += (2 * PI**2 * arb(m)**4 * (u * arb(9) / 2).exp()
              - 3 * PI * arb(mm) * (u * arb(5) / 2).exp()) * (-PI * arb(mm) * (2 * u).exp()).exp()
    s = 2 * s
    # m-tail: |tail| <= 4 pi^2 e^{4.5 Re u} * 1e-295  (valid for Re e^{2u} >= 1)
    err = (4 * PI**2 * ((u.real + u.real.rad()) * arb(9) / 2).exp() * MTAIL)
    err_hi = err.mid() + err.rad()
    return s + acb(arb(0, err_hi), 0)

def integrand(s, c, theta, k):
    return (s - c)**k * (theta * s).exp() * s.exp() * phi_ball(s.exp())

def u_star(nu):
    x = max(0.05, 0.5 * math.log(max(nu, 2.0) / math.pi))
    for _ in range(60):
        x -= (2 * x + math.log(x) - math.log(nu / math.pi)) / (2 + 1 / x)
        if x <= 0:
            x = 1e-4
    return x

def window_params(theta):
    us = u_star(theta / 2.0)
    c = math.log(us)
    sig = 1.0 / math.sqrt(theta * (2 * us + 1))
    s_lo = min(-1.0, c - 30 * sig)
    # adaptive minimal s_hi: extend just far enough that the elementary right
    # tail bound is valid (pi E > B+1 for k<=4); keeps acb_calc out of the
    # deep double-exponential collapse zone.
    margin = 0.15
    while margin < 2.5:
        s_hi = c + 25 * sig + margin
        if right_tail_bound(theta, c, s_hi, 4) is not None:
            umax = math.exp(s_hi)
            E = math.exp(2 * umax)
            B = theta + 5 + 4.5 * umax
            if math.pi * E > 2 * B + 2:
                break
        margin += 0.05
    return c, sig, s_lo, s_hi

def right_tail_bound(theta, c, s_hi, k):
    """int_{s_hi}^inf (s-c)^k e^{(theta+1)s} Phi(e^s) ds
       <= e^A / (pi E - B) with E=e^{2 u_max}, u_max=e^{s_hi},
       A=(theta+1)s_hi + k log(1+|s_hi-c|) + log(4 pi^2 C0) + 4.5 u_max - pi E,
       valid when pi E > B+1, B=(theta+1)+k+4.5*u_max  (uses e^{e^t growth} >= linear)."""
    umax = arb(s_hi).exp()
    E = (2 * umax).exp()
    B = arb(theta + 1 + k) + arb(9) / 2 * umax
    ok = (PI * E > B + 1)
    if not ok:
        return None
    C0 = arb("0.05")  # sum m^4 e^{-pi m^2} < 0.0433; use 0.05
    A = arb(theta + 1) * arb(s_hi) + arb(k) * arb(1 + abs(s_hi - c)).log() \
        + (4 * PI**2 * C0).log() + arb(9) / 2 * umax - PI * E
    return A.exp() / (PI * E - B)

def left_tail_bound(theta, c, s_lo, k):
    """int_{-inf}^{s_lo} (c-s)^k e^{(theta+1)s} Phi(e^s) ds with Phi <= 1.71 e^{4.5 e^{s_lo}}:
       closed form e^{beta a} sum_j (k!/(k-j)!) (c-a)^{k-j} / beta^{j+1}."""
    beta = arb(theta + 1)
    a = arb(s_lo)
    pref = arb("1.71") * (arb(9) / 2 * arb(s_lo).exp()).exp()
    tot = arb(0)
    fac = arb(1)
    for j in range(0, k + 1):
        term = fac * (arb(c) - a)**(k - j) / beta**(j + 1)
        tot += term
        fac *= (k - j)
    return pref * (beta * a).exp() * tot

def rig_integral(theta, c, k, s_lo, s_hi, half=None, sig=None):
    """Rigorous integral of (s-c)^k e^{theta s} V(s) over [s_lo, s_hi] (or a half).

    Manual piecewise splitting near the peak: acb_calc's error hull inflates on
    double-exponential integrands when a single call spans the whole range
    (complex excursions hit exploding e^{-pi e^{2 e^s}}); sigma-sized pieces
    keep the complex disks tame and the enclosures tight."""
    ta = arb(theta)
    ca = arb(c)
    lo, hi = s_lo, s_hi
    if half == "pos":
        lo = max(0.0, s_lo)
    elif half == "neg":
        hi = min(0.0, s_hi)
        if hi <= lo:
            return arb(0)
    if sig is None:
        sig = 0.05
    # Architecture: acb_calc rigorous integration ONLY on the core
    # [c-12sig, c+12sig] (well-conditioned there); outside the core the true
    # contribution is negligible (Gaussian factor e^{-72} at 12 sigma), so
    # direct interval-hull bounds |int| <= sup|f| * width suffice and are
    # deterministic -- acb_calc is never exposed to the double-exponential
    # collapse zone where its ellipse bounds explode.
    core_lo = max(lo, c - 12 * sig)
    core_hi = min(hi, c + 12 * sig)
    f = lambda s, _: integrand(s, ca, ta, k)
    tot = arb(0)
    if core_hi > core_lo:
        ncore = 24
        for i in range(ncore):
            a_ = core_lo + (core_hi - core_lo) * i / ncore
            b_ = core_lo + (core_hi - core_lo) * (i + 1) / ncore
            tot += acb.integral(f, arb(a_), arb(b_), rel_tol=arb(2)**(-80)).real
    # outside-core: march in sigma/2 steps with interval-hull bounds
    def hull_add(a_, b_):
        sball = arb((a_ + b_) / 2, (b_ - a_) / 2)
        v = integrand(acb(sball), ca, ta, k)
        R = (abs(v.real.mid()) + v.real.rad()) * arb(b_ - a_)
        # pure-arb symmetric enclosure [-R, R]; no float conversion (overflow-safe)
        return arb(0).union(R).union(-R)
    step = max(sig / 2, 1e-6)
    x = core_lo
    while x > lo:
        a_ = max(lo, x - step)
        tot += hull_add(a_, x)
        x = a_
    x = core_hi
    while x < hi:
        b_ = min(hi, x + step)
        tot += hull_add(x, b_)
        x = b_
    return tot

def kappa3_ball(theta):
    c, sig, s_lo, s_hi = window_params(theta)
    rt = [right_tail_bound(theta, c, s_hi, k) for k in range(4)]
    lt = [left_tail_bound(theta, c, s_lo, k) for k in range(4)]
    if any(r is None for r in rt):
        return None, c
    I = []
    for k in range(4):
        v = rig_integral(theta, c, k, s_lo, s_hi, sig=sig)
        v = v.union(v + rt[k] + lt[k]).union(v - rt[k] - lt[k])
        I.append(v)
    a1 = I[1] / I[0]
    a2 = I[2] / I[0]
    a3 = I[3] / I[0]
    k3 = a3 - 3 * a2 * a1 + 2 * a1**3
    return k3, c

def mu4_box_bound(ta_, tb_, c, s_lo, s_hi, sig=None):
    """sup over theta in [ta,tb] of mu_4 <= N/D via monotone split at s=0."""
    Np = rig_integral(tb_, c, 4, s_lo, s_hi, half="pos", sig=sig)
    Nn = rig_integral(ta_, c, 4, s_lo, s_hi, half="neg", sig=sig)
    Dp = rig_integral(ta_, c, 0, s_lo, s_hi, half="pos", sig=sig)
    Dn = rig_integral(tb_, c, 0, s_lo, s_hi, half="neg", sig=sig)
    rtN = right_tail_bound(tb_, c, s_hi, 4)
    ltN = left_tail_bound(tb_, c, s_lo, 4)
    N = Np + Nn + rtN + ltN
    D = Dp + Dn
    return N / D

def certify_window(theta_min, theta_max, h0=1.0, progress_path=None):
    t = theta_min
    segments = []
    t0 = time.time()
    plog = open(progress_path, "a") if progress_path else None
    while t < theta_max:
        k3, c = kappa3_ball(t)
        if k3 is None:
            return {"ok": False, "fail_at": t, "reason": "tail condition"}
        _, sig, s_lo, s_hi = window_params(t)
        h = h0
        placed = False
        for _ in range(8):
            tb_ = min(t + h, theta_max)
            mu4b = mu4_box_bound(t, tb_, c, s_lo, s_hi, sig=sig)
            lip = 3 * mu4b * arb(tb_ - t)
            floor = -1 / arb(tb_)**2
            lower_k3 = k3 - lip
            if lower_k3 > floor:
                margin = (lower_k3 - floor) * arb(tb_)**2
                seg = {
                    "theta": [t, tb_],
                    "kappa3_mid": float(k3.mid()),
                    "kappa3_t2": float((k3 * arb(t)**2).mid()),
                    "lip": float(lip.mid() + lip.rad()),
                    "margin_rel": float(margin.mid()),
                }
                segments.append(seg)
                if plog:
                    plog.write(json.dumps(seg) + "\n")
                    plog.flush()
                t = tb_
                placed = True
                break
            h /= 2
        if not placed:
            return {"ok": False, "fail_at": t, "reason": "no step placed", "segments": len(segments)}
    return {"ok": True, "segments": segments, "count": len(segments), "secs": round(time.time() - t0, 1)}

if __name__ == "__main__":
    tmin = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 160.0
    h0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    ppath = sys.argv[4] if len(sys.argv) > 4 else None
    res = certify_window(tmin, tmax, h0, ppath)
    if res.get("ok"):
        segs = res["segments"]
        print(json.dumps({
            "ok": True, "window": [tmin, tmax], "count": res["count"], "secs": res["secs"],
            "first": segs[0], "last": segs[-1],
            "min_margin_rel": min(s["margin_rel"] for s in segs),
        }, indent=1))
    else:
        print(json.dumps(res, indent=1))
