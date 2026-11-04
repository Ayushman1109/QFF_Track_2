import sys, io, contextlib, math, random, statistics, json, os
_h = os.path.dirname(os.path.abspath(__file__))
_t = os.path.join(_h, "..", "testing")
sys.path[:0] = [os.path.join(_t, "mock"), _t, os.path.join(_h, ".."), _h]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    import test2
import my_solver as M
import subgroup_structure as S_sub
from silentgroup import Group, Oracle, make_instance, load_public

EPS = (math.pi * M.NOISE_SIGMA) ** 2
tri = lambda a: a * (a + 1) // 2


def model(L, sims=20000, rng=random.Random(1)):
    m_o = 2 * 61 + M.ORDER_EXTRA_BITS


    P_o = math.exp(-0.68 * EPS * tri(61 + L)) * 0.97
    shots_ne, shots_e = 6, 15
    E_order = 6 + tri(m_o) / P_o + (1 / P_o - 1) * shots_ne * 61 + shots_e * L
    chk = M.JOINT_CHECKS * max(1, L - 3)
    C_g, C_x = tri(L + M.G_EXTRA_BITS) + chk, tri(L + M.EXTRA_BITS) + chk
    P_acc = math.exp(-EPS * tri(L + 1)) * 0.97
    tot = 0
    for _ in range(sims):
        todo, cost = 3, 0
        while todo:
            cost += C_g + todo * C_x
            if rng.random() < P_acc:
                todo -= sum(rng.random() < P_acc for _ in range(todo))
        tot += cost
    E_logs = tot / sims + 3 * shots_e * (L + 2)
    floor = 6 + tri(m_o) + shots_e * L + C_g + 3 * C_x + 3 * shots_e * (L + 2)
    return E_order, E_logs, floor

def random_prime(bits, rng):
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if M._is_probable_prime(p):
            return p

def measure(L, seeds):
    rng = random.Random(L)
    o_q, l_q = [], []
    for s in range(seeds):
        p = random_prime(L, rng)
        g = (rng.randrange(1, p),)
        ex = {X: rng.randrange(p) for X in "ABC"}
        inst = test2.custom([p], g, ex, {}, p)
        orc = Oracle(inst, noise=1.0, seed=10_000 + 97 * L + s)
        S = M.Solver(orc)
        n = S.find_order(); q_order = orc.queries_used
        S.solve_logs(n); o_q.append(q_order); l_q.append(orc.queries_used - q_order)
    return o_q, l_q

Ls = [8, 16, 24, 32, 40, 48, 56, 60]
rows = []
for L in Ls:
    o_q, l_q = measure(L, 60)
    Eo, El, fl = model(L)
    tq = [a + b for a, b in zip(o_q, l_q)]
    meas_mean = statistics.mean(tq)
    rel_err = (Eo + El - meas_mean) / meas_mean
    rows.append(dict(L=L, meas=meas_mean, sem=statistics.stdev(tq) / math.sqrt(len(tq)),
                     meas_order=statistics.mean(o_q), meas_logs=statistics.mean(l_q),
                     model=Eo + El, model_order=Eo, model_logs=El, floor=fl, rel_err=rel_err))
    print(f"L={L:2d}  measured {meas_mean:7.0f} ± {rows[-1]['sem']:4.0f}   model {Eo+El:7.0f} ({rel_err:+.1%})   "
          f"(order {statistics.mean(o_q):6.0f} vs {Eo:6.0f}, logs {statistics.mean(l_q):6.0f} vs {El:6.0f})")


noise_levels = [0, 0.5, 1, 1.5, 2, 2.5, 3]
sweep = []
for nz in noise_levels:
    sc, ok, left = [], 0, []
    for kind in test2.kinds:
        for seed in range(20):
            inst = test2.gen(kind, seed); orc = Oracle(inst, noise=nz, seed=seed + 5000)
            res = M.solve(orc)
            good = res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"]
            ok += good; left.append(orc.budget_left() / 80000)
            sc.append((1 + 0.5 * orc.budget_left() / 80000) if good else (0.3 if res["order"] == inst.answer["order"] else 0))
    sweep.append(dict(noise=nz, score=statistics.mean(sc), success=ok / len(sc), left=statistics.mean(left)))
    print(f"noise {nz:3.1f}x: success {ok/len(sc):.4f}  mean score {statistics.mean(sc):.4f}  budget left {statistics.mean(left):.4f}")


public_bench = {}
for pub_name in ("public_1", "public_2"):
    inst = load_public(pub_name)
    sc, qs, ok = [], [], 0
    for seed in range(50):
        orc = Oracle(inst, noise=1.0, seed=seed + 100)
        res = M.solve(orc)
        good = (res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"])
        ok += good; qs.append(orc.queries_used)
        sc.append((1 + 0.5 * orc.budget_left() / 80000) if good else (0.3 if res["order"] == inst.answer["order"] else 0))
    public_bench[pub_name] = dict(
        success=ok / 50,
        mean_score=statistics.mean(sc),
        mean_queries=statistics.mean(qs),
        max_queries=max(qs),
        mean_budget_left=1 - statistics.mean(qs) / 80000,
    )
    print(f"{pub_name}: {ok}/50 correct, score={statistics.mean(sc):.4f}, mean_q={statistics.mean(qs):.0f}, max_q={max(qs)}")


structure_rows = []
for name in ("public_1", "public_2"):
    inst = load_public(name)
    truth = S_sub.true_subgroup_structure(inst)
    qs, matches_80k, matches_unlim = [], 0, 0
    for s in range(15):
        orc = Oracle(inst, noise=1.0, budget=80000, seed=42 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        matches_80k += (got["invariants"] == truth)
        qs.append(orc.queries_used)
    structure_rows.append(dict(
        case=name, invariants=truth, rank=len([x for x in truth if x > 1]),
        rate_80k=matches_80k / 15, mean_queries=statistics.mean(qs), max_queries=max(qs)
    ))
    print(f"Subgroup {name:24s}: 80k recovery={matches_80k}/15 ({matches_80k/15:.1%}), mean_q={statistics.mean(qs):.0f}, truth={truth}")

for kind in test2.kinds:
    inst0 = test2.gen(kind, 1)
    truth0 = S_sub.true_subgroup_structure(inst0)
    qs, matches_80k = [], 0
    for s in range(15):

        inst = test2.gen(kind, 2 * s + 1 if kind in ("3^37 elem. abelian", "g = e") else s)
        truth = S_sub.true_subgroup_structure(inst)
        orc = Oracle(inst, noise=1.0, budget=80000, seed=42 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        matches_80k += (got["invariants"] == truth)
        qs.append(orc.queries_used)
    structure_rows.append(dict(
        case=kind, invariants=truth0, rank=len([x for x in truth0 if x > 1]),
        rate_80k=matches_80k / 15, mean_queries=statistics.mean(qs), max_queries=max(qs)
    ))
    print(f"Subgroup {kind:24s}: 80k recovery={matches_80k}/15 ({matches_80k/15:.1%}), mean_q={statistics.mean(qs):.0f}, truth={truth0}")


budgets = [15000, 25000, 35000, 45000, 55000, 65000, 80000, 100000]
budget_curves = {
    "Genuine cyclic (public_1, Zp×Z2)": [],
    "On-grid rank-2 forgery (public_2, Zn×Zn)": [],
    "Single off-grid forgery (Zp×Z2 forged A)": [],
    "Multi 58-bit forgery (big forged ABC)": [],
}


for B in budgets:

    ok_g = 0
    for s in range(15):
        inst = load_public("public_1") if s % 2 == 0 else test2.gen("Zp x Z2 genuine", s)
        truth = S_sub.true_subgroup_structure(inst)
        orc = Oracle(inst, noise=1.0, budget=B, seed=100 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        ok_g += (got["invariants"] == truth)
    budget_curves["Genuine cyclic (public_1, Zp×Z2)"].append(ok_g / 15)


    ok_on = 0
    for s in range(15):
        inst = load_public("public_2") if s % 2 == 0 else test2.gen("ZnxZn on-grid forgery", s)
        truth = S_sub.true_subgroup_structure(inst)
        orc = Oracle(inst, noise=1.0, budget=B, seed=200 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        ok_on += (got["invariants"] == truth)
    budget_curves["On-grid rank-2 forgery (public_2, Zn×Zn)"].append(ok_on / 15)


    ok_off1 = 0
    for s in range(15):
        inst = test2.gen("Zp x Z2, forged A", s)
        truth = S_sub.true_subgroup_structure(inst)
        orc = Oracle(inst, noise=1.0, budget=B, seed=300 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        ok_off1 += (got["invariants"] == truth)
    budget_curves["Single off-grid forgery (Zp×Z2 forged A)"].append(ok_off1 / 15)


    ok_big = 0
    for s in range(15):
        inst = test2.gen("big forged ABC", s)
        truth = S_sub.true_subgroup_structure(inst)
        orc = Oracle(inst, noise=1.0, budget=B, seed=400 + s)
        got = S_sub.find_subgroup_structure(orc, target_shots=6)
        ok_big += (got["invariants"] == truth)
    budget_curves["Multi 58-bit forgery (big forged ABC)"].append(ok_big / 15)


shot_counts = [1, 2, 3, 4, 5, 6, 7, 8]
shot_study = []
for T in shot_counts:

    ok_p2, ok_p3 = 0, 0
    trials = 40
    for s in range(trials):

        inst2 = test2.gen("ZnxZn on-grid forgery", 2 * s)
        truth2 = S_sub.true_subgroup_structure(inst2)
        orc2 = Oracle(inst2, noise=0.0, budget=None, seed=500 + s)

        rng2 = random.Random(500 + s)
        mat2 = []
        for _ in range(T):
            chi = [rng2.randrange(d) for d in inst2.moduli]
            row = [sum((c * xi) % d for c, xi, d in zip(chi, inst2.elements[el], inst2.moduli)) % (2**29) for el in ("g", "B")]
            mat2.append(row)
        inv2 = S_sub._p_adic_snf_invariants(mat2, 2, 29)
        ok_p2 += (inv2 == [2**29, 2**29])


        inst3 = test2.gen("3^37 elem. abelian", 2 * s + 1)
        rng3 = random.Random(600 + s)
        mat3 = []
        for _ in range(T):
            chi = [rng3.randrange(3) for _ in inst3.moduli]
            row = [sum(c * xi for c, xi in zip(chi, inst3.elements[el])) % 3 for el in ("g", "B")]
            mat3.append(row)
        inv3 = S_sub._p_adic_snf_invariants(mat3, 3, 1)
        ok_p3 += (inv3 == [3, 3])

    shot_study.append(dict(
        T=T,
        meas_p2=ok_p2 / trials,
        theo_p2=S_sub.theoretical_span_prob(T, 2, 2),
        meas_p3=ok_p3 / trials,
        theo_p3=S_sub.theoretical_span_prob(T, 2, 3),
    ))

with open(os.path.join(_h, "analysis_results.json"), "w") as f:
    json.dump(dict(
        cost=rows,
        noise=sweep,
        public=public_bench,
        subgroup_structure=structure_rows,
        subgroup_budget_curves=dict(budgets=budgets, curves=budget_curves),
        subgroup_shot_study=shot_study,
    ), f, indent=1)


fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
x = [r["L"] for r in rows]
ax[0].errorbar(x, [r["meas"] for r in rows], yerr=[2 * r["sem"] for r in rows], fmt="o", color="#2563eb", label="measured total (±2 s.e.)")
ax[0].plot(x, [r["model"] for r in rows], "-", color="#2563eb", lw=1.8, alpha=.8, label="model total")
ax[0].plot(x, [r["meas_order"] for r in rows], "s", color="#ea580c", ms=4.5, label="measured: order finding")
ax[0].plot(x, [r["model_order"] for r in rows], "--", color="#ea580c", lw=1.5, alpha=.8, label="model: order finding")
ax[0].plot(x, [r["meas_logs"] for r in rows], "^", color="#16a34a", ms=4.5, label="measured: discrete logs")
ax[0].plot(x, [r["model_logs"] for r in rows], "--", color="#16a34a", lw=1.5, alpha=.8, label="model: discrete logs")
ax[0].plot(x, [r["floor"] for r in rows], ":", color="#64748b", lw=1.4, label="noise-free floor")
ax[0].axhline(80000, color="#dc2626", lw=1.0); ax[0].text(46, 76500, "budget 80k", color="#dc2626", fontsize=8.5, fontweight="bold")
ax[0].set_xlabel("L = bit length of ord(g)  (cyclic, prime order)")
ax[0].set_ylabel("queries per meeting")
ax[0].set_title("Query cost: analytical model vs measured", fontweight="bold")
ax[0].legend(fontsize=7.5, loc="upper left")
ax[0].set_ylim(0, 85000)
ax[0].grid(True, alpha=0.2)

nz = [s["noise"] for s in sweep]
ax[1].plot(nz, [s["success"] for s in sweep], "o-", color="#4f46e5", lw=2.0, label="fraction fully correct")
ax[1].plot(nz, [s["score"] / 1.5 for s in sweep], "s-", color="#0284c7", lw=1.8, label="mean score / max (1.5)")
ax[1].plot(nz, [s["left"] for s in sweep], "^-", color="#ea580c", lw=1.8, label="mean budget left")
ax[1].axvline(1, color="#64748b", lw=1.0, ls="--"); ax[1].text(1.04, 0.08, "stated noise (1.0×)", fontsize=8.5, color="#475569")
ax[1].set_xlabel("noise multiplier (× 0.002·√cost turns)")
ax[1].set_ylabel("fraction / normalized metric")
ax[1].set_ylim(0, 1.05)
ax[1].set_title("Degradation with noise (8 group shapes × 20 seeds)", fontweight="bold")
ax[1].legend(fontsize=8, loc="lower left")
ax[1].grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(_h, "open_track_analysis.png"), dpi=180)
plt.close(fig)


fig2, ax2 = plt.subplots(1, 2, figsize=(11.5, 4.4))
b_k = [b / 1000 for b in budgets]
colors = ["#16a34a", "#2563eb", "#ea580c", "#dc2626"]
markers = ["o", "s", "^", "D"]
for (label, vals), col, mk in zip(budget_curves.items(), colors, markers):
    ax2[0].plot(b_k, [v * 100 for v in vals], f"{mk}-", color=col, lw=1.9, ms=5, label=label)
ax2[0].axvline(80, color="#475569", lw=1.2, ls="--")
ax2[0].text(81.5, 78, "80k limit", fontsize=8.5, color="#334155", fontweight="bold")
ax2[0].set_xlabel("query budget (thousands of queries)")
ax2[0].set_ylabel("exact structure recovery (%)")
ax2[0].set_title("Subgroup structure recovery vs. query budget", fontweight="bold")
ax2[0].set_ylim(-3, 105)
ax2[0].legend(fontsize=7.5, loc="lower right")
ax2[0].grid(True, alpha=0.2)

Ts = [r["T"] for r in shot_study]
ax2[1].plot(Ts, [r["theo_p2"] * 100 for r in shot_study], "--", color="#2563eb", lw=1.8, label="theory rank-2 (p=2): ∏(1 − 2⁻ᵏ)")
ax2[1].plot(Ts, [r["meas_p2"] * 100 for r in shot_study], "o", color="#2563eb", ms=6, label="simulated rank-2 (p=2, 40 trials reused)")
ax2[1].plot(Ts, [r["theo_p3"] * 100 for r in shot_study], "--", color="#ea580c", lw=1.8, label="theory rank-2 (p=3): ∏(1 − 3⁻ᵏ)")
ax2[1].plot(Ts, [r["meas_p3"] * 100 for r in shot_study], "s", color="#ea580c", ms=6, label="simulated rank-2 (p=3, 40 trials reused)")
ax2[1].axvline(6, color="#16a34a", lw=1.0, ls=":")
ax2[1].text(6.1, 20, "T=6 target\n(>95% span)", fontsize=8, color="#15803d")
ax2[1].set_xlabel("clean joint character shots T")
ax2[1].set_ylabel("full rank-2 span recovery (%)")
ax2[1].set_title("Dual-group span probability vs. shots T", fontweight="bold")
ax2[1].set_ylim(-3, 105)
ax2[1].legend(fontsize=7.5, loc="lower right")
ax2[1].grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(_h, "subgroup_structure_analysis.png"), dpi=180)
plt.close(fig2)
print("saved open_track_analysis.png, subgroup_structure_analysis.png, and analysis_results.json")
