import sys, time, random, math, statistics
import os; _h = os.path.dirname(os.path.abspath(__file__)); _t = os.path.join(_h, "..", "testing"); sys.path[:0] = [os.path.join(_t, "mock"), _t, os.path.join(_h, ".."), _h]
from silentgroup import Group, Oracle, make_instance
import my_solver

P59 = 576460752303423433

def custom(moduli, g, exps, forged_elems, n):
    pw = lambda e: tuple(e * x % d for x, d in zip(g, moduli))
    els = {"g": g}
    for X in "ABC":
        els[X] = forged_elems[X] if X in forged_elems else pw(exps[X])
    secret = "none" if forged_elems else exps["A"] * exps["B"] * exps["C"] % n
    return Group(moduli, els, {"order": n, "secret": secret})

def gen(kind, seed):
    rng = random.Random(seed)
    if kind == "cyclic 2^59":
        return make_instance([2**59], seed)
    if kind == "Zp x Z2, forged A":
        n = P59; g = (rng.randrange(1, P59), 0)
        ex = {X: rng.randrange(n) for X in "BC"}
        return custom([P59, 2], g, ex, {"A": (rng.randrange(P59), 1)}, n)
    if kind == "Zp x Z2 genuine":
        n = P59; g = (rng.randrange(1, P59), 0)
        return custom([P59, 2], g, {X: rng.randrange(n) for X in "ABC"}, {}, n)
    if kind == "3^37 elem. abelian":
        mods = [3] * 37
        g = tuple([1] + [rng.randrange(3) for _ in range(36)])
        ex = {X: rng.randrange(3) for X in "ABC"}
        if seed % 2:
            f = list(g); f[1] = (f[1] + 1) % 3
            return custom(mods, g, ex, {"B": tuple(f)}, 3)
        return custom(mods, g, ex, {}, 3)
    if kind == "g = e":
        mods = [2**20, 3**10]
        g = (0, 0)
        if seed % 2:
            return custom(mods, g, {"A": 0, "B": 0, "C": 0}, {"C": (5, 0)}, 1)
        return custom(mods, g, {"A": 0, "B": 0, "C": 0}, {}, 1)
    if kind == "big forged ABC":
        n = 2**40; g = (2**18 * rng.randrange(1, 2**18) * 2 + 2**18, 0)
        g = (2**18, 0);
        return custom([2**58, 4], g, {}, {"A": (1, 0), "B": (3, 1), "C": (2**18, 1)}, 2**40)
    if kind == "ZnxZn on-grid forgery":
        n = 2**29 * 3 if seed % 2 else 999999937
        g = (1, 0)
        ex = {X: rng.randrange(n) for X in "ABC"}
        return custom([n, n], g, ex, {"B": (ex["B"], rng.randrange(1, n))}, n)
    if kind == "Z_n1 x Z_n2 coprime-ish":
        return make_instance([1000003 * 7**4, 999983 * 2**3], seed)

kinds = ["cyclic 2^59", "Zp x Z2, forged A", "Zp x Z2 genuine", "3^37 elem. abelian",
         "g = e", "big forged ABC", "ZnxZn on-grid forgery", "Z_n1 x Z_n2 coprime-ish"]
if __name__ == '__main__':
    T0 = time.time(); scores = []
    for kind in kinds:
        ok = 0; lefts = []
        for seed in range(int(sys.argv[1]) if len(sys.argv) > 1 else 10):
            inst = gen(kind, seed)
            orc = Oracle(inst, noise=float(sys.argv[2]) if len(sys.argv) > 2 else 1.0, seed=seed + 1000)
            res = my_solver.solve(orc)
            good = res == inst.answer or (res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"])
            sc = (1 + 0.5 * orc.budget_left() / 80000) if good else (0.3 if res["order"] == inst.answer["order"] else 0)
            scores.append(sc); ok += good; lefts.append(orc.budget_left() / 80000)
            if not good: print("   BAD", kind, seed, res, inst.answer)
        print(f"{kind:26s} {ok}/{len(lefts)} correct, mean budget left {statistics.mean(lefts):.2f}")
    print(f"mean per-meeting score {statistics.mean(scores):.3f} (max 1.5), time {time.time()-T0:.0f}s")
