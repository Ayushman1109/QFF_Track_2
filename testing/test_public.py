import os, sys, time, statistics
_h = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(_h, "mock"), os.path.join(_h, ".."), _h]

from silentgroup import Oracle, load_public
import my_solver

repeats = int(sys.argv[1]) if len(sys.argv) > 1 else 30
noise = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0

all_scores = []
for name in ("public_1", "public_2"):
    inst = load_public(name)
    scores, queries, lefts, ok = [], [], [], 0
    t0 = time.time()
    for seed in range(repeats):
        orc = Oracle(inst, noise=noise, budget=80000, seed=seed + 100)
        res = my_solver.solve(orc)
        good = (res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"])
        pts = (1 + 0.5 * orc.budget_left() / 80000) if good else (0.3 if res["order"] == inst.answer["order"] else 0.0)
        ok += good
        scores.append(pts)
        all_scores.append(pts)
        queries.append(orc.queries_used)
        lefts.append(orc.budget_left() / 80000)
        if not good:
            print(f"  BAD {name} seed={seed}: got={res} expected={inst.answer}")
    print(f"{name:10s}: {ok}/{repeats} correct, mean score {statistics.mean(scores):.4f} "
          f"(mean queries {statistics.mean(queries):.0f}, max {max(queries)}, "
          f"mean budget left {statistics.mean(lefts):.2%}, {time.time()-t0:.1f}s)")
print(f"Overall public mean score: {statistics.mean(all_scores):.4f} / 1.5000")
