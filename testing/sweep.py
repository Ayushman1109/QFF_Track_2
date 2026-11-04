import sys, io, contextlib, statistics, importlib
import os; _h = os.path.dirname(os.path.abspath(__file__)); _t = os.path.join(_h, "..", "testing"); sys.path[:0] = [os.path.join(_t, "mock"), _t, os.path.join(_h, ".."), _h]
oc, jc, seeds, noise = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
sys.argv = ["x", "0"]
with contextlib.redirect_stdout(io.StringIO()):
    import test2
import my_solver as M
M.ORDER_CHECKS, M.JOINT_CHECKS = oc, jc
from silentgroup import Oracle
scores, bad = [], 0
for kind in test2.kinds:
    for seed in range(seeds):
        inst = test2.gen(kind, seed); orc = Oracle(inst, noise=noise, seed=seed + 1000)
        res = M.solve(orc)
        good = res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"]
        bad += not good
        scores.append((1 + 0.5 * orc.budget_left() / 80000) if good else (0.3 if res["order"] == inst.answer["order"] else 0))
print(f"order_checks={oc} joint_checks={jc} noise={noise}: mean score {statistics.mean(scores):.4f}, wrong {bad}/{len(scores)}")
