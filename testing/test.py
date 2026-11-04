import sys, time, random
T0 = time.time()
import os; _h = os.path.dirname(os.path.abspath(__file__)); _t = os.path.join(_h, "..", "testing"); sys.path[:0] = [os.path.join(_t, "mock"), _t, os.path.join(_h, ".."), _h]
from silentgroup import make_instance, Oracle
import my_solver

P61 = (1 << 61) - 1
P59 = 576460752303423433
cases = [
    ("cyclic 2^20 x 6", [2**20, 6], ""),
    ("prime ~2^59",     [P59], ""),
    ("2^30 x 2^29",     [2**30, 2**29], ""),
    ("smooth ~2^58",    [2**10*3**8*5**6*7**5*11**2], ""),
    ("mixed",           [1000003*999983, 12, 8], ""),
    ("forged C small",  [2**20, 6], "C"),
    ("forged AB",       [2**16, 2**16], "AB"),
    ("tiny Z_6 x Z_6",  [6, 6], ""),
]
assert my_solver._is_probable_prime(P59)
tot = ok = 0
for name, mods, forged in cases:
    for seed in range(25):
        inst = make_instance(mods, seed, forged)
        orc = Oracle(inst, noise=1.0, seed=seed + 100)
        t = time.time()
        res = my_solver.solve(orc)
        good = res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"]
        tot += 1; ok += good
        print(f"{name:16s} s{seed} {'OK ' if good else 'BAD'} used={orc.queries_used:6d} "
              f"left={orc.budget_left()/80000:.2f} t={time.time()-t:5.1f}s  got={res} want={inst.answer}" if not good else
              f"{name:16s} s{seed} OK  used={orc.queries_used:6d} left={orc.budget_left()/80000:.2f} t={time.time()-t:5.1f}s")
print(ok, "/", tot, "total time", round(time.time()-T0), "s")
