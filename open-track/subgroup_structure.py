import itertools
import json
import math
import os
import sys
from fractions import Fraction

_h = os.path.dirname(os.path.abspath(__file__))
_t = os.path.join(_h, "..", "testing")
sys.path[:0] = [os.path.join(_t, "mock"), _t, os.path.join(_h, ".."), _h]

import my_solver as M
from silentgroup import Group, Oracle, load_public


def theoretical_span_prob(T, r, p):
    if r <= 0:
        return 1.0
    if T < r:
        return 0.0
    prob = 1.0
    for i in range(r):
        prob *= 1.0 - (p ** -(T - i))
    return prob


def _p_adic_snf_invariants(mat, p, max_e, filter_outliers=False):
    Q = p ** max_e
    A = [[x % Q for x in row] for row in mat]
    R, C = len(A), len(A[0])
    if R == 0 or C == 0:
        return []
    invs = []
    r_idx = 0
    for c_idx in range(min(R, C)):
        best_v, best_pos = max_e, None
        for i in range(r_idx, R):
            for j in range(c_idx, C):
                val = A[i][j] % Q
                if val == 0:
                    continue
                v = 0
                while val % p == 0:
                    val //= p
                    v += 1
                if v < best_v:
                    best_v, best_pos = v, (i, j)
        if best_pos is None:
            break
        pi, pj = best_pos
        A[r_idx], A[pi] = A[pi], A[r_idx]
        for i in range(R):
            A[i][c_idx], A[i][pj] = A[i][pj], A[i][c_idx]
        pv = p ** best_v
        u = (A[r_idx][c_idx] // pv) % Q
        u_inv = pow(u, -1, Q)
        for j in range(c_idx, C):
            A[r_idx][j] = (A[r_idx][j] * u_inv) % Q

        qq = p ** (max_e - best_v)
        dep_cols = {}
        if filter_outliers and max_e == 1 and qq > 1000 and R >= 4:
            for j in range(c_idx + 1, C):
                ratios = {}
                for i in range(r_idx, R):
                    v_piv = A[i][c_idx] % Q
                    v_col = A[i][j] % Q
                    if v_piv != 0 and v_piv % pv == 0 and (v_piv // pv) % p != 0 and v_col % pv == 0:
                        r_val = ((v_col // pv) * pow((v_piv // pv) % qq, -1, qq)) % qq
                        ratios[r_val] = ratios.get(r_val, 0) + 1
                if ratios:
                    top_r, top_c = max(ratios.items(), key=lambda t: t[1])
                    if top_c >= 2:
                        dep_cols[j] = top_r
        for i in range(R):
            if i != r_idx and A[i][c_idx] % Q != 0:
                factor = (A[i][c_idx] // pv) % Q
                for j in range(c_idx, C):
                    A[i][j] = (A[i][j] - factor * A[r_idx][j]) % Q
            for j in dep_cols:
                A[i][j] = 0
        invs.append(p ** (max_e - best_v))
        r_idx += 1
    return sorted([x for x in invs if x > 1], reverse=True)


def true_subgroup_structure(inst):
    moduli = inst.moduli
    names = [k for k in ("g", "A", "B", "C") if k in inst.elements]
    ords = [math.lcm(*[d // math.gcd(vi, d) for vi, d in zip(inst.elements[k], moduli)]) for k in names]
    M_exp = math.lcm(*ords)
    if M_exp == 1:
        return [1]
    primes = M._prime_factors(M_exp)
    prime_invs = {}
    for p in primes:
        e_mods = []
        for d in moduli:
            e = 0
            while d % (p ** (e + 1)) == 0:
                e += 1
            e_mods.append(e)
        max_e = max(e_mods)
        if max_e == 0:
            continue
        Q = p ** max_e
        mat = []
        for k in names:
            row = []
            for vi, d, e in zip(inst.elements[k], moduli, e_mods):
                if e == 0:
                    row.append(0)
                else:
                    row.append(((vi % (p ** e)) * (p ** (max_e - e))) % Q)
            mat.append(row)
        prime_invs[p] = _p_adic_snf_invariants(mat, p, max_e)
    max_len = max((len(v) for v in prime_invs.values()), default=0)
    if max_len == 0:
        return [1]
    factors = [1] * max_len
    for p, invs in prime_invs.items():
        for idx, val in enumerate(invs):
            factors[idx] *= val
    return factors


def _minimize_elem_order(solver, elem, c):
    if c > 1 and solver.is_identity([(elem, 1)]):
        return 1
    for p in sorted(M._prime_factors(c), reverse=True):
        e = 0
        tmp = c
        while tmp % p == 0:
            tmp //= p
            e += 1
        if e > 1 and solver.is_identity([(elem, c // (p ** e))]):
            c //= p ** e
            continue
        while c % p == 0:
            rem = 0
            t2 = c
            while t2 % p == 0:
                t2 //= p
                rem += 1
            step = 1
            while step * 2 <= rem:
                step *= 2
            if step > 1 and solver.is_identity([(elem, c // (p ** step))]):
                c //= p ** step
            elif solver.is_identity([(elem, c // p)]):
                c //= p
            else:
                break
    return c


def _find_element_order(solver, elem, known_mult=1):
    if solver.is_identity([(elem, 1)]):
        return 1
    if known_mult > 1:
        for mult in (1, 2, 3, 4):
            cand = known_mult * mult
            if cand <= M.N_MAX and solver.is_identity([(elem, cand)]):
                return _minimize_elem_order(solver, elem, cand)
    base_R = known_mult if (known_mult > 1 and known_mult <= M.N_MAX // 2) else 1
    R = base_R
    tested = {1, known_mult, 2 * known_mult, 3 * known_mult, 4 * known_mult} if known_mult > 1 else {1}
    divisors = [base_R] if base_R > 1 else []
    stalls = 0
    while True:
        bound = M.N_MAX // R
        if bound < 2:
            R, stalls = 1, 0
            continue
        m = 2 * bound.bit_length() + M.ORDER_EXTRA_BITS
        y, _, ok = solver.sample_phases([(elem, m, R, M.ORDER_EXTRA_BITS, 0)], keep=100)[elem]
        frac = Fraction(y, 1 << m).limit_denominator(bound)
        r, num = frac.denominator, frac.numerator
        if r >= 4096 and (r & (r - 1)) == 0 and (num % r in (1, r - 1)):
            r = 1
        c = R * r
        cands = {c}
        for d in divisors:
            l = c * d // math.gcd(c, d)
            if l <= M.N_MAX:
                cands.add(l)
        if c > 1 and c not in divisors:
            divisors = (divisors + [c])[-8:]
        for cand in sorted(cands, reverse=True):
            if cand in tested:
                continue
            tested.add(cand)
            if solver.is_identity([(elem, cand)]):
                return _minimize_elem_order(solver, elem, cand)
        if r == 1:
            if R > 1:
                for mult in (2, 3):
                    cand = R * mult
                    if cand <= M.N_MAX and cand not in tested:
                        tested.add(cand)
                        if solver.is_identity([(elem, cand)]):
                            return _minimize_elem_order(solver, elem, cand)
            stalls += 1
            if stalls >= 3:
                R, stalls = 1, 0
            continue
        stalls = 0
        R = c if c <= M.N_MAX // 2 else 1


def find_subgroup_structure(oracle, target_shots=6):
    S = M.Solver(oracle)
    n_g = 1
    try:
        n_g = S.find_order()
    except M.OutOfBudget:
        n_g = getattr(S, "_best_order_guess", 1) or 1
        return {"order": n_g, "secret": 0, "orders": {S.g: n_g}, "invariants": [n_g], "rank": 1 if n_g > 1 else 0}

    if n_g == 1:
        ords = {S.g: 1}
        forged = []
        try:
            for X in S.X:
                if S.is_identity([(X, 1)]):
                    ords[X] = 1
                else:
                    ords[X] = _find_element_order(S, X, known_mult=1)
                    forged.append(X)
        except M.OutOfBudget:
            pass
        if not forged:
            return {"order": 1, "secret": 0, "orders": ords, "invariants": [1], "rank": 0}
        inv = [math.lcm(*[ords.get(X, 1) for X in forged])]
        return {"order": 1, "secret": "none", "orders": ords, "invariants": inv, "rank": 1}


    pe_g = {}
    for p in M._prime_factors(n_g):
        q = 1
        while n_g % (q * p) == 0:
            q *= p
        pe_g[p] = q

    m = n_g.bit_length() + M.EXTRA_BITS
    mg = n_g.bit_length() + M.G_EXTRA_BITS
    votes = {X: {p: {} for p in pe_g} for X in S.X}
    partial = {X: {p: {} for p in pe_g} for X in S.X}
    tested = {X: {} for X in S.X}
    failed = {X: 0 for X in S.X}
    tainted = {X: set() for X in S.X}
    inconsistent = {X: 0 for X in S.X}
    clean_pairs = {X: [] for X in S.X}
    order_divides_ng = {}
    logs = {}
    forged_set = set()
    saved_shots = []

    def to_grid(y, bits, tol, mod_n):
        v = Fraction(y * mod_n, 1 << bits)
        k = round(v)
        return k % mod_n, abs(v - k) <= tol

    def check_div_ng(X):
        if X not in order_divides_ng:
            order_divides_ng[X] = S.is_identity([(X, n_g)])
        return order_divides_ng[X]

    rounds = 0
    try:
        while len(logs) + len(forged_set) < len(S.X):
            todo = [X for X in S.X if X not in logs and X not in forged_set]
            spec = [(S.g, mg, 1, M.G_EXTRA_BITS, M.JOINT_CHECKS)] + [
                (X, m, 1, M.EXTRA_BITS, M.JOINT_CHECKS) for X in todo
            ]
            res = S.sample_phases(spec, keep=200)
            rounds += 1
            j, grid_g = to_grid(res[S.g][0], mg, M.G_CLOSE_TOL, n_g)
            ok_g = grid_g and res[S.g][2]
            shot_record = {S.g: j} if ok_g else {}

            for X in list(todo):
                s, grid_x = to_grid(res[X][0], m, M.CLOSE_TOL, n_g)
                chk_x = res[X][2]
                if not grid_x and chk_x:
                    if not check_div_ng(X):
                        forged_set.add(X)
                        continue
                ok_x = grid_x and chk_x
                if ok_g and ok_x:
                    shot_record[X] = s
                    g_jn = math.gcd(j, n_g)
                    if s % g_jn != 0:
                        if not check_div_ng(X):
                            forged_set.add(X)
                            continue
                        inconsistent[X] += 1
                        if inconsistent[X] >= 3:
                            forged_set.add(X)
                            continue
                    else:
                        clean_pairs[X].append((j, s))
                        for p, q in pe_g.items():
                            if j % p:
                                a_q = s * pow(j, -1, q) % q
                                entry = votes[X][p].setdefault(a_q, [0, 0])
                            else:
                                v, jj = 0, j
                                while jj % p == 0 and p ** v < q:
                                    jj //= p
                                    v += 1
                                pv = p ** v
                                qq = q // pv
                                a_q = (s // pv) * pow(jj, -1, qq) % qq if qq > 1 else 0
                                entry = partial[X][p].setdefault((qq, a_q), [0, 0])
                            entry[0] += 1
                            entry[1] = rounds

                if X in forged_set:
                    continue
                verdict = S._try_candidates(
                    X, n_g, pe_g, votes[X], partial[X], tested[X],
                    failed, tainted[X], check_div_ng
                )
                if verdict == "forged":
                    forged_set.add(X)
                elif verdict is not None:
                    logs[X] = verdict
                else:
                    pairs = clean_pairs[X]
                    if len(pairs) >= 3:
                        clique = []
                        for j_c, s_c in pairs:
                            if all((s_c * j_o - s_o * j_c) % n_g != 0 for j_o, s_o in clique):
                                clique.append((j_c, s_c))
                        if len(clique) >= 4 or (len(clique) >= 3 and tested[X]):
                            forged_set.add(X)

            if ok_g:
                saved_shots.append(shot_record)
    except M.OutOfBudget:

        for X in S.X:
            if X not in logs and X not in forged_set:
                if failed[X] > 0 or inconsistent[X] > 0:
                    forged_set.add(X)


    if not forged_set:
        secret = 1
        for X in S.X:
            secret = (secret * logs.get(X, 0)) % n_g
        ords = {S.g: n_g}
        for X in S.X:
            a_x = logs.get(X, 1)
            ords[X] = n_g // math.gcd(n_g, a_x) if a_x else 1
        return {"order": n_g, "secret": secret, "orders": ords, "invariants": [n_g], "rank": 1}


    ords = {S.g: n_g}
    for X, a_x in logs.items():
        ords[X] = n_g // math.gcd(n_g, a_x) if a_x else 1

    forged_list = sorted(forged_set)
    try:
        for X in forged_list:
            if check_div_ng(X):

                ords[X] = _minimize_elem_order(S, X, n_g)
            else:
                ords[X] = _find_element_order(S, X, known_mult=n_g)
    except M.OutOfBudget:
        for X in forged_list:
            if X not in ords:
                ords[X] = n_g

    gens = ([S.g] if n_g > 1 else []) + [X for X in forged_list if ords.get(X, 1) > 1]
    if len(gens) <= 1:
        inv = [ords[gens[0]]] if gens else [1]
        return {"order": n_g, "secret": "none", "orders": ords, "invariants": inv, "rank": len(gens)}

    exp_H = math.lcm(*[ords[el] for el in gens])
    primes = M._prime_factors(exp_H)
    pe = {}
    for p in primes:
        e = 0
        while exp_H % (p ** (e + 1)) == 0:
            e += 1
        pe[p] = (p ** e, e)

    prime_invs = {}
    ambiguous_primes = []
    for p, (q, max_e) in pe.items():
        active_gens = [el for el in gens if ords[el] % p == 0]
        if len(active_gens) <= 1:
            prime_invs[p] = [q]
        else:
            ambiguous_primes.append(p)

    if ambiguous_primes:

        rows_mod_M = []
        if exp_H == n_g:
            for sh in saved_shots:
                if all(el in sh for el in gens):
                    rows_mod_M.append([sh[el] % exp_H for el in gens])

        spec = []
        for el in gens:
            n_el = ords[el]
            eb = 9 if n_el > (1 << 40) else 7
            m_el = n_el.bit_length() + eb
            spec.append((el, m_el, 1, eb, 3))

        attempts = 0
        try:
            while len(rows_mod_M) < target_shots and attempts < target_shots * 3:
                attempts += 1
                res = S.sample_phases(spec, keep=0)
                row = []
                clean = True
                for el in gens:
                    n_el = ords[el]
                    y, m_el, ok_chk = res[el]
                    v = Fraction(y * n_el, 1 << m_el)
                    k = round(v)
                    if not (ok_chk and abs(v - k) <= 0.03):
                        clean = False
                        break
                    row.append(((k % n_el) * (exp_H // n_el)) % exp_H)
                if clean:
                    rows_mod_M.append(row)
        except M.OutOfBudget:
            pass

        for p in ambiguous_primes:
            q, max_e = pe[p]
            if not rows_mod_M:
                prime_invs[p] = [q]
                continue
            mat_p = [[x % q for x in r] for r in rows_mod_M]
            prime_invs[p] = _p_adic_snf_invariants(mat_p, p, max_e, filter_outliers=True)

    max_len = max((len(v) for v in prime_invs.values()), default=0)
    factors = [1] * max_len if max_len else [1]
    for p, invs in prime_invs.items():
        for idx, val in enumerate(invs):
            factors[idx] *= val
    return {
        "order": n_g,
        "secret": "none",
        "orders": ords,
        "invariants": factors,
        "rank": len([x for x in factors if x > 1]),
    }


if __name__ == "__main__":
    import test2

    print("=== Subgroup Structure Recovery (Strict 80,000 Query Budget) ===")
    for name in ("public_1", "public_2"):
        inst = load_public(name)
        orc = Oracle(inst, noise=1.0, budget=80000, seed=42)
        truth = true_subgroup_structure(inst)
        got = find_subgroup_structure(orc, target_shots=6)
        print(f"{name:26s}: invariants={got['invariants']}  truth={truth}  "
              f"match={got['invariants']==truth}  queries={orc.queries_used}/80000")

    for kind in test2.kinds:
        inst = test2.gen(kind, 1)
        orc = Oracle(inst, noise=1.0, budget=80000, seed=42)
        truth = true_subgroup_structure(inst)
        got = find_subgroup_structure(orc, target_shots=6)
        print(f"{kind:26s}: invariants={got['invariants']}  truth={truth}  "
              f"match={got['invariants']==truth}  queries={orc.queries_used}/80000")
