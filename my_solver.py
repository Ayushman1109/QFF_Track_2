import itertools
import math
import random
import time
from fractions import Fraction

from qiskit import QuantumCircuit
from silentgroup import oracle_gate


N_MAX = 1 << 60
ORDER_EXTRA_BITS = 6
EXTRA_BITS = 6
CLOSE_TOL = 0.08
G_EXTRA_BITS = 7
G_CLOSE_TOL = 0.03
CORR_DEPTH = 12
NOISE_SIGMA = 0.002
ID_MISS = 5e-5
ID_FALSE = 1e-5
ID_MAX_SHOTS = 50
RESERVE = 250
SMALL_ORDER = 1 << 32
ORDER_CHECKS = 0
JOINT_CHECKS = 3
CHECK_SLACK = 3
TIME_GUARD = 540.0
MAX_RETRIES = 8

DEBUG = False


class OutOfBudget(Exception):
    pass


def _add_phase_estimation(qc, element, m, offset, mult=1):
    for k in range(m, 0, -1):
        qc.reset(0)
        qc.h(0)
        qc.append(oracle_gate(element, mult << (k - 1)), [0])
        for j in range(1, min(CORR_DEPTH, m - k) + 1):
            with qc.if_test((qc.clbits[offset + k - 1 + j], 1)):
                qc.p(-2 * math.pi / (1 << (j + 1)), 0)
        qc.h(0)
        qc.measure(0, offset + k - 1)


def _add_checks(qc, element, m, offset, mult, ts, check_offset):
    for c, t in enumerate(ts):
        qc.reset(0)
        qc.h(0)
        qc.append(oracle_gate(element, mult * t), [0])
        for i in range(1, m + 1):
            num = (t << (m - i)) % (1 << m)
            if num:
                with qc.if_test((qc.clbits[offset + i - 1], 1)):
                    qc.p(-2 * math.pi * num / (1 << m), 0)
        qc.h(0)
        qc.measure(0, check_offset + c)


def _bits_from_key(key, i):
    if not isinstance(key, str):
        key = str(key)
    s = key.replace(" ", "").replace("_", "")
    if s.startswith(("0x", "0X")):
        try:
            s = bin(int(s, 16))[2:]
        except Exception:
            pass
    idx = len(s) - 1 - i
    if idx < 0 or idx >= len(s):
        return 0
    return 1 if s[idx] == "1" else 0


def _cf_extract(y, m, bound, min_jump=1500, max_ulp=3.0, allow_last=True):
    a, b = y, 1 << m
    p0, p1, q0, q1 = 0, 1, 1, 0
    convs = []
    while b:
        quot = a // b
        a, b = b, a - quot * b
        p0, p1 = p1, quot * p1 + p0
        q0, q1 = q1, quot * q1 + q0
        if q1 > bound:
            break
        next_q = (a // b) if b else (1 << m)
        ulp = abs(y * q1 - p1 * (1 << m)) / q1
        if q1 > 1 and not (q1 >= 4096 and (q1 & (q1 - 1)) == 0 and (p1 % q1 in (1, q1 - 1))):
            convs.append((q1, p1, next_q, ulp))
    if not convs:
        return (1, 0, 1 << m, 0.0) if allow_last else None
    both = [c for c in convs if c[2] >= min_jump and c[3] <= max_ulp]
    if both:
        return max(both, key=lambda c: c[2])
    if allow_last:
        huge = [c for c in convs if c[2] >= 10000]
        if huge:
            return max(huge, key=lambda c: c[2])
        last = convs[-1]
        if last[3] <= 1.8:
            return last
    return None


def _cf_with_repair(y, m, bound):
    if y == 0:
        return 1, True, True
    mj = 1200 if bound > (1 << 20) else 20
    direct = _cf_extract(y, m, bound, min_jump=mj, max_ulp=3.5, allow_last=True)
    if direct is not None:
        return direct[0], True, True
    if bound <= (1 << 20):
        return 1, False, False
    found = {}
    for d in (-1, 1, -2, 2, -3, 3):
        for bit in range(1, m):
            y2 = y + d * (1 << bit)
            if 0 <= y2 < (1 << m):
                res = _cf_extract(y2, m, bound, min_jump=1800, max_ulp=2.2, allow_last=False)
                if res is not None and res[0] > 100:
                    q_val = res[0]
                    score = (res[2], -res[3])
                    if q_val not in found or score > found[q_val][0]:
                        found[q_val] = (score, res)
    if found:
        ranked = [item[1] for item in sorted(found.values(), key=lambda x: x[0], reverse=True)]
        best = ranked[0]
        extra = [r[0] for r in ranked[1:3] if r[2] >= 10000 and r[3] <= 0.5]
        return (best[0], extra), True, (best[2] >= 10_000_000)
    return 1, False, False


class Solver:
    def __init__(self, oracle):
        self.oracle = oracle
        self.t0 = time.monotonic()
        try:
            names = list(getattr(oracle, "elements", None) or ["g", "A", "B", "C"])
        except Exception:
            names = ["g", "A", "B", "C"]
        self.g = "g" if "g" in names else names[0]
        others = [x for x in names if x != self.g]
        self.X = [x for x in ("A", "B", "C") if x in others] or others[:3]
        self.rng = random.Random(0x5EED)

    def _left(self):
        try:
            left = self.oracle.budget_left()
            return float("inf") if left is None else left
        except Exception:
            return float("inf")

    def _cost(self, qc, shots):
        for _ in range(3):
            try:
                return self.oracle.cost(qc, shots)
            except Exception:
                pass
        return 0

    def _run(self, qc, shots, keep=0):
        if time.monotonic() - self.t0 > TIME_GUARD:
            raise OutOfBudget()
        price = self._cost(qc, shots)
        if price > self._left() - keep:
            raise OutOfBudget()
        last_exc = None
        for _ in range(MAX_RETRIES):
            if time.monotonic() - self.t0 > TIME_GUARD:
                raise OutOfBudget()
            if price > self._left() - keep:
                raise OutOfBudget()
            try:
                res = self.oracle.run(qc, shots)
                if isinstance(res, dict) and res:
                    return res
            except Exception as exc:
                last_exc = exc
                msg = (type(exc).__name__ + " " + str(exc)).lower()
                if "budget" in msg:
                    raise OutOfBudget() from exc
        raise OutOfBudget() from last_exc

    def _build_phase_circuit(self, spec):
        total = sum(s[1] + s[4] for s in spec)
        qc = QuantumCircuit(1, total)
        off, layout = 0, []
        for elem, m, mult, margin, nck in spec:
            _add_phase_estimation(qc, elem, m, off, mult)
            hi = max(2, m - margin - CHECK_SLACK)
            ts = [self.rng.randrange(1 << (hi - 1), 1 << hi) | 1 for _ in range(nck)]
            _add_checks(qc, elem, m, off, mult, ts, off + m)
            layout.append((elem, m, off, nck))
            off += m + nck
        return qc, layout

    def sample_phases(self, spec, keep=RESERVE):
        qc, layout = self._build_phase_circuit(spec)
        counts = self._run(qc, 1, keep)
        key = next(iter(counts), "0") if isinstance(counts, dict) else "0"
        out = {}
        for elem, m, o, nck in layout:
            y = 0
            for k in range(1, m + 1):
                y |= _bits_from_key(key, o + k - 1) << (m - k)
            ok = all(_bits_from_key(key, o + m + c) == 0 for c in range(nck))
            out[elem] = (y, m, ok)
        return out

    def is_identity(self, terms, keep=0, id_miss=ID_MISS):
        terms = [(e, int(p)) for e, p in terms if int(p) != 0]
        if not terms:
            return True
        cost = sum(max(1, abs(p).bit_length()) for _, p in terms)
        p0 = min(0.05, 1.5 * (math.pi * NOISE_SIGMA) ** 2 * cost)
        w1 = math.log(0.5 / p0)
        w0 = math.log(0.5 / (1 - p0))
        upper = math.log((1 - id_miss) / ID_FALSE)
        lower = math.log(id_miss / (1 - ID_FALSE))
        qc = QuantumCircuit(1, 1)
        qc.h(0)
        for e, p in terms:
            qc.append(oracle_gate(e, p), [0])
        qc.h(0)
        qc.measure(0, 0)
        llr, shots = 0.0, 0
        while shots < ID_MAX_SHOTS:
            need = math.ceil((llr - lower) / -w0)
            batch = 2 if shots == 0 else max(2, min(need, 5))
            counts = self._run(qc, batch, keep)
            ones = sum(c for k, c in counts.items() if _bits_from_key(k, 0) == 1)
            shots += batch
            llr += ones * w1 + (batch - ones) * w0
            if llr >= upper:
                return False
            if llr <= lower:
                return True
        return llr < 0

    def find_order(self):
        self._pre_samples = []
        if self.is_identity([(self.g, 1)]):
            return 1
        if self.is_identity([(self.g, 2)]):
            return 2
        R, tested, divisors, stalls, d_first = 1, {1, 2}, [], 0, None
        rec_tried = set()
        self._best_order_guess = 1
        first_full_shot = True
        while True:
            full_bound = N_MAX // R
            if full_bound < 2:
                R, stalls = 1, 0
                continue
            if R > 1 and R not in rec_tried and full_bound > (1 << 9):
                bound = 1 << 9
                rec_tried.add(R)
                o_eb = 4
            else:
                bound = full_bound
                o_eb = 4 if R == 1 else 6
            m = 2 * (bound - 1).bit_length() + o_eb
            if first_full_shot and R == 1 and len(self.X) > 0:
                first_full_shot = False
                mx = 64
                spec = [(self.g, m, 1, o_eb, 0)] + [
                    (X, mx, 1, 6, 0) for X in self.X
                ]
                res = self.sample_phases(spec, keep=100)
                self._pre_samples.append(res)
                y, _, ok = res[self.g]
            else:
                res = self.sample_phases(
                    [(self.g, m, R, o_eb, 0)], keep=100
                )
                if R == 1:
                    self._pre_samples.append(res)
                y, _, ok = res[self.g]
            r_out, ok_cf, is_direct = _cf_with_repair(y, m, bound)
            if not ok_cf:
                if R > 1:
                    stalls += 1
                    if stalls >= 2:
                        R, stalls = 1, 0
                continue
            if isinstance(r_out, tuple):
                r, extra_rs = r_out
            else:
                r, extra_rs = r_out, []
            c = R * r
            cands = {c} | {R * er for er in extra_rs if R * er <= N_MAX}
            for d in divisors:
                l = c * d // math.gcd(c, d)
                if l <= N_MAX:
                    cands.add(l)
            if is_direct and c > 1 and c not in divisors:
                divisors = (divisors + [c])[-8:]
                self._best_order_guess = max(cands)
                if d_first is None and R == 1 and ok_cf:
                    d_first = c
            for cand in sorted(cands, reverse=True):
                if cand in tested or (not ok and cand > SMALL_ORDER):
                    continue
                tested.add(cand)
                if self.is_identity([(self.g, cand)], id_miss=2e-3):
                    ep = _prime_factors(cand // d_first) if (d_first and cand % d_first == 0) else None
                    return self._minimise_order(cand, extra_primes=ep, false_divs=tested - {cand})
            if r == 1:
                if R > 1:
                    for mult in (2, 3):
                        cand = R * mult
                        if cand <= N_MAX and cand not in tested:
                            tested.add(cand)
                            if self.is_identity([(self.g, cand)], id_miss=2e-3):
                                return self._minimise_order(cand, extra_primes={2, 3}, false_divs=tested - {cand})
                stalls += 1
                if stalls >= 2:
                    R, stalls = 1, 0
                continue
            stalls = 0
            if is_direct:
                R = c if c <= N_MAX // 2 else 1


    def _minimise_order(self, c, extra_primes=None, false_divs=()):
        self._best_order_guess = c
        if extra_primes is not None:
            primes = sorted(set(extra_primes) | set(_prime_factors(c)), reverse=True)
        else:
            primes = sorted(_prime_factors(c), reverse=True)
        for p in primes:
            if c % p != 0:
                continue
            sub = c // p
            if any(d % sub == 0 for d in false_divs):
                continue
            if not self.is_identity([(self.g, c // p)]):
                continue
            c //= p
            self._best_order_guess = c
            while c % p == 0:
                if self.is_identity([(self.g, c // p)]):
                    c //= p
                    self._best_order_guess = c
                else:
                    break
        self._best_order_guess = c
        return c

    def solve_logs(self, n):
        if n == 1:
            for X in self.X:
                if not self.is_identity([(X, 1)]):
                    return None
            return {X: 0 for X in self.X}

        pe = {}
        for p in _prime_factors(n):
            q = 1
            while n % (q * p) == 0:
                q *= p
            pe[p] = q

        n_odd = n
        while n_odd % 2 == 0:
            n_odd //= 2
        pseudo_pow2 = any(min(n % (1 << k), (1 << k) - (n % (1 << k))) * (1 << 20) < (1 << k) for k in range(28, n.bit_length()))
        can_grid_repair = (n_odd >= 2048) and not pseudo_pow2
        nck = 2 if can_grid_repair else 3

        eb = 6 if can_grid_repair else 7
        geb = 7 if can_grid_repair else 8
        m = n.bit_length() + eb
        mg = n.bit_length() + geb
        votes = {X: {p: {} for p in pe} for X in self.X}
        partial = {X: {p: {} for p in pe} for X in self.X}
        base_mod = {X: None for X in self.X}
        tested = {X: {} for X in self.X}
        failed = {X: 0 for X in self.X}
        tainted = {X: set() for X in self.X}
        inconsistent = {X: 0 for X in self.X}
        refuted = {X: 0 for X in self.X}
        clean_pairs = {X: [] for X in self.X}
        order_checked = set()
        logs, rounds = {}, 0

        self._log_state = (pe, votes, partial, logs, failed, inconsistent, clean_pairs)

        def to_grid_strict(y, bits, tol):
            v = Fraction(y * n, 1 << bits)
            k = round(v)
            return k % n, abs(v - k) <= tol

        def try_repair_grid_list(y, bits):
            if not can_grid_repair:
                return []
            margin = bits - n.bit_length()
            if margin < 5:
                return []
            tol_rep = min(0.008, 1.4 / (1 << margin))
            v0 = Fraction(y * n, 1 << bits)
            k0 = round(v0)
            err0 = float(abs(v0 - k0))
            if err0 <= min(0.03, 2.0 / (1 << margin)):
                return [(k0 % n, err0)]
            found = {}
            for d in (-1, 1, -2, 2, -3, 3):
                pen = 0 if abs(d) <= 2 else 1
                for b in range(1, bits):
                    y2 = y + d * (1 << b)
                    if 0 <= y2 < (1 << bits):
                        v2 = Fraction(y2 * n, 1 << bits)
                        k2 = round(v2)
                        err = float(abs(v2 - k2))
                        if err <= tol_rep:
                            val = k2 % n
                            key = (pen, err)
                            if val not in found or key < found[val]:
                                found[val] = key
            ranked = sorted(found.items(), key=lambda t: t[1])
            return [(val, info[1]) for val, info in ranked[:3]]

        def try_repair_grid(y, bits):
            lst = try_repair_grid_list(y, bits)
            return lst[0] if lst else (None, 1.0)

        def try_guided_lift(y, bits, j_val, n_sub, a0):
            d = n // n_sub
            if math.gcd(j_val, d) != 1:
                return None
            best = None
            for t in range(d):
                cand = (a0 + t * n_sub) % n
                s_pred = (cand * j_val) % n

                y_pred = (s_pred * (1 << bits)) / n
                diff = abs(y - y_pred) % (1 << bits)
                diff = min(diff, (1 << bits) - diff)
                if best is None or diff < best[0]:
                    best = (diff, cand)
            if best is not None and best[0] <= 2.5:
                return best[1]
            return None

        def verify_order_divides_n(X):
            if X in order_checked:
                return True
            order_checked.add(X)
            return self.is_identity([(X, n)])

        def max_clique_size(pairs):
            clique = []
            for j_c, s_c in pairs:
                if all((s_c * j_o - s_o * j_c) % n != 0 for j_o, s_o in clique):
                    clique.append((j_c, s_c))
            return len(clique)

        pre_queue = list(getattr(self, "_pre_samples", []))

        while len(logs) < len(self.X):
            todo = [X for X in self.X if X not in logs]
            suspects = [X for X in todo if refuted[X] >= 1 or failed[X] >= 1 or inconsistent[X] >= 1 or max_clique_size(clean_pairs[X]) >= 2]
            active_todo = list(suspects) if (suspects and not pre_queue) else list(todo)

            is_pre_round = False
            if pre_queue:
                is_pre_round = True
                res = pre_queue.pop(0)
                active_todo = [X for X in active_todo if X in res]
                if not active_todo:
                    continue
            else:
                eb_cur = (eb + 1) if len(active_todo) == 1 else eb
                geb_cur = (geb + 1) if len(active_todo) == 1 else geb
                m_cur_target = n.bit_length() + eb_cur
                mg_cur_target = n.bit_length() + geb_cur
                use_short = (rounds == 1 and can_grid_repair and all(base_mod[X] is not None for X in active_todo) and not any(refuted[X] >= 1 or inconsistent[X] >= 1 or failed[X] >= 1 for X in self.X))
                while len(active_todo) > 1:
                    spec_try = [(self.g, mg_cur_target, 1, geb_cur, nck)] + [
                        (X, 15 if use_short else m_cur_target, 1, 4 if use_short else eb_cur, 0 if use_short else nck)
                        for X in active_todo
                    ]
                    qc_try, _ = self._build_phase_circuit(spec_try)
                    if self._cost(qc_try, 1) <= self._left() - RESERVE:
                        break
                    active_todo.pop()

                spec = [(self.g, mg_cur_target, 1, geb_cur, nck)] + [
                    (X, 15 if use_short else m_cur_target, 1, 4 if use_short else eb_cur, 0 if use_short else nck)
                    for X in active_todo
                ]
                res = self.sample_phases(spec, keep=RESERVE)

            rounds += 1
            yg, mg_cur, chkg = res[self.g]
            hi_margin_g = False
            if mg_cur > n.bit_length() + 14:
                mg_use = n.bit_length() + 14
                sh = mg_cur - mg_use
                yg_top = (yg + (1 << (sh - 1))) >> sh
                j, grid_g = to_grid_strict(yg_top, mg_use, 1.8 / (1 << 14))
                j_opt, err_g = (j, 0.0) if grid_g else try_repair_grid(yg_top, mg_use)
                ok_g = grid_g and (n % 2 == 1)
                hi_margin_g = (grid_g or (j_opt is not None and err_g <= 2e-4)) and (n % 2 == 1) and can_grid_repair
            else:
                j, grid_g = to_grid_strict(yg, mg_cur, G_CLOSE_TOL)
                ok_g = grid_g and chkg
                j_opt, _ = (j, 0.0) if ok_g else try_repair_grid(yg, mg_cur)

            grid = {}
            opt_s = {}
            opt_s_list = {}
            hi_margin_x = {}
            for X in active_todo:
                yx, mx_cur, chkx = res[X]
                if mx_cur < n.bit_length() and base_mod[X] is not None:

                    grid[X] = (0, False)
                    opt_s[X] = None
                    opt_s_list[X] = []
                    hi_margin_x[X] = False
                    continue
                if mx_cur >= n.bit_length() + 14:
                    mx_use = n.bit_length() + 14
                    shx = mx_cur - mx_use
                    yx_top = (yx + (1 << (shx - 1))) >> shx if shx > 0 else yx
                    s, grid_x = to_grid_strict(yx_top, mx_use, 1.8 / (1 << 14))
                    rep_list = [(s, 0.0)] if grid_x else try_repair_grid_list(yx_top, mx_use)
                    ok_x = grid_x and (n % 2 == 1)
                    hi_margin_x[X] = ok_x and can_grid_repair
                    opt_s[X] = rep_list[0][0] if rep_list else None
                    opt_s_list[X] = [v for v, _ in rep_list]
                else:
                    hi_margin_x[X] = False
                    tol_strict = 0.02 if is_pre_round else CLOSE_TOL
                    s, grid_x = to_grid_strict(yx, mx_cur, tol_strict)
                    if not grid_x and chkx and not is_pre_round:
                        if not verify_order_divides_n(X):
                            return None
                    ok_x = grid_x and (chkx and not is_pre_round)
                    rep_list = [(s, 0.0)] if grid_x else try_repair_grid_list(yx, mx_cur)
                    opt_s[X] = rep_list[0][0] if rep_list else None
                    opt_s_list[X] = [v for v, _ in rep_list]
                grid[X] = (s, ok_x)

            if ok_g:
                g_jn = math.gcd(j, n)
                for X in active_todo:
                    s, ok_x = grid[X]
                    if not ok_x:

                        if is_pre_round and 1 < g_jn <= 64 and opt_s[X] is not None and opt_s[X] % g_jn == 0:
                            n_sub = n // g_jn
                            inv_j_sub = pow(j // g_jn, -1, n_sub)
                            base_mod[X] = (n_sub, ((opt_s[X] // g_jn) * inv_j_sub) % n_sub)
                        continue
                    clean_pairs[X].append((j, s))
                    if s % g_jn != 0:
                        inconsistent[X] += 1
                        refuted[X] += 1
                        if hi_margin_g and hi_margin_x.get(X, False):
                            return None
                        if refuted[X] >= (2 if can_grid_repair else 4):
                            return None
                        if can_grid_repair and max_clique_size(clean_pairs[X]) >= 2:
                            return None
                        continue
                    if 1 < g_jn <= 64 and base_mod[X] is None:
                        n_sub = n // g_jn
                        inv_j_sub = pow(j // g_jn, -1, n_sub)
                        base_mod[X] = (n_sub, ((s // g_jn) * inv_j_sub) % n_sub)
                    for p, q in pe.items():
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


            if is_pre_round and can_grid_repair and (n % 2 == 1) and j_opt is not None and math.gcd(j_opt, n) > 1 and (any(res[X][1] >= n.bit_length() + 14 for X in active_todo) or math.gcd(j_opt, n) <= 64):
                d_lf = math.gcd(j_opt, n)
                n_sub_lf = n // d_lf
                inv_j_lf = pow(j_opt // d_lf, -1, n_sub_lf) if n_sub_lf > 1 else 0
                lf_cands = {}
                for X in active_todo:
                    s_lf = opt_s.get(X)
                    if s_lf is not None:
                        if s_lf % d_lf != 0:
                            if hi_margin_g and hi_margin_x.get(X, False):
                                return None
                            refuted[X] += 1
                        else:
                            a0_lf = ((s_lf // d_lf) * inv_j_lf) % n_sub_lf
                            lf_cands[X] = (d_lf * a0_lf) % n
                if lf_cands and not any(refuted[X] >= 1 for X in active_todo):
                    for X, Ex in lf_cands.items():
                        if not self.is_identity([(X, d_lf), (self.g, -Ex)], id_miss=0.10):
                            if hi_margin_g and hi_margin_x.get(X, False):
                                if not self.is_identity([(X, d_lf), (self.g, -Ex)]):
                                    return None
                            refuted[X] += 1
                            break


            opt_cands = {}
            if j_opt is not None:
                d_opt = math.gcd(j_opt, n)
                if d_opt == 1:
                    inv_j_opt = pow(j_opt, -1, n)
                    for X in active_todo:
                        if base_mod[X] is not None and res[X][1] < n.bit_length():
                            n_sub, a0 = base_mod[X]
                            yx, mx_cur, _ = res[X]
                            gl = try_guided_lift(yx, mx_cur, j_opt, n_sub, a0)
                            if gl is not None:
                                opt_cands[X] = [gl]
                        elif opt_s_list.get(X):
                            opt_cands[X] = [((sv * inv_j_opt) % n) for sv in opt_s_list[X]]
                        elif base_mod[X] is not None:
                            n_sub, a0 = base_mod[X]
                            yx, mx_cur, _ = res[X]
                            gl = try_guided_lift(yx, mx_cur, j_opt, n_sub, a0)
                            if gl is not None:
                                opt_cands[X] = [gl]

            batch_ready = {X: clist[0] for X, clist in opt_cands.items() if clist and clist[0] not in tested[X]}
            if False and can_grid_repair and len(batch_ready) >= 2:
                multi_x = [X for X in batch_ready if len(opt_cands[X]) > 1]
                batch_variants = [dict(batch_ready)]
                if len(multi_x) == 1:
                    xm = multi_x[0]
                    for alt_c in opt_cands[xm][1:3]:
                        v_map = dict(batch_ready)
                        v_map[xm] = alt_c
                        batch_variants.append(v_map)
                for b_map in batch_variants:
                    weights = [1, 1, 1][:len(b_map)]
                    terms = [(X, w) for X, w in zip(b_map.keys(), weights)]
                    g_pow = -sum(w * b_map[X] for X, w in zip(b_map.keys(), weights)) % n
                    terms.append((self.g, g_pow))
                    if self.is_identity(terms, id_miss=5e-3):
                        w2 = [1, -1, 1][:len(b_map)]
                        t2 = [(X, w) for X, w in zip(b_map.keys(), w2)]
                        g_pow2 = -sum(w * b_map[X] for X, w in zip(b_map.keys(), w2)) % n
                        t2.append((self.g, g_pow2))
                        if self.is_identity(t2, id_miss=5e-3):
                            for X, cand in b_map.items():
                                tested[X][cand] = True
                                logs[X] = cand
                            active_todo = [X for X in active_todo if X not in logs]
                            break

            new_suspect = False
            rem_opt = [X for X in active_todo if X in opt_cands and any(c not in tested[X] for c in opt_cands[X])]
            if len(rem_opt) == 1:
                X = rem_opt[0]
                is_strict_pair = ok_g and grid[X][1]
                solved_x = False
                for cand in opt_cands[X]:
                    if cand in tested[X]:
                        continue
                    ok_full = self.is_identity([(X, 1), (self.g, -cand)], id_miss=(1e-2 if can_grid_repair else ID_MISS))
                    tested[X][cand] = ok_full
                    if ok_full:
                        logs[X] = cand
                        active_todo.remove(X)
                        solved_x = True
                        break
                if not solved_x and is_strict_pair:
                    if hi_margin_g and hi_margin_x.get(X, False):
                        return None
                    refuted[X] += 1
                    if refuted[X] >= (2 if can_grid_repair else 3):
                        return None
                    new_suspect = True
            elif len(rem_opt) >= 2:
                screen_passed = []
                for X in rem_opt:
                    is_strict_pair = ok_g and grid[X][1]
                    solved_quick = False
                    for cand in opt_cands[X]:
                        if cand in tested[X]:
                            continue
                        ok_quick = self.is_identity([(X, 1), (self.g, -cand)], id_miss=0.10)
                        if ok_quick:
                            screen_passed.append((X, cand, is_strict_pair))
                            solved_quick = True
                            break
                        else:
                            tested[X][cand] = False
                    if not solved_quick and is_strict_pair:
                        if hi_margin_g and hi_margin_x.get(X, False):
                            if not self.is_identity([(X, 1), (self.g, -opt_cands[X][0])]):
                                return None
                        refuted[X] += 1
                        if refuted[X] >= (2 if can_grid_repair else 3):
                            return None
                        new_suspect = True
                        break
                for X, cand, is_strict_pair in screen_passed:
                    if new_suspect and is_strict_pair:

                        continue
                    ok_full = self.is_identity([(X, 1), (self.g, -cand)], id_miss=(1e-2 if can_grid_repair else ID_MISS))
                    tested[X][cand] = ok_full
                    if ok_full:
                        logs[X] = cand
                        if X in active_todo:
                            active_todo.remove(X)
                    elif is_strict_pair:
                        if hi_margin_g and hi_margin_x.get(X, False):
                            return None
                        refuted[X] += 1
                        if refuted[X] >= (2 if can_grid_repair else 3):
                            return None
                        new_suspect = True
                        break
            if new_suspect:
                continue

            for X in active_todo:
                tied_lifts = 1
                for p_k, q_k in pe.items():
                    if not votes[X][p_k]:
                        if partial[X][p_k]:
                            best_qq = max(qq for (qq, _) in partial[X][p_k])
                            tied_lifts *= (q_k // best_qq)
                        else:
                            tied_lifts = 999
                verdict = None if tied_lifts > 2 else self._try_candidates(
                    X, n, pe, votes[X], partial[X], tested[X],
                    failed, tainted[X], verify_order_divides_n,
                    max_lifts=16, can_grid_repair=can_grid_repair,
                    inconsistent_cnt=inconsistent[X],
                )
                if verdict == "forged":
                    return None
                if verdict is not None:
                    logs[X] = verdict
                    continue
                clq = max_clique_size(clean_pairs[X])
                if clq >= 4 or (clq >= 3 and tested[X]):
                    return None
                if can_grid_repair and clq >= 2 and (refuted[X] >= 1 or failed[X] + inconsistent[X] >= 1 or (n % 2 == 1 and len(clean_pairs[X]) >= 2)):
                    return None
        return logs


    def _build_candidates(self, pe, votes, partial, tested, tainted, max_lifts=4):
        options, consensus, lifted = [], True, False
        for p, q in pe.items():
            ranked = []
            for a, (c, last) in votes[p].items():
                supp = c + sum(pc for (qq, pa), (pc, _) in partial[p].items() if qq > 1 and a % qq == pa)
                is_tainted = (p, a) in tainted and supp < 2
                ranked.append((a, [supp, not is_tainted, last, c]))
            ranked.sort(key=lambda t: (-t[1][1], -t[1][0], -t[1][3], -t[1][2]))

            best_partial_lifts = None
            if partial[p]:
                scored_p = []
                for (qq, pa), (pc, plast) in partial[p].items():
                    if q // qq <= max_lifts:
                        supp = pc + sum(fc for fa, (fc, _) in votes[p].items() if qq > 1 and fa % qq == pa)
                        supp += sum(
                            oc for (oqq, opa), (oc, _) in partial[p].items()
                            if oqq != qq and qq > 1 and oqq > 1 and (pa % min(qq, oqq) == opa % min(qq, oqq))
                        )
                        scored_p.append(((qq, pa), (supp, qq, plast)))
                if scored_p:
                    scored_p.sort(key=lambda t: (-t[1][0], -t[1][1], -t[1][2]))
                    best_partial_lifts = scored_p[0]

            if ranked:
                if ranked[0][1][0] < 2:
                    consensus = False
                if not ranked[0][1][1] and best_partial_lifts is not None:
                    (qq, pa), (supp, _, _) = best_partial_lifts
                    lifted = True
                    cand_list = [(pa + qq * t, supp) for t in range(q // qq)]
                    cand_list.extend([(a, info[0]) for a, info in ranked[:2]])
                    options.append(cand_list)
                    continue
                cand_list = [(a, info[0]) for a, info in ranked[:2]]
                if len(pe) == 1 and all(a in tested for a, _ in cand_list) and best_partial_lifts is not None:
                    (qq, pa), (supp, _, _) = best_partial_lifts
                    lifted = True
                    for t in range(q // qq):
                        lift_a = pa + qq * t
                        if lift_a not in tested:
                            cand_list.append((lift_a, supp))
                options.append(cand_list)
                continue

            consensus = False
            if best_partial_lifts is None:
                return []
            (qq, a), (c, _, _) = best_partial_lifts
            lifted = True
            options.append([(a + qq * t, c) for t in range(q // qq)])

        primes = list(pe)
        cands = []
        for idx, combo in enumerate(itertools.product(*options)):
            if idx >= 64:
                break
            r, mod = 0, 1
            for p, (a, _) in zip(primes, combo):
                r, mod = _crt(r, mod, a, pe[p])
            cands.append((idx, r, combo, consensus, lifted))
        return cands

    def _try_candidates(self, X, n, pe, votes, partial, tested, failed,
                        tainted=None, verify_order_divides_n=None,
                        per_round=2, max_lifts=4, can_grid_repair=False,
                        inconsistent_cnt=0):
        if tainted is None:
            tainted = set()
        for p in pe:
            if votes[p]:
                raw_ranked = sorted(votes[p].items(), key=lambda t: (-t[1][0], -t[1][1]))
                if len(raw_ranked) >= (3 if can_grid_repair else 4) and raw_ranked[0][1][0] < 2:
                    return "forged"
        cands = self._build_candidates(pe, votes, partial, tested, tainted, max_lifts)
        if not cands:
            return None
        consensus = cands[0][3]
        lifted = cands[0][4]
        primes = list(pe)
        limit_tries = 4 if (lifted and len(cands) <= 4) else per_round
        tries = 0
        for idx, cand, combo, _, _ in cands:
            if cand in tested:
                if idx == 0 and consensus and not tested[cand]:
                    return "forged"
                continue
            if tries >= limit_tries:
                break
            tries += 1
            ok = self.is_identity([(X, 1), (self.g, -cand)])
            tested[cand] = ok
            if ok:
                return cand
            if idx == 0 and consensus:
                return "forged"
            if idx == 0 and not lifted:
                used_tainted = any((p, a) in tainted and c < 2 for p, (a, c) in zip(primes, combo))
                for p, (a, c) in zip(primes, combo):
                    if c < 2:
                        tainted.add((p, a))
                if not used_tainted:
                    failed[X] += 1
                    if failed[X] + inconsistent_cnt >= (2 if can_grid_repair else 3):
                        return "forged"
        return None

    def _fallback_secret(self):
        if not self.n or not hasattr(self, "_log_state"):
            return 0
        pe, votes, partial, logs, failed, inconsistent, clean_pairs = self._log_state
        for X in self.X:
            if X in logs:
                continue
            if failed.get(X, 0) >= 2 or inconsistent.get(X, 0) >= 2:
                return "none"
            pairs = clean_pairs.get(X, [])
            if len(pairs) >= 2:
                clique = []
                for j_c, s_c in pairs:
                    if all((s_c * j_o - s_o * j_c) % self.n != 0 for j_o, s_o in clique):
                        clique.append((j_c, s_c))
                if len(clique) >= 3:
                    return "none"
            for p in pe:
                if len(votes[X][p]) >= 3 and max((v[0] for v in votes[X][p].values()), default=0) < 2:
                    return "none"
        guess_logs = dict(logs)
        for X in self.X:
            if X in guess_logs:
                continue
            r, mod = 0, 1
            for p, q in pe.items():
                if votes[X][p]:
                    a = max(votes[X][p].items(), key=lambda t: (t[1][0], t[1][1]))[0]
                elif partial[X][p]:
                    (qq, a), _ = max(partial[X][p].items(), key=lambda t: (t[0][0], t[1][0], t[1][1]))
                else:
                    return "none" if any(failed.values()) else 0
                r, mod = _crt(r, mod, a, q)
            guess_logs[X] = r
        secret = 1
        for X in self.X:
            secret = (secret * guess_logs[X]) % self.n
        return secret

    def solve(self):
        self.n = None
        self._best_order_guess = 1
        try:
            self.n = self.find_order()
            logs = self.solve_logs(self.n)
            if logs is None:
                return {"order": int(self.n), "secret": "none"}
            secret = 1
            for X in self.X:
                secret = (secret * logs[X]) % self.n
            return {"order": int(self.n), "secret": int(secret)}
        except Exception:
            order_ans = int(self.n or getattr(self, "_best_order_guess", 1) or 1)
            if order_ans < 1:
                order_ans = 1
            try:
                sec = self._fallback_secret() if self.n else 0
            except Exception:
                sec = 0
            return {
                "order": order_ans,
                "secret": sec if sec == "none" else int(sec) % order_ans,
            }


def _crt(r1, m1, r2, m2):
    g = math.gcd(m1, m2)
    if (r1 - r2) % g:
        return None
    l = m1 // g * m2
    t = ((r2 - r1) // g) * pow(m1 // g, -1, m2 // g) % (m2 // g) if m2 // g > 1 else 0
    return ((r1 + m1 * t) % l, l)


def _is_probable_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n):
    if n % 2 == 0:
        return 2
    rng = random.Random(n)
    while True:
        y, c, mm = rng.randrange(1, n), rng.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(mm, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += mm
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def _prime_factors(n):
    out, small = set(), 2
    while small < 1000 and n > 1:
        while n % small == 0:
            out.add(small)
            n //= small
        small += 1
    stack = [n] if n > 1 else []
    while stack:
        x = stack.pop()
        if x == 1:
            continue
        if _is_probable_prime(x):
            out.add(x)
            continue
        f = _pollard_rho(x)
        stack += [f, x // f]
    return out


def solve(oracle) -> dict:
    try:
        return Solver(oracle).solve()
    except Exception:
        return {"order": 1, "secret": 0}
