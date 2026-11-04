import cmath, json, math, os, random
from fractions import Fraction

class _OracleGate:
    def __init__(self, element, power): self.element, self.power = element, int(power)

def oracle_gate(element, power=1): return _OracleGate(element, power)

def _call_cost(k): return max(1, abs(k).bit_length())

def _circuit_cost(ops):
    c = 0
    for op in ops:
        if op[0] == "oracle": c += _call_cost(op[2].power)
        elif op[0] == "if": c += _circuit_cost(op[3])
    return c

class Group:
    def __init__(self, moduli, elements, answer, name="mock"):
        self.moduli, self.elements, self.answer, self.name = tuple(moduli), elements, answer, name

def make_instance(moduli, seed, forged=""):
    rng = random.Random(seed)
    rnd = lambda: tuple(rng.randrange(d) for d in moduli)
    g = rnd()
    n = math.lcm(*[d // math.gcd(x, d) for x, d in zip(g, moduli)])
    pw = lambda e: tuple(e * x % d for x, d in zip(g, moduli))
    els, ex = {"g": g}, {}
    for name in "ABC":
        if name in forged:
            while True:
                cand = rnd()

                if n > 10**7 or all(pw(e) != cand for e in range(n)):
                    break
            els[name] = cand
        else:
            ex[name] = rng.randrange(n); els[name] = pw(ex[name])
    secret = "none" if forged else ex["A"] * ex["B"] * ex["C"] % n
    return Group(moduli, els, {"order": n, "secret": secret})

def MockGroup(moduli, *, g=None, a=None, b=None, c=None, forgery=False, seed=None):
    rng = random.Random(seed)
    moduli = tuple(moduli)
    if g is None:
        g = tuple(rng.randrange(d) for d in moduli)
    else:
        g = tuple(g)
    ord_g = math.lcm(*[d // math.gcd(x, d) for x, d in zip(g, moduli)])
    if a is None: a = rng.randrange(ord_g)
    if b is None: b = rng.randrange(ord_g)
    if c is None: c = rng.randrange(ord_g)
    pw = lambda e: tuple((e * x) % d for x, d in zip(g, moduli))
    A, B, C = pw(a), pw(b), pw(c)
    secret = (a * b * c) % ord_g
    if forgery:
        A = tuple(rng.randrange(d) for d in moduli)
        secret = "none"
    return Group(moduli, {"g": g, "A": A, "B": B, "C": C}, {"order": ord_g, "secret": secret})

def load_public(name):
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "instances", f"{name}.json")
    raw = json.load(open(path))
    els = {k: tuple(v) for k, v in raw["elements"].items()}
    return Group(raw["moduli"], els, raw["answer"], name=raw.get("name", name))

class Oracle:
    def __init__(self, inst, noise=1.0, budget=80000, seed=None):
        self.inst, self.noise, self.budget = inst, noise, budget
        self.queries_used = 0; self.rng = random.Random(seed)
        self.elements = list(inst.elements)
    def cost(self, qc, shots): return _circuit_cost(qc.ops) * shots
    def budget_left(self):
        return None if self.budget is None else max(0, self.budget - self.queries_used)
    def run(self, qc, shots=1):
        assert qc.num_qubits <= 16
        c = self.cost(qc, shots)
        if self.budget is not None and c > self.budget_left():
            raise RuntimeError("out of budget")
        self.queries_used += c
        counts = {}
        for _ in range(shots):
            chi = [self.rng.randrange(d) for d in self.inst.moduli]; self.last_chi = chi
            bits = [0] * qc.num_clbits
            self._exec(qc.ops, [1+0j, 0j], bits, chi)
            key = "".join(str(b) for b in reversed(bits))
            counts[key] = counts.get(key, 0) + 1
        return counts
    def _phase(self, elem, k, chi):
        x = self.inst.elements[elem]
        f = sum(Fraction(c * xi, d) for c, xi, d in zip(chi, x, self.inst.moduli)) * k
        return float(f - math.floor(f))
    def _exec(self, ops, s, bits, chi):
        r2 = 1 / math.sqrt(2)
        for op in ops:
            t = op[0]
            if t == "h": s[0], s[1] = r2*(s[0]+s[1]), r2*(s[0]-s[1])
            elif t == "p": s[1] *= cmath.exp(1j * op[2])
            elif t == "reset": s[0], s[1] = 1+0j, 0j
            elif t == "oracle":
                k = op[2].power
                ph = self._phase(op[2].element, k, chi)
                ph += self.rng.gauss(0, 0.002 * math.sqrt(_call_cost(k))) * self.noise
                s[1] *= cmath.exp(2j * math.pi * ph)
            elif t == "measure":
                p1 = abs(s[1])**2
                b = 1 if self.rng.random() < p1 else 0
                s[0], s[1] = (0j, 1+0j) if b else (1+0j, 0j)
                bits[op[2]] = b
            elif t == "if":
                if bits[op[1]] == op[2]: self._exec(op[3], s, bits, chi)
