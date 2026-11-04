class _Clbit:
    def __init__(self, i): self.index = i

class _IfCtx:
    def __init__(self, qc, cond): self.qc, self.cond, self.block = qc, cond, []
    def __enter__(self):
        self.qc._stack.append(self.block); return self
    def __exit__(self, *a):
        self.qc._stack.pop()
        self.qc._stack[-1].append(("if", self.cond[0].index, self.cond[1], self.block))

class QuantumCircuit:
    def __init__(self, nq, nc):
        self.num_qubits, self.num_clbits = nq, nc
        self.clbits = [_Clbit(i) for i in range(nc)]
        self.ops = []; self._stack = [self.ops]
    def _add(self, op): self._stack[-1].append(op)
    def h(self, q): self._add(("h", q))
    def p(self, t, q): self._add(("p", q, t))
    def reset(self, q): self._add(("reset", q))
    def measure(self, q, c): self._add(("measure", q, c))
    def append(self, gate, qargs): self._add(("oracle", qargs[0], gate))
    def if_test(self, cond): return _IfCtx(self, cond)
