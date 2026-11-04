# The Silent Group — Qiskit Fall Fest 2026 (Track 2 Submission)

Complete submission package for **Track 2: The Silent Group** ([ProAltro/qff-hackathon](https://github.com/ProAltro/qff-hackathon)), including the 1-qubit semiclassical quantum solver (`my_solver.py`), the Open Track extension on **Abelian Subgroup Structure Recovery** (`open-track/`), the 16-slide presentation deck (`The_Silent_Group_deck.pdf`), and the local benchmark suite (`testing/`).

---

## Package Layout

```text
silent-group/
├── my_solver.py                  # Auto-graded Track 2 solver: solve(oracle) -> {"order": int, "secret": int | "none"}
├── The_Silent_Group_deck.pdf     # 16-slide presentation (13 main slides + 3 technical appendix slides)
├── README.md                     # Documentation and reproduction commands
├── open-track/
│   ├── subgroup_structure.py     # Main Open Track extension: find_subgroup_structure(oracle) within 80k budget
│   ├── analysis.py               # Supporting characterization: 80k budget sweep, cost model, & noise sweep
│   ├── analysis_results.json     # Cached empirical metrics from full benchmark runs
│   ├── subgroup_structure_analysis.png
│   └── open_track_analysis.png
└── testing/
    ├── test_public.py            # Multi-trial benchmark on public_1 and public_2
    ├── sweep.py                  # 8-category (640-meeting) stress test across group families
    ├── test2.py                  # Wrapper for sweep.py across 8 group categories
    ├── test.py                   # 200-trial quick regression suite
    └── mock/                     # Local mock oracle for offline batch sweeps
```

---

## 1. Core Solver (`my_solver.py`) — Auto Score (60%)

### Key Algorithmic & Engineering Features (Directly from `my_solver.py`)
1. **1-Qubit Semiclassical Phase Estimation (`QuantumCircuit(1, total)`):**
   Every circuit uses **1 control qubit** (`qc.reset(0) -> qc.h(0) -> oracle_gate(X, mult << (k-1)) -> classical if_test phase feedback up to depth CORR_DEPTH = 12 -> qc.h(0) -> qc.measure(0, bit)`), while the group register remains implicit inside the oracle. Keeping `CORR_DEPTH = 12` caps circuit feedback complexity so each meeting runs in **0.12s–0.15s** mean wall-clock time on the official Qiskit simulator.
2. **Joint Shot-1 Pre-Seeding (`find_order` $\to$ `solve_logs`):**
   Because the hidden group register starts in $\lvert e \rangle = \frac{1}{\sqrt{|G|}}\sum_{\chi \in \widehat{G}} \lvert \chi \rangle$ and is never reset within a single `oracle.run(qc, shots=1)` call, Shot 1 of `find_order()` estimates $\chi_1(g)$ to 124 bits ($2 \times 60 + 4$, costing 7,750 queries for $g$) and **simultaneously** estimates $\chi_1(A), \chi_1(B), \chi_1(C)$ to 64 bits on the exact same random character $\chi_1$. When Shot 1 recovers $n = \text{ord}(g)$ and $\gcd(j_1, n) = 1$, `solve_logs(n)` extracts $a, b, c \pmod n$ directly from Shot 1's pre-seeded phases with **zero extra QPE circuits** (~16.1k–17.8k total queries, **1.388–1.399 / 1.500** score).
3. **Partial-Quotient Jump + Single-Bit-Flip Convergent Recovery (`_cf_with_repair`) & Full Prime-Factor Order Stripping (`_minimise_order`):**
   Selects continued-fraction convergents $j/r$ of $y / 2^{124}$ using partial-quotient jumps ($a_{i+1} \ge 1200$) and tight Legendre ULP bounds, repairs single bit flips $y \pm 2^b$ when $a_{i+1} \ge 1800$, and strips **every prime factor** of the candidate multiple $c$ in logarithmic steps while skipping divisors of already-refuted exponents (`false_divs`). Stripping all prime factors of $c$ ensures that even on small orders where a noisy sample produces a random multiple of $n$, spurious odd prime factors are completely removed.
4. **Sequential Identity Test (`is_identity`), Congruence Forgery Screen, & Two-Stage Individual Verification:**
   - **Wald's SPRT Identity Test:** `H -> oracle_gate(X, k) -> H -> measure` stops as soon as `1`s refute identity ($X^k \ne e$ in ~2–4.5 shots) or confirms $X^k = e$ in ~15 shots ($\alpha = 10^{-5}, \beta = 5 \times 10^{-5}$).
   - **3-Anchor In-Circuit Parity Checks (`_add_checks`):** Appends 3 parity checks on orders with small odd part or long binary `0`/`1` runs (`pseudo_pow2`), reducing accepted bit-flip errors on $n = 2^{59}$ from **6.80% to 0.75%**.
   - **Lift-Free Congruence Forgery Screen ($X^d = g^{d a_0 \bmod n}$):** Even when $d = \gcd(j, n) > 1$ is large, tests $X^d \cdot g^{-(d a_0 \bmod n)} = e$ directly in Round 1 to refute forgeries (such as `B` on `public_2`) without enumerating $d$ lifts.
   - **Two-Stage Individual Element Verification ($X \cdot g^{-a_X} = e$):** Screens each candidate with a 2-shot quick test (10% miss rate, `id_miss = 0.10`) before running the full SPRT identity confirmation on each broadcast individually—guaranteeing soundness against multiple forgeries that share an order-2 or order-3 defect.
5. **Fail-Safe Reliability Wrapper (`solve`):**
   Top-level `try / except Exception` fallback (preserving `KeyboardInterrupt` and `SystemExit`), up to 8 free retries per `oracle.run` call on refused/transient runs, tolerant bitstring key parsing, and a 540-second wall-clock guard (against the 600s limit).

### Verified Performance
- **Official `silentgroup.check` (`--repeat 20`, run from official repo root outside `testing/`):**
  - `public_1` ($n \approx 2^{56}$, genuine): **20 / 20 correct**, mean score **1.385 / 1.500** (`18,499` mean queries; typical 1-shot run **1.393 pts / 17,079 queries**), mean wall-clock **0.13s / run**.
  - `public_2` ($n \approx 2^{35}$, forged): **20 / 20 correct**, mean score **1.387 / 1.500** (`18,129` mean queries; typical 1-shot run **1.399 pts / 16,117 queries**), mean wall-clock **0.11s / run**.
  - Overall official checker accuracy: **40 / 40 (100%)**, overall mean score **1.386 / 1.500**.
- **15% Random `oracle.run` Refusal Injection Test:** **20 / 20 correct, 0 crashes**.
- **Extended 50-Trial Public Benchmark (`testing/test_public.py 50 1.0`):**
  - `public_1`: **50 / 50 (100%)**, mean score **1.3816**, mean queries **18,950** (max `26,718`).
  - `public_2`: **50 / 50 (100%)**, mean score **1.3882**, mean queries **17,892** (max `26,590`).
  - **Combined public mean score: 1.3849 / 1.5000 (~1.385 / 1.500, 18,421 mean queries).**
- **8-Category 640-Meeting Sweep (`testing/sweep.py 0 3 80 1.0`):** **640 / 640 (100%, 0 wrong)** across all 8 group families, mean score **1.3655**.
- **200-Meeting Regression Suite (`testing/test.py`):** **200 / 200 (100%)** correct.

### How to Run on the Official Simulator
Run from the root of the official `ProAltro/qff-hackathon` repository (do **not** run inside `testing/`, where `testing/mock/silentgroup` would shadow the installed `silentgroup` package):
```bash
python3 -m silentgroup.check /path/to/silent-group/my_solver.py --repeat 20
```

---

## 2. Open Track Extension (`open-track/`) — Open Track (20%)

### Primary Extension: Abelian Subgroup Structure Recovery (`open-track/subgroup_structure.py`)
Recovers the full invariant factor decomposition
$$\langle g, A, B, C \rangle \cong \mathbb{Z}_{n_1} \times \mathbb{Z}_{n_2} \times \dots \times \mathbb{Z}_{n_r}, \qquad n_{i+1} \mid n_i$$
**strictly inside the 80,000-query competition budget** by:
1. **Reusing Joint QPE Trajectories:** Reuses the clean 4-element phase vectors $(\chi(g), \chi(A), \chi(B), \chi(C))$ already collected during `solve(oracle)`, requiring **zero additional quantum circuits** when all elements are genuine (`public_1`) and only a handful of supplemental rounds on on-grid forgeries (`public_2`).
2. **Pontryagin Dual Sampling & Leave-One-Out $p$-adic Smith Normal Form:** Every character of $H = \langle g, A, B, C \rangle$ extends to $|G|/|H|$ characters of $G$, so each joint 1-qubit QPE shot samples $\widehat{H} \cong H$ uniformly. Scaling phases by $M = \text{lcm}(\text{ord}(g), \text{ord}(A), \text{ord}(B), \text{ord}(C))$ and applying leave-one-out $p$-adic Smith Normal Form over each prime power $p^e \parallel M$ yields the exact invariant factors.

### Worked Verification on `public_2` (Appendix A, Slide 14)
On `public_2`, $\text{ord}(g) = 37,199,287,125 = 3^3 \cdot 5^3 \cdot 7^2 \cdot 11^3 \cdot 13^2$ and $B$ is an on-grid forgery with $\text{ord}(B) = 1,771,394,625$. Computing the $p$-adic SNF of $(g, B)$ modulo each prime power $p^e \parallel M$ gives:
- $p = 3$ ($p^e = 27$): `[27, 9]`
- $p = 5$ ($p^e = 125$): `[125, 5]`
- $p = 7$ ($p^e = 49$): `[49, 7]`
- $p = 11$ ($p^e = 1331$): `[1331, 11]`
- $p = 13$ ($p^e = 169$): `[169, 169]`
Combining invariant columns gives $n_1 = 27 \times 125 \times 49 \times 1331 \times 169 = 37,199,287,125$ and $n_2 = 9 \times 5 \times 7 \times 11 \times 169 = 585,585$, so $\langle g, A, B, C \rangle \cong \mathbb{Z}_{37,199,287,125} \times \mathbb{Z}_{585,585}$.

### Supporting Characterization (`open-track/analysis.py` & Appendix Slides 14–16)
- **Exact-Recovery Rate vs Query Budget Curve (15 runs per family, $\pm 13\%$ binomial SE):** 8 of 10 group families reach **100% (15 / 15)** exact structure recovery inside the 80,000-query limit, while documenting the two exact failure boundaries: rank-2 even-order span tail (`Zn×Zn forged B`: **93.3% / 14 of 15** at 80k, matching the $T=6, p=2$ theoretical $4.6\%$ span tail) and multi-forged 58-bit off-grid orders (`big forged ABC`: **53.3% / 8 of 15** at 80k $\to$ **86.7%** at 100k).
- **Analytical Query-Cost Model:** Verified across $L \in [8, 60]$ with **3.9% RMS error** (range `-5.8%` to `+5.9%`).
- **Noise Scaling Sweep:** **100% correct from $0\times$ to $1.0\times$ stated noise (160 meetings per level)** and **98.8%** (158 / 160) at $1.5\times$ and $2.0\times$ noise (plus **640 / 640** on the 8-category sweep at $1.0\times$).

```bash
python3 open-track/subgroup_structure.py
python3 open-track/analysis.py
```
