# Gated stages in driver.py (gates.py) — Nathan A. Meeks & Grok. [HEUR-gate] / [NUM] candidate, unrefereed; not RH.
Both gates are search/diagnostic aids; acceptance is always Lemma CP on the whole window. Every firing → logs/gates.jsonl (tag HEUR-gate) and the ledger `notes` column.

## 0/0 gate (C0 removable singularity; Dedekind cut, Nathan A. Meeks)
- The enclosure itself is in v9p (C0_FIX.md): when |h|=|p−p0|<0.01, p0∈{1/4,3/4}, Ψ is computed from the sinc form as a proven cut L≤Ψ≤U; outside, the original quotient.
- driver: `gates.c0_gate` checks per window whether p(t)=√(t/2π) passes within 0.01 of n+1/4 or n+3/4 and logs a firing (gate C0_cut, t of p0). It does not alter the run.

## Close-pair gate (centre-out, Nathan A. Meeks) — chunked reuse is the DEFAULT
- Fires only after the normal v9 ladder + MAXR=200 rung ends UNPROVEN with Turing − certified = 2 or 4, Na>0, brackets = certified count. Disable: GATE_CLOSEPAIR=0.
- No full-window rerun. The original certified brackets are kept. Bisection: at the midpoint of the current [L,R] run a small certified chunk (~320 t) whose Lemma CP gives
  exact N at its ends (cut points snapped off original brackets); compare with original bracket counts left/middle/right; recurse into the deficient part. Stops when [L,R]≤640
  or a chunk itself holds the whole deficit.
- Only that chunk [L,R] is rerun (delta 0.02 → 0.005 → 0.0005, MAXR=400); it must be PROVEN with N(L), N(R) equal to the bracketing exact values. Its brackets replace the
  original ones lying inside [L,R]. Accept iff spliced count = Turing count of the window, brackets sorted & non-overlapping, recovered = deficit, and the 40-t overlap with the
  previous window agrees (ovB ≤ ovN ≤ ovB+1). Otherwise FAILED.md and stop (no skip).
- Rigor note: chunk-boundary counts ARE certified — each is v9p's own Lemma CP (Turing bounds at chunk ends that coincide because sign changes reach them), so N(L),N(R) are exact.
  Cut points never fall inside an original bracket, so no zero is double-counted. If a probe chunk fails to prove, the gate fails (no fallback to a full rerun, per Nathan).
- If the deficit is split between parts, the gate reruns the current [L,R] as one chunk (smallest certified sub-interval it can justify).

## Tests (offline, scratch copies scratch/gatetest_8758, scratch/gatetest_10202; real ledger untouched)
- 8758 from pre-fix checkpoint: ladder fails 44169/44171 → gate: chunk u[8302,8827], +2 at t≈3000090783911.63 (matches WINDOW8758.md) → PROVEN, ov 173/173. C0 gate fired (p=3/4 near t≈3.00009075e12 within |h|<0.01).
- 10202: ladder fails 44167/44169 → chunk u[7542,7866], +2 at t≈3000105627072.83 → PROVEN, ov 170/170.

## OPEN-3 patch (2026-10-30, Nathan A. Meeks & Grok)
ZDUMP endpoints now %.17g (exact certified doubles); no brackets dumped if any certified-sign disagreement. Gate reuse of pre-patch brackets is backed by Arb endpoint re-certification (spotcheck/fixwin.py). [NUM]

## Independent-check gates for OPEN-1 (G-cut → G-Turing → G-RS) — Nathan A. Meeks & Grok, approved by Nathan 2026-10-30
- Code: spotcheck/indgates.py, wired in by gates.ind_gate in driver.py. Disable with GATE_IND=0; IND_PROCS sets the number of worker processes (default 6, nice, cpus 0-5).
- Triggers: the close-pair gate fired (block = the spliced chunk only; a failed close-pair → whole window); v9 refine_disagreements>0; overlap mismatch; idx listed in ind_blocks.txt (Nathan names the blocks).
- Order is by measured cost, cheapest first (see the test numbers below): G-cut on known-uncertain points, then G-Turing (fails fast on a mismatch), then G-RS (with G-cut run on any Arb-uncertain endpoint).
- Resolve rule (Rule G5, WP §G.5): independent N(lo) and N(hi) exact (and equal v9's Na, Nb for whole windows; chunk blocks: bracket count = N(hi)−N(lo)); G-RS certifies exactly Nb−Na disjoint opposite-sign brackets; 0 conflicts; 0 unresolved cuts. No single gate certifies a block by itself: G-cut and G-RS give signs only, and G-Turing gives the count only. A failure writes FAILED.md and stops the run, with no gap skipped.
- Journal: one row per triggered block in fixes_ledger.csv (gate=independent:…, result with cost per gate) and in logs/gates.jsonl (gate=independent).

### OPEN-1 gate tests and install, 2026-10-30 (MST) [NUM], candidate, unrefereed; not RH (Nathan A. Meeks & Grok)
- Offline scratch tests in chunk-block mode, IND_PROCS=2 (spotcheck/ind_test2.jsonl). G-Turing mesh now includes known bracket endpoints within ±25 of each block end. This is search-only and fixes an ambiguous N(lo) next to a close pair (try 1: 10202 N(lo)∈[..129,..130] → FAIL, correctly refused).
  - 10202 [3000105627062,3000105627386]: RESOLVED. G-Turing 208 s; G-RS 222 s, 1388/1388; 0 G-cut.
  - 8758 [3000090783750,3000090784070]: RESOLVED. G-Turing 159 s; G-RS 224 s, 1370/1370.
  - 34157 [3000351881920,3000351882242]: RESOLVED. G-Turing 202 s; G-RS 143 s, 1378/1378.
- Live, whole window (old trigger): 34158 RESOLVED. G-Turing 63 s; G-RS 2509 s, 44169/44169 over 51569 endpoints; 6 procs.
- Measured cost order: G-cut ≈0.2–0.5 s per uncertain point (none needed so far); G-Turing ≈1–3.5 min per block (2 ends); G-RS ≈0.17 s per endpoint (chunk ≈2.5–4 min, whole window ≈42 min).
- Installed 14:50 MST in gates.py, driver.py and spotcheck/indgates.py. Backups are in scratch/*.bak_1440. The changes take effect at the next run_batch, and the running batch keeps its old in-memory code.

## Trigger rule (t90, Oct 30 2026)
Error-only: v9-PROVEN windows (incl. MAXR-only) get no gates; v9_sign_uncertain triggers only if the window is not proven; close-pair windows get independent checks on the repaired chunk only. OPEN-1b resolved (Arb ball arithmetic; trust basis = Arb correctness, not formally verified) — see WHITEPAPER G.5.

### G.5.t92 Optimized error-only fixer (version f) — Nathan A. Meeks & Grok. [NUM] candidate, unrefereed; not RH.
The resolve rule is unchanged: a window that had an error counts only when Lemma CP holds on the whole window (v9) AND, on the repaired block, independent G-Turing N at both ends, G-RS signs at every bracket endpoint, and G-cut all agree. No independent check was dropped. Changes, each with its soundness argument:
1. **Fast path to the close-pair gate** (driver.py, v9run.sh `V9_FIRST_ONLY`): if attempt 1 is short by 2 or 4 and every bracket it found is certified, go straight to the close-pair gate (Nathan's center-out method). If the gate does not splice, the full v9 ladder and the MAXR=200 rung run exactly as before. *Sound:* the ladder is search-only, and acceptance still needs full-window Lemma CP.
2. **Smaller chunk** (gates.py `CH` 320 → 40, env `GATE_CH`): the bisection on exact chunk-boundary counts continues down to ≤80 wide (6–66 wide on tests). *Sound:* every boundary N is v9-certified, the splice is accepted only if it closes Lemma CP on the window, and the independent gates re-derive N at the block ends.
3. **Adaptive G-Turing** (indgates.gturing_fast): try h=12 with mesh 0.1 first, then fall back to h=24 with mesh 0.05. *Sound:* N is accepted only when the certified Trudgian lower and upper bounds coincide; the mesh is search-only.
4. **G-RS on the small block, 7 workers.** Every endpoint in the block is still checked (no sampling), because the resolve rule needs it.
Not applied: (a) lowering Arb precision to 96 bits, because acb.zeta costs about the same at 64, 96 and 128 bits (≈0.16 s), so it gains nothing; (b) reusing v9's N at chunk ends instead of G-Turing, because that would remove the independence the rule requires; window-end G-Turing values aren't cached from previous windows. G-Turing is now the floor (≈30–50 s), and the remaining speedup is the Platt multi-evaluation proposal below.
Measured on scratch copies under load (pair brackets removed, then re-found): 10202 380 s → 40 s (gates), 8758 47 s, 27053 77 s; end-to-end driver on 10202: 93 s (8 s window + close-pair, 49 s G-Turing, 35 s G-RS), versus about 534 s before.

## t96: Platt fast path for G-Turing + G-RS (Nathan A. Meeks & Grok) — [NUM], candidate, unrefereed; not RH
- Wrapper: `/workspace/opt/flint3/platt_iso` (C, `platt_wrap.c`), linked against the FLINT 3.6.0 already bundled with python-flint 0.9.0 (no build needed; system python-flint untouched). Python side: `spotcheck/plattgate.py`; wired into `indgates.resolve` as step 2' (env `PLATT=1` default, `PLATT_THREADS=4`, `PLATT_MAXZ=4000`).
- API: `acb_dirichlet_platt_local_hardy_z_zeros(res, n, len, prec)` (FLINT docs, acb_dirichlet "Riemann zeta function zeros (Platt's method)"; D. J. Platt, "Isolating some non-trivial zeros of zeta", Math. Comp. 86 (2017) 2449–2467). Returns k ≤ len certified balls for zeros n..n+k−1 of Z(t) (labels certified by Turing's method inside FLINT; refinement by Illinois on Platt's interpolant with rigorous error bounds). Endpoints are passed out exactly (mantissa·2^exp).
- Assumptions: documented range 10^4 ≤ n ≤ 10^23 (we are at n≈1.236e13). Grid/interpolation tuning is FLINT's "heuristic context": it affects only success and speed, not soundness. No RH assumption found in the code path used. Soundness rests on FLINT's implementation of Platt's bounds (not formally verified; same trust class as Arb).
- Accept rule (G-Platt): request zeros Na−2..Nb+3; require all returned, ball(Na) < lo < ball(Na+1), ball(Nb) < hi < ball(Nb+1) (⇒ N(lo)=Na, N(hi)=Nb independently), and a one-to-one sorted containment of the Nb−Na Platt balls in the v9 brackets of [lo,hi]. Anything else → refusal → the existing Arb G-Turing → G-RS path runs unchanged. Resolve rule unchanged.
- Validation (spotcheck/platt_validate.jsonl): 10202 chunk 1388/1388 matched, N ends = Arb G-Turing; 27053 chunk 1370/1370, N ends = Arb G-Turing (…510295 / …511665); 8758 live spliced brackets 1370/1370, N ends = Arb. 8758 scratch pre-fix copy (1369 brackets, one straddling hi) correctly refused. Negative tests: one bracket removed → refused; wrong N(lo) → refused.
- Cost: ~50 s per 320-wide chunk (4 threads, under batch load) replacing G-Turing+G-RS ≈ 84 s; per-error ≈ 93 s → ≈ 59 s.
