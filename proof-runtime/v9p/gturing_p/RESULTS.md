# turing_fast — Nathan A. Meeks & Grok. [NUM] candidate, unrefereed; not RH. Needs independent audit + 4 replays before adoption.
Box: 8 cores, load 14–16 throughout (live run), no GPU (`nvidia-smi` not found). Timings are noisy.

## A. G-Turing-P (pme.c + gturing_p.py)
Mesh signs over [a-12,a+12] from FLINT `acb_dirichlet_platt_multieval` (A=16, B=4096, FLINT's own 1e7–5e22 params; prec 128),
i.e. Platt DFT/Gaussian-interpolation path: different code path from rsbatch.c (G-RS) and from acb_zeta. Sign of Platt's scaled
Lambda = sign Z (positive factor). Same-sign |value| dips refined with Arb acb_zeta points (search only). Trudgian bound as in indgates.
Accept N(a) only if certified lower=upper bound; any uncertified mesh sign -> refuse -> existing Arb gturing_fast (fallback).
| end | Arb gturing_fast (current, 1 proc) | G-Turing-P | N |
| 10202 lo 3000105619600 | 284.8 s, 1208 Arb evals | 84.6 s, 0 Arb | 12362855067193 = ledger Na |
| 10202 hi | – | 59.1 s, 0 Arb | 12362855111362 = Nb |
| 8758 lo/hi | – | 87.7 / 78.4 s (12 Arb) | 12362791534078 / ...578249 = ledger |
| 27053 lo/hi | – | 85.8 / 78.0 s (8 Arb) | 12363596478367 / ...522535 = ledger |
| 38498 lo/hi | – | 85.9 / 80.7 s (22 Arb) | 12364100036980 / ...081149 = ledger |
8/8 ends certified, 8/8 equal ledger Na/Nb (and window Turing counts 44171/44169/44168/44169). ~3.4–4.8x less CPU per end vs measured current.
First try A=8 (mesh 0.125) refused 7/8 (bounds did not meet) — correctly refused, no wrong N. FLINT multieval threads gave no gain on the loaded box.
Caveat: one pme run is ~60–85 s wall single-core; it wins on CPU (what a saturated pipeline pays), not latency, unless the two ends run on 2 cores.

## Split proof note (Lemma CP-split)
Cuts a=c0<c1<...<cm=b, each with Z(ck) certified nonzero (sign ball excludes 0) and N(ck) certified exactly (Turing lower=upper).
N(b)-N(a) = sum_k [N(c_{k+1})-N(c_k)] (telescoping), and [N(c_{k+1})-N(c_k)] = # zeros in (c_k,c_{k+1}] exactly.
If piece k has s_k certified sign changes of Z with s_k = N(c_{k+1})-N(c_k), then (sign changes are distinct zeros on the line, each
counted >=1) every zero in that piece is simple and on the line (Lemma CP on the piece). Nonzero Z at cuts => no zero on a cut, no double count.
Hence piece-wise certification <=> whole-range certification, total = sum. Interior N(ck) is computed once and shared by both neighbours.
Each piece is certified individually; cut placement is search only (Lemma CP-S). Status: proof sketch, needs audit.

## B. Framework overhead (64 Arb Z evals, 4 workers)
serial ~10.9 s; C OpenMP 2.7–3.3 s; Python Pool 2.6–6.3 s; Ray cold default 22.7 s (5.0 s init), Ray cold coarse tuned 12.6 s (3.2 s init);
Ray warm fine 2.8 s (= Pool). Causes: ray.init + cold worker spawn + per-worker `import flint`; num_cpus default=8 on a box already at load 15
(oversubscription); one-task-per-point scheduling; object store spilling to disk if /dev/shm small (here 7.9 G, fine).
Fixes: long-lived warm workers (actors, init once per run), num_cpus=cores minus live-run use, coarse tasks (one sub-range per task),
OMP/OPENBLAS/FLINT_NUM_THREADS=1, object_store_memory small and _temp_dir on /dev/shm (no /tmp/ray disk spill), pin with taskset.
Verdict: on one box, a plain C OpenMP or multiprocessing pool is as fast; Ray only pays off multi-node.

## C. Independent audit of G-Turing-P (turing_audit/, Nathan A. Meeks & Grok) — verdict SOUND [NUM], OPEN-P CLOSED
Audited copies in ../turing_audit (pme binary sha256 c82df083… = original). pme links FLINT 3.6.0 (rh-research/tools/flintc), source fetched (flint-3.6.0); diffs vs 3.1.3 in platt_* are cosmetic.
1. OPEN-P [NUM]: `acb_dirichlet_platt_multieval` (src/acb_dirichlet/platt_multieval.c) assumes no prior verification. Every truncation or
   discretisation error is added as an Arb ball radius: Lemma A5 :582, A7 :599, 3.2 :617, B1 :621, B2 :627, A9 :636, A11 :646 (Platt 2017, App. A/B).
   These are unconditional analytic bounds with no RH, known N(T), or Turing input. When a lemma's hypothesis fails, the lemma returns ±inf
   (A5:43, A11:197, A7:112, A9:149, B1:23, c_bound:211 sigma+1/2<=t0), so the sign is uncertified and G-Turing-P refuses. RH, Gram blocks, and
   Turing's bound live only in the zero-finding layer (platt_local_hardy_z_zeros.c, platt_i_bound.c `acb_dirichlet_turing_method_bound`), and pme never calls them.
   The parameter fits (local_hardy_z_zeros.c, range [1e7,5e22]) only tune speed. Docs (acb_dirichlet.rst:776-806): output is e^{πt/4}Λ(t) = (positive)·Z(t), grid T-B/2+i/A, which matches pme's u=(i-N/2)/A.
2. Turing rule re-derived [NUM]: N(t)=θ/π+1+S(t) a.e. Trudgian (Math. Comp. 80 (2011)) gives |∫_{t1}^{t2}S| <= 2.067+0.059 log t2 for t2>t1>168π.
   It is used with t2=a+h (upper bound) and t2=a (lower bound), which is correct. Each certified sign-change cell contributes h-r (or l+h), an underestimate, so the bound is conservative.
   Several zeros in one cell count once, which is also conservative. θ antiderivative checked by differentiation; the remainder after 1/(48t) is about 7/(5760t^3), well below the (b-a)/a^3 used.
   Floor/ceil via float is safe: N<2^53, integers are exact doubles, and rounding is monotone. Acceptance requires the certified lower and upper bounds to give the same integer, so the mesh is search-only. Lemma CP-split re-derived (telescoping +
   Z(c_k)≠0 certified at cuts + s_k=ΔN_k ⇒ simple, on-line, no double count): correct.
3. Independence [NUM]: gturing_p imports only flint/subprocess. Mesh comes from Platt DFT (pme). Refinement uses acb_zeta+lgamma. rsbatch_t100.c uses
   double RS + arb/acb_hypgeom only, no platt and no acb_zeta. Neither reads library, brackets, or finder output. (Note: the live G-Turing does read finder brackets as search hints.)
4. Replays [NUM]: 8/8 ends certified, 8/8 = ledger Na/Nb (8758, 10202, 27053, 38498), 0 wrong, 0 refusals.
5. Cost, own-process CPU via getrusage (perf not installed; cycles = CPU-s × 2.4 GHz nominal, so [HEUR]); 1 core, box load ~13–16:
   | gate (per end, window 10202 unless noted) | CPU-s/end | CPU-s/window (2 ends) | cycles/zero (44,169 zeros) |
   | G-Turing current, full mesh (1208 acb_zeta evals, 0.149 s each) | 179.9 | ~360 | ~1.96e7 |
   | G-Turing current, live mode HYB_TMESH=1 (meshless with 214 finder hints failed, then mesh; 1509 evals) | 237.3 | ~475 | ~2.6e7 |
   | G-Turing-P (pme + Arb refinement, 8 ends) | 51.0–56.1 (mean ~53.6) | ~107 | ~5.8e6 |
   | finder clean window (ledger sec_setup+nufft+points, wall in 1-thread worker) [HEUR] | – | ~0.27 | ~1.5e4 |
   | G-RS rsbatch close-pair chunk (t100 replay, 4 threads, wall) [HEUR] | – | 3.1–14.9 | – (1-core full-window meter still running) |
   Puzzle [HEUR]: 285 s was wall time for 1208 evals on a loaded box. CPU time is 180 s (0.15 s per eval). The live ~17 s is wall time over the 8-worker pool,
   and clean windows pass the meshless config (~216 evals ≈ 32 CPU-s). On close-pair windows like 10202 the meshless config fails, and the gate costs ~237 CPU-s per end.
   G-Turing-P: 3.4× less CPU than full mesh, 4.4× less than the live path on close-pair windows. It would lose on clean windows, where the meshless config wins (~32 vs ~54 CPU-s).
   Recommendation: order it as meshless Arb first, then G-Turing-P, then the Arb mesh as the last fallback.
Verdict: SOUND (no fix needed). OPEN-P CLOSED for platt_multieval. Remaining trust is in the FLINT implementation of Platt's lemmas (published, tested code).

### C.5b Full per-gate cost breakdown, own process CPU (getrusage, 1 core, nice 19), window 10202, 44,169 zeros [NUM; cycles = CPU-s × 2.4 GHz nominal, HEUR]
| stage | CPU-s per window | cycles/zero (÷44,169) | notes |
| finder fftzi_v9, attempt 1 (δ=0.24) | 1.56 (10202 P=131072); 1.62 / 2.14 (8758 P=131072 / 262144) | ~8.5e4 | wall 2–152 s on loaded box; CPU is what counts |
| G-RS rsbatch_t100, 1 thread, close-pair chunk u[7542,7866] | 20.09 (1388 brackets, 2773 pts, 0 uncertified) | ~1.1e6 (3.5e7 per checked zero) | only runs on gated chunk |
| G-Turing current (acb_zeta), 2 ends | 360 full mesh / 475 live path | ~2.0e7 / ~2.6e7 | from C.5 |
| G-Turing-P, 2 ends | ~107 | ~5.8e6 | from C.5 |
On a gated window the Turing gate dominates: about 96% of CPU with the current gate, about 83% with G-Turing-P. (Earlier full-window meter_grs run was too slow and was stopped; replaced by this chunk-only run.) Nathan A. Meeks & Grok.
