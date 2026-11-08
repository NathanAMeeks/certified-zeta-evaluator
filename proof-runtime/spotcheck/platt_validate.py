# [NUM] validate Platt fast path vs Arb results on scratch copies (Nathan A. Meeks & Grok). Candidate, unrefereed; not RH.
import json, sys, time, indgates, plattgate
from fractions import Fraction as F
R = '/workspace/rh_extend/'
cases = [  # idx, bracket file, T, chunk lo, hi, Na (Arb G-Turing, ind_test2.jsonl) or None
 (8758, R+'scratch/gatetest_8758/brackets/w0008758_T3000090775200.txt.xz', 3000090775200, 3000090783750, 3000090784070, 12362791570331, 12362791571701),
 (10202, R+'scratch/gatetest_10202/brackets/w0010202_T3000105619520.txt.xz', 3000105619520, 3000105627062, 3000105627386, 12362855099129, 12362855100517),
 (27053, R+'brackets/w0027053_T3000278847800.txt.xz', 3000278847800, 3000278855340, 3000278855660, None, None)]
for idx, fn, T, lo, hi, Na, Nb in cases[int(sys.argv[1]):int(sys.argv[2])]:
    allb = indgates.read_br(fn); br = [x for x in allb if T + F(x[0]) >= lo and T + F(x[1]) <= hi]
    src = 'Arb G-Turing'
    if Na is None:  # 27053: v9 ledger N(window lo) + brackets below chunk; Arb G-Turing run below for agreement
        Na = 12363596478367 + sum(1 for x in allb if T + F(x[1]) < lo and T + F(x[0]) >= 3000278847880); Nb = Na + len(br); src = 'v9 ledger+brackets'
    s = time.time(); g = plattgate.gate(T, lo, hi, br, Na, Nb); tp = time.time() - s
    out = dict(idx=idx, block=[lo, hi], brackets=len(br), Na=Na, Nb=Nb, Na_src=src, platt=g, platt_sec=round(tp, 1))
    if src != 'Arb G-Turing':
        from multiprocessing import Pool
        ext = [indgates._exact(T, u) for ab in allb for u in ab if abs(T + float(u) - lo) < 25 or abs(T + float(u) - hi) < 25]
        with Pool(4) as p: a = indgates.gturing_fast(lo, p, ext); b = indgates.gturing_fast(hi, p, ext)
        out['arb_turing'] = [a[0], b[0]]
    print(json.dumps(out), flush=True)
