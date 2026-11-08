# [HEUR-gate] Gated stages for driver.py — Nathan A. Meeks & Grok. [NUM] candidate, unrefereed; not RH. See GATES.md.
# Close-pair gate (method: Nathan A. Meeks, centre-out): chunked reuse of the ORIGINAL certified brackets; no full-window rerun.
import json, math, os, subprocess, datetime
def now(): return datetime.datetime.now().astimezone().isoformat(timespec='seconds')
FLCOLS = ['idx','window_lo','window_hi','gate','trigger','chunk_interval','recovered_brackets','method_ref','lemma_ref_gturing_p','lemma_ref_lemma_h','result','timestamp']   # t106: +lemma refs (Research option A)
import fcntl, threading
_LK = threading.Lock()
def _locked(fn, line):
    with _LK, open(fn + '.lock', 'w') as L:   # t98: serialize log appends across parallel workers
        fcntl.flock(L, fcntl.LOCK_EX)
        with open(fn, 'a') as f: f.write(line); f.flush()
LEM_GTP = 'G-Turing-P (v9p/gturing_p, SOURCE.sha256; WP G.5.t96 Lemma G5-P labels)'
LEM_H = 'LEMMA_EXT.md Lemma H (Arb pair certificate: sgn Z(a)=sgn Z(b)!=sgn Z(c) => 2 sign changes)'
def flog(H, **k):
    """[HEUR] t106 RECORD-ONLY fixes ledger (Nathan A. Meeks & Grok; Research option A). Appends one row per gate firing to {H}/fixes_ledger.csv
    with O_APPEND + fsync. Never returns a value, never raises, never alters a verdict: any error is logged to {H}/logs/fixes_ledger.err.
    FIXLEDGER=0 turns it off. A scratch run (gates.py or H under scratch/) refuses to write a live ledger; a live run writes it (t106d)."""
    if os.environ.get('FIXLEDGER', '1') != '1': return
    try:
        import csv, io; fn = os.environ.get('FIXLEDGER_PATH', f'{H}/fixes_ledger.csv'); rf = os.path.realpath(fn)
        # t106d guard fix: refuse a live-ledger write only from a SCRATCH run (this gates.py or H under scratch/); live runs may write it
        _S = '/workspace/rh_extend/scratch/'; _scr = os.path.realpath(os.path.abspath(__file__)).startswith(_S) or os.path.realpath(H).startswith(_S)
        if _scr and rf.startswith('/workspace/rh_extend/') and not rf.startswith(_S): raise PermissionError(f'refusing live ledger {rf} from scratch run')
        # t106d F1 (Audit): a live run may append ONLY to the live fixes ledger (symlinks resolved); any other FIXLEDGER_PATH is refused
        if not _scr and rf != os.path.realpath('/workspace/rh_extend/fixes_ledger.csv'): raise PermissionError(f'live run: refusing ledger path {rf} (only rh_extend/fixes_ledger.csv)')
        # t106d F2 (Audit re-check F2a): width comes from the EXISTING file's header (10 = pre-t106, 12 = FLCOLS); FIXLEDGER_NCOLS (10|12, default 12)
        # applies only when the file is new. An explicit FIXLEDGER_NCOLS that disagrees with the header, an unknown header, or any row whose width
        # differs from the header is refused and logged (never a mixed-width ledger).
        _env_nc = os.environ.get('FIXLEDGER_NCOLS'); _nc = _env_nc or '12'
        if _nc not in ('10', '12'): raise ValueError(f'FIXLEDGER_NCOLS={_nc}')
        C10 = [c for c in FLCOLS if not c.startswith('lemma_ref_')]
        k.setdefault('lemma_ref_gturing_p', LEM_GTP if 'G-Turing-P' in str(k.get('tiers', '')) or 'Platt' in str(k.get('gate', '')) else 'n/a')
        k.setdefault('lemma_ref_lemma_h', LEM_H if 'Lemma H' in str(k.get('method_ref', '')) + str(k.get('trigger', '')) else 'n/a')
        _r = [k.get(c, 'unknown') if k.get(c) not in (None, '') else 'unknown' for c in FLCOLS[:-1]] + [now()]
        with _LK, open(fn + '.lock', 'w') as L:
            fcntl.flock(L, fcntl.LOCK_EX)   # header read + append under the same lock (no create race, no swap race)
            new = not os.path.exists(fn) or os.path.getsize(fn) == 0
            if new: hdr = FLCOLS if _nc == '12' else C10
            else:
                with open(fn, newline='') as f0: hdr = next(csv.reader(f0), [])   # re-read every write (one small read; no stale cache after a swap)
                if hdr not in (FLCOLS, C10): raise ValueError(f'unknown ledger header ({len(hdr)} cols) in {rf}')
                if _env_nc and int(_env_nc) != len(hdr): raise ValueError(f'FIXLEDGER_NCOLS={_env_nc} but existing ledger header has {len(hdr)} cols')
            row = [v for c, v in zip(FLCOLS, _r) if c in hdr]
            if len(row) != len(hdr): raise ValueError(f'row width {len(row)} != header width {len(hdr)}; refused')
            buf = io.StringIO(); w = csv.writer(buf)
            if new: w.writerow(hdr)
            w.writerow(row)
            fd = os.open(fn, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
            try: os.write(fd, buf.getvalue().encode()); os.fsync(fd)
            finally: os.close(fd)
    except Exception as e:
        try: open(f'{H}/logs/fixes_ledger.err', 'a').write(f'{now()} fixes_ledger write failed (non-fatal, record-only): {type(e).__name__}: {e}\n')
        except Exception: pass
def _recovered(T, new, orig):
    # t106: recovered brackets = brackets in the gate's RETURNED list that were not in its input list (from the actual return value)
    if new is None: return 'none (gate returned None)'
    o = {(float(a), float(b)) for a, b in (orig or [])}; r = [(float(a), float(b)) for a, b in new if (float(a), float(b)) not in o]
    return f'{len(r)}: ' + ';'.join(f'{T+a:.4f}..{T+b:.4f}' for a, b in r)
def glog(H, **k):
    if k.get('gate') == 'close_pair' and 'T' in k:
        T = k['T']; ch = k.get('chunk')
        flog(H, idx=k.get('idx'), window_lo=T + 80, window_hi=T + 10400, gate='close_pair_chunked',
             trigger=f"CP short by {k.get('deficit')}", chunk_interval=(f"u[{ch[0]},{ch[1]}] = t[{T+ch[0]},{T+ch[1]}]" if ch else None),
             recovered_brackets=_recovered(T, k.get('returned'), k.get('orig')),
             method_ref=('Lemma H pair certificate (Arb signs); ' if 'Lemma H' in str(k.get('method', '')) else '') + 'WP §G.3 (Thm G3, chunked splice)', result=k.get('result'))
    k.pop('returned', None); k.pop('orig', None)
    _locked(f'{H}/logs/gates.jsonl', json.dumps(dict(ts=now(), tag='HEUR-gate', **k)) + '\n')
def c0_gate(H, i, T, A, B):
    """0/0 gate: report whether [T+A,T+B] has |p-p0|<0.01 for p0 in {1/4,3/4} (v9p psi_cut fires there; it is in-path in the evaluator)."""
    hits = []
    for t in (T + A, T + B):
        pass
    pa, pb = math.sqrt((T + A) / (2 * math.pi)), math.sqrt((T + B) / (2 * math.pi))
    for n in range(int(pa) - 1, int(pb) + 2):
        for p0 in (0.25, 0.75):
            c = n + p0
            if pa < c + 0.01 and pb > c - 0.01: hits.append(dict(p0=p0, t_center=round(2 * math.pi * c * c, 1)))
    if hits:
        glog(H, gate='C0_cut', idx=i, T=T, hits=hits, action='v9p Dedekind-cut enclosure (|h|<0.01) in path')
        flog(H, idx=i, window_lo=T + A, window_hi=T + B, gate='c0_cut', trigger=';'.join(f"|h|<0.01 near p={h['p0']} (t~{h['t_center']})" for h in hits),
             chunk_interval='n/a (in-path)', recovered_brackets='n/a', method_ref='WP §G.1 (Lemma G1)', result='fired; window status in ledger.csv row idx')
    return bool(hits)
def _run(V9BIN, T, delta, b, env, zd, A=80):
    need = 1.2 * (b + 140) / delta; P = 1024
    while P < need: P *= 2
    p = subprocess.run(['taskset', '-c', env.get('CPU', '7'), V9BIN, str(T), str(delta), str(P), str(A), str(b), '60', '/dev/null', 'dbl'],
                       env=dict(env, ZDUMP=zd, MAXR='400'), capture_output=True, text=True)
    try: d = json.loads(p.stdout.strip().splitlines()[-1])
    except Exception: d = {}
    br = [tuple(map(float, l.split()[:2])) for l in open(zd)] if os.path.exists(zd) else []
    ok = d.get('status') == 'PROVEN' and d.get('found') == len(br) and int(d.get('Nb_upper', 0)) - int(d.get('Na_lower', 0)) == len(br)
    return ok, int(d.get('Na_lower', 0) or 0), int(d.get('Nb_upper', 0) or 0), br, P
def chunk(V9BIN, T, lo, hi, env, zd, deltas=(0.24, 0.05, 0.02, 0.005, 0.0005)):
    """certify chunk u in [lo,hi] (integers, relative to T); returns exact N(T+lo), N(T+hi) and brackets (in window u) or None."""
    for dl in deltas:
        ok, na, nb, br, P = _run(V9BIN, T + lo - 80, dl, hi - lo + 80, env, zd)
        if ok: return dict(Na=na, Nb=nb, br=[(a + lo - 80, b + lo - 80) for a, b in br], delta=dl, P=P)
    return None
def snap(x, br):
    """move integer cut x upward until no original bracket contains it."""
    x = int(math.ceil(x))
    for a, b in br:
        if a < x < b: x = int(math.ceil(b)) + 1
    for a, b in br:
        if a < x < b: return snap(x, br)
    return x
def cnt(br, lo, hi): return sum(1 for a, b in br if a >= lo and b <= hi)
def close_pair(H, V9BIN, i, T, A, B, Na, Nb, br, env, scr, CH=int(os.environ.get("GATE_CH", "40"))):  # t98b multi-chunk (Nathan A. Meeks & Grok)
    """br: original certified brackets (u). Recursive centre-out bisection on exact (v9-certified) chunk-boundary counts;
    every part with a positive deficit is bisected further, down to <= 2*CH wide; each deficient chunk is rerun finely and
    spliced. Never reruns the whole window (Nathan's rule). Accept only if every chunk's N matches its boundary N, the splice
    is disjoint, and the total equals the window's Turing count (Lemma CP). Returns (spliced brackets, info) or (None, info)."""
    os.makedirs(scr, exist_ok=True); br = sorted(br); tur = Nb - Na; deficit = tur - len(br)
    info = dict(gate='close_pair', idx=i, T=T, deficit=deficit, method='Nathan A. Meeks centre-out; multi-chunk reuse (t98b)', probes=[])
    if deficit not in (2, 4): info['result'] = 'not_eligible'; return None, info
    k = [0]; segs = []
    def loc(L, R, NL, NR, depth=0):
        d = (NR - NL) - cnt(br, L, R)
        if d < 0: raise ValueError('negative_deficit')
        if d == 0: return
        if R - L <= 2 * CH or depth > 40: segs.append((L, R, NL, NR)); return
        m = snap((L + R) / 2, br); hi = snap(m + CH, br)
        if not (L < m < hi < R): segs.append((L, R, NL, NR)); return
        c = chunk(V9BIN, T, m, hi, env, f'{scr}/probe{k[0]}.z'); k[0] += 1
        if c is None: raise ValueError(f'probe_unproven at u={m}')
        info['probes'].append(dict(m=m, hi=hi, dL=(c['Na'] - NL) - cnt(br, L, m), dM=(c['Nb'] - c['Na']) - cnt(br, m, hi), dR=(NR - c['Nb']) - cnt(br, hi, R)))
        loc(L, m, NL, c['Na'], depth + 1); loc(m, hi, c['Na'], c['Nb'], depth + 1); loc(hi, R, c['Nb'], NR, depth + 1)
    try: loc(A, B, Na, Nb)
    except ValueError as e: info['result'] = str(e); return None, info
    new = list(br); rec = 0; chunks = []; recs = []
    for j, (L, R, NL, NR) in enumerate(segs):
        if R - L > 4 * CH + 200: info['result'] = f'segment_too_wide [{L},{R}] (no whole-window rerun)'; return None, info
        c = chunk(V9BIN, T, L, R, env, f'{scr}/final{j}.z', deltas=(0.02, 0.005, 0.0005))
        if c is None or c['Na'] != NL or c['Nb'] != NR: info['result'] = f'final_chunk_unproven_or_N_mismatch [{L},{R}]'; return None, info
        old = [x for x in new if x[0] >= L and x[1] <= R]
        if any(x[0] < L < x[1] or x[0] < R < x[1] for x in new): info['result'] = 'cut_inside_bracket'; return None, info
        rec += len(c['br']) - len(old); new = sorted([x for x in new if not (x[0] >= L and x[1] <= R)] + c['br'])
        chunks.append([L, R]); recs += [x for x in c['br'] if not any(o[0] < x[1] and x[0] < o[1] for o in old)]
    ok = len(new) == tur and all(new[j][1] <= new[j + 1][0] for j in range(len(new) - 1)) and rec == deficit
    info.update(chunk=[min(c[0] for c in chunks), max(c[1] for c in chunks)] if chunks else None, chunks=chunks, chunk_len=sum(c[1] - c[0] for c in chunks),
                recovered=rec, recovered_brackets=recs[:8], result='SPLICED' if ok else 'splice_check_failed')
    return (new if ok else None), info
# ---- t98c: Arb pair certificate (Lemma H), adopted from unified prototype (Nathan A. Meeks & Grok). Finder proposes, Arb disposes; fallback = multi-chunk close_pair.
def _dump(V9BIN, T0, delta, P, env):
    import numpy as np
    import threading; f = f'/tmp/gz_{os.getpid()}_{threading.get_ident()}.txt'
    subprocess.run([V9BIN, str(T0), str(delta), str(P), '80', '200', '60', '/dev/null', 'dbl'], env=dict(env, GZDUMP=f, AHCACHE=env.get('AHCACHE', '')), capture_output=True)
    d = np.loadtxt(f); os.remove(f); return d[:, 0], d[:, 1]
def pair_cert(V9BIN, T, A, B, Na, Nb, br, env, K=int(os.environ.get('PC_K', '16'))):
    """Heuristic proposes (uncertified fftzi grid dips, fine local grid), Arb disposes (indgates.sgn, independent of fftzi).
    Accepts only: Arb-certified sgn Z(a)=sgn Z(b)!=sgn Z(c), a<c<b inside one gap of br, and the spliced set closes
    Lemma CP on the whole window (len == Nb-Na, pairwise disjoint). Returns (new brackets, info) or (None, info)."""
    import numpy as np, sys
    sp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spotcheck')
    if sp not in sys.path: sys.path.insert(0, sp)
    import indgates
    from decimal import Decimal as D, getcontext; getcontext().prec = 60
    br = sorted(br); deficit = (Nb - Na) - len(br); info = dict(gate='pair_cert', deficit=deficit, tried=0, arb=0)
    if deficit != 2: info['result'] = 'not_eligible'; return None, info
    u, z = _dump(V9BIN, T, 0.02, 524288, env); a = np.abs(z)
    k = np.where((np.sign(z[:-2]) == np.sign(z[1:-1])) & (np.sign(z[1:-1]) == np.sign(z[2:])) & (a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]))[0] + 1
    k = k[(u[k] > A + 1) & (u[k] < B - 1)]
    from numpy.lib.stride_tricks import sliding_window_view as sw
    M = np.pad(sw(a, 31).max(1), 15, mode='edge'); k = k[np.argsort(a[k] / M[k])]
    mids = np.array([(x + y) / 2 for x, y in br]); lefts = np.array([x for x, y in br]); rights = np.array([y for x, y in br])
    for kk in k[:K]:
        c0 = u[kk]; j = np.searchsorted(rights, c0)   # gap (rights[j-1], lefts[j])
        if j < len(br) and lefts[j] < c0: continue    # dip inside an old bracket: skip (needs 3-zero logic)
        gl = rights[j - 1] if j else A; gr = lefts[j] if j < len(br) else B
        info['tried'] += 1
        T0 = int(T + c0) - 1; du, dz = _dump(V9BIN, T0, 1e-4, 32768, env); du = du + (T0 - T)
        s0 = np.sign(z[kk]); flip = np.where((np.sign(dz) == -s0) & (np.abs(du - c0) < 0.03))[0]
        if not len(flip): continue
        ci = flip[np.argmin(-np.abs(dz[flip]))]; lo_i, hi_i = flip.min() - 1, flip.max() + 1
        if lo_i < 0 or hi_i >= len(du): continue
        pa, pc, pb = du[lo_i], du[ci], du[hi_i]
        if not (gl < pa < pc < pb < gr): continue
        ex = lambda x: str(D(T) + D(repr(float(x))))
        fs = indgates.pc_fast_signs(T, [repr(float(x)) for x in (pa, pc, pb)]) if os.environ.get('PC_FAST', '0') == '1' else None   # t106c S1 (off by default)
        if fs is not None:
            sa, sc, sb = (s_ if s_ else indgates.sgn(ex(x)) for s_, x in zip(fs, (pa, pc, pb))); info['fast'] = info.get('fast', 0) + sum(1 for s_ in fs if s_); info['arb'] += sum(1 for s_ in fs if not s_)
        else: sa, sc, sb = (indgates.sgn(ex(x)) for x in (pa, pc, pb)); info['arb'] += 3
        if sa and sc and sb and sa == sb and sc != sa:
            new = sorted(br + [(float(repr(pa) and pa), float(pc)), (float(pc), float(pb))])
            ok = len(new) == Nb - Na and all(new[i][1] <= new[i + 1][0] for i in range(len(new) - 1))
            info.update(result='CERT' if ok else 'splice_check_failed', pair=[pa, pc, pb], chunk=[snap(pa - 20, br), snap(pb + 20, br)], recovered=2)
            if ok: return new, info
    info.setdefault('result', 'no_candidate'); return None, info
def nstep_gate(H, i, T, A, B, Nw, tN):
    flog(H, idx=i, window_lo=T + A, window_hi=T + B, gate='nstep', trigger=f'window crosses t_N={tN:.1f} (N {Nw}->{Nw+1})',
         chunk_interval='n/a (in-path)', recovered_brackets='n/a', method_ref='WP §G.2 (Lemma G2)', result='fired; window status in ledger.csv row idx')
# ---- Independent-check gates for OPEN-1 (G-cut, G-Turing, G-RS) — Nathan A. Meeks & Grok. [NUM] candidate, unrefereed; not RH. GATES.md, WP §G.5
def ind_named(H):
    fn = f'{H}/ind_blocks.txt'   # blocks Nathan names: one window idx per line
    return {int(x) for x in open(fn).read().split() if x.isdigit()} if os.path.exists(fn) else set()
SHARED_POOL = None; SCHED = None   # t102 (Nathan A. Meeks & Grok): set by driver; scheduling only, no math change
def ind_gate(H, i, T, lo, hi, zd, Na, Nb, trigger, block=None, blocks=None):
    if SCHED is not None and not getattr(SCHED.tl, 'held', False):
        with SCHED.slot(i): return ind_gate(H, i, T, lo, hi, zd, Na, Nb, trigger, block=block, blocks=blocks)
    if blocks and len(blocks) > 1:   # t98b: each repaired chunk independently (Rule G5 per chunk); all must resolve
        logs = []
        for bl in blocks:
            ok, lg = ind_gate(H, i, T, lo, hi, zd, Na, Nb, trigger, block=tuple(bl)); logs.append(lg)
            if not ok: return False, dict(steps=[s for l in logs for s in l['steps']], need=sum(l.get('need') or 0 for l in logs))
        return True, dict(steps=[s for l in logs for s in l['steps']], need=sum(l.get('need') or 0 for l in logs))
    if blocks: block = tuple(blocks[0])
    if block:   # t99 root-cause fix (Nathan A. Meeks & Grok): the chunk snap (pb+20) could overrun the window [lo,hi]; zeros past hi
        # belong to the next window and have no brackets here, so G-Turing's count exceeded the bracket count (38498: 162 vs 137).
        # Clip the referee block to the window. Sound: Lemma CP-S, cut positions are search only.
        block = (max(block[0], lo), min(block[1], hi))
    # block=(blo,bhi) absolute: the spliced close-pair chunk -> independent Lemma CP on that chunk only (Na=Nb=None: N from G-Turing)
    hint = None
    if block and Na:   # t98 H2 (Nathan A. Meeks & Grok): v9 window count -> chunk count hint (tested by Platt, never trusted)
        bb = [tuple(map(float, l.split()[:2])) for l in open(zd)]
        ha = Na + sum(1 for a, b in bb if T + b <= block[0]); hint = (ha, ha + sum(1 for a, b in bb if T + a >= block[0] and T + b <= block[1]))
    if block: lo, hi, Na, Nb = block[0], block[1], None, None
    sys_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spotcheck')
    import sys
    if sys_path not in sys.path: sys.path.insert(0, sys_path)
    import indgates
    res, log = indgates.resolve(T, lo, hi, indgates.read_br(zd), Na, Nb, trigger, procs=int(os.environ.get('IND_PROCS', '7')), hint=hint, shared_pool=SHARED_POOL)
    gl = ' -> '.join(f"{s['gate']}({s.get('sec', 'refused' if s.get('refused') else '?')}s)" for s in log['steps'])
    glog(H, gate='independent', idx=i, T=T, result=res, **log)
    flog(H, tiers=[t for st_ in log['steps'] for t in (st_.get('tiers') or [])], idx=i, window_lo=lo, window_hi=hi, gate='independent:' + '+'.join(s['gate'] for s in log['steps']), trigger=trigger,
         chunk_interval=(f'[{lo},{hi}]' if block else 'whole window'), recovered_brackets='n/a', method_ref='WP §G.5 (Rule G5); GATES.md',
         result=f"{res}{' CONDITIONAL(OPEN-P)' if any((s.get('gate') == 'G-Platt' and not s.get('refused')) or (s.get('escalation') or {}).get('platt') == 'RESOLVED CONDITIONAL(OPEN-P)' for s in log['steps']) else ''} need={log.get('need')} cost={log.get('total_sec')}s [{gl}]")
    return res == 'RESOLVED', log
