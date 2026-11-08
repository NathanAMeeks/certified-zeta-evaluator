#!/usr/bin/env python3
# [NUM] candidate, unrefereed. Contiguous certified zeta-zero verification above 3e12 with the v9 evaluator
# (v9 = v7 speed + v8 gated self-correcting loop, directed by Nathan A. Meeks). Credit: Nathan A. Meeks & Grok.
# Extends a verified height only; does NOT prove RH.
# t102: 1 orchestrator core (main/commit/xz) + 7 v9 workers; gate map calls go to one shared pinned pool, head window first (Nathan A. Meeks & Grok).
# t98: parallel independent windows on WORKER_CPUS; ordered atomic commit of the contiguous PROVEN prefix only (Nathan A. Meeks & Grok).
import json, math, os, subprocess, sys, time, datetime
os.environ['HYB_VC'] = '0'   # t98c: referee fully independent (no G-Turing/G-RS value sharing); TMESH stays (search-only)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import gates
V9BIN = os.environ.get('V9BIN', '/workspace/rh_extend/v9p/fftzi_v9')
H = os.environ.get('RHX_HOME', os.path.dirname(os.path.abspath(__file__)))
V9 = os.environ.get('V9RUN', '/workspace/rh_extend/v9p/v9/v9run.sh')
CK, LED, FAIL = os.environ.get('CKPT', f'{H}/checkpoint.json'), f'{H}/ledger.csv', os.environ.get('FAILMD', f'{H}/FAILED.md')
LED = os.environ.get('LEDGER', LED)
A = 80                       # fftzi arg a: window starts at T+A
P, B = 65536, 10400           # fftzi P and b (SPEED2 v9 benchmark config, 60/60 PROVEN): window [T+80, T+10400], width 10320
OVL = 40                     # overlap between consecutive windows (t-units)
STEP = (B - A) - OVL         # 10280
NWIN = int(os.environ.get('BATCH_WINDOWS', '1300'))
T0 = 3_000_000_000_000 - A   # first window starts exactly at t = 3e12
CPU = os.environ.get('CPU', '7'); XZCPUS = os.environ.get('XZCPUS', '0-5'); XZPAR = int(os.environ.get('XZPAR', '3'))
COLS = 'idx,T_arg,lo,hi,zeros_turing,cert_sign_changes,Na,Nb,status,delta,P_used,attempts,sec_v9,wall_s,overlap_N,overlap_brackets,timestamp,notes'
def now(): return datetime.datetime.now().astimezone().isoformat(timespec='seconds')
def load():
    if os.path.exists(CK): return json.load(open(CK))
    return dict(next_T=T0, next_idx=0, N_start=None, N_last_hi=None, covered_lo=3e12, covered_hi=3e12, windows=0, zeros=0,
                failures=0, last_run=None, last_result=None, params=dict(A=A, P=P, B=B, OVL=OVL, STEP=STEP))
def save(c):
    tmp = CK + '.tmp'
    with open(tmp, 'w') as f: json.dump(c, f, indent=1); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, CK)
def ledger(row):
    new = not os.path.exists(LED)
    with open(LED, 'a') as f:
        if new: f.write(COLS + '\n')
        f.write(','.join(str(x) for x in row) + '\n'); f.flush(); os.fsync(f.fileno())
def window(T, i, env):
    """[t98] one window: all work independent of the previous window (Nathan A. Meeks & Grok)."""
    lo, hi = T + A, T + B; info = None
    braw = f'{H}/brackets/w{i:07d}_T{T}.txt'
    with open(braw, 'w') as f: f.write(f'# [NUM] candidate, unrefereed. Certified sign-change brackets of Z(t), t = {T} + u; window [{lo}, {hi}]; v9p (v9 + C0 removable-singularity fix); Nathan A. Meeks & Grok\n')
    zd = braw + '.z'; s = time.time(); notes = []
    if gates.c0_gate(H, i, T, A, B): notes.append('HEUR-gate:C0_cut(Dedekind cut, Nathan A. Meeks)')
    # cheap N-step pre-check: does [T, T+B+Lw] cross t_N = 2pi(N+1)^2?  v9p handles it in-path (NSTEP_FIX.md); log it and
    # fail fast (no MAXR rung) if such a window still fails, since that would mean the N-step path itself is broken.
    Nw = math.isqrt(int(T / (2 * math.pi))); tN = 2 * math.pi * (Nw + 1) ** 2; nstep = T <= tN <= T + B + 60 + 1
    if nstep: open(f'{H}/logs/nstep.log', 'a').write(f'{now()} idx={i} T={T} crosses t_N={tN:.1f} (N {Nw}->{Nw+1})\n')
    if nstep: gates.nstep_gate(H, i, T, A, B, Nw, tN)
    fast = os.environ.get('GATE_FASTCP', '1') == '1' and not nstep
    p = subprocess.run([V9, str(T), str(P), str(B)], env=dict(env, ZDUMP=zd, **({'V9_FIRST_ONLY': '1'} if fast else {})), capture_output=True, text=True)
    att = 1 + p.stderr.count('escalate'); err = p.stderr
    if fast and p.returncode != 0:   # t92 fast path (Nathan A. Meeks & Grok): attempt 1 short by 2/4 with all found brackets certified -> close-pair gate directly
        try: d1 = json.loads(p.stdout.strip().splitlines()[-1])
        except Exception: d1 = {}
        nb1 = sum(1 for _ in open(zd)) if os.path.exists(zd) else -1
        dfc = int(d1.get('Nb_upper', 0) or 0) - int(d1.get('Na_lower', 0) or 0) - int(d1.get('found', -1))
        if str(d1.get('status', '')).startswith('UNPROVEN') and int(d1.get('Na_lower', 0) or 0) > 0 and nb1 == d1.get('found') and dfc in (2, 4):
            br1 = [tuple(map(float, l.split()[:2])) for l in open(zd)]
            new1, info1 = (gates.pair_cert(V9BIN, T, A, B, int(d1['Na_lower']), int(d1['Nb_upper']), br1, env) if dfc == 2 and os.environ.get('PAIR_CERT', '1') == '1' else (None, {}))  # t98c Lemma H first
            if new1 is not None: gates.glog(H, idx=i, T=T, **info1); info1.update(gate='close_pair', idx=i, T=T, deficit=2, method='Lemma H pair certificate (Arb signs)', result='SPLICED', chunks=[info1['chunk']])
            else: new1, info1 = gates.close_pair(H, V9BIN, i, T, A, B, int(d1['Na_lower']), int(d1['Nb_upper']), br1, env, f'{H}/scratch/gate_w{i:07d}')
            gates.glog(H, fastpath=True, **info1)
            if new1 is not None:
                with open(zd, 'w') as g: g.writelines(f'{a!r} {b!r}\n' for a, b in new1)
                d1.update(status='PROVEN', found=len(new1)); info = info1
                p = subprocess.CompletedProcess(p.args, 0, stdout=json.dumps(d1) + '\n', stderr='')
                notes.append(f"HEUR-gate:close_pair(centre-out, Nathan A. Meeks) chunk u[{info['chunk'][0]},{info['chunk'][1]}] +{info['recovered']} fastpath")
        if p.returncode != 0:   # not resolved by the fast path: full v9 ladder exactly as before
            err += '--- fastpath miss: full ladder ---\n'
            p = subprocess.run([V9, str(T), str(P), str(B)], env=dict(env, ZDUMP=zd), capture_output=True, text=True)
            att += 1 + p.stderr.count('escalate'); err += p.stderr
    if p.returncode != 0 and not nstep:   # extra driver rung (search-only): rerun the whole v9 ladder with refine round cap MAXR=200 instead of 20
        p = subprocess.run([V9, str(T), str(P), str(B)], env=dict(env, ZDUMP=zd, MAXR='200'), capture_output=True, text=True)
        att += 1 + p.stderr.count('escalate'); err += '--- MAXR=200 rung ---\n' + p.stderr
    w = time.time() - s
    if err.strip():
        open(f'{H}/logs/w{i:07d}.err', 'w').write(err)
    d = {}
    try: d = json.loads(p.stdout.strip().splitlines()[-1])
    except Exception: pass
    st = d.get('status', 'NO_OUTPUT'); found = d.get('found', -1)
    Na, Nb = int(d.get('Na_lower', 0) or 0), int(d.get('Nb_upper', 0) or 0); tur = Nb - Na
    nbr = sum(1 for _ in open(zd)) if os.path.exists(zd) else -1
    if (p.returncode != 0 and st.startswith('UNPROVEN') and Na > 0 and nbr == found and tur - found in (2, 4)
            and os.environ.get('GATE_CLOSEPAIR', '1') == '1'):
        br0 = [tuple(map(float, l.split()[:2])) for l in open(zd)]
        new, info = gates.close_pair(H, V9BIN, i, T, A, B, Na, Nb, br0, env, f'{H}/scratch/gate_w{i:07d}')
        gates.glog(H, **info)
        if new is not None:
            with open(zd, 'w') as g: g.writelines(f'{a:.9f} {b:.9f}\n' for a, b in new)
            found = nbr = len(new); st = 'PROVEN'; p = subprocess.CompletedProcess(p.args, 0)
            notes.append(f"HEUR-gate:close_pair(centre-out, Nathan A. Meeks) chunk u[{info['chunk'][0]},{info['chunk'][1]}] +{info['recovered']}")
        else:
            notes.append(f"HEUR-gate:close_pair FAILED ({info.get('result')})")
            if os.environ.get('T104_DIP', '1') == '1':   # t104 (Nathan A. Meeks & Grok, t132s656): window end count not pinned (e.g. close pair just
                # beyond hi): whole-window referee with N from G-Turing (+ t104 dip refinement), never from v9; all brackets G-RS certified
                tok, tlog = gates.ind_gate(H, i, T, lo, hi, zd, None, None, 'turing_end_unpinned')
                tN = [x for x in tlog['steps'] if x.get('gate') == 'G-Turing']
                if tok and tN and tN[-1]['Na'][0] is not None and tN[-1]['Nb'][0] is not None and tN[-1]['Nb'][0] - tN[-1]['Na'][0] == nbr == found:
                    Na, Nb = tN[-1]['Na'][0], tN[-1]['Nb'][0]; tur = Nb - Na; st = 'PROVEN'; p = subprocess.CompletedProcess(p.args, 0)
                    d.update(Na_lower=Na, Nb_upper=Nb, status='PROVEN')
                    notes.append(f'IND-gate:t104 turing_end_unpinned RESOLVED N=[{Na},{Nb}] (G-Turing, not v9)')
                else: notes.append('IND-gate:t104 turing_end_unpinned FAILED')
    pre = None   # t98: run the independent gate in the worker (it needs only this window); commit re-runs it only if overlap fails
    okp = p.returncode == 0 and st == 'PROVEN' and found == tur and nbr == found and Na > 0
    tp = []; bk = None
    if any(x.startswith('HEUR-gate:close_pair(') for x in notes): tp.append('close_pair_chunk'); bk = [(T + a_, T + b_) for a_, b_ in info.get('chunks') or [info['chunk']]]
    elif any('close_pair' in x for x in notes): tp.append('gate_ladder')
    if i in gates.ind_named(H): tp.append('named_by_Nathan')
    if tp and okp and os.environ.get('GATE_IND', '1') == '1' and os.path.exists(zd):
        pre = ('+'.join(tp),) + tuple(gates.ind_gate(H, i, T, lo, hi, zd, Na, Nb, '+'.join(tp), blocks=bk if tp == ['close_pair_chunk'] else None))
    return dict(T=T, i=i, lo=lo, hi=hi, p=p, d=d, st=st, found=found, Na=Na, Nb=Nb, tur=tur, nbr=nbr, notes=notes, w=time.time() - s, att=att, braw=braw, zd=zd, info=info, pre=pre)
def main():
    c = load(); c['params'] = dict(A=A, P=P, B=B, OVL=OVL, STEP=STEP)
    if os.path.exists(FAIL):
        line = f"[NUM] batch HALTED: open failure recorded in {FAIL} (window T={c['next_T']}); covered to {c['covered_hi']:.9e}; fix/resolve then move FAILED.md aside"
        c['last_run'] = now(); c['last_result'] = line; save(c); print(line); return 3
    env = dict(os.environ, AHCACHE=f'{H}/ahcache.bin', CPU=CPU)
    xz = []; t_start = time.time(); h0 = c['covered_hi']; z0 = c['zeros']; done = 0; rc = 0; failmsg = ''
    jl = open(f'{H}/logs/windows.jsonl', 'a')
    from concurrent.futures import ThreadPoolExecutor
    import queue
    cpus = [x for x in os.environ.get('WORKER_CPUS', '1,2,3,4,5,6').split(',') if x]; cq = queue.Queue(); [cq.put(x) for x in cpus]
    def job(T, i):
        cpu = cq.get()
        try: return window(T, i, dict(env, CPU=cpu))
        finally: cq.put(cpu)
    # t98: pre-extend the shared AH cache once (sequential) for the batch's top window, so parallel workers only read it
    # (ah_load is checksummed: a torn read can only force a recompute, never a wrong value)
    Te = c['next_T'] + (NWIN - 1) * STEP
    subprocess.run(['taskset', '-c', CPU, V9BIN, str(Te), '0.24', str(P), str(A), str(B), '60', '/dev/null', 'dbl'], env=dict(env, ZDUMP='/dev/null'), capture_output=True)
    # t102 (Nathan A. Meeks & Grok): orchestrator core + shared pinned gate pool + head-first gate scheduler (scheduling only)
    if os.environ.get('T102', '1') == '1':
        import multiprocessing as mp, threading, contextlib
        gcpus = sorted({int(x) for x in os.environ.get('GATE_CPUS', ','.join(cpus)).split(',') if x})
        def _pin(cs):
            os.sched_setaffinity(0, cs)
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(gates.__file__)), 'spotcheck')); import indgates  # before fork: children resolve pickled refs
        gp = mp.get_context('fork').Pool(int(os.environ.get('GATE_PROCS', str(len(gcpus)))), initializer=_pin, initargs=(gcpus,))
        class Sched:
            # at most K gated windows use the shared pool at once; the lowest waiting idx is admitted first;
            # the window that blocks the in-order checkpoint (head) is always admitted at once
            def __init__(s, K): s.K = K; s.cv = threading.Condition(); s.run = set(); s.wait = set(); s.tl = threading.local()
            @contextlib.contextmanager
            def slot(s, i):
                with s.cv:
                    s.wait.add(i)
                    while not (i == st8['head'] or (len(s.run) < s.K and i == min(s.wait))): s.cv.wait(0.5)
                    s.wait.discard(i); s.run.add(i)
                s.tl.held = True
                try: yield
                finally:
                    s.tl.held = False
                    with s.cv: s.run.discard(i); s.cv.notify_all()
        gates.SHARED_POOL = gp; gates.SCHED = Sched(int(os.environ.get('GATE_K', '2')))
    ex = ThreadPoolExecutor(len(cpus)); pend = []; st8 = dict(T=c['next_T'], i=c['next_idx'], n=0, head=c['next_idx'])
    def fill():
        while len(pend) < int(os.environ.get('LOOKAHEAD', '48')) * len(cpus) and st8['n'] < NWIN:   # consecutive windows in idx order
            pend.append(ex.submit(job, st8['T'], st8['i'])); st8['T'] += STEP; st8['i'] += 1; st8['n'] += 1
    fill()   # t101 (Nathan A. Meeks & Grok): lookahead 4->48 per CPU; a ~40 s gated window no longer stalls all 8 workers (ordered commit unchanged)
    while pend:
        R = pend.pop(0).result(); pre = R['pre']
        T, i, lo, hi, p, d, st, found, Na, Nb, tur, nbr, notes, w, att, braw, zd, info = (R[k] for k in ('T','i','lo','hi','p','d','st','found','Na','Nb','tur','nbr','notes','w','att','braw','zd','info'))
        assert T == c['next_T'] and i == c['next_idx']   # ordered commit: contiguous prefix only
        ovN = ovB = ''
        jl.write(json.dumps(dict(idx=i, rc=p.returncode, **d)) + '\n'); jl.flush()
        ok = p.returncode == 0 and st == 'PROVEN' and found == tur and nbr == found and Na > 0
        if ok and c['N_last_hi'] is not None:      # contiguity with previous window: N(prev hi) >= N(this lo), overlap counts consistent
            ovN = c['N_last_hi'] - Na
            ovB = sum(1 for l in open(zd) if float(l.split()[1]) <= A + OVL)
            if not (0 <= ovN and ovB <= ovN <= ovB + 1): ok = False; st = st + '_OVERLAP_MISMATCH'
        trig = []                                   # independent-check gates (OPEN-1), Nathan A. Meeks & Grok
        blk = None   # MAXR=200-only windows are v9-standard PROVEN (no uncertain sign) -> not flagged (GATES.md rule G5)
        if any(x.startswith('HEUR-gate:close_pair(') for x in notes): trig.append('close_pair_chunk'); blk = [(T + a_, T + b_) for a_, b_ in info.get('chunks') or [info['chunk']]]
        elif any('close_pair' in x for x in notes): trig.append('gate_ladder')
        if not ok and int(d.get('refine_disagreements', 0) or 0) > 0: trig.append('v9_sign_uncertain')   # error-only (Nathan, t90): a v9-PROVEN window gets no gates
        if 'OVERLAP_MISMATCH' in st: trig.append('overlap_mismatch')
        if i in gates.ind_named(H): trig.append('named_by_Nathan')
        if trig and (ok or 'overlap_mismatch' in trig) and os.environ.get('GATE_IND', '1') == '1' and os.path.exists(zd):
            if pre is not None and pre[0] == '+'.join(trig): iok, ilog = pre[1], pre[2]   # identical call already made in the worker
            else: iok, ilog = gates.ind_gate(H, i, T, lo, hi, zd, Na, Nb, '+'.join(trig), blocks=blk if trig == ['close_pair_chunk'] else None)
            notes.append(f"IND-gate:{'RESOLVED' if iok else 'FAILED'}({'+'.join(s['gate'] for s in ilog['steps'])})")
            if not iok: ok = False; st = st + '_IND_GATE_FAILED'
        row = [i, T, lo, hi, tur, found, Na, Nb, st, d.get('v9_delta', ''), d.get('v9_P', ''), att, d.get('sec_all_attempts', ''), round(w, 3), ovN, ovB, now(), ';'.join(notes)]
        ledger(row)
        if not ok:
            failmsg = f'window idx={i} T_arg={T} [{lo},{hi}] status={st} found={found} turing={tur} brackets={nbr} rc={p.returncode} attempts={att}'
            with open(FAIL, 'a') as f:
                f.write(f'# OPEN FAILURE {now()}\n{failmsg}\nAfter full v9 retry ladder + MAXR=200 rung + close-pair gate (if eligible) + independent gates (if triggered); notes: {';'.join(notes)}. Checkpoint left at this window (not skipped). stderr: logs/w{i:07d}.err, json: logs/windows.jsonl idx={i}\n\n')
            c['failures'] += 1; rc = 2
            for f_ in pend: f_.cancel()
            for f_ in pend:   # in-flight windows finish, then are discarded (never committed past the failure)
                if f_.cancelled(): continue
                r_ = f_.result()
                for fn_ in (r_['braw'], r_['zd']):
                    if os.path.exists(fn_): os.remove(fn_)
            pend = []; break
        with open(braw, 'a') as f, open(zd) as g: f.write(g.read())
        os.remove(zd)
        while len([q for q in xz if q.poll() is None]) >= XZPAR: time.sleep(0.05)
        xz = [q for q in xz if q.poll() is None] + [subprocess.Popen(['taskset', '-c', XZCPUS, 'nice', 'xz', '-T1', '-6', '-f', braw])]
        if c['N_start'] is None: c['N_start'] = Na
        c['N_last_hi'] = Nb; c['zeros'] = Nb - c['N_start']; c['covered_hi'] = float(hi); c['windows'] += 1
        c['next_T'] = T + STEP; c['next_idx'] = i + 1; save(c); done += 1; st8['head'] = i + 1
        fill()
    ex.shutdown(wait=True)
    if gates.SHARED_POOL is not None: gates.SHARED_POOL.close(); gates.SHARED_POOL.join()
    for q in xz: q.wait()
    el = time.time() - t_start
    line = (f"[NUM] batch {'OK' if rc == 0 else 'FAILED'}: T {h0:.9e}->{c['covered_hi']:.9e}, +{c['zeros'] - z0} zeros, {done} windows PROVEN, "
            f"total covered to {c['covered_hi']:.9e} ({c['zeros']} zeros since 3e12); failures {1 if rc else 0}; {el/60:.1f} min"
            + (f' | FAIL {failmsg}' if rc else '') + ' (candidate, unrefereed; not RH)')
    c['last_run'] = now(); c['last_result'] = line; save(c)
    open(f'{H}/logs/runs.log', 'a').write(c['last_run'] + ' ' + line + '\n')
    print(line); return rc
if __name__ == '__main__':
    try: sys.exit(main())
    except SystemExit: raise
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'[NUM] batch ERROR: driver crashed ({type(e).__name__}: {e}); checkpoint unchanged past last PROVEN window; see logs/driver.err'); sys.exit(5)
