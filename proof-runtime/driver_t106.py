#!/usr/bin/env python3
# [NUM] candidate, unrefereed. Contiguous certified zeta-zero verification above 3e12 with the v9 evaluator
# (v9 = v7 speed + v8 gated self-correcting loop, directed by Nathan A. Meeks). Credit: Nathan A. Meeks & Grok.
# Extends a verified height only; does NOT prove RH.
# t106 [NUM] candidate, unrefereed (Nathan A. Meeks & Grok): S3=1 -> 4-window parallel referee lanes (Creative S3, refine/hybrid/RESULTS.md);
#   explicit contiguous-prefix commit guard + commit log; G-RS batch evaluator selectable via GRS_BIN in spotcheck/indgates.py. Copy only, not deployed.
# t103 (Nathan A. Meeks & Grok): checkpoint every CK_EVERY windows with ledger_rows + truncate-on-restart; xz batched per CK_EVERY windows, raw .txt swept on restart.
# t102: 1 orchestrator core (main/commit/xz) + 7 v9 workers; gate map calls go to one shared pinned pool, head window first (Nathan A. Meeks & Grok).
# t98: parallel independent windows on WORKER_CPUS; ordered atomic commit of the contiguous PROVEN prefix only (Nathan A. Meeks & Grok).
import shutil, json, math, os, subprocess, sys, time, datetime
os.environ['HYB_VC'] = '0'   # t98c: referee fully independent (no G-Turing/G-RS value sharing); TMESH stays (search-only)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import gates
XZL = '-' + os.environ.get('XZLEVEL', '1')   # t106c (Nathan A. Meeks & Grok; Audit v9_vs_t106b #1): xz -6 -> -1; brackets_sha is on the raw bytes, unchanged
V9BIN = os.environ.get('V9BIN', '/workspace/rh_extend/v9p/fftzi_v9')
H = os.environ.get('RHX_HOME', os.path.dirname(os.path.abspath(__file__)))
V9 = os.environ.get('V9RUN', '/workspace/rh_extend/v9p/v9/v9run.sh')
CK, LED, FAIL = os.environ.get('CKPT', f'{H}/checkpoint.json'), f'{H}/ledger.csv', os.environ.get('FAILMD', f'{H}/FAILED.md')
LED = os.environ.get('LEDGER', LED)
A = 80                       # fftzi arg a: window starts at T+A
P, B = 65536, 10400           # fftzi P and b (SPEED2 v9 benchmark config, 60/60 PROVEN): window [T+80, T+10400], width 10320
B = int(os.environ.get('WIN_LEN', str(B)))   # t106c (Nathan A. Meeks & Grok; scratch testing only): window length switch, default 10400 = live
if B != 10400:   # fftzi_v9 needs P a power of 2 (P=98304 -> certified=0); smallest 2^k with 0.24*P >= 1.2*(B+140) (v9 ladder 'need'); gives 65536 at 10400
    P = 1
    while 0.24 * P < (12 * (B + 140) + 9) // 10: P *= 2
OVL = 40                     # overlap between consecutive windows (t-units)
STEP = (B - A) - OVL         # 10280
NWIN = int(os.environ.get('BATCH_WINDOWS', '1300'))
T0 = 3_000_000_000_000 - A   # first window starts exactly at t = 3e12
CPU = os.environ.get('CPU', '7'); XZCPUS = os.environ.get('XZCPUS', '0-5'); XZPAR = int(os.environ.get('XZPAR', '3'))
XZSTAG = os.environ.get('XZSTAG', '0') == '1'   # t106: t105's staggered low-priority xz (ported, Nathan A. Meeks & Grok)
XZPAR5 = int(os.environ.get('XZPAR5', '2')); XZCHUNK = int(os.environ.get('XZCHUNK', '8'))
XZCMD = ['taskset', '-c', XZCPUS, 'nice', '-n', '19'] + (['ionice', '-c3'] if shutil.which('ionice') else []) + ['xz', '-T1', XZL, '-f']
def xz_tick(xz, q, drain=False):
    # launch at most one small xz per call; never block the commit loop unless draining
    xz = [p for p in xz if p.poll() is None or (p.wait() and open(f'{H}/logs/xz.log', 'a').write(f'{now()} xz rc={p.returncode} (logged, batch not failed)\n') and False)]
    while q and len(xz) < XZPAR5:
        xz.append(subprocess.Popen(XZCMD + q[:XZCHUNK])); del q[:XZCHUNK]
        if not drain: break
    return xz
COLS = 'idx,T_arg,lo,hi,zeros_turing,cert_sign_changes,Na,Nb,status,delta,P_used,attempts,sec_v9,wall_s,overlap_N,overlap_brackets,timestamp,brackets_sha,notes'   # t106 (e): +brackets_sha before notes (old rows: 18 cols, no sha)
def now(): return datetime.datetime.now().astimezone().isoformat(timespec='seconds')
def load():
    if os.path.exists(CK): return json.load(open(CK))
    return dict(next_T=T0, next_idx=0, N_start=None, N_last_hi=None, covered_lo=3e12, covered_hi=3e12, windows=0, zeros=0,
                failures=0, last_run=None, last_result=None, ledger_rows=0, clog_bytes=0, params=dict(A=A, P=P, B=B, OVL=OVL, STEP=STEP))
CLOG = f'{H}/logs/commit_t106.jsonl'
def save(c):
    # t106 B2 (re-audit): clog_bytes is NOT stamped here; the commit loop sets it to the offset just after the last COMMITTED entry
    tmp = CK + '.tmp'
    with open(tmp, 'w') as f: json.dump(c, f, indent=1); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, CK)
def ledger(row):
    new = not os.path.exists(LED)
    with open(LED, 'a') as f:
        if new: f.write(COLS + '\n')
        f.write(','.join(str(x) for x in row) + '\n'); f.flush(); os.fsync(f.fileno())
def nrows(fn):
    if not os.path.exists(fn): return 0
    with open(fn, 'rb') as f: return sum(1 for _ in f)
def _legacy_offsets(c):
    # t106 B2: checkpoint without ledger_rows/clog_bytes (legacy ck or older t106): the trusted prefix is the rows/entries with idx < next_idx,
    # in file order, up to the first one at or past next_idx (everything after was written after the checkpoint and will be redone).
    lr = 0
    if os.path.exists(LED):
        with open(LED, 'rb') as f:
            for l in f:
                x = l.split(b',', 1)[0]
                if x.isdigit() and int(x) >= c['next_idx']: break
                lr += 1
    cb = 0
    if os.path.exists(CLOG):
        with open(CLOG, 'rb') as f:
            for l in f:
                try: ok = l.endswith(b'\n') and json.loads(l)['commit'] < c['next_idx']
                except Exception: ok = False
                if not ok: break
                cb += len(l)
    return lr, cb
def t103_restart(c):
    if 'ledger_rows' not in c or 'clog_bytes' not in c:
        lr, cb = _legacy_offsets(c); c['ledger_rows'] = c.get('ledger_rows') if isinstance(c.get('ledger_rows'), int) else lr; c['clog_bytes'] = cb   # stored clog_bytes of a legacy ck is untrusted (old save() stamped file size)
        open(f'{H}/logs/t103.log', 'a').write(f'{now()} t106 legacy checkpoint: trusted ledger_rows={c["ledger_rows"]} clog_bytes={c["clog_bytes"]}\n')
    # t106 B2: drop commit-log entries written after the last checkpoint (they will be redone)
    c['clog_bytes'] = min(c['clog_bytes'], _legacy_offsets(c)[1])   # never trust past the prefix of entries with commit < next_idx
    if os.path.exists(CLOG) and os.path.getsize(CLOG) > c['clog_bytes']:
        with open(CLOG, 'r+b') as f: f.truncate(c['clog_bytes']); f.flush(); os.fsync(f.fileno())
        open(f'{H}/logs/t103.log', 'a').write(f'{now()} t106 truncated commit log to {c["clog_bytes"]} bytes\n')
    # t103: ledger rows past the checkpoint's ledger_rows belong to windows the checkpoint has not committed (crash between saves);
    # they will be redone, so cut them (kept in ledger.csv.t103cut). Sweep committed raw .txt -> xz.
    n = nrows(LED)
    if n > c['ledger_rows']:
        with open(LED, 'rb') as f: keep = b''.join(l for _, l in zip(range(c['ledger_rows']), f))
        open(LED + '.t103cut', 'ab').write(open(LED, 'rb').read()[len(keep):])
        if c['ledger_rows'] <= 1: os.remove(LED); c['ledger_rows'] = 0   # header-only/empty -> fresh file (header rewritten by ledger())
        else:
            with open(LED, 'r+b') as f: f.truncate(len(keep)); f.flush(); os.fsync(f.fileno())
        open(f'{H}/logs/t103.log', 'a').write(f'{now()} truncated ledger {n}->{c["ledger_rows"]} rows (cut rows kept in ledger.csv.t103cut)\n')
    import glob
    for t_ in glob.glob(f'{H}/brackets/w*_T*.txt.tmp*'): os.remove(t_)   # t106 (f): stray temp files of an interrupted atomic write
    raw = [f for f in glob.glob(f'{H}/brackets/w*_T*.txt') if int(os.path.basename(f)[1:8]) < c['next_idx']]
    for k in range(0, len(raw), 200):
        r = subprocess.run(['taskset', '-c', XZCPUS, 'nice', 'xz', '-T1', XZL, '-f'] + raw[k:k + 200])
        if r.returncode: open(f'{H}/logs/xz.log', 'a').write(f'{now()} sweep xz rc={r.returncode} files={len(raw[k:k+200])}\n')
import hashlib, re
def atomic_write(fn, data):
    # t106 (f): temp file in the same dir, fsync, rename, fsync the dir -> a reader/crash sees the old file or the complete new one, never a torn file
    d = os.path.dirname(fn); tmp = f'{fn}.tmp{os.getpid()}'
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        mv = memoryview(data)
        while mv: n = os.write(fd, mv); mv = mv[n:]
        os.fsync(fd)
    finally: os.close(fd)
    os.replace(tmp, fn); dfd = os.open(d, os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(dfd)
    finally: os.close(dfd)
def brackets_sha_of(fn):
    # sha256 of the (uncompressed) bracket file; .xz is decompressed. None if absent.
    import lzma
    for p, op in ((fn, open), (fn + '.xz', lzma.open)):
        if os.path.exists(p):
            h = hashlib.sha256()
            with op(p, 'rb') as f:
                for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
            return h.hexdigest()
    return None
def _is_sha(x): return len(x) == 64 and all(ch in '0123456789abcdef' for ch in x)
def verify_brackets(c, last_k=int(os.environ.get('VERIFY_LAST', '64')), full=False, sample=int(os.environ.get('VERIFY_SAMPLE', '64'))):
    # t106b (Nathan A. Meeks & Grok; Audit speed REPORT): incremental (e) check. Checked: every row at/after c['verified_rows'] (rows committed since
    # the last clean verification), every still-uncompressed .txt, the last last_k rows, and a random sample of older rows. full=True: every row.
    # Returns (mismatches, rows_checked, total_sha_rows). Old rows without a sha are accepted (unchanged).
    if not os.path.exists(LED): return [], 0, 0
    rows = []; vr = c.get('verified_rows', 0) if not full else 0; n = 0
    with open(LED) as f:
        for l in f:
            n += 1; x = l.rstrip('\n').split(',')
            if x[0].isdigit() and len(x) > 18 and x[8] == 'PROVEN' and _is_sha(x[17]) and int(x[0]) < c['next_idx']: rows.append((n, int(x[0]), int(x[1]), x[17]))
    old = [k for k, r in enumerate(rows) if r[0] <= vr and k < len(rows) - last_k]
    samp = set()
    if not full and old:   # t106b (Audit nit): rotating sample, cursor persisted in the checkpoint (every old row checked once per len(old)/sample restarts)
        cur = c.get('verify_cursor', 0) % len(old); samp = {old[(cur + q) % len(old)] for q in range(min(sample, len(old)))}; c['verify_cursor'] = (cur + sample) % len(old)
    bad = []; seen = set(); m = 0
    for k, (ln, i, T, want) in enumerate(rows):
        fn = f'{H}/brackets/w{i:07d}_T{T}.txt'
        if not (full or ln > vr or k >= len(rows) - last_k or k in samp or os.path.exists(fn)) or (i, T) in seen: continue
        seen.add((i, T)); got = brackets_sha_of(fn); m += 1
        if got != want: bad.append((i, want, got))
    return bad, m, len(rows)
class _R:
    def __init__(s, rc, out, err): s.returncode, s.stdout, s.stderr, s.args = rc, out, err, ['v9direct']
def _g6(x):   # awk 'print a+b' (OFMT %.6g)
    return '%.6g' % x
def v9direct(T, P, B, env):
    """t106c (Nathan A. Meeks & Grok; Audit v9_vs_t106b #4): v9run.sh ported to Python, same ladder, same fftzi_v9 argv, same stdout/stderr/rc
    contract; no bash/grep/awk/sed spawns. Ladder is search-only; only a PROVEN fftzi_v9 result is accepted (unchanged)."""
    T, P, B = int(T), int(P), int(B); need = (12 * (B + 140) + 9) // 10; spent = 0.0; err = []; o = ''
    e = dict(env); e.setdefault('AHCACHE', os.path.join(os.path.dirname(V9BIN), 'ahcache.bin')); cpu = e.get('CPU', '7')
    def run(dl, Pn):
        nonlocal spent, o
        r = subprocess.run(['taskset', '-c', cpu, V9BIN, str(T), dl, str(Pn), '80', str(B), '60', '/dev/null', 'dbl'], env=e, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        o = r.stdout.rstrip('\n'); m = re.search(r'"sec_total":([0-9.]*)', o); spent = float(_g6(spent + float(m.group(1) if m and m.group(1) else 0)))
        if '"status":"PROVEN"' in o:
            return _R(0, '\n'.join(re.sub(r'}$', f',"v9_delta":{dl},"v9_P":{Pn},"sec_all_attempts":{_g6(spent)}}}', l) for l in o.split('\n')) + '\n', ''.join(err))
        err.append(f'v9run: T={T} delta={dl} P={Pn} not proven, escalate\n'); return None
    def grow(Pn, dl):
        while Pn * float(dl) < need: Pn *= 2
        return Pn
    r = run('0.24', P)
    if r: return r
    if e.get('V9_FIRST_ONLY', '0') == '1': return _R(1, o + '\n', ''.join(err))
    Pg = grow(P, '0.24')
    if Pg != P:
        r = run('0.24', Pg)
        if r: return r
    for dl in ('0.2', '0.16', '0.12'):
        r = run(dl, grow(P, dl))
        if r: return r
    return _R(1, o + '\n', ''.join(err))
def V9CALL(T, P, B, env):
    if os.environ.get('V9DIRECT', '1') == '1': return v9direct(T, P, B, env)
    return subprocess.run([V9, str(T), str(P), str(B)], env=env, capture_output=True, text=True)
def window(T, i, env):
    """[t98] one window: all work independent of the previous window (Nathan A. Meeks & Grok)."""
    lo, hi = T + A, T + B; info = None
    braw = f'{H}/brackets/w{i:07d}_T{T}.txt'
    hdr_line = (f'# [NUM] candidate, unrefereed. Certified sign-change brackets of Z(t), t = {T} + u; window [{lo}, {hi}]; v9p (v9 + C0 removable-singularity fix); Nathan A. Meeks & Grok\n')   # t106 (f): header kept in memory; file written atomically at commit
    zd = braw + '.z'; s = time.time(); notes = []
    if gates.c0_gate(H, i, T, A, B): notes.append('HEUR-gate:C0_cut(Dedekind cut, Nathan A. Meeks)')
    # cheap N-step pre-check: does [T, T+B+Lw] cross t_N = 2pi(N+1)^2?  v9p handles it in-path (NSTEP_FIX.md); log it and
    # fail fast (no MAXR rung) if such a window still fails, since that would mean the N-step path itself is broken.
    Nw = math.isqrt(int(T / (2 * math.pi))); tN = 2 * math.pi * (Nw + 1) ** 2; nstep = T <= tN <= T + B + 60 + 1
    if nstep: open(f'{H}/logs/nstep.log', 'a').write(f'{now()} idx={i} T={T} crosses t_N={tN:.1f} (N {Nw}->{Nw+1})\n')
    if nstep: gates.nstep_gate(H, i, T, A, B, Nw, tN)
    fast = os.environ.get('GATE_FASTCP', '1') == '1' and not nstep
    p = V9CALL(T, P, B, dict(env, ZDUMP=zd, **({'V9_FIRST_ONLY': '1'} if fast else {})))
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
            gates.glog(H, fastpath=True, returned=new1, orig=br1, **info1)
            if new1 is not None:
                with open(zd, 'w') as g: g.writelines(f'{a!r} {b!r}\n' for a, b in new1)
                d1.update(status='PROVEN', found=len(new1)); info = info1
                p = subprocess.CompletedProcess(p.args, 0, stdout=json.dumps(d1) + '\n', stderr='')
                notes.append(f"HEUR-gate:close_pair(centre-out, Nathan A. Meeks) chunk u[{info['chunk'][0]},{info['chunk'][1]}] +{info['recovered']} fastpath")
        if p.returncode != 0:   # not resolved by the fast path: full v9 ladder exactly as before
            err += '--- fastpath miss: full ladder ---\n'
            p = V9CALL(T, P, B, dict(env, ZDUMP=zd))
            att += 1 + p.stderr.count('escalate'); err += p.stderr
    if p.returncode != 0 and not nstep:   # extra driver rung (search-only): rerun the whole v9 ladder with refine round cap MAXR=200 instead of 20
        p = V9CALL(T, P, B, dict(env, ZDUMP=zd, MAXR='200'))
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
        gates.glog(H, returned=new, orig=br0, **info)
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
    return dict(T=T, i=i, lo=lo, hi=hi, p=p, d=d, st=st, found=found, Na=Na, Nb=Nb, tur=tur, nbr=nbr, notes=notes, w=time.time() - s, att=att, braw=braw, zd=zd, info=info, pre=pre, hdr=hdr_line)
def _runlock():
    # t106 B3 (re-audit, Nathan A. Meeks & Grok): the driver takes ITS OWN locks with fresh, non-inheritable fds; inherited fds are never reused.
    #  (1) {H}/driver.lock (separate from run_batch.sh's run.lock -> no deadlock under run_batch.sh), LOCK_EX|LOCK_NB, then fstat ino == stat(path) ino;
    #  (2) the home directory itself (cannot be unlinked while it holds files), LOCK_EX|LOCK_NB.  Either busy/mismatch -> None (exit 4).
    import fcntl
    fn = f'{H}/driver.lock'; fds = []
    try:
        dfd = os.open(H, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC); fds.append(dfd); fcntl.flock(dfd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fd = os.open(fn, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o644); fds.append(fd); fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if os.fstat(fd).st_ino != os.stat(fn).st_ino: raise OSError('driver.lock replaced/deleted between open and flock')
        global _LOCKFD; _LOCKFD = (fd, fn); return fd
    except OSError:
        for x in fds: os.close(x)
        return None
_LOCKFD = None
def _lock_still_ours():
    # re-checked before every commit: driver.lock must still be the inode we locked
    try: return _LOCKFD is not None and os.fstat(_LOCKFD[0]).st_ino == os.stat(_LOCKFD[1]).st_ino
    except OSError: return False
def main():
    if _runlock() is None:
        print(f'[NUM] batch SKIPPED: another driver holds {H}/driver.lock or {H} (t106 driver lock)'); return 4
    c = load(); c['params'] = dict(A=A, P=P, B=B, OVL=OVL, STEP=STEP)
    if os.path.exists(FAIL):
        line = f"[NUM] batch HALTED: open failure recorded in {FAIL} (window T={c['next_T']}); covered to {c['covered_hi']:.9e}; fix/resolve then move FAILED.md aside"
        c['last_run'] = now(); c['last_result'] = line; save(c); print(line); return 3
    _t0 = time.time(); CKE = int(os.environ.get('CK_EVERY', '1')); t103_restart(c); _t_restart = time.time() - _t0
    save(c)   # t106 B2: persist the truncation point at once, before any work
    _ts = time.time(); crashed = c.get('clean') is not True   # t106b: no clean-shutdown marker (crash/kill/legacy ck) -> full sweep
    bad = []; vinfo = 'off'
    if os.environ.get('VERIFY_BRACKETS', '1') == '1':
        stale = time.time() - c.get('last_full_verify', 0) > 86400   # t106b: forced full sweep at least once per 24 h
        full = crashed or stale or os.environ.get('VERIFY_FULL', '0') == '1'   # VERIFY_FULL=1 REQUIRED before any export / referee hand-off
        bad, m_, nt_ = verify_brackets(c, full=full)
        if bad and not full: bad, m_, nt_ = verify_brackets(c, full=True); full = True   # mismatch -> full sweep before halting
        vinfo = f"{'full' if full else 'incremental'} checked={m_}/{nt_} bad={len(bad)}"
    open(f'{H}/logs/t103.log', 'a').write(f'{now()} t106b startup verify {vinfo} crash_flag={crashed} sec={time.time() - _ts:.2f}\n')
    if not bad:
        c['verified_rows'] = c['ledger_rows']
        if vinfo.startswith('full'): c['last_full_verify'] = time.time(); c['last_full_verify_iso'] = now()
    c['clean'] = False; save(c)   # t106b: unclean marker while running; set True only on a normal return
    if bad:   # t106 (e): refuse on any committed bracket file whose sha differs from the ledger (forged / corrupted / missing)
        msg = '; '.join(f'idx {i}: ledger {w_[:12]} file {(g_ or "MISSING")[:12]}' for i, w_, g_ in bad[:10])
        with open(FAIL, 'a') as f: f.write(f'# OPEN FAILURE {now()}\nt106 brackets_sha mismatch on resume ({len(bad)} files): {msg}\nNothing run; checkpoint unchanged.\n\n')
        print(f'[NUM] batch HALTED: brackets_sha mismatch on resume ({len(bad)}): {msg}'); return 3
    env = dict(os.environ, AHCACHE=f'{H}/ahcache.bin', CPU=CPU)
    xz = []; t_start = time.time(); h0 = c['covered_hi']; z0 = c['zeros']; done = 0; rc = 0; failmsg = ''
    jl = open(f'{H}/logs/windows.jsonl', 'a'); XZB = []; XZQ = []; FIN = []; cl = open(CLOG, 'a')
    LIDX = set(); LAST = None   # t106 B4 (re-audit): committed (status PROVEN) idx in the ledger, and the last committed (idx, hi)
    if os.path.exists(LED):
        for l_ in open(LED):
            x_ = l_.split(',')
            if x_[0].isdigit() and len(x_) > 8 and x_[8] == 'PROVEN': LIDX.add(int(x_[0])); LAST = (int(x_[0]), int(x_[3]), int(x_[1]))
    from concurrent.futures import ThreadPoolExecutor
    import queue
    cpus = [x for x in os.environ.get('WORKER_CPUS', '1,2,3,4,5,6').split(',') if x]; cq = queue.Queue(); [cq.put(x) for x in cpus]
    def job(T, i):
        cpu = cq.get()
        try:
            r = window(T, i, dict(env, CPU=cpu)); FIN.append(i); return r   # t106: record real finish order (may differ from idx order)
        finally: cq.put(cpu)
    # t98: pre-extend the shared AH cache once (sequential) for the batch's top window, so parallel workers only read it
    # (ah_load is checksummed: a torn read can only force a recompute, never a wrong value)
    _t1 = time.time(); Te = c['next_T'] + (NWIN - 1) * STEP
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
        open(f'{H}/logs/t103.log', 'a').write(f'{now()} t106b startup phases: restart={_t_restart:.2f}s ah_preextend+pool={time.time() - _t1:.2f}s\n')
        indgates.grs_reset_timeouts()
        gates.SHARED_POOL = gp; gates.SCHED = Sched(int(os.environ['GATE_K']) if 'GATE_K' in os.environ else (int(os.environ.get('S3_K', '4')) if os.environ.get('S3', '0') == '1' else 2))   # t106 S3: 4 gated windows in parallel unless GATE_K is set explicitly (GATE_K wins)
    ex = ThreadPoolExecutor(len(cpus)); pend = []; st8 = dict(T=c['next_T'], i=c['next_idx'], n=0, head=c['next_idx'])
    def fill():
        while len(pend) < int(os.environ.get('LOOKAHEAD', '48')) * len(cpus) and st8['n'] < NWIN:   # consecutive windows in idx order
            pend.append(ex.submit(job, st8['T'], st8['i'])); st8['T'] += STEP; st8['i'] += 1; st8['n'] += 1
    fill()   # t101 (Nathan A. Meeks & Grok): lookahead 4->48 per CPU; a ~40 s gated window no longer stalls all 8 workers (ordered commit unchanged)
    while pend:
        R = pend.pop(0).result(); pre = R['pre']
        T, i, lo, hi, p, d, st, found, Na, Nb, tur, nbr, notes, w, att, braw, zd, info = (R[k] for k in ('T','i','lo','hi','p','d','st','found','Na','Nb','tur','nbr','notes','w','att','braw','zd','info'))
        # (t103 assert replaced by the explicit t106 B4 guard below, which also survives python -O)
        # t106 B4 (every window, incl. the first of a run): adjacency to the committed prefix, idx == checkpoint, idx never already in ledger
        gerr = []
        # expectation from the LEDGER (not the checkpoint the job was generated from): i = last committed idx + 1, T = T0 + i*STEP, lo = last hi - OVL
        if LAST is not None:
            if i != LAST[0] + 1: gerr.append(f'idx {i} != last ledger idx {LAST[0]} + 1')
            if T != LAST[2] + STEP: gerr.append(f'T {T} != last ledger T {LAST[2]} + STEP')   # (live grid is not T0+i*STEP: idx/T offset from history)
            if lo != LAST[1] - OVL: gerr.append(f'lo {lo} != last ledger hi {LAST[1]} - OVL')
        else:   # empty ledger (fresh home / first ever window): checkpoint and grid only
            if i != c['next_idx'] or T != c['next_T']: gerr.append(f'idx/T {i}/{T} != checkpoint {c["next_idx"]}/{c["next_T"]}')
        if c['windows'] > 0 and float(lo) != c['covered_hi'] - OVL: gerr.append(f"lo {lo} != covered_hi-OVL {c['covered_hi'] - OVL}")
        if i in LIDX: gerr.append(f'idx {i} already committed in ledger')
        if not _lock_still_ours(): gerr.append('driver.lock no longer ours (deleted/replaced)')
        if gerr:
            msg = f't106 gap/dup guard REFUSED window idx={i}: ' + '; '.join(gerr)
            with open(FAIL, 'a') as f: f.write(f'# OPEN FAILURE {now()}\n{msg}\nNothing committed for this window; checkpoint unchanged.\n\n')
            raise RuntimeError(msg)
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
        bsha = ''
        if ok:   # t106 (e)+(f): final bracket file = header + certified brackets, written atomically BEFORE the ledger row; its sha goes in the row
            with open(zd, 'rb') as g: data = R.get('hdr', '#\n').encode() + g.read()
            atomic_write(braw, data); bsha = hashlib.sha256(data).hexdigest()
        row = [i, T, lo, hi, tur, found, Na, Nb, st, d.get('v9_delta', ''), d.get('v9_P', ''), att, d.get('sec_all_attempts', ''), round(w, 3), ovN, ovB, now(), bsha, ';'.join(notes)]
        hdr = not os.path.exists(LED); ledger(row); c['ledger_rows'] = c.get('ledger_rows', 0) + 1 + hdr   # rows = file lines incl. header
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
        os.remove(zd)
        XZB.append(braw)
        if c['N_start'] is None: c['N_start'] = Na
        c['N_last_hi'] = Nb; c['zeros'] = Nb - c['N_start']; c['covered_hi'] = float(hi); c['windows'] += 1
        c['next_T'] = T + STEP; c['next_idx'] = i + 1; done += 1
        cl.write(json.dumps(dict(commit=i, finished_so_far=FIN[:], t=round(time.time(), 3))) + '\n'); cl.flush(); os.fsync(cl.fileno()); c['clog_bytes'] = cl.tell(); LIDX.add(i); LAST = (i, hi, T)   # t106 commit log (B2: fsync)
        if done % CKE == 0 or not pend:
            save(c)
            if XZSTAG: XZQ += XZB; XZB = []
            if XZB:   # xz only committed windows' files, one process per checkpoint group
                while len([q for q in xz if q.poll() is None]) >= XZPAR: time.sleep(0.05)
                xz = [q for q in xz if q.poll() is None] + [subprocess.Popen(['taskset', '-c', XZCPUS, 'nice', 'xz', '-T1', XZL, '-f'] + XZB)]; XZB = []
        if XZSTAG: xz = xz_tick(xz, XZQ)
        st8['head'] = i + 1
        fill()
    ex.shutdown(wait=True)
    if gates.SHARED_POOL is not None: gates.SHARED_POOL.close(); gates.SHARED_POOL.join()
    save(c)
    if XZSTAG:
        XZQ += XZB; XZB = []
        while XZQ or any(p.poll() is None for p in xz): xz = xz_tick(xz, XZQ, drain=True); time.sleep(0.05)
    if XZB: xz.append(subprocess.Popen(['taskset', '-c', XZCPUS, 'nice', 'xz', '-T1', XZL, '-f'] + XZB)); XZB = []
    for q in xz:
        if q.wait(): open(f'{H}/logs/xz.log', 'a').write(f'{now()} xz rc={q.returncode} (logged, batch not failed)\n')
    el = time.time() - t_start
    line = (f"[NUM] batch {'OK' if rc == 0 else 'FAILED'}: T {h0:.9e}->{c['covered_hi']:.9e}, +{c['zeros'] - z0} zeros, {done} windows PROVEN, "
            f"total covered to {c['covered_hi']:.9e} ({c['zeros']} zeros since 3e12); failures {1 if rc else 0}; {el/60:.1f} min"
            + (f' | FAIL {failmsg}' if rc else '') + ' (candidate, unrefereed; not RH)')
    c['last_run'] = now(); c['last_result'] = line; c['clean'] = (rc == 0); save(c)   # t106b clean-shutdown marker
    open(f'{H}/logs/runs.log', 'a').write(c['last_run'] + ' ' + line + '\n')
    print(line); return rc
if __name__ == '__main__':
    try: sys.exit(main())
    except SystemExit: raise
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f'[NUM] batch ERROR: driver crashed ({type(e).__name__}: {e}); checkpoint unchanged past last PROVEN window; see logs/driver.err'); sys.exit(5)
