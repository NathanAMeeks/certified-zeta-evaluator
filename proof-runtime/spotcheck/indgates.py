# [NUM] Independent-check gates for OPEN-1 (Nathan A. Meeks & Grok). Candidate, unrefereed; not RH.
# Shares NO code with v9/fftzi: Z(t) from Arb (python-flint acb.zeta, rigorous ball arithmetic), theta from Arb lgamma.
#  G-cut    : Dedekind-cut bisection (method: Nathan A. Meeks) on any endpoint whose Arb sign is uncertain.
#  G-Turing : independent N(t) at block ends, Turing's method with Trudgian's bound |int_{t1}^{t2} S| <= 2.067+0.059 log t2
#             (Trudgian, Math. Comp. 80 (2011) 2259-2279, t2>t1>=168pi).
#  G-RS     : independent direct Z at every bracket endpoint of the block.
import os, sys, json, lzma, math, time
from multiprocessing import Pool
from flint import arb, acb, ctx
TB = (arb('2.067'), arb('0.059'))
def Z(t, prec=64):
    old = ctx.prec; ctx.prec = prec
    try:
        t = arb(t) if not isinstance(t, arb) else +t
        th = acb(arb('0.25'), t / 2).lgamma().imag - t / 2 * arb.pi().log()
        return (acb(0, th).exp() * acb(arb('0.5'), t).zeta()).real
    finally: ctx.prec = old
def sgn(t):
    """certified sign of Z at exact decimal point t (string), raising precision; 0 = uncertain."""
    for p in (128, 256, 512):
        z = Z(t, p)
        if not z.contains(0): return 1 if z > 0 else -1
    return 0
def gcut(p, left, right, maxit=60):
    """[G-cut] Dedekind cut (method: Nathan A. Meeks). Endpoint p (exact decimal string) has an uncertain Arb sign. Replace it by
    the nearest point p' = p +/- eps (eps = 1e-12 * 2^k) strictly between its neighbours (left, right) in the sorted endpoint list
    whose sign IS certified. Order is preserved, so brackets stay disjoint; returns (p', sign) or None."""
    from decimal import Decimal as D, getcontext
    getcontext().prec = 60; P, L, R = D(p), D(left), D(right); e = D('1e-12')
    for _ in range(maxit):
        for q in (P - e, P + e):
            if L < q < R:
                s = sgn(str(q))
                if s: return str(q), s
        e *= 2
    return None
def _sg(x): return sgn(x)
def theta_int(a, b):
    """certified enclosure of int_a^b theta(t) dt via Stirling antiderivative; |remainder| <= 1/t^3 (t>=1e3) (Gabcke 1979)."""
    a, b = arb(a), arb(b); tp = 2 * arb.pi()
    F = lambda t: t * t / 4 * (t / tp).log() - 3 * t * t / 8 - arb.pi() * t / 8 + t.log() / 48
    r = (b - a) / (a * a * a)
    return F(b) - F(a) + arb(0, r.upper())
def _sv(x):
    for p in (128, 256, 512):
        z = Z(x, p)
        if not z.contains(0): return float(z.mid())
    return 0.0
import threading as _th; _TL = _th.local()
def _vc(): return _TL.__dict__.setdefault('d', {})   # t98: per-thread (parallel workers)
# cleared per resolve() call (t98) so values are shared only within one block's G-Turing+G-RS   # [HYB] (Nathan A. Meeks & Grok): Arb-certified values shared by G-Turing and G-RS; identical sign logic (_sv/sgn)
def _cmap(M, xs):
    need = [x for x in dict.fromkeys(xs) if x not in _vc()]
    _vc().update(zip(need, M(_sv, need)))
    return [_vc()[x] for x in xs]
def _gt_vals(xs):
    """t106b GT_BATCH (Nathan A. Meeks & Grok; candidate, needs Audit): G-Turing values from the pinned certified batch tiers instead of Arb _sv.
    Sign comes only from a tier's certified sign (|mid| > rad, Gabcke t>200 guard); the magnitude is the tier's mid (search-only, as with _sv's mid).
    Any point a tier leaves at 0 goes to Arb _sv exactly as before. Points are passed as u = x - T0 (T0 integer), same contract as G-RS."""
    from decimal import Decimal as D
    if not xs: return []
    T0 = int(D(min(xs, key=D))) - 1; us = [str(D(x) - T0) for x in xs]
    sg = [0] * len(xs); md = [0.0] * len(xs); ok = sorted((k for k, u in enumerate(us) if not _t_le_200(T0, u)), key=lambda k: D(us[k]))   # sorted: _rsbatch1 groups by |u-c0|<=900
    for name, b in _grs_tiers():
        miss = [k for k in ok if sg[k] == 0]
        if not miss: break
        mm = [0.0] * len(miss); sub = _rsbatch1(b, T0, [us[k] for k in miss], mm)
        for k, x, m in zip(miss, sub, mm):
            if x: sg[k] = x; md[k] = abs(m) * x
    out = [md[k] if sg[k] else None for k in range(len(xs))]
    need = [k for k, v in enumerate(out) if v is None]
    if need:
        fb = list(map(_sv, [xs[k] for k in need]))
        for k, v in zip(need, fb): out[k] = v
    _TL.__dict__['gt_batch'] = _TL.__dict__.get('gt_batch', 0) + len(xs) - len(need); _TL.__dict__['gt_arb'] = _TL.__dict__.get('gt_arb', 0) + len(need)
    return out
def gturing(a, h=24, step=0.05, pool=None, rounds=4, extra=(), meshless=False):
    """[G-Turing] independent exact N(a) via Turing/Trudgian on [a-h,a] and [a,a+h]; returns (N or None, lo, up, evals).
    Mesh is refined (search only) at same-sign local |Z| dips so close pairs are seen; soundness never depends on the search."""
    from decimal import Decimal as D
    a = int(a); n = int(round(2 * h / step)); ctx.prec = 256
    pts = [D(a - h) + D(str(step)) * k for k in range(n + 1)] if not meshless else [D(a - h), D(a + h)]
    pts += [D(x) for x in extra if a - h < D(x) < a + h]   # search-only: known bracket endpoints (catch close pairs); soundness unaffected
    pts = sorted(set(pts))
    M = pool.map if pool else map
    C = (lambda xs: _cmap(M, xs)) if os.environ.get('HYB_VC', '0') == '1' else (lambda xs: list(M(_sv, xs)))
    if os.environ.get('GT_BATCH', '1') == '1' and os.environ.get('GRS_T100', '1') == '1' and _grs_tiers(): C = _gt_vals   # t106b
    V = dict(zip(pts, C([str(x) for x in pts]))); ev = len(pts)
    for _ in range(rounds):
        P = sorted(V); add = []
        for k in range(1, len(P) - 1):
            v0, v1, v2 = V[P[k - 1]], V[P[k]], V[P[k + 1]]
            if v0 * v1 > 0 and v1 * v2 > 0 and abs(v1) < abs(v0) and abs(v1) < abs(v2):
                add += [(P[k - 1] + P[k]) / 2, (P[k] + P[k + 1]) / 2]
        add = [x for x in dict.fromkeys(add) if x not in V]
        if not add: break
        V.update(zip(add, C([str(x) for x in add]))); ev += len(add)
    P = sorted(V); S = [V[x] for x in P]
    cells = [(float(P[k] - a), float(P[k + 1] - a)) for k in range(len(P) - 1) if S[k] and S[k + 1] and (S[k] > 0) != (S[k + 1] > 0)]
    A = arb(a); Hh = arb(h); bnd_r = TB[0] + TB[1] * (A + Hh).log(); bnd_l = TB[0] + TB[1] * A.log()
    intL = sum(h - r for l, r in cells if l >= 0)          # right: N(t) >= N(a)+L(t)
    # t106b (Audit N2): the arb(0, 1e-8) slack below covers the float/format error of intL/intR (cells as floats, printed to 9 decimals: <=5.1e-10),
    # NOT any evaluator shift. NEVER reduce it below 1e-9.
    up = (bnd_r + theta_int(A, A + Hh) / arb.pi() + Hh - arb(f'{intL:.9f}') - arb(0, 1e-8)) / Hh
    intR = sum(l + h for l, r in cells if r <= 0)          # left:  N(t) <= N(a)-R(t)
    lo = (theta_int(A - Hh, A) / arb.pi() + Hh + arb(f'{intR:.9f}') - arb(0, 1e-8) - bnd_l) / Hh
    U = int(math.floor(float(up.upper()))); Lw = int(math.ceil(float(lo.lower())))
    return (U if U == Lw else None), Lw, U, ev
def gturing_fast(a, pool, ext):
    """[t92, Nathan A. Meeks & Grok] Try the cheap Turing/Trudgian configuration first (h=12, step 0.1), then the old one
    (h=24, step 0.05). Sound: a configuration is accepted only if its certified lower and upper bounds for N(a) coincide; the mesh
    is search-only, so a coarse mesh can only make the bounds fail to meet (-> fallback), never give a wrong N."""
    ev = 0
    cf = (((12, 0.1, True),) if os.environ.get('HYB_TMESH', '0') == '1' else ()) + (('P', 0, 0),) + ((12, 0.1, False), (24, 0.05, False))
    for h, st, ml in cf:   # [HYB] meshless first: Arb signs at known endpoints only (mesh is search-only; bounds decide)
        s0 = time.time()
        if h == 'P':   # t102 middle tier (Nathan A. Meeks & Grok): G-Turing-P (Platt multieval mesh, RESULTS.md §C); accepted only if its bounds coincide
            if os.environ.get('GTURING_P', '1') != '1' or not os.path.exists(GTP + '/pme'): continue
            try: r = pool.apply(_gtp, (int(a),)) if pool else _gtp(int(a))
            except Exception: r = (None, None, None, 0)
            _tier().append(('G-Turing-P', round(time.time() - s0, 2), r[0] is not None))
        else:
            r = gturing(a, h=h, step=st, pool=pool, extra=ext, meshless=ml)
            _tier().append((f"{'meshless' if ml else 'mesh'}(h={h},step={st})", round(time.time() - s0, 2), r[0] is not None))
        ev += r[3]
        if r[0] is not None: return r[0], r[1], r[2], ev
    if os.environ.get('T104_DIP', '1') == '1':   # t104 (Nathan A. Meeks & Grok, t132s656): end count not pinned -> deeper same-sign dip refinement
        # on [a-h, a+h] (covers [hi, hi+h] and [lo-h, lo]); extra certified Arb points are search only; N accepted only if bounds coincide
        for h, st in ((12, 0.1), (24, 0.05)):
            s0 = time.time(); r = gturing(a, h=h, step=st, pool=pool, extra=ext, rounds=int(os.environ.get('T104_ROUNDS', '10'))); ev += r[3]
            _tier().append((f't104-dip(h={h},step={st})', round(time.time() - s0, 2), r[0] is not None))
            if r[0] is not None: return r[0], r[1], r[2], ev
    return None, r[1], r[2], ev
GTP = os.environ.get('GTP_DIR', '/workspace/rh_extend/v9p/gturing_p')
def _tier(): return _TL.__dict__.setdefault('tier', [])
def _gtp(a):
    if GTP not in sys.path: sys.path.insert(0, GTP)
    import gtp_wrap; return gtp_wrap.N(a)
def _ep(args):
    T, u = args; return sgn(_exact(T, u))
def _exact(T, u):
    # exact decimal T+u as a string (u is the dumped decimal of a double)
    from decimal import Decimal, getcontext
    getcontext().prec = 60; return str(Decimal(T) + Decimal(u))
RSB = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rsbatch_t100')   # [NUM] t100 G-RS tier 1 binary (WP G.5.t100)
RSB = os.environ.get('GRS_BIN', RSB)   # t106 [NUM] candidate, unrefereed (Nathan A. Meeks & Grok): batch G-RS evaluator selectable (rs_fast/rs_fast512/rs_tb_fix, same stdin/stdout contract); sign 0 -> unchanged Arb sgn ladder
# t106 B1 (Nathan A. Meeks & Grok; Audit t106): pin the G-RS batch binary and libflint by sha256 once at import. Allowlist =
# grs_stage/BUILD_SHA256.txt bin/ lines + the live t100 binary (rsbatch_t100). Any mismatch/missing -> RSB=None (all endpoints go to Arb), logged.
import hashlib as _hl
def _sha(p):
    h = _hl.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
# t106 B1 re-audit (Nathan A. Meeks & Grok): allowlist is PINNED IN CODE (no env override). Values = grs_stage/BUILD_SHA256.txt bin/ lines as of
# 2026-11-01 19:23 MST + live rsbatch_t100. The binary and its three libraries are COPIED into a private mkdtemp dir (0700, files 0500); the COPIES are
# hashed and only the copies are executed, with a minimal environment (no inherited LD_* / GLIBC_TUNABLES; LD_LIBRARY_PATH = the private dir only).
GRS_ALLOW = {'01e9ad3b012db45986e75a66ad865584389027a58b9c6cb7e274e136e3df03a3',   # rs_fast (x86-64-v4, prefer-256)
             '6f27558e2a5b6969d28765469ff0963854cac6c3f7d2a6ceaff1ab39ccea3678',   # rs_fast512
             '0a113f1ee25955aacf44eb98a06e9aa5338d27c39cf50cc1db6b290c2ac47082',   # rs_tb_fix
             '4dc5a03485c1b58f92469b695b5b91c0079ad40e738cee56a0da11deed7b554c'}   # live spotcheck/rsbatch_t100 (t100, WP G.5.t100)
LIB_PIN = {'libflint.so.24': ('/workspace/rh-research/tools/flintc/root/usr/lib/x86_64-linux-gnu/libflint.so.24', 'd4d65f1045ba48bba16dd796c1b2003c4b599b5c8741f42f61fed47db55d8e2d'),
           'libmpfr.so.6': ('/lib/x86_64-linux-gnu/libmpfr.so.6', '4e3ccf5b442fc91054b0153140fe62939009f98ca1367f45de62376d5b18d28c'),
           'libgmp.so.10': ('/lib/x86_64-linux-gnu/libgmp.so.10', 'd15e1585c6c267296459570c8ac217ca9123044bcd4e896d2f686b32333445d2')}
RSB_PIN = {}; RSB_ENV = None; _PIN_DIR = None; _PINNED = {}
def _pinlog(msg):
    sys.stderr.write('[t106 G-RS pin] ' + msg + '\n')
    try:
        d = os.path.join(os.environ.get('RHX_HOME', '/tmp'), 'logs')
        if os.path.isdir(d): open(os.path.join(d, 'grs_pin.log'), 'a').write(time.strftime('%Y-%m-%dT%H:%M:%S%z ') + msg + '\n')
    except Exception: pass
def _pin_dir():
    # one private dir per process: 0700 mkdtemp; pinned library copies (0500) checked once; minimal env points LD_LIBRARY_PATH here only
    global _PIN_DIR, RSB_ENV
    if _PIN_DIR is not None: return _PIN_DIR or None
    import tempfile, shutil, atexit
    try:
        d = tempfile.mkdtemp(prefix='grs_pin_'); os.chmod(d, 0o700); atexit.register(shutil.rmtree, d, True)
        for name, (src, want) in LIB_PIN.items():
            lc = os.path.join(d, name); shutil.copyfile(os.path.realpath(src), lc); os.chmod(lc, 0o500); lh = _sha(lc); RSB_PIN[name] = lh
            if lh != want: _pinlog(f'REFUSED {name} sha {lh} != pinned {want}; G-RS -> Arb only'); _PIN_DIR = ''; return None
        RSB_ENV = {'PATH': '/usr/bin:/bin', 'LD_LIBRARY_PATH': d, 'LC_ALL': 'C'}   # minimal env: no inherited LD_* / GLIBC_TUNABLES
        _PIN_DIR = d; return d
    except Exception as e:
        _pinlog(f'REFUSED pin dir: {type(e).__name__}: {e}; G-RS -> Arb only'); _PIN_DIR = ''; return None
def _pin_rsb(p):
    """private 0500 copy of binary p in the pin dir, hashed AFTER copying; returns the copy's path (the only thing ever executed) or None."""
    import shutil
    if p in _PINNED: return _PINNED[p]
    r = None
    try:
        d = _pin_dir()
        if d and p and os.path.exists(p):
            cp = os.path.join(d, f'tier{len(_PINNED)}_' + os.path.basename(p)); shutil.copyfile(p, cp); os.chmod(cp, 0o500); h = _sha(cp)
            RSB_PIN[os.path.basename(p)] = h
            if h in GRS_ALLOW: r = cp
            else: _pinlog(f'REFUSED {p} sha {h} not in pinned allowlist; this tier skipped')
    except Exception as e: _pinlog(f'REFUSED {p}: {type(e).__name__}: {e}; this tier skipped')
    _PINNED[p] = r; return r
RSB = _pin_rsb(RSB)
# [NUM] G-RS tiers (merged from grs_stage/integration/indgates.py.patch, t106; Nathan A. Meeks & Grok). GRS_TIER1=1: rs_tb_fix, then on its sign-0 points
# rs_fast512 (or GRS_FAST=rs_fast; both need avx512f), then rsbatch_t100, then Arb (_ep). GRS_TIER1=0 (default): [RSB] only (= t100 / GRS_BIN).
# EVERY tier binary, rsbatch_t100 included, goes through _pin_rsb (private copy + sha allowlist + minimal env). Each tier is independently certified;
# a 0 only moves a point DOWN a tier.
GRS_DIR = os.environ.get('GRS_BIN_DIR', '/workspace/rh_extend/scratch/grs_stage/bin')
T100_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rsbatch_t100')
def _avx512():
    try: return ' avx512f' in open('/proc/cpuinfo').read()
    except Exception: return False
_TIERS = None
def _grs_tiers():
    global _TIERS
    if _TIERS is None:
        if os.environ.get('GRS_TIER1', '0') != '1': src = [os.environ.get('GRS_BIN', T100_BIN)]
        else: src = [os.path.join(GRS_DIR, 'rs_tb_fix')] + ([os.path.join(GRS_DIR, os.environ.get('GRS_FAST', 'rs_fast512'))] if _avx512() else []) + [T100_BIN]
        if os.environ.get('PC_FAST', '0') == '1': src = [p for p in src if os.path.basename(p) != os.path.basename(_pc_bin())]   # t106c S1: referee (G-RS/G-Turing) never uses the pair_cert evaluator
        _TIERS = [(os.path.basename(p), c) for p in src for c in [_pin_rsb(p)] if c]
    return _TIERS
def _pc_bin():   # t106c S1: pair_cert sign evaluator = rs_fast512 (avx512f) else rs_fast; disjoint from the G-RS tiers when PC_FAST=1
    return os.path.join(GRS_DIR, os.environ.get('PC_BIN', 'rs_fast512' if _avx512() else 'rs_fast'))
def pc_fast_signs(T, us):
    """t106c S1 (Nathan A. Meeks & Grok; candidate, needs Audit): certified signs from the pinned pair_cert evaluator (sign 0 iff |Z| within its proven
    radius -> caller uses Arb sgn). Same u-as-double contract as G-RS (point = T + double(u), u = repr of the finder's double)."""
    b = _pin_rsb(_pc_bin())
    if not b or any(_t_le_200(T, u) for u in us): return None
    o = sorted(range(len(us)), key=lambda k: float(us[k])); r = _rsbatch1(b, T, [us[k] for k in o]); out = [0] * len(us)
    for k, v in zip(o, r): out[k] = v
    return out
_grs_tiers()   # pin all tiers once at import
def _grs_line(l):
    x = l.split()
    if len(x) != 3: raise ValueError(l)
    v = int(x[0]); float(x[1]); float(x[2])
    if v not in (1, 0, -1): raise ValueError(l)
    return v
def _t100(): return _TL.__dict__.setdefault('t100', {})   # t106: per-THREAD telemetry (S3 runs gated windows in parallel threads; a shared dict raced -> KeyError)
def _t_le_200(T, u):
    # t106 (c): Gabcke R0 bound needs t>=200; any point with t=T+u <= 200 is never sent to ANY batch tier (sign 0 -> Arb). Exact decimal compare.
    from decimal import Decimal, getcontext
    getcontext().prec = 60; return Decimal(T) + Decimal(u) <= 200
def _rsbatch(T, pts):
    """tiered: run each pinned binary only on points still 0 after the previous tier; t<=200 points skip all tiers."""
    sg = [0] * len(pts); _t100()['tiers'] = tl = []
    ok = [k for k, u in enumerate(pts) if not _t_le_200(T, u)]
    if len(ok) < len(pts): tl.append(('t<=200 -> Arb', len(pts) - len(ok), 0, 1.0))
    for name, b in _grs_tiers():
        miss = [k for k in ok if sg[k] == 0]
        if not miss: break
        sub = _rsbatch1(b, T, [pts[k] for k in miss])
        for k, x in zip(miss, sub): sg[k] = x
        got = sum(1 for x in sub if x)   # tier telemetry: (binary, points in, certified, refusal rate)
        tl.append((name, len(miss), got, round(1 - got / len(miss), 4)))
    return sg
GRS_TO_MIN = float(os.environ.get('GRS_TO_MIN', '30')); GRS_TO_PER_PT = float(os.environ.get('GRS_TO_PER_PT', '0.05')); GRS_TO_MAX = int(os.environ.get('GRS_TO_MAX', '2')); _GRS_TO = {}; _GRS_SKIP = {}; GRS_TO_PROBE = int(os.environ.get('GRS_TO_PROBE', '20'))
def grs_reset_timeouts(): _GRS_TO.clear(); _GRS_SKIP.clear()   # t106b (Audit): timeout count is per batch; driver calls this at batch start
def _rsbatch1(RSB, T, pts, mids=None):
    """signs from rsbatch (independent of acb_zeta/G-Turing and of v9); 0 = not certified. Points grouped so |u-c0|<=900."""
    import subprocess as sp
    sg = [0] * len(pts); k = 0
    while k < len(pts):
        c0 = int(float(pts[k])) + 900; j = k
        while j < len(pts) and float(pts[j]) - c0 <= 900: j += 1
        try:
            if _GRS_TO.get(RSB, 0) >= GRS_TO_MAX:   # t106c (Audit follow-up): dropped after GRS_TO_MAX timeouts IN A ROW; probed again every GRS_TO_PROBE groups
                _GRS_SKIP[RSB] = _GRS_SKIP.get(RSB, 0) + 1
                if _GRS_SKIP[RSB] % GRS_TO_PROBE: k = j; continue
            to = max(float(os.environ.get('GRS_TO_MIN', GRS_TO_MIN)), float(os.environ.get('GRS_TO_PER_PT', GRS_TO_PER_PT)) * (j - k))   # t106c: env read per call   # t106b (Nathan A. Meeks & Grok; Audit speed REPORT): scaled timeout, was fixed 600 s
            pr = sp.Popen([RSB, os.environ.get('RSB_THREADS', '1')], env=RSB_ENV, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, text=True, start_new_session=True)
            try: out, _ = pr.communicate(f'{T} {c0}\n' + '\n'.join(pts[k:j]) + '\n', timeout=to)
            except sp.TimeoutExpired:
                try: os.killpg(pr.pid, 9)   # kill child + its process group; points stay 0 -> Arb (certification unchanged)
                except Exception: pass
                try: pr.communicate(timeout=5)
                except Exception: pass
                _GRS_TO[RSB] = _GRS_TO.get(RSB, 0) + 1
                _pinlog(f'GRS_TIMEOUT tier={os.path.basename(RSB)} T={T} pts={j - k} timeout={to:.2f}s n={_GRS_TO[RSB]}' + ('; tier DROPPED (probed again later)' if _GRS_TO[RSB] >= GRS_TO_MAX else ''))
                _t100().setdefault('timeouts', []).append((os.path.basename(RSB), j - k, round(to, 2)))
                k = j; continue
            ls = [l for l in out.split('\n') if l.strip()]
            if pr.returncode == 0 and len(ls) == j - k:
                vs = [_grs_line(l) for l in ls]
                if _GRS_TO.get(RSB, 0):   # t106c: a completed call resets the consecutive-timeout count (tier is back)
                    if _GRS_TO[RSB] >= GRS_TO_MAX: _pinlog(f'GRS_TIER_BACK tier={os.path.basename(RSB)} after successful probe')
                    _GRS_TO[RSB] = 0   # t106 B1: any unparsable line raises -> whole group stays 0 (Arb)
                for m, v in enumerate(vs): sg[k + m] = v
                if mids is not None:
                    for m, l in enumerate(ls): mids[k + m] = float(l.split()[1])
        except Exception: pass
        k = j
    return sg
def grs(T, br, pool=None):
    """[G-RS] independent Z sign at every bracket endpoint; G-cut on uncertain ones. Returns dict."""
    pts = sorted({u for ab in br for u in ab}, key=float)
    if os.environ.get('HYB_VC', '0') == '1':   # [HYB] reuse cached Arb values; _sv!=0 iff sgn certified, same sign
        vs = _cmap(pool.map if pool else map, [_exact(T, u) for u in pts]); sg = [(1 if v > 0 else -1) if v else 0 for v in vs]
    elif os.environ.get('GRS_T100', '1') == '1' and _grs_tiers():   # t100 tier 1 (Nathan A. Meeks & Grok): certified batch RS; Arb on any sign 0
        sg = _rsbatch(T, pts); miss = [k for k, x in enumerate(sg) if x == 0]
        if miss:
            fb = list(pool.map(_ep, [(T, pts[k]) for k in miss], chunksize=16) if pool else [_ep((T, pts[k])) for k in miss])
            for k, x in zip(miss, fb): sg[k] = x
        _t100()['fallback'] = len(miss)
    else: sg = list(pool.map(_ep, [(T, u) for u in pts], chunksize=64) if pool else [_ep((T, u)) for u in pts])
    cut = unres = 0
    for k, u in enumerate(pts):
        if sg[k]: continue
        cut += 1
        L = _exact(T, pts[k - 1]) if k else _exact(T, str(float(u) - 1e-3))
        R = _exact(T, pts[k + 1]) if k + 1 < len(pts) else _exact(T, str(float(u) + 1e-3))
        r = gcut(_exact(T, u), L, R)
        if r is None: unres += 1
        else: sg[k] = r[1]
    esc = None
    if unres and os.environ.get('GRS_ESC', '1') == '1':   # t99 escalation (Nathan A. Meeks & Grok), ONLY on still-unresolved endpoints
        esc = dict(n=unres, prec_fixed=0, bisect_fixed=0, left=0); s0 = time.time()
        for k, u in enumerate(pts):
            if sg[k]: continue
            x = _exact(T, u)
            for pr in (1024, 2048):                       # (1) higher Arb precision (128/256/512 already tried in sgn)
                z = Z(x, pr)
                if not z.contains(0): sg[k] = 1 if z > 0 else -1; esc['prec_fixed'] += 1; break
            if sg[k]: continue
            L = _exact(T, pts[k - 1]) if k else _exact(T, str(float(u) - 1e-3))   # (2) local referee bisection: midpoints of
            R = _exact(T, pts[k + 1]) if k + 1 < len(pts) else _exact(T, str(float(u) + 1e-3))   # [L,p],[p,R] (never v9 values)
            from decimal import Decimal as D
            for q in [(D(L) + D(x)) / 2, (D(x) + D(R)) / 2, (3 * D(L) + D(x)) / 4, (D(x) + 3 * D(R)) / 4]:
                sq = sgn(str(q))
                if sq: sg[k] = sq; esc['bisect_fixed'] += 1; break   # replacement point keeps order (strictly between neighbours)
            if not sg[k]: esc['left'] += 1
        unres = esc['left']; esc['sec'] = round(time.time() - s0, 2)
    S = dict(zip(pts, sg)); ok = bad = 0
    for a, b in br:
        if S[a] and S[b]:
            if S[a] != S[b]: ok += 1
            else: bad += 1
    r = dict(endpoints=len(pts), brackets=len(br), certified=ok, conflict=bad, gcut_used=cut, gcut_unresolved=unres)
    if 'fallback' in _t100(): r['rsbatch_fallback'] = _t100().pop('fallback')
    if 'tiers' in _t100(): r['grs_tiers'] = _t100().pop('tiers')
    if 'timeouts' in _t100(): r['grs_timeouts'] = _t100().pop('timeouts')
    if esc: r['escalation'] = esc
    return r
def read_br(fn):
    op = lzma.open if fn.endswith('.xz') else open
    with op(fn, 'rt') as f: return [tuple(l.split()[:2]) for l in f if l[0] != '#' and l.strip()]
def resolve(T, lo, hi, br, Na, Nb, trigger, procs=7, uncertain=(), hint=None, shared_pool=None):
    """Run gates in measured-cost order; stop at the earliest point soundness allows. Rule: GATES.md / WP G.5."""
    log = dict(trigger=trigger, steps=[]); t0 = time.time(); _vc().clear()
    ext = [_exact(T, u) for ab in br for u in ab if abs(T + float(u) - lo) < 25 or abs(T + float(u) - hi) < 25]
    br = [x for x in br if T + float(x[0]) >= lo and T + float(x[1]) <= hi]
    import contextlib   # t102 (Nathan A. Meeks & Grok): shared pinned orchestrator pool if given; same map calls, same order -> same values
    with (contextlib.nullcontext(shared_pool) if shared_pool is not None else Pool(procs)) as pool:
        if uncertain:   # 1. G-cut (cheapest) on given uncertain points (p, left, right); inside G-RS it runs on every Arb-uncertain endpoint
            s = time.time(); r = [gcut(_exact(T, p), _exact(T, l), _exact(T, rr)) for p, l, rr in uncertain]
            log['steps'].append(dict(gate='G-cut', sec=round(time.time() - s, 2), fixed=sum(x is not None for x in r), n=len(r)))
            if any(x is None for x in r): return 'FAIL', log
        pNa, pNb = (Na, Nb) if Na is not None else (hint or (None, None))   # t98 H2: v9 count is a hint only; Platt re-certifies N(lo),N(hi) or refuses
        if os.environ.get('PLATT', '1') == '1' and pNa is not None and len(br) >= int(os.environ.get('PLATT_MINZ', '300')):   # 2'. [t96/t98] Platt fast path; Arb for small blocks
            try:
                import plattgate; pg = plattgate.gate(T, lo, hi, br, pNa, pNb)
            except Exception as e: pg = None
            if pg is not None:
                log['steps'].append(pg); need = pNb - pNa
                res = 'RESOLVED' if (len(br) == need and pg['zeros'] == need and pg['brackets_matched'] == need) else 'FAIL'
                if res == 'RESOLVED':
                    log.update(need=need, total_sec=round(time.time() - t0, 1)); return res, log
            log['steps'].append(dict(gate='G-Platt', refused=True, note='fallback to Arb path'))
        _tier().clear(); s = time.time(); na = gturing_fast(lo, pool, ext); nb = gturing_fast(hi, pool, ext)   # 2. G-Turing (adaptive h/step; t92)
        log['steps'].append(dict(gate='G-Turing', sec=round(time.time() - s, 2), Na=na[:3], Nb=nb[:3], evals=na[3] + nb[3], tiers=list(_tier()))); _tier().clear()
        if na[0] is None or nb[0] is None or (Na is not None and (na[0] != Na or nb[0] != Nb)): return 'FAIL', log   # fail fast, before the expensive gate
        s = time.time(); g = grs(T, br, pool=pool)                                       # 3. G-RS (+G-cut inside)
        log['steps'].append(dict(gate='G-RS', sec=round(time.time() - s, 2), **g))
    need = nb[0] - na[0]
    if len(br) != need: g['coverage'] = dict(brackets=len(br), need=need, note='bracket count != G-Turing count on block (t99: check block lies inside window)')
    if (g.get('escalation') and g['gcut_unresolved'] and g['conflict'] == 0 and len(br) == need):   # (3) residue -> G-Platt with
        try:                                                                                          # G-Turing's (not v9's) N
            import plattgate; pg = plattgate.gate(T, lo, hi, br, na[0], nb[0])
        except Exception: pg = None
        g['escalation']['platt'] = 'RESOLVED CONDITIONAL(OPEN-P)' if pg else 'refused'
        if pg: g['certified'] = need; g['gcut_unresolved'] = 0
    res = 'RESOLVED' if (len(br) == need and g['certified'] >= need and g['conflict'] == 0 and g['gcut_unresolved'] == 0) else 'FAIL'
    log.update(need=need, total_sec=round(time.time() - t0, 1)); return res, log
if __name__ == '__main__':
    # usage: indgates.py T lo hi Na Nb bracketfile trigger [procs]
    T, lo, hi, Na, Nb = map(int, sys.argv[1:6])
    r, log = resolve(T, lo, hi, read_br(sys.argv[6]), Na, Nb, sys.argv[7], int(sys.argv[8]) if len(sys.argv) > 8 else 6)
    print(json.dumps(dict(result=r, **log)))
