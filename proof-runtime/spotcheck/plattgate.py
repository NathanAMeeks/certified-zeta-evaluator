# [NUM] Platt fast path for G-Turing + G-RS (Nathan A. Meeks & Grok). Candidate, unrefereed; not RH.
# Calls /workspace/opt/flint3/platt_iso (C wrapper over FLINT 3.6 acb_dirichlet_platt_isolate_local_hardy_z_zeros,
# Platt, Math. Comp. 86 (2017) 2449-2467). FLINT returns k<=len disjoint balls (refined, mode r); ball i is certified (Turing's
# method, rigorous Platt error bounds) to contain exactly one zero, the (n+i)-th, counted with multiplicity.
# Tuning parameters are chosen heuristically by FLINT; that affects only success/speed, never soundness. Valid 1e4<=n<=1e23.
# Any refusal/short return/mismatch -> None, and the caller falls back to the Arb path.
import subprocess, os, time
from fractions import Fraction as F
BIN = '/workspace/opt/flint3/platt_iso'
def isolate(n, ln, prec=256, threads=int(os.environ.get('PLATT_THREADS', '4')), timeout=1800):
    try: out = subprocess.run([BIN, str(n), str(ln), str(prec), str(threads), 'r'], capture_output=True, text=True, timeout=timeout).stdout.split('\n')
    except Exception: return None
    k = int(out[0]); r = []
    for l in out[1:1 + k]:
        ma, ea, mb, eb = map(int, l.split())
        r.append((F(ma) * F(2) ** ea, F(mb) * F(2) ** eb))
    return r
def gate(T, lo, hi, br, Na, Nb, margin=3):
    """Returns (dict) on success or None (refuse -> Arb fallback). Success means: the Platt intervals for zeros
    Na-margin+1 .. Nb+margin certify N(lo)=Na and N(hi)=Nb independently, and each v9 bracket in [lo,hi] contains exactly
    one Platt zero interval (so every bracket's sign change is an isolated, independently certified zero)."""
    if Na is None or Nb is None or Nb - Na > int(os.environ.get('PLATT_MAXZ', '4000')): return None   # big blocks -> Arb path
    s = time.time(); n0 = Na - margin + 1; ln = Nb - Na + 2 * margin
    for pr in [int(x) for x in os.environ.get('PLATT_PRECS', '128,256').split(',')]:   # t98 H1 (Nathan A. Meeks & Grok): precision only affects success, never soundness
        r = _check(T, lo, hi, br, Na, Nb, n0, ln, isolate(n0, ln, prec=pr))
        if r is not None: r.update(sec=round(time.time() - s, 2), prec=pr); return r
    return None
def _check(T, lo, hi, br, Na, Nb, n0, ln, Z):
    if Z is None or len(Z) != ln: return None
    lo, hi = F(lo), F(hi)
    j = Na - n0                                       # list position of zero Na
    if not (Z[j][1] < lo < Z[j + 1][0]): return None    # N(lo) = Na certified
    m = Nb - n0
    if not (Z[m][1] < hi < Z[m + 1][0]): return None    # N(hi) = Nb certified
    inner = Z[j + 1:m + 1]; B = sorted((F(T) + F(a), F(T) + F(b)) for a, b in br)
    if len(B) != len(inner): return None
    for (a, b), (za, zb) in zip(B, inner):
        if not (a < za and zb < b): return None       # one-to-one containment (both lists sorted, disjoint)
    return dict(gate='G-Platt', conditional='OPEN-P: FLINT zero labels may assume earlier verification; Arb remains certifying fallback', Na=Na, Nb=Nb, zeros=len(inner), brackets_matched=len(B))
