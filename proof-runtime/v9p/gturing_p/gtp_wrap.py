# [NUM] t102 G-Turing-P middle tier wrapper (Nathan A. Meeks & Grok; G-Turing-P by CMather, audited in RESULTS.md §C). Candidate, unrefereed; not RH.
# Source gturing_p.py/pme are byte-identical copies (SOURCE.sha256); this wrapper only points signs() at the absolute pme path (same args).
import os, sys, subprocess
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); import gturing_p as _g
def _signs(a, h):
    o = subprocess.run([f'{D}/pme', str(a), str(h), '1', '128', '16', '4096'], capture_output=True, text=True, check=True)
    return [(float(l.split()[0]), int(l.split()[1]), float(l.split()[2])) for l in o.stdout.split('\n') if l.strip()]
_g.signs = _signs
def N(a): return _g.N(a)   # (N or None, L, U, n_points); N only when certified bounds coincide
