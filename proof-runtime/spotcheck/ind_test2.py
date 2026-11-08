# [NUM] offline chunk-block test of OPEN-1 gates (Nathan A. Meeks & Grok); candidate, unrefereed; not RH
import glob, json, os, shutil, sys, time
sys.path.insert(0, '..'); import gates, indgates as g
os.makedirs('../scratch/indtest', exist_ok=True)
for i, T, c0, c1 in [(10202, 3000105619520, 7542, 7866), (8758, 3000090775200, 8550, 8870), (34157, 3000351876920, 5000, 5320)]:
    f = sorted(glob.glob(f'../brackets/w{i:07d}_*'))[-1]; sc = f'../scratch/indtest/{os.path.basename(f)}'; shutil.copy(f, sc)
    br = g.read_br(sc); fb = sorted((float(a), float(b)) for a, b in br)
    L, R = gates.snap(c0, fb), gates.snap(c1, fb)
    res, log = g.resolve(T, T + L, T + R, br, None, None, f'offline chunk test idx={i}', procs=int(os.environ.get('P', '2')))
    out = dict(idx=i, block=[T + L, T + R], result=res, **log); print(json.dumps(out), flush=True)
    open('ind_test2.jsonl', 'a').write(json.dumps(out) + '\n')
