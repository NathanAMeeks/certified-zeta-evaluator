# [NUM] offline test of independent gates on scratch copies (Nathan A. Meeks & Grok); candidate, unrefereed
import csv, glob, json, shutil, os, sys, indgates as g
os.makedirs('../scratch/indtest', exist_ok=True)
rows = {}
for r in csv.reader(open('../ledger.csv')):
    if r[0] != 'idx' and r[8] == 'PROVEN': rows[int(r[0])] = r
for i in map(int, sys.argv[1:]):
    r = rows[i]; T, lo, hi, Na, Nb = int(r[1]), int(r[2]), int(r[3]), int(r[6]), int(r[7])
    f = glob.glob(f'../brackets/w{i:07d}_*')[0]; sc = f'../scratch/indtest/{os.path.basename(f)}'; shutil.copy(f, sc)
    res, log = g.resolve(T, lo, hi, g.read_br(sc), Na, Nb, f'offline test idx={i}')
    out = dict(idx=i, result=res, **log); print(json.dumps(out), flush=True)
    open('ind_test.jsonl', 'a').write(json.dumps(out) + '\n')
