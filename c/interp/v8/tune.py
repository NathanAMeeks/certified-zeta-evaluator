# v8 gated self-tuner (search-only knobs). Proof logic untouched: every run must be PROVEN (Lemma CP) or the knob is rejected.
import subprocess,json,os,statistics as st,sys
import os as _o; D=_o.environ.get('V8DIR', _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))); BIN=sys.argv[1] if len(sys.argv)>1 else D+'/fftzi_v8'
W=[(t,65536,10400) for t in (2999999999980,3000000019980,3000000039980,3000000010300,9999999999980,10000000010300,10000000020620,3000000051580,3000000061900,3000000072220)]+[(3000000051580,131072,31040),(9999999999980,131072,31040)]
def run(env,reps=3,delta='0.24'):
  out=[]
  for r in range(reps):
    for T,P,B in W:
      e=dict(os.environ,AHCACHE=D+'/ahcache.bin',**env)
      p=subprocess.run(['taskset','-c','7',BIN,str(T),delta,str(P),'80',str(B),'60','/dev/null','dbl'],env=e,capture_output=True,text=True,cwd=D)
      d=json.loads(p.stdout); out.append((d['ms_per_zero'],d['status'],d['evals_per_zero'],d['sec_setup'],d['found'],d['Na_lower'],d['Nb_upper']))
  return out
def summ(o):
  m=[x[0] for x in o]; return dict(min=min(m),med=st.median(m),max=max(m),proven=sum(x[1]=='PROVEN' for x in o),n=len(o),evals=st.median([x[2] for x in o]),setup=st.median([x[3] for x in o]))
if __name__=='__main__':
  base=run({}); b=summ(base)
  # noise: spread of per-rep medians of baseline
  reps=[st.median([x[0] for x in base[i*12:(i+1)*12]]) for i in range(3)]; noise=max(reps)-min(reps)
  print('BASE',b,'noise',round(noise,5),flush=True)
  cands=[('KGATE','0.2'),('KGATE','0.3'),('KT5','0.4'),('KT5','0.6'),('KFLAT','0.05'),('KFLAT','0.2'),('NDIV','6'),('NDIV','10'),('DIPTAU','2.5'),('DIPTAU','3.5'),('MAXR','4'),('DELTA','0.22'),('DELTA','0.23')]
  kept={}
  for k,v in cands:
    o=run({},delta=v) if k=='DELTA' else run({k:v}); s=summ(o)
    ok=s['proven']==s['n'] and s['med']<b['med']-noise and all(x[4:]==y[4:] for x,y in zip(o,base))
    print('CAND',k,v,s,'KEEP' if ok else 'REVERT',flush=True)
    if ok: kept[k]=(v,s['med'])
  print('KEPT',kept)
