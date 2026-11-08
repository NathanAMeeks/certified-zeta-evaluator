# v7 vs v8 vs v9 interleaved, same session, pinned CPU, 3 reps. Cache pre-extended (covers <=5.2e13) and read-only for all versions.
import subprocess,json,os,time,statistics as st,sys
D=os.environ.get('V9DIR',os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); C=os.environ['BENCHCACHE']; CPU=os.environ.get('CPU','7')
W=[(t,65536,10400,'bench') for t in (2999999999980,3000000019980,3000000039980,3000000010300,9999999999980,10000000010300,10000000020620,3000000051580,3000000061900,3000000072220)]
W+=[(3000000051580,131072,31040,'chain3e12'),(9999999999980,131072,31040,'chain1e13')]
W+=[(t,65536,10400,'high') for t in (4000000000000,6000000000000,8000000000000,15000000000000,20000000000000,20000000020800,30000000000000,50000000000000)]
def one(tag,T,P,B):
  e=dict(os.environ,AHCACHE=C+'.'+tag,CPU=CPU,V8DIR=D,V9DIR=D)
  cmd={'v7':['taskset','-c',CPU,D+'/fftzi_v7',str(T),'0.24',str(P),'80',str(B),'60','/dev/null','dbl'],'v8':[D+'/v8/v8run.sh',str(T),str(P),str(B)],'v9':[D+'/v9/v9run.sh',str(T),str(P),str(B)]}[tag]
  t0=time.time(); o=subprocess.run(cmd,env=e,capture_output=True,text=True); w=time.time()-t0
  d=json.loads(o.stdout.strip().splitlines()[-1]); d['wall']=w; d['tries']=o.stderr.count('not proven')+o.stderr.count('escalate'); return d
res=[]
for r in range(3):
  for T,P,B,g in W:
    for tag in ('v7','v8','v9'):
      d=one(tag,T,P,B); res.append(dict(tag=tag,g=g,T=T,P=P,rep=r,ms=d['ms_per_zero'],st=d['status'],found=d['found'],turing=int(d['Nb_upper'])-int(d['Na_lower']),wall=d['wall'],tries=d['tries'],sec=d['sec_total'],secall=d.get('sec_all_attempts',d['sec_total']),Pused=d.get('v9_P',d['P'] if 'P' in d else P),dl=d.get('v9_delta',d.get('v8_delta',0.24))))
      print(tag,T,P,r,d['status'],d['ms_per_zero'],round(d['wall'],3),flush=True)
json.dump(res,open(os.environ.get('OUT','bench.json'),'w'))
