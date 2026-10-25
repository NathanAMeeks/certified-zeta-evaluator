import subprocess,json,os,statistics as st
D='/workspace/meeks_test/refine/speed2'
W=[(t,65536,10400) for t in (2999999999980,3000000019980,3000000039980,3000000010300,9999999999980,10000000010300,10000000020620,3000000051580,3000000061900,3000000072220)]+[(3000000051580,131072,31040),(9999999999980,131072,31040)]
NW=[(t,65536,10400) for t in (4000000000000,6000000000000,8000000000000,15000000000000,20000000000000,20000000020800)]
def one(tag,T,P,B):
  e=dict(os.environ,AHCACHE=D+'/ahcache.bin')
  cmd=['taskset','-c','7',D+'/fftzi_v7',str(T),'0.24',str(P),'80',str(B),'60','/dev/null','dbl'] if tag=='v7' else [D+'/v8/v8run.sh',str(T),str(P),str(B)]
  d=json.loads(subprocess.run(cmd,env=e,capture_output=True,text=True).stdout.strip().splitlines()[-1]); return d
res=[]
for r in range(3):
  for grp,L in (('bench',W),('new',NW)):
    for T,P,B in L:
      for tag in ('v7','v8'):
        d=one(tag,T,P,B); res.append(dict(tag=tag,grp=grp,T=T,P=P,rep=r,ms=d['ms_per_zero'],st=d['status'],ev=d['evals_per_zero'],setup=d['sec_setup'],tot=d['sec_total'],found=d['found'],Na=d['Na_lower'],Nb=d['Nb_upper'],dl=d.get('v8_delta',0.24)))
json.dump(res,open('final.json','w'))
for grp in ('bench','new'):
  for tag in ('v7','v8'):
    x=[r for r in res if r['tag']==tag and r['grp']==grp]; m=[r['ms'] for r in x]
    print(grp,tag,'min %.4f med %.4f max %.4f setup_med %.3f evals_med %.3f PROVEN %d/%d'%(min(m),st.median(m),max(m),st.median([r['setup'] for r in x]),st.median([r['ev'] for r in x]),sum(r['st']=='PROVEN' for r in x),len(x)))
for r in res:
  if r['rep']==0 and r['grp']=='new': print(r['tag'],r['T'],r['st'],r['found'],r['Na'],r['Nb'],r['ms'],r['dl'])
