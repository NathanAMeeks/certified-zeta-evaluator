exec(open('fftarb.py').read())
import sys
T,dens,P,a,b,Lw={'a':(1e6,20.,131072,500.,6000.,250.),'a2':(1e6,80.,131072,300.,1300.,250.),'d':(3e12,250.,131072,60.,460.,50.),'d2':(3e12,60.,131072,300.,1800.,100.),'e':(1e13,250.,131072,60.,460.,50.),'c':(267653395647,250.,131072,60.,460.,50.)}[sys.argv[1]]
delta=1/dens
rf,uf,sgf=batch(T,0.0,delta,P)          # float path (fftb)
t1=time.time(); ra,Z=batch_arb(T,delta,P); ta=time.time()-t1
sg=np.array([0 if (z is None or z.contains(0)) else (1 if z>0 else -1) for z in Z])
# referee: observed float-vs-Arb error vs float bound (minus Gabcke)
exec(open('fftb.py').read().split('def batch')[0])
zf_mid=None
# recompute float Z values (batch returns signs only) -> rerun quickly with values
src=open('fftb.py').read().replace("return dict(N=N","return Z,B,dict(N=N")
exec(src)
Zf,Bf,_,_,_=batch(T,0.0,delta,P)
mid=np.array([float(z.mid()) if z is not None else np.nan for z in Z]); rad=np.array([float(z.rad()) if z is not None else np.nan for z in Z])
gb=float(gab(arb(T)).upper())
obs=np.abs(Zf-mid); ok=~np.isnan(mid)
viol=int((obs[ok]>(Bf[ok]-gb)).sum())
fb=np.where(sg==0)[0]; N=ra['N']; lnA=[arb(n).log() for n in range(1,N+1)]; unc=0; t2=time.time()
for j in fb:
    zz=Zdirect(arb(T)+arb(delta)*int(j),N,lnA)
    if zz is None or zz.contains(0): unc+=1
    else: sg[j]=1 if zz>0 else -1
tfb=time.time()-t2
t3=time.time(); [Zdirect(arb(T)+arb(delta)*int(j),N,lnA) for j in range(0,P,P//20)]; tdir=(time.time()-t3)/20
tu=turing(T,delta,sg,a,b,Lw)
res=dict(T=T,P=P,N=N,L=ra['L'],float_sec=rf['batch_sec'],arb_sec=ta,fallback_sec=tfb,arb_cost_vs_float=ta/rf['batch_sec'],
  certified_frac=1-len(fb)/P,uncertified=unc,referee_max_obs_err=float(np.nanmax(obs)),referee_min_slack_ratio=float(np.nanmin((Bf[ok]-gb)/np.maximum(obs[ok],1e-300))),
  referee_violations=viol,arb_median_rad=float(np.nanmedian(rad)),
  speedup_vs_arb_direct=tdir*P/(ta+tfb),turing=tu,label='NUM')
if T>1e11:
    z3=np.array([float(x) for x in open('zeros3').read().split('\n') if x.strip() and x.strip()[0].isdigit() and ' ' not in x.strip()])
    res['table_zeros_in_window']=int(((z3>a)&(z3<=b)).sum())
print(res,flush=True); json.dump(res,open(f'fftarb_{sys.argv[1]}.json','w'),indent=1,default=str)
