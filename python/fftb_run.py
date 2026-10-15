exec(open('fftb.py').read())
tr={round(d['T']):d for d in json.load(open('tune_results2.json'))}
z3=np.array([float(x) for x in open('zeros3').read().split('\n') if x.strip() and x.strip()[0].isdigit() and ' ' not in x.strip()])
out=[]
for T,dens,P in [(1e6,20.,10000),(1e6,20.,131072),(1e8,50.,10000),(1e8,50.,131072),(267653395647,250.,10000),(267653395647,250.,131072)]:
    r,u,sg=batch(T,0.0,1/dens,P)
    c2=tr[round(T)]['taylor_sec']/10000*P      # cert2 auto_HK cost scaled linearly (same density)
    rng=np.random.default_rng(1); smp=rng.choice(P,100,replace=False); Ta=arb(T)
    N=r['N']; lnA=[arb(n).log() for n in range(1,N+1)]
    t1=time.time(); dis=0; chk=0
    for j in smp:
        z=Zdirect(Ta+arb(float(u[j])),N,lnA)
        if z is not None and not z.contains(0) and sg[j]!=0: chk+=1; dis+=int((1 if z>0 else -1)!=sg[j])
    tdir=(time.time()-t1)/100
    r.update(T=T,span=P/dens,cert2_sec_est=c2,speedup_vs_cert2=c2/(r['batch_sec']+r['fallback_sec']),
             arb_direct_sec_est=tdir*P,speedup_vs_direct=tdir*P/(r['batch_sec']+r['fallback_sec']),spot_checked=chk,spot_disagree=dis,label='NUM*')
    if T>1e11:
        uu=u[sg!=0]; r['table_zeros']=int(((z3>=uu[0])&(z3<=uu[-1])).sum()); r['zeros_below_table_start']=int(r['signchanges']-r['table_zeros'])
    print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()},flush=True); out.append(r)
json.dump(out,open('fftb_results.json','w'),indent=1)
