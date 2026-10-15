import sys,math,json,numpy as np
src=open('cert2.py').read().split('z1=np.loadtxt')[0]; exec(src)
from flint import ctx
def Kfor(T,H):
    N=int(math.sqrt(T/(2*math.pi))); R=math.log(N)/2; S0=2*math.sqrt(N); tol=0.1*0.127*(T/(2*math.pi))**-0.75
    K=1
    while math.log(S0)+K*math.log(R*H)-math.lgamma(K+1)+R*H>math.log(tol): K+=1
    return N,K
def tune(T,density):   # cost/pt ~ N*K/(2H*density) + c*K  (c ~ per-point overhead in units of term-ops)
    best=None
    for H in np.exp(np.linspace(np.log(0.05),np.log(32),60)):
        N,K=Kfor(T,H); cost=(N*K+5500)/(2*H*density)+10*K
        if best is None or cost<best[0]: best=(cost,float(H),K)
    prec=int(64+math.log2(T)+math.log2(max(N,2)))
    return best[1],best[2],prec
out=[]
for T,span in [(1e4,800.),(1e5,600.),(1e6,500.),(1e8,200.),(267653395647,40.)]:
    dens=10000/span; H,K,prec=tune(T,dens)
    for HH in [H]:
        ctx.prec=prec if HH==H else 110
        r=run(T,span,10000,HH,None); r.update(tuned=(HH==H),prec=ctx.prec,density=dens); out.append(r)
json.dump(out,open('tune_results2.json','w'),indent=1)
