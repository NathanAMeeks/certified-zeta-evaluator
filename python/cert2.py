# [NUM] fully-Arb certified Taylor-reuse Z(t); only external constant: Gabcke 0.127
import numpy as np,time,json,gzip,math
from flint import acb,arb,acb_series,ctx
ctx.prec=110; PI=arb.pi()
def gab(t): return arb(0.127)*t**arb(-0.75)   # Gabcke 1979: |R0|<=0.127 t^{-3/4}, t>=200
def theta_ser(c,H,K=4):   # theta(c+v) coeffs on point c, and ball enclosure of 3rd coeff on [c-H,c+H]
    ctx.cap=K
    f=lambda t:(acb_series([acb(0.25,t/2),acb(0,0.5)]).lgamma())
    s=f(c).coeffs(); sb=f(arb(c,H)).coeffs()
    co=[x.imag for x in s]; co[0]-=c/2*PI.log(); co[1]-=PI.log()/2
    return co[:3], abs(sb[3].imag)   # |theta'''/6| bound on ball
def C0(t,N):
    a=(t/(2*PI)).sqrt(); p=a-N; den=(2*PI*p).cos()
    if abs(den)<arb(1e-3): return None
    return (-1)**(N-1)*a**arb(-0.5)*(2*PI*(p*p-p-arb(1)/16)).cos()/den
def Zdirect(t,N,lnA):  # Arb direct RS + Gabcke, returns ball
    th=theta_ser(t,0)[0][0]; s=acb(0)
    for n in range(1,N+1): s+=acb(0,th-t*lnA[n-1]).exp()/arb(n).sqrt()
    c=C0(t,N)
    if c is None: return None
    v=2*s.real+c; return arb(v.mid(),v.rad()+gab(t).upper())
def run(T,span,npts,H,tab):
    T=arb(T); N=int((T/(2*PI)).sqrt().floor().unique_fmpz()); lnA=[arb(n).log() for n in range(1,N+1)]
    S0=sum(arb(n).rsqrt() for n in range(1,N+1))
    u=np.linspace(0,span,npts); out=[None]*npts
    cents=np.arange(H,span+2*H,2*H); t0=time.time()
    lc=(lnA[0]+lnA[-1])/2; w=[x-lc for x in lnA]; R=max(abs(x) for x in w)
    K=1
    while S0*(R*H)**K/math.factorial(K)*(R*H).exp()>gab(T)*0.1: K+=1
    for uc in cents:
        idx=[i for i in range(npts) if abs(u[i]-uc)<=H*(1+1e-12) and out[i] is None]
        if not idx: continue
        c=T+arb(float(uc))
        P=[acb(0)]*K
        for n in range(N):
            b=acb(0,-c*lnA[n]).exp()*arb(n+1).rsqrt(); wn=w[n]
            for k in range(K): P[k]+=b; b*=wn
        thc,th3=theta_ser(c,H)
        for i in idx:
            v=arb(float(u[i]))-arb(float(uc)); av=abs(v)
            S=acb(0); fac=acb(1)
            for k in range(K): S+=P[k]*fac; fac=fac*acb(0,-v)/(k+1)
            trunc=S0*(R*av)**K/math.factorial(K)*(R*av).exp()
            S=S+acb(arb(0,trunc.upper()),arb(0,trunc.upper()))
            th=thc[0]+thc[1]*v+thc[2]*v*v; th=arb(th.mid(),th.rad()+(th3*av**3).upper())
            z=2*(acb(0,th-v*lc).exp()*S).real
            cc=C0(c+v,N)
            out[i]=None if cc is None else arb((z+cc).mid(),(z+cc).rad()+gab(c+v).upper())
    tt=time.time()-t0
    sg=np.array([0 if (x is None or x.contains(0)) else (1 if x>0 else -1) for x in out]); fb=np.where(sg==0)[0]
    # certified direct cost per point (timed on 3 points) and actual fallbacks
    t1=time.time(); [Zdirect(T+arb(float(u[j])),N,lnA) for j in range(3)]; td=(time.time()-t1)/3
    t1=time.time(); unc=0
    for j in fb:
        r=Zdirect(T+arb(float(u[j])),N,lnA)
        if r is None or r.contains(0): unc+=1
        else: sg[j]=1 if r>0 else -1
    tf=time.time()-t1
    s=sg[sg!=0]; uu=u[sg!=0]; sc=int((np.diff(s)!=0).sum())
    r=dict(T=float(T.mid()),H=H,K=K,N=N,points=npts,certified_frac=1-len(fb)/npts,fallbacks=len(fb),uncertified=unc,
      taylor_sec=tt,fallback_sec=tf,arb_direct_sec_per_pt=td,speedup_vs_arb_direct=npts*td/(tt+tf),
      certified_signchanges=sc,label='NUM')
    if tab is not None: r['table_zeros']=int(((tab>=uu[0])&(tab<=uu[-1])).sum())
    print(r,flush=True); return r

def auto_HK(T,density):  # tuned rule: min (N*K+5500)/(2H*density)+10K s.t. Taylor rem <= 0.1*Gabcke; prec=64+log2(T)+log2(N)
    import math
    N=int(math.sqrt(T/(2*math.pi))); R=math.log(N)/2; S0=2*math.sqrt(N); lt=math.log(0.1*0.127*T**-0.75)
    best=None
    for i in range(60):
        H=0.05*(640**(i/59)); K=1
        while math.log(S0)+K*math.log(R*H)-math.lgamma(K+1)+R*H>lt: K+=1
        c=(N*K+5500)/(2*H*density)+10*K
        if best is None or c<best[0]: best=(c,H,K)
    return best[1],best[2],int(64+math.log2(T)+math.log2(max(N,2)))
z1=np.loadtxt(gzip.open('zeros1.gz'))
z3=np.array([float(x) for x in open('zeros3').read().split('\n') if x.strip() and x.strip()[0].isdigit() and ' ' not in x.strip()])
import sys; which=sys.argv[1]
cfg={'a':(1e4,800.,10000,2,z1-1e4),'b':(1e6,500.,10000,2,None),'c':(267653395647,40.,10000,2,z3)}[which]
json.dump(run(*cfg),open(f'cert2_{which}.json','w'),indent=1)
