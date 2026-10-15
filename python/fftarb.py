# [NUM] fully Arb NUFFT (acb.dft) + Turing/Trudgian completeness
exec(open('fftb.py').read())
from flint import acb,arb,ctx
def batch_arb(T,delta,P,Hth=1.0):
    ctx.prec=96; t0=time.time(); Ta=arb(T); N=int((Ta/(2*PI)).sqrt().floor().unique_fmpz())
    M=1<<int(math.ceil(math.log2(2*P))); d=arb(delta); half=P//2; TWO_PI=2*PI
    G=[[acb(0)]*M for _ in range(1)]  # placeholder
    S0=arb(0); fs=[];ms=[];es=[]
    for n in range(1,N+1):
        l=arb(n).log(); w=arb(n).rsqrt(); S0+=w
        x=d*l/TWO_PI; x=x-x.floor()
        # f_n = n^{-1/2} e^{-i T ln n} e^{-2 pi i half x}
        ph=Ta*l/TWO_PI+half*x; ph=ph-ph.floor()
        fs.append(w*acb(0,-TWO_PI*ph).exp())
        mm=int(round(float((x*M).mid()))); ms.append(mm%M); es.append(x*M-mm)
    r=arb(half)/M; L=1
    while float((S0*(PI*r)**L/arb(math.factorial(L))*(PI*r).exp()).upper())>1e-4*float(gab(Ta).mid()): L+=1
    trunc=S0*(PI*r)**L/arb(math.factorial(L))*(PI*r).exp()
    pp=list(range(-half,P-half)); F=[acb(0)]*P
    cur=list(fs)
    for l in range(L):
        G=[acb(0)]*M
        for n in range(N): G[ms[n]]+=cur[n]
        Gh=acb.dft(G)  # sum_m G[m] e^{-2pi i q m/M}
        if l==0: sgn=None
        for i,q in enumerate(pp):
            F[i]+=(acb(0,-2*PI*q/M)**l)/math.factorial(l)*Gh[q%M]
        cur=[cur[n]*es[n] for n in range(N)]
    tr=arb(0,trunc.upper())
    Z=[None]*P; nb=int(math.ceil(P*delta/(2*Hth)))
    for b in range(nb):
        lo=int(math.ceil(b*2*Hth/delta)); hi=min(P,int(math.ceil((b+1)*2*Hth/delta)))
        if lo>=hi: continue
        cu=arb((b*2+1)*Hth); c=Ta+cu; thc,th3=theta_ser(c,Hth)
        for p in range(lo,hi):
            v=d*p-cu; th=thc[0]+thc[1]*v+thc[2]*v*v; th=arb(th.mid(),th.rad()+(th3*abs(v)**3).upper())
            Fp=F[p]+acb(tr,tr)
            z=2*(acb(0,th).exp()*Fp).real; cc=C0(c+v,N)
            Z[p]=None if cc is None else arb((z+cc).mid(),(z+cc).rad()+gab(c+v).upper())
    return dict(N=N,M=M,L=L,sec=time.time()-t0),Z
def turing(T,delta,sg,a,b,Lw):
    # Trudgian 2011: |int_{t1}^{t2} S| <= 2.067+0.059 log t2 (t2>t1>168 pi). Rigorous N bounds via integral form.
    ctx.prec=96; u=np.arange(len(sg))*delta; Ta=arb(T)
    idx=np.where(sg!=0)[0]; s=sg[idx]; ch=np.where(np.diff(s)!=0)[0]
    lo_e=[(u[idx[j]],u[idx[j+1]]) for j in ch]   # zero in (s_j,e_j)
    def thint(x0,x1):  # int theta(T+x)/pi dx over [x0,x1], Taylor about mid + remainder
        h=(x1-x0)/2; c=Ta+arb((x0+x1)/2); co,b3=theta_ser(c,h)
        I=2*co[0]*h+2*co[2]*h**3/3; R=b3*2*arb(h)**4/4
        return (I+arb(0,R.upper()))/PI
    C=lambda t2: arb(2.067)+arb(0.059)*(Ta+arb(t2)).log()
    # upper bound on N(T+b): int_b^{b+L} [N(b)+n_b(t)-theta/pi-1] <= C
    I_n=sum(arb(b+Lw-e) for s_,e in lo_e if s_>=b and e<=b+Lw)
    Ub=1+(C(b+Lw)-I_n+thint(b,b+Lw))/Lw
    I_m=sum(arb(s_-(a-Lw)) for s_,e in lo_e if s_>=a-Lw and e<=a)
    La=1+(-C(a)+I_m+thint(a-Lw,a))/Lw
    found=sum(1 for s_,e in lo_e if s_>=a and e<=b)
    up=int(Ub.upper().floor().unique_fmpz()) if hasattr(Ub.upper(),'floor') else math.floor(float(Ub.upper()))
    lw=math.ceil(float(La.lower()))
    return dict(a=T+a,b=T+b,L=Lw,found=found,N_b_upper=float(Ub.upper()),N_a_lower=float(La.lower()),
       max_possible=math.floor(float(Ub.upper()))-lw,proven=(math.floor(float(Ub.upper()))-lw==found),
       Na_int=lw, Nb_int=math.floor(float(Ub.upper())))
