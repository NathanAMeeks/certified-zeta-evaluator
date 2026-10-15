# [NUM*] Odlyzko-Schonhage-style batched certified Z on a uniform grid.
# F(u_p)=sum_n a_n e^{-i u_p ln n} = nonuniform DFT; Taylor-NUFFT (L terms) + FFT, rigorous float error bound.
import numpy as np,math,time,json,sys
from flint import acb,arb,acb_series,ctx
ctx.prec=120; PI=arb.pi(); U=2.0**-53
exec(open('cert2.py').read().split('def run(')[0].split('ctx.prec=110; PI=arb.pi()')[1])  # theta_ser, C0, gab, Zdirect
def batch(T,u0,delta,P,Hth=1.0):
    t0=time.time(); Ta=arb(T); N=int((Ta/(2*PI)).sqrt().floor().unique_fmpz())
    base=Ta+arb(u0)
    lnA=[arb(n).log() for n in range(1,N+1)]
    # f_n = n^{-1/2} e^{-i base ln n}: phase reduced mod 2pi in Arb, then to float (|err|<=2u*2pi)
    ph=np.array([float(((base*l)/(2*PI)).frac_mid() if hasattr(arb,'frac_mid') else 0) for l in lnA]) if False else None
    tp=[]
    for l in lnA:
        x=base*l/(2*PI); x=x-x.floor(); tp.append(float(x.mid()))
    lnn=np.log(np.arange(1,N+1)); w=np.arange(1,N+1)**-0.5
    f=w*np.exp(-2j*np.pi*np.array(tp)); S0=w.sum()
    M=1<<int(math.ceil(math.log2(2*P)))
    x=(delta*lnn/(2*np.pi))%1.0                       # e^{-i p delta ln n}=e^{-2pi i p x_n}
    f=f*np.exp(-2j*np.pi*(P//2)*x)                    # shift p=p'+P/2, |p'|<=P/2
    mi=np.rint(x*M).astype(np.int64); e=x*M-mi; mi%=M  # |e|<=1/2
    r=(P/2)/M                                          # |2pi p' d| <= pi*r <= pi/4
    tol=1e-4*float(gab(Ta).mid()); L=1
    while S0*(math.pi*r)**L/math.factorial(L)*math.exp(math.pi*r)>tol: L+=1
    pp=np.arange(P)-P//2; F=np.zeros(P,complex); el=np.ones(N); G1=0.0
    for l in range(L):
        G=np.bincount(mi,weights=(f*el).real,minlength=M)+1j*np.bincount(mi,weights=(f*el).imag,minlength=M)
        Gh=np.fft.fft(G)[pp%M]
        F+=(-2j*np.pi*pp/M)**l/math.factorial(l)*Gh
        el=el*e
    trunc=S0*(math.pi*r)**L/math.factorial(L)*math.exp(math.pi*r)
    # float error bound (Higham): phase (tp,x,shift,e) ~ (8+2pi*P/2*... ) handled as relative err per term
    phase_err=(4*math.pi*U*(1+P//2)+8*U)      # rounding in e^{-2pi i tp}, e^{-2pi i (P/2) x}
    xerr=delta*lnn[-1]/(2*np.pi)*4*U*M        # error in e (scaled d) from x rounding
    fft_err=5*U*math.log2(M)*math.sqrt(M)
    rnd=S0*math.exp(math.pi*r)*(phase_err+(N+L+4)*U+fft_err+math.pi*r*xerr*2)+S0*2*math.pi*(P//2)*abs(delta*lnn[-1])*U*0
    mainerr=2*(trunc+rnd)
    # theta Taylor per block (Arb coeffs), evaluated in float
    u=u0+delta*np.arange(P); Z=np.empty(P); B=np.empty(P)
    nb=int(math.ceil(P*delta/(2*Hth)))
    for b in range(nb):
        lo=b*2*Hth/delta; idx=np.arange(int(math.ceil(lo)),min(P,int(math.ceil((b+1)*2*Hth/delta))))
        if len(idx)==0: continue
        cu=u0+(b*2+1)*Hth; c=Ta+arb(cu)
        thc,th3=theta_ser(c,Hth)
        t0c=thc[0]/(2*PI); t0c=float(((t0c-t0c.floor())*2*PI).mid())
        v=u[idx]-cu; th=t0c+float(thc[1].mid())*v+float(thc[2].mid())*v*v
        therr=float(th3.upper())*np.abs(v)**3+float(thc[1].rad()+thc[2].rad())*4+16*U*(7+abs(float(thc[1].mid()))*Hth)
        Z[idx]=2*np.real(np.exp(1j*th)*F[idx]); B[idx]=mainerr+2*S0*(therr+2*U)
    # C0 + Gabcke, float with bound
    t=T+u; a=np.sqrt(t/(2*np.pi)); p=a-N; den=np.cos(2*np.pi*p); bad=np.abs(den)<1e-3
    c0=(-1)**(N-1)*a**-0.5*np.cos(2*np.pi*(p*p-p-1/16))/np.where(bad,1,den)
    Z+=c0; B+=1e3*a**-0.5*64*U*(1+a*4*np.pi)+float(gab(Ta).upper())
    tt=time.time()-t0
    sg=np.where((np.abs(Z)>B)&~bad,np.sign(Z),0)
    fb=np.where(sg==0)[0]; t1=time.time(); unc=0
    for j in fb:
        rr=Zdirect(Ta+arb(float(u[j])),N,lnA)
        if rr is None or rr.contains(0): unc+=1
        else: sg[j]=1 if rr>0 else -1
    tf=time.time()-t1
    return dict(N=N,P=P,M=M,L=L,batch_sec=tt,fallback_sec=tf,certified_frac=1-len(fb)/P,fallbacks=len(fb),uncertified=unc,
        max_bound=float(B.max()),med_bound=float(np.median(B)),signchanges=int((np.diff(sg[sg!=0])!=0).sum())),u,sg
