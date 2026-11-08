# G-Turing-P (Nathan A. Meeks & Grok) [NUM candidate]. Certified N(a) by Turing/Trudgian, mesh signs from FLINT Platt multieval (pme.c).
# Independent of rsbatch.c (G-RS) and of acb_zeta. Mesh is search-only; N accepted only if certified bounds coincide, else None -> Arb fallback.
import subprocess, math, sys, time
from flint import arb, ctx
TB=(arb('2.067'),arb('0.059'))
def theta_int(a,b):
    a,b=arb(a),arb(b); tp=2*arb.pi()
    F=lambda t: t*t/4*(t/tp).log()-3*t*t/8-arb.pi()*t/8+t.log()/48
    return F(b)-F(a)+arb(0,((b-a)/(a*a*a)).upper())
def signs(a,h):
    o=subprocess.run(['./pme',str(a),str(h),'1','128','16','4096'],capture_output=True,text=True,check=True)
    return [(float(l.split()[0]),int(l.split()[1]),float(l.split()[2])) for l in o.stdout.split('\n') if l.strip()]
def arbZ(a,u):
    from flint import acb
    from decimal import Decimal as D, getcontext
    getcontext().prec=60; t0=str(D(a)+D(repr(u)))
    for p in (128,256,512):
        old=ctx.prec; ctx.prec=p
        t=arb(t0); th=acb(arb('0.25'),t/2).lgamma().imag-t/2*arb.pi().log()
        z=(acb(0,th).exp()*acb(arb('0.5'),t).zeta()).real; ctx.prec=old
        if not z.contains(0): return (1 if z>0 else -1, abs(float(z.mid())))
    return (0,0.0)
def N(a,h=12,rounds=4):
    ctx.prec=256; a=int(a); S=signs(a,h)
    if any(x[1]==0 for x in S): return None,None,None,len(S)   # any uncertified mesh point -> refuse (conservative)
    V={u:(s,v) for u,s,v in S}; ev=0
    for _ in range(rounds):   # search-only refinement at same-sign |value| dips; new points' signs certified by Arb acb_zeta (not rsbatch)
        P=sorted(V); add=[]
        for k in range(1,len(P)-1):
            (s0,v0),(s1,v1),(s2,v2)=V[P[k-1]],V[P[k]],V[P[k+1]]
            if s0==s1==s2 and abs(v1)<abs(v0) and abs(v1)<abs(v2): add+=[(P[k-1]+P[k])/2,(P[k]+P[k+1])/2]
        add=[x for x in dict.fromkeys(add) if x not in V]
        if not add: break
        for x in add:
            z=arbZ(a,x); V[x]=(z,abs(zv)) if False else (z[0],z[1]); ev+=1
        if any(V[x][0]==0 for x in add): return None,None,None,len(S)+ev
    P=sorted(V); S=[(u,)+V[u] for u in P]
    cells=[(S[k][0],S[k+1][0]) for k in range(len(S)-1) if S[k][1]!=S[k+1][1]]
    A=arb(a);H=arb(h)
    intL=sum(h-r for l,r in cells if l>=0); intR=sum(l+h for l,r in cells if r<=0)
    up=(TB[0]+TB[1]*(A+H).log()+theta_int(A,A+H)/arb.pi()+H-arb(f'{intL:.9f}')-arb(0,1e-8))/H
    lo=(theta_int(A-H,A)/arb.pi()+H+arb(f'{intR:.9f}')-arb(0,1e-8)-TB[0]-TB[1]*A.log())/H
    U=math.floor(float(up.upper())); L=math.ceil(float(lo.lower()))
    ctx.prec=256
    return (U if U==L else None),L,U,len(S)
if __name__=='__main__':
    a=int(sys.argv[1]); s=time.time(); r=N(a); print(a,r,round(time.time()-s,1),'s',flush=True)
