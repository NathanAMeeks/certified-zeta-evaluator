# [LEM C0FIX] validation: psi_c0 enclosure vs mpmath (original quotient at 80 digits; limit at p0)
import mpmath as mp, subprocess, numpy as np
mp.mp.dps=110
def truth(p):
    p=mp.mpf(p); 
    for p0 in (mp.mpf(1)/4, mp.mpf(3)/4):
        if p==p0: return mp.limit(lambda q: mp.cos(2*mp.pi*(q*q-q-mp.mpf(1)/16))/mp.cos(2*mp.pi*q), p0)
    return mp.cos(2*mp.pi*(p*p-p-mp.mpf(1)/16))/mp.cos(2*mp.pi*p)
pts=[]
for p0 in (0.25,0.75):
    pts+=list(p0+np.linspace(-0.012,0.012,4001))
    pts+=[p0]+[p0+s*h for s in (1,-1) for h in (1e-16,1e-13,1e-10,1e-7,1e-5,1.6e-4,1e-3,0.01*(1-1e-12),0.01,0.01*(1+1e-12))]
    pts+=list(p0+np.linspace(-2e-4,2e-4,2001))
pts+=list(np.linspace(0.001,0.999,999))
inp='\n'.join(repr(float(x)) for x in pts)
tot=bad=0; wmax={0.25:0,0.75:0}; okfail=0
for prec in (96,256):
  for rad in (0.0,1e-12,1e-9):
    out=subprocess.run(['./t_psi',str(prec),str(rad)],input=inp,capture_output=True,text=True).stdout.split('\n')
    for line in out:
        if not line: continue
        x,ok,m,r=line.split(); x=float(x); ok=int(ok); r=mp.mpf(r); m=mp.mpf(m)
        if not ok: okfail+=1; continue
        for xx in ((x,) if rad==0 else (x-rad,x,x+rad)):
            v=truth(mp.mpf(x) if xx==x else mp.mpf(x)+(mp.mpf(rad) if xx>x else -mp.mpf(rad)))
            tot+=1
            if abs(v-m)>r: bad+=1; print('VIOLATION',prec,rad,x,m,r,v)
        for p0 in wmax:
            if abs(x-p0)<=0.01 and rad==0 and prec==96: wmax[p0]=max(wmax[p0],2*float(r))
print(f'checks={tot} violations={bad} not_ok={okfail} max_width_prec96_pointballs: p0=1/4 {wmax[0.25]:.3g}, p0=3/4 {wmax[0.75]:.3g}')
