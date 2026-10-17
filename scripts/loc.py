import numpy as np,math
T=3e12;delta=0.01;P=1<<20;a,b=60,10380
sg=np.fromfile('sg0.bin',dtype=np.int8).astype(int)
u=np.arange(P)*delta
ch=np.where(np.diff(sg)!=0)[0]   # all certified, so sign change between p and p+1
zc=(u[ch]+u[ch+1])/2
zc=zc[(zc>a)&(zc<b)]
L=0.5*math.log(T/(2*math.pi))
th=lambda x:(x*L+x*x/(4*T))/math.pi   # (theta(T+x)-theta(T))/pi approx
grid=np.arange(a,b,0.05)
cnt=np.searchsorted(zc,grid)
d=cnt-(th(grid)-th(a))
# moving mean over 20 units
w=400; m=np.convolve(d,np.ones(w)/w,mode='valid'); g2=grid[w//2:w//2+len(m)]
print('d start/end mean',d[:2000].mean(),d[-2000:].mean())
# find step: best split maximizing mean difference
cs=np.cumsum(d); n=len(d); best=None
for k in range(2000,n-2000,50):
    l=cs[k]/k; r=(cs[-1]-cs[k])/(n-k); 
    if best is None or l-r>best[0]: best=(l-r,grid[k])
print('step',best)
# also local: gaps between consecutive zeros much larger than mean spacing
sp=np.diff(zc); mean=1/(L/math.pi)
idx=np.argsort(-sp)[:15]
for i in idx: print('gap',zc[i],zc[i+1],sp[i]/mean)
