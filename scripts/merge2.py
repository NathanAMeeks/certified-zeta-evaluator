import numpy as np
P=1<<20; a,b=60.,10380.
s0=np.fromfile('sg0.bin',dtype=np.int8).astype(int); u0=np.arange(P)*0.01
s1=np.fromfile('sg2.bin',dtype=np.int8).astype(int); u1=5476+np.arange(len(s1))*1e-3
# consistency on coincident points
k=np.arange(547600,547600+205); j=(k-547600)*10
print('overlap pts',len(k),'disagree',int((s0[k]!=s1[j]).sum()))
u=np.concatenate([u0,u1]); s=np.concatenate([s0,s1]); o=np.argsort(u,kind='stable'); u,s=u[o],s[o]
m=s!=0; u,s=u[m],s[m]
ch=np.where(np.diff(s)!=0)[0]; zs,ze=u[ch],u[ch+1]
found=int(((zs>=a)&(ze<=b)).sum())
print('found_refined',found,'max_possible',12362403064087-12362403019918, 'proven',found==12362403064087-12362403019918)
print('extra certified points',len(s1),'= %.2f%% of window'%(100*len(s1)/P))
