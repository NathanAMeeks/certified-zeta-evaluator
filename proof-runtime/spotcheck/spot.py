# [NUM] independent Arb (python-flint) spot-check of v9p ledger; Nathan A. Meeks & Grok. candidate, unrefereed.
import csv,random,lzma,glob,json,math
from flint import arb,acb,ctx
ctx.prec=128; random.seed(20261030)
def theta(t): return acb(arb('0.25'),t/2).lgamma().imag - t/2*arb.pi().log()
def Z(t): return (acb(0,theta(t)).exp()*acb(arb('0.5'),t).zeta()).real
rows={}
for r in csv.DictReader(open('../ledger.csv')):
    if r['status']=='PROVEN': rows[int(r['idx'])]=r
ids=sorted(rows); samp=sorted(random.sample(ids,40))
res=dict(windows=0,count_ok=0,turing_ok=0,ends=0,ends_ok=0,ends_unc=0,bad=[])
for i in samp:
    r=rows[i]; T=int(r['T_arg']); lo,hi=int(r['lo']),int(r['hi']); Na,Nb=int(r['Na']),int(r['Nb'])
    f=glob.glob(f'../brackets/w{i:07d}_T*.txt.xz')
    if not f: res['bad'].append([i,'nofile']); continue
    br=[tuple(l.split()) for l in lzma.open(f[0],'rt') if l[0]!='#']
    res['windows']+=1
    inn=[b for b in br if T+float(b[0])>=lo and T+float(b[1])<=hi]
    if len(inn)==Nb-Na==int(r['cert_sign_changes']): res['count_ok']+=1
    else: res['bad'].append([i,'count',len(inn),Nb-Na])
    S=0.112*math.log(hi)+0.278*math.log(math.log(hi))+2.51
    ok=all(abs(N-(theta(arb(x))/arb.pi()+1)).upper()<S for N,x in ((Na,lo),(Nb,hi)))
    res['turing_ok']+=ok
    if not ok: res['bad'].append([i,'turing'])
    for a,b in random.sample(br,min(5,len(br))):
        za,zb=Z(T+arb(a)),Z(T+arb(b)); res['ends']+=2
        if za.contains(0) or zb.contains(0): res['ends_unc']+=1; res['bad'].append([i,'unc',a,b]); continue
        if (za>0)!=(zb>0): res['ends_ok']+=2
        else: res['bad'].append([i,'samesign',a,b])
res['sample']=samp; json.dump(res,open('spot.json','w'))
print({k:v for k,v in res.items() if k!='sample'})
