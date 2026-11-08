# [NUM] Arb re-certification of fixed windows 8758,10202 (Nathan A. Meeks & Grok); candidate, unrefereed
import lzma,glob,random,json
from flint import arb,acb,ctx
ctx.prec=128; random.seed(7)
def Z(t):
    th=acb(arb('0.25'),t/2).lgamma().imag - t/2*arb.pi().log()
    return (acb(0,th).exp()*acb(arb('0.5'),t).zeta()).real
out={}
for i,c,T,lo,hi,need in ((8758,8711.633,3000090775200,3000090775280,3000090785600,44171),(10202,7552.83,3000105619520,3000105619600,3000105629920,44169)):
    br=[tuple(l.split()) for l in lzma.open(glob.glob(f'../brackets/w{i:07d}_*')[0],'rt') if l[0]!='#']
    inn=[b for b in br if T+float(b[0])>=lo and T+float(b[1])<=hi]
    narrow=[b for b in br if abs(float(b[0])-c)<2]
    print(i,len(narrow),flush=True); pick=list(dict.fromkeys(narrow+random.sample(br,100)))
    ok=unc=bad=0
    for a,b in pick:
        za,zb=Z(T+arb(a)),Z(T+arb(b))
        if za.contains(0) or zb.contains(0): unc+=1
        elif (za>0)!=(zb>0): ok+=1
        else: bad+=1
    out[i]=dict(in_window=len(inn),need=need,checked=len(pick),narrow=len(narrow),ok=ok,unc=unc,bad=bad)
json.dump(out,open('fixwin.json','w')); print(out)
