#define _GNU_SOURCE
/* [NUM] rsbatch: certified batch Riemann-Siegel evaluator for the G-RS referee stage. Nathan A. Meeks & Grok. Candidate, unrefereed; not RH.
   Independent of G-Turing (which uses Arb acb_zeta) and of the finder (fftzi). Reads nothing from the library.
   Z(t) = 2 Re(e^{i theta(t)} sum_{n<=N} n^{-1/2} n^{-it}) + (-1)^{N-1} tau^{-1/4} C0(p) + R,  |R| <= 0.127 tau^{-3/4} (Gabcke 1979, t>=200),
   tau=t/2pi, N=floor(sqrt tau), p=frac(sqrt tau).  Batch trick: n^{-it} = c_n * e^{-i v log n}, c_n = n^{-1/2-i t0} (Arb, once per block),
   t = t0 + v, |v| small.  Inner sum in double with an explicit a-priori error bound (see ERR below); everything else in Arb.
   AUDIT FIX: libm sincos replaced by rigorous rsc() (rsc.h), no libm ulp assumption.
   stdin: T c0 (integers), then lines "u" (decimal offsets from T). stdout per point: sign mid rad   (sign 0 = not certified -> caller falls back). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rsc.h"
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#define PR 192
static void arb_sub_d(arb_t z,const arb_t a,double d,slong p){arb_t q;arb_init(q);arb_set_d(q,d);arb_sub(z,a,q,p);arb_clear(q);}
static double *cr,*ci,*ln_,*ce,*ca;
int main(int argc,char**argv){
  int nth=argc>1?atoi(argv[1]):3; char Ts[64],cs[64]; if(scanf("%63s %63s",Ts,cs)!=2) return 1;
  int cap=1<<16,np=0; char **us=malloc(cap*sizeof(char*)); char buf[128];
  while(scanf("%127s",buf)==1){ if(np==cap){cap*=2; us=realloc(us,cap*sizeof(char*)); if(!us) return 2;} us[np++]=strdup(buf); } /* t100 audit: fixed fixed-cap overflow */
  arb_t t0,tp,x,y; arb_init(t0);arb_init(tp);arb_init(x);arb_init(y);
  arb_set_str(t0,Ts,PR); arb_set_str(x,cs,PR); arb_add(t0,t0,x,PR); arb_const_pi(tp,PR); arb_mul_2exp_si(tp,tp,1);
  /* Nmax from largest t: use t0+1000 (points must lie within |v|<1000) */
  arb_add_ui(x,t0,1000,PR); arb_div(x,x,tp,PR); arb_sqrt(x,x,PR); fmpz_t f; fmpz_init(f); arb_floor(x,x,PR); arb_get_unique_fmpz(f,x);
  long Nmax=fmpz_get_si(f)+1;
  cr=malloc(Nmax*8);ci=malloc(Nmax*8);ln_=malloc(Nmax*8);ce=malloc(Nmax*8);ca=malloc(Nmax*8);
  #pragma omp parallel num_threads(nth)
  { arb_t L,a,s,c,r,e; arb_init(L);arb_init(a);arb_init(s);arb_init(c);arb_init(r);arb_init(e); mag_t m; mag_init(m);
    #pragma omp for schedule(static)
    for(long n=1;n<Nmax;n++){
      arb_log_ui(L,n,PR); arb_mul(a,L,t0,PR); arb_neg(a,a); arb_sin_cos(s,c,a,PR); arb_rsqrt_ui(r,n,PR); arb_mul(s,s,r,PR); arb_mul(c,c,r,PR);
      double dr=arf_get_d(arb_midref(c),ARF_RND_NEAR), di=arf_get_d(arb_midref(s),ARF_RND_NEAR);
      arb_sub_d(e,c,dr,PR); arb_get_mag(m,e); double er=mag_get_d(m); arb_sub_d(e,s,di,PR); arb_get_mag(m,e); double ei=mag_get_d(m);
      cr[n]=dr;ci[n]=di; ce[n]=(er+ei)*1.0001; ca[n]=(fabs(dr)+fabs(di))*1.0001;   /* |c-c_d| <= ce, |re|+|im| <= ca */
      ln_[n]=arf_get_d(arb_midref(L),ARF_RND_NEAR); }
    arb_clear(L);arb_clear(a);arb_clear(s);arb_clear(c);arb_clear(r);arb_clear(e); mag_clear(m); }
  int *sg=malloc(np*sizeof(int)); double *md=malloc(np*8),*rd=malloc(np*8);
  const double u53=ldexp(1,-53);
  #pragma omp parallel num_threads(nth)
  { arb_t t,v,w,tau,sq,p,th,S,Zr,g,h; acb_t z,lg; fmpz_t N; fmpz_init(N);
    arb_init(t);arb_init(v);arb_init(w);arb_init(tau);arb_init(sq);arb_init(p);arb_init(th);arb_init(S);arb_init(Zr);arb_init(g);arb_init(h); acb_init(z);acb_init(lg);
    mag_t m; mag_init(m);
    #pragma omp for schedule(dynamic,1)
    for(int k=0;k<np;k++){
      sg[k]=0; md[k]=0; rd[k]=INFINITY;
      arb_set_str(v,us[k],PR); arb_set_str(w,cs,PR); arb_sub(v,v,w,PR);            /* v = u - c exact-ish ball */
      arb_add(t,t0,v,PR); arb_div(tau,t,tp,PR); arb_sqrt(sq,tau,PR); arb_floor(p,sq,PR);
      if(!arb_get_unique_fmpz(N,p)) continue; long Nn=fmpz_get_si(N); if(Nn>=Nmax||fabs(arf_get_d(arb_midref(v),ARF_RND_NEAR))>1000) continue;
      double vd=arf_get_d(arb_midref(v),ARF_RND_NEAR); arb_sub_d(w,v,vd,PR); arb_get_mag(m,w); double ev=mag_get_d(m);
      /* inner sum, blocked (512) recursive summation */
      double Sr=0,Si=0,Ab=0,Eb=0; int okk=1; long B=512, nb=0;
      for(long n0=1;n0<=Nn;n0+=B){ double br=0,bi=0; long n1=n0+B-1<Nn?n0+B-1:Nn; nb++;
        for(long n=n0;n<=n1;n++){ double ph=vd*ln_[n],s,c; if(!rsc(ph,&s,&c)) okk=0;
          br+= cr[n]*c + ci[n]*s; bi+= ci[n]*c - cr[n]*s;          /* c_n * e^{-i ph} */
          double dph=ln_[n]*(ev+3.0*fabs(vd)*u53)*1.01;              /* |phi - fl(vd*l_d)| (l_d err <= ulp/2) */
          Eb+= ce[n] + ca[n]*(dph+ 2*RSC_E + 6*u53);                   /* coeff err + phase err + sincos (4ulp each) + products */
          Ab+= ca[n]; }
        Sr+=br; Si+=bi; }
      if(!okk) continue; double ERR = (Eb + (double)(B+nb+4)*u53*1.01*Ab*2.0)*1.001;    /* ERR: per-term + recursive-sum gamma bounds, both components */
      acb_set_d_d(z,Sr,Si); mag_set_d(m,ERR); arb_add_error_mag(acb_realref(z),m); arb_add_error_mag(acb_imagref(z),m);
      /* theta = Im lgamma(1/4 + i t/2) - t/2 log pi */
      acb_set_d(lg,0.25); arb_mul_2exp_si(acb_imagref(lg),t,-1); acb_lgamma(lg,lg,PR); arb_const_pi(g,PR); arb_log(g,g,PR); arb_mul(g,g,t,PR); arb_mul_2exp_si(g,g,-1);
      arb_sub(th,acb_imagref(lg),g,PR);
      arb_sin_cos(g,h,th,PR); /* Re(e^{i th} z) = cos*zr - sin*zi */
      arb_mul(Zr,h,acb_realref(z),PR); arb_mul(g,g,acb_imagref(z),PR); arb_sub(Zr,Zr,g,PR); arb_mul_2exp_si(Zr,Zr,1);
      /* C0 term */
      arb_sub(p,sq,p,PR); arb_const_pi(g,PR); arb_mul_2exp_si(g,g,1); arb_mul(h,g,p,PR); arb_cos(h,h,PR);
      if(arf_cmpabs_2exp_si(arb_midref(h),-10)<0) continue;          /* |cos 2pi p| < ~1e-3: refuse (fallback) */
      arb_mul(S,p,p,PR); arb_sub(S,S,p,PR); arb_set_d(w,0.0625); arb_sub(S,S,w,PR); arb_mul(S,S,g,PR); arb_cos(S,S,PR); arb_div(S,S,h,PR);
      arb_root_ui(w,tau,4,PR); arb_div(S,S,w,PR); if(Nn%2==0) arb_neg(S,S); arb_add(Zr,Zr,S,PR);
      /* Gabcke |R0| <= 0.127 tau^{-3/4} */
      arb_set_d(g,-0.75); arb_pow(g,tau,g,PR); arb_set_ui(w,127); arb_div_ui(w,w,1000,PR); arb_mul(g,g,w,PR); arb_get_mag(m,g); arb_add_error_mag(Zr,m);
      md[k]=arf_get_d(arb_midref(Zr),ARF_RND_NEAR); rd[k]=mag_get_d(arb_radref(Zr));
      if(!arb_contains_zero(Zr)) sg[k]=arb_is_positive(Zr)?1:-1; }
  }
  for(int k=0;k<np;k++) printf("%d %.17g %.3g\n",sg[k],md[k],rd[k]);
  return 0; }
