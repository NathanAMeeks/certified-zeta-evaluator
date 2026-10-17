/* [NUM] batched certified Z(t) on uniform grid u_p=p*delta, p<P, t=T+u (Taylor-NUFFT + acb_dft) + Turing/Trudgian.
   usage: fftz T delta P a b Lw [signfile]   (T integer-valued double ok) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/acb_dft.h>
#include <complex.h>
#undef I
static int USE_DBL=0;
/* iterative radix-2 DIT FFT, sign -1 */
static void fft_d(double complex*a,slong M,const double complex*W){
  for(slong i=1,j=0;i<M;i++){slong bit=M>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j){double complex t=a[i];a[i]=a[j];a[j]=t;}}
  for(slong len=2;len<=M;len<<=1){slong st=M/len;for(slong i=0;i<M;i+=len)for(slong k=0;k<len/2;k++){double complex u=a[i+k],v=a[i+k+len/2]*W[k*st];a[i+k]=u+v;a[i+k+len/2]=u-v;}}}

#define PREC 96
static void arb_sub_d(arb_t r,const arb_t a,double x,slong p){arb_t t;arb_init(t);arb_set_d(t,x);arb_sub(r,a,t,p);arb_clear(t);}
static void arb_add_d(arb_t r,const arb_t a,double x,slong p){arb_t t;arb_init(t);arb_set_d(t,x);arb_add(r,a,t,p);arb_clear(t);}
static void arb_mul_d(arb_t r,const arb_t a,double x,slong p){arb_t t;arb_init(t);arb_set_d(t,x);arb_mul(r,a,t,p);arb_clear(t);}
static void arb_div_d(arb_t r,const arb_t a,double x,slong p){arb_t t;arb_init(t);arb_set_d(t,x);arb_div(r,a,t,p);arb_clear(t);}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
/* theta series coeffs at c (radius rad): co[0..3] real, theta(c+v)=sum co[k] v^k */
static void theta_ser(arb_t *co,const arb_t c,double rad,slong prec){
  acb_poly_t h,r; acb_t z; arb_t t,lp; acb_poly_init(h);acb_poly_init(r);acb_init(z);arb_init(t);arb_init(lp);
  arb_set(t,c); if(rad>0) arb_add_error_2exp_si(t,0), mag_set_d(arb_radref(t),rad);
  arb_mul_2exp_si(t,t,-1); acb_set_d(z,0.25); arb_set(acb_imagref(z),t);
  acb_poly_set_coeff_acb(h,0,z); acb_set_d_d(z,0,0.5); acb_poly_set_coeff_acb(h,1,z);
  acb_poly_lgamma_series(r,h,4,prec);
  arb_const_pi(lp,prec); arb_log(lp,lp,prec);
  for(int k=0;k<4;k++){acb_poly_get_coeff_acb(z,r,k); arb_set(co[k],acb_imagref(z));}
  arb_set(t,c); if(rad>0) mag_set_d(arb_radref(t),rad);
  arb_mul(t,t,lp,prec); arb_mul_2exp_si(t,t,-1); arb_sub(co[0],co[0],t,prec);
  arb_mul_2exp_si(t,lp,-1); arb_sub(co[1],co[1],t,prec);
  acb_poly_clear(h);acb_poly_clear(r);acb_clear(z);arb_clear(t);arb_clear(lp);
}
static void gabcke(arb_t g,const arb_t t,slong prec){ arb_t x; arb_init(x); arb_const_pi(x,prec); arb_mul_2exp_si(x,x,1);
  arb_div(x,t,x,prec); arb_set_d(g,-0.75); arb_pow(x,x,g,prec); arb_set_d(g,0.127); arb_mul(g,g,x,prec); arb_clear(x);}
/* C0 term; returns 0 if |cos 2pi p|<1e-3 */
static int c0term(arb_t res,const arb_t t,slong N,slong prec){
  arb_t a,p,pi2,num,den,x; arb_init(a);arb_init(p);arb_init(pi2);arb_init(num);arb_init(den);arb_init(x);
  arb_const_pi(pi2,prec); arb_mul_2exp_si(pi2,pi2,1); arb_div(a,t,pi2,prec); arb_sqrt(a,a,prec); arb_sub_si(p,a,N,prec);
  arb_mul(den,pi2,p,prec); arb_cos(den,den,prec); int ok=1;
  arb_abs(x,den); if(!(arf_cmpabs_d(arb_midref(x),1e-3)>0 && !arb_contains_zero(den))) ok=0;
  arb_mul(x,p,p,prec); arb_sub(x,x,p,prec); arb_set_d(num,1.0/16); arb_sub(x,x,num,prec); arb_mul(x,x,pi2,prec); arb_cos(num,x,prec);
  arb_div(res,num,den,prec); arb_rsqrt(x,a,prec); arb_mul(res,res,x,prec); if((N-1)&1) arb_neg(res,res);
  arb_clear(a);arb_clear(p);arb_clear(pi2);arb_clear(num);arb_clear(den);arb_clear(x); return ok;}
/* certified signs of Z(T+p*delta), p<P, into sg (0 = not certified) */
static void zgrid(double Td,double delta,slong P,signed char*sg,slong*Nout,int*Lout,slong*ncout,double*tm){
  slong prec=PREC; double t0=now();
  arb_t T,d,pi,tp,l,w,x,ph,tmp,S0,g; arb_init(T);arb_init(d);arb_init(pi);arb_init(tp);arb_init(l);arb_init(w);arb_init(x);arb_init(ph);arb_init(tmp);arb_init(S0);arb_init(g);
  arb_set_d(T,Td); arb_set_d(d,delta); arb_const_pi(pi,prec); arb_mul_2exp_si(tp,pi,1);
  slong N=(slong)floor(sqrt(Td/(2*M_PI))); slong half=P/2, M=1; while(M<2*P) M<<=1;
  acb_ptr cur=_acb_vec_init(N), G=_acb_vec_init(M), Gh=_acb_vec_init(M), F=_acb_vec_init(P), fac=_acb_vec_init(P);
  arb_ptr e=_arb_vec_init(N); slong *mi=malloc(N*sizeof(slong));
  acb_t z; acb_init(z); arb_zero(S0); double complex *a=NULL,*W=NULL;
  for(slong n=1;n<=N;n++){
    arb_log_ui(l,n,prec); arb_set_ui(w,n); arb_rsqrt(w,w,prec); arb_add(S0,S0,w,prec);
    arb_mul(x,d,l,prec); arb_div(x,x,tp,prec); arb_floor(tmp,x,prec); arb_sub(x,x,tmp,prec);
    arb_mul(ph,T,l,prec); arb_div(ph,ph,tp,prec); arb_mul_si(tmp,x,half,prec); arb_add(ph,ph,tmp,prec);
    arb_floor(tmp,ph,prec); arb_sub(ph,ph,tmp,prec); arb_mul(ph,ph,tp,prec); arb_neg(ph,ph);
    arb_sin_cos(acb_imagref(z),acb_realref(z),ph,prec); acb_mul_arb(cur+n-1,z,w,prec);
    arb_mul_si(tmp,x,M,prec); slong m=(slong)llround(arf_get_d(arb_midref(tmp),ARF_RND_NEAR));
    arb_sub_si(e+n-1,tmp,m,prec); mi[n-1]=((m%M)+M)%M;
  }
  /* L: S0 (pi r)^L/L! e^{pi r} <= 1e-4 Gabcke, r=half/M */
  arb_t pr,trunc,fl; arb_init(pr);arb_init(trunc);arb_init(fl);
  arb_set_si(pr,half); arb_div_si(pr,pr,M,prec); arb_mul(pr,pr,pi,prec); gabcke(g,T,prec);
  int Lt=1; for(;;Lt++){ arb_pow_ui(trunc,pr,Lt,prec); arb_fac_ui(fl,Lt,prec); arb_div(trunc,trunc,fl,prec); arb_exp(tmp,pr,prec);
     arb_mul(trunc,trunc,tmp,prec); arb_mul(trunc,trunc,S0,prec); if(arf_get_d(arb_midref(trunc),ARF_RND_UP)<=1e-4*arf_get_d(arb_midref(g),ARF_RND_DOWN)) break;}
  double t1=now();
  for(slong i=0;i<P;i++) acb_one(fac+i);
  acb_dft_pre_t pre; acb_dft_precomp_init(pre,M,prec);
  for(int k=0;k<Lt;k++){
    _acb_vec_zero(G,M); for(slong n=0;n<N;n++) acb_add(G+mi[n],G+mi[n],cur+n,prec);
    if(USE_DBL){ /* rigorous: Higham Thm 24.2 bound + input rounding */
      if(!a){a=malloc(M*sizeof(double complex));W=malloc(M/2*sizeof(double complex));
        arb_t ang,sn,cs; arb_init(ang);arb_init(sn);arb_init(cs);
        for(slong k=0;k<M/2;k++){arb_mul_si(ang,tp,-k,prec);arb_div_si(ang,ang,M,prec);arb_sin_cos(sn,cs,ang,prec);
          W[k]=arf_get_d(arb_midref(cs),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(sn),ARF_RND_NEAR);} arb_clear(ang);arb_clear(sn);arb_clear(cs);}
      double dG1=0,n2=0,uu=ldexp(1,-53);
      for(slong m=0;m<M;m++){double re=arf_get_d(arb_midref(acb_realref(G+m)),ARF_RND_NEAR),im=arf_get_d(arb_midref(acb_imagref(G+m)),ARF_RND_NEAR);
        a[m]=re+_Complex_I*im; dG1+=mag_get_d(arb_radref(acb_realref(G+m)))+mag_get_d(arb_radref(acb_imagref(G+m)))+2*uu*(fabs(re)+fabs(im)); n2+=re*re+im*im;}
      dG1*=1+1e-12; n2=sqrt(n2)*(1+1e-12)+dG1;
      fft_d(a,M,W);
      double lg=log2((double)M), mu=2*uu, g4=4*uu/(1-4*uu), eta=mu+g4*(sqrt(2.0)+mu);
      double err=(lg*eta/(1-lg*eta))*sqrt((double)M)*n2*(1+1e-9)+dG1+ (sqrt((double)M)*n2)*4*uu; /* + rounding of outputs */
      for(slong m=0;m<M;m++){arb_set_d(acb_realref(Gh+m),creal(a[m]));arb_set_d(acb_imagref(Gh+m),cimag(a[m]));
        arb_add_error_2exp_si(acb_realref(Gh+m),-2000); mag_set_d(arb_radref(acb_realref(Gh+m)),err); mag_set_d(arb_radref(acb_imagref(Gh+m)),err);}
    } else acb_dft_precomp(Gh,G,pre,prec);
    for(slong i=0;i<P;i++){ slong q=i-half; acb_addmul(F+i,fac+i,Gh+((q%M)+M)%M,prec);
       /* fac *= (-2 pi i q/M)/(k+1) */
       arb_mul_si(tmp,tp,-q,prec); arb_div_si(tmp,tmp,M*(k+1),prec); acb_mul_onei(z,fac+i); acb_mul_arb(fac+i,z,tmp,prec);}
    for(slong n=0;n<N;n++) acb_mul_arb(cur+n,cur+n,e+n,prec);
  }
  acb_dft_precomp_clear(pre);
  double t2=now();
   arb_t co[4],th3,v,th,c,cc,gg,zr; for(int k=0;k<4;k++)arb_init(co[k]); arb_init(th3);arb_init(v);arb_init(th);arb_init(c);arb_init(cc);arb_init(gg);arb_init(zr);
  double H=1.0; slong nb=(slong)ceil(P*delta/(2*H)); slong ncert=0;
  for(slong b=0;b<nb;b++){
    slong lo=(slong)ceil(b*2*H/delta), hi=(slong)ceil((b+1)*2*H/delta); if(hi>P)hi=P; if(lo>=hi)continue;
    double cu=(2*b+1)*H; arb_set_d(c,cu); arb_add(c,c,T,prec);
    theta_ser(co,c,0,prec); arb_t cb[4]; for(int k=0;k<4;k++)arb_init(cb[k]); theta_ser(cb,c,H,prec); arb_abs(th3,cb[3]);
    for(slong p=lo;p<hi;p++){
      arb_mul_si(v,d,p,prec); arb_sub_d(v,v,cu,prec);   /* v = p*delta - cu (exact delta) */
      arb_mul(th,co[2],v,prec); arb_add(th,th,co[1],prec); arb_mul(th,th,v,prec); arb_add(th,th,co[0],prec);
      arb_abs(tmp,v); arb_pow_ui(tmp,tmp,3,prec); arb_mul(tmp,tmp,th3,prec); arb_add_error(th,tmp);
      acb_set(z,F+p); arb_add_error(acb_realref(z),trunc); arb_add_error(acb_imagref(z),trunc);
      acb_t e1; acb_init(e1); arb_sin_cos(acb_imagref(e1),acb_realref(e1),th,prec); acb_mul(z,z,e1,prec); acb_clear(e1);
      arb_mul_2exp_si(zr,acb_realref(z),1);
      arb_add(tmp,c,v,prec); int ok=c0term(cc,tmp,N,prec); arb_add(zr,zr,cc,prec); gabcke(gg,tmp,prec); arb_add_error(zr,gg);
      sg[p]= (!ok||arb_contains_zero(zr))?0:(arb_is_positive(zr)?1:-1); if(sg[p])ncert++;
    }
    for(int k=0;k<4;k++)arb_clear(cb[k]);
  }
  double t3=now(); tm[0]=t1-t0; tm[1]=t2-t1; tm[2]=t3-t2; *Nout=N; *Lout=Lt; *ncout=ncert;
  if(a){free(a);free(W);}
  for(int k=0;k<4;k++)arb_clear(co[k]); arb_clear(th3);arb_clear(v);arb_clear(th);arb_clear(c);arb_clear(cc);arb_clear(gg);arb_clear(zr);
  arb_clear(pr);arb_clear(trunc);arb_clear(fl); acb_clear(z);
  _acb_vec_clear(cur,N);_acb_vec_clear(G,M);_acb_vec_clear(Gh,M);_acb_vec_clear(F,P);_acb_vec_clear(fac,P);_arb_vec_clear(e,N);free(mi);
  arb_clear(T);arb_clear(d);arb_clear(pi);arb_clear(tp);arb_clear(l);arb_clear(w);arb_clear(x);arb_clear(ph);arb_clear(tmp);arb_clear(S0);arb_clear(g);
}
