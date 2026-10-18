/* Referee harness for docs/DP_CHECK.md: fftzsp.c with main replaced. Runs the all-Arb and double paths on the same grid and compares F point by point.
   usage: dpcheck T delta P   (P a power of two); prints one JSON line. */
/* Original fftzsp header: batched certified Z(t) on uniform grid u_p=p*delta, p<P, t=T+u (Taylor-NUFFT + acb_dft) + Turing/Trudgian. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/acb_dft.h>
#include <complex.h>
#undef I
static acb_ptr G_FOUT=NULL; static double G_EPH=0;
static int USE_DBL=0; static slong G_AN=0; static int MF=2;
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
static void gabcke(arb_t g,const arb_t t,slong prec){ /* Gabcke 1979: |R0| <= 0.127 t^{-3/4}, t>=200 */
  arb_t x; arb_init(x); arb_set_d(g,-0.75); arb_pow(x,t,g,prec); arb_set_d(g,0.127); arb_mul(g,g,x,prec); arb_clear(x);}
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
typedef struct{int init; slong N,M,P; double delta; arb_ptr l,w,x,e; slong*mi; arb_t S0;} zcache;
static void zgrid_c(zcache*C,double Td,double delta,slong P,signed char*sg,slong*Nout,int*Lout,slong*ncout,double*tm){
  slong prec=PREC; double t0=now();
  arb_t T,d,pi,tp,l,w,x,ph,tmp,S0,g; arb_init(T);arb_init(d);arb_init(pi);arb_init(tp);arb_init(l);arb_init(w);arb_init(x);arb_init(ph);arb_init(tmp);arb_init(S0);arb_init(g);
  arb_set_d(T,Td); arb_set_d(d,delta); arb_const_pi(pi,prec); arb_mul_2exp_si(tp,pi,1);
  slong N=(slong)floor(sqrt(Td/(2*M_PI))); slong half=P/2, M=1; while(M<MF*P) M<<=1;
  slong NC=USE_DBL?1:N; acb_ptr cur=_acb_vec_init(NC), G=_acb_vec_init(M), Gh=_acb_vec_init(M), F=_acb_vec_init(P), fac=_acb_vec_init(P);
  int skipall=(USE_DBL&&G_AN>=N); slong NE=skipall?1:N; arb_ptr e=_arb_vec_init(NE); slong *mi=malloc(N*sizeof(slong));
  acb_t z; acb_init(z); arb_zero(S0); double complex *a=NULL,*W=NULL;
  if(C&&C->init&&(C->N!=N||C->M!=M||C->P!=P||C->delta!=delta)){_arb_vec_clear(C->l,C->N);_arb_vec_clear(C->w,C->N);_arb_vec_clear(C->x,C->N);_arb_vec_clear(C->e,C->N);free(C->mi);arb_clear(C->S0);C->init=0;}
  int fresh=(!C||!C->init);
  if(skipall){fresh=0; arb_set_d(S0,2*sqrt((double)N)+1);}
  arb_ptr L_=skipall?NULL:fresh?_arb_vec_init(N):C->l, W_=skipall?NULL:fresh?_arb_vec_init(N):C->w, X_=skipall?NULL:fresh?_arb_vec_init(N):C->x;
  if(fresh){ for(slong n=1;n<=N;n++){
    arb_log_ui(L_+n-1,n,prec); arb_set_ui(W_+n-1,n); arb_rsqrt(W_+n-1,W_+n-1,prec); arb_add(S0,S0,W_+n-1,prec);
    arb_mul(x,d,L_+n-1,prec); arb_div(x,x,tp,prec); arb_floor(tmp,x,prec); arb_sub(x,x,tmp,prec); arb_set(X_+n-1,x);
    arb_mul_si(tmp,x,M,prec); slong m=(slong)llround(arf_get_d(arb_midref(tmp),ARF_RND_NEAR));
    arb_sub_si(e+n-1,tmp,m,prec); mi[n-1]=((m%M)+M)%M; }
    if(C){C->init=1;C->N=N;C->M=M;C->P=P;C->delta=delta;C->l=L_;C->w=W_;C->x=X_;C->e=_arb_vec_init(N);_arb_vec_set(C->e,e,N);C->mi=malloc(N*sizeof(slong));memcpy(C->mi,mi,N*sizeof(slong));arb_init(C->S0);arb_set(C->S0,S0);}
  } else if(!skipall){ if(!USE_DBL){_arb_vec_set(e,C->e,N); memcpy(mi,C->mi,N*sizeof(slong));} arb_set(S0,C->S0); }
  double S0d=0,uu=ldexp(1,-53),EPH=0,rho=0,gc=0,*ed=NULL; double complex*cd=NULL; slong*md=NULL;
  if(USE_DBL){ static double complex TB[1025]; static int tbi=0;
    if(!tbi){arb_t an,sn,cs; arb_init(an);arb_init(sn);arb_init(cs);
      for(int j=0;j<=1024;j++){arb_mul_si(an,tp,-j,prec);arb_div_si(an,an,1024,prec);arb_sin_cos(sn,cs,an,prec);
        TB[j]=arf_get_d(arb_midref(cs),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(sn),ARF_RND_NEAR);} tbi=1; arb_clear(an);arb_clear(sn);arb_clear(cs);}
    cd=malloc(N*sizeof(double complex)); ed=malloc(N*sizeof(double)); md=malloc(N*sizeof(slong));
    double amax=log((double)N+1)/(2*M_PI)*(1+1e-12), dpsi=0, dX=0; const double TWOPI=6.283185307179586;
    static double *AH=NULL,*AL=NULL; slong AN=G_AN; static double ARAD=0;
    if(AN<N){ arb_t an,r2; arb_init(an); arb_init(r2); AH=realloc(AH,N*sizeof(double)); AL=realloc(AL,N*sizeof(double)); ARAD=0;
      for(slong n=1;n<=N;n++){ arb_div(an,L_+n-1,tp,prec); double h=arf_get_d(arb_midref(an),ARF_RND_NEAR);
        arb_sub_d(r2,an,h,prec); double l=arf_get_d(arb_midref(r2),ARF_RND_NEAR); AH[n-1]=h; AL[n-1]=l;
        double rr=mag_get_d(arb_radref(r2))+fabs(l)*uu; if(rr>ARAD)ARAD=rr; } AN=N; G_AN=N; arb_clear(an); arb_clear(r2);}
    if(!(Td==floor(Td)&&fabs(Td)<9e15)) { fprintf(stderr,"T must be integer\n"); exit(1);}
    dpsi=fabs(Td)*ARAD*1.000001+12*uu; double hf=(double)half;
    for(slong n=1;n<=N;n++){ double ah=AH[n-1]; double xr=delta*ah, X=xr-floor(xr);
      double xm=X*(double)M, m=nearbyint(xm); ed[n-1]=xm-m; md[n-1]=((((slong)m)%M)+M)%M;
      double p=Td*ah, pe=fma(Td,ah,-p); double f=(p-floor(p))+pe; f+=Td*AL[n-1]; {double q=hf*X; q-=floor(q); f+=q;} f-=floor(f);
      double jj=nearbyint(f*1024.0), r=f-jj/1024.0, x=r*TWOPI, x2=x*x;
      double c=1-x2*(0.5-x2*(1.0/24-x2/720)), sn=x*(1-x2*(1.0/6-x2*(1.0/120-x2/5040)));
      double w=1.0/sqrt((double)n); S0d+=w;
      cd[n-1]=w*(TB[(int)jj]*(c-_Complex_I*sn)); }
    dX=2.5*uu*delta*amax+1e-300;
    EPH=TWOPI*1.0000001*S0d*1.000001*(dpsi*1.000001+(double)P*dX); /* phase-model error, |F| abs */
    if(getenv("DBGPH"))fprintf(stderr,"dpsi=%g S0d=%g dX=%g\n",dpsi,S0d,dX);
    if(!(dpsi<1e-12)) EPH=1e300; /* unsafe ball -> nothing certified */
    rho=16*uu; /* per-term relative error of cd (table+Taylor+products+sqrt), see SPEED.md */
    int*cnt=calloc(M,sizeof(int)); int cm=0; for(slong n=0;n<N;n++){int v=++cnt[md[n]]; if(v>cm)cm=v;} free(cnt);
    gc=cm*uu/(1-cm*uu);
  } else
  for(slong n=1;n<=N;n++){
    arb_mul(ph,T,L_+n-1,prec); arb_div(ph,ph,tp,prec); arb_mul_si(tmp,X_+n-1,half,prec); arb_add(ph,ph,tmp,prec);
    arb_floor(tmp,ph,prec); arb_sub(ph,ph,tmp,prec); arb_mul(ph,ph,tp,prec); arb_neg(ph,ph);
    arb_sin_cos(acb_imagref(z),acb_realref(z),ph,prec); acb_mul_arb(cur+n-1,z,W_+n-1,prec);
  }
  if(!C&&!skipall){_arb_vec_clear(L_,N);_arb_vec_clear(W_,N);_arb_vec_clear(X_,N);}
  /* L: S0 (pi r)^L/L! e^{pi r} <= 1e-4 Gabcke, r=half/M */
  arb_t pr,trunc,fl; arb_init(pr);arb_init(trunc);arb_init(fl);
  arb_set_si(pr,half); arb_div_si(pr,pr,M,prec); arb_mul(pr,pr,pi,prec); gabcke(g,T,prec);
  int Lt=1; for(;;Lt++){ arb_pow_ui(trunc,pr,Lt,prec); arb_fac_ui(fl,Lt,prec); arb_div(trunc,trunc,fl,prec); arb_exp(tmp,pr,prec);
     arb_mul(trunc,trunc,tmp,prec); arb_mul(trunc,trunc,S0,prec); if(arf_get_d(arb_midref(trunc),ARF_RND_UP)<=1e-4*arf_get_d(arb_midref(g),ARF_RND_DOWN)) break;}
  double t1=now();
  for(slong i=0;i<P;i++) acb_one(fac+i);
  double EF=0,ESUM=0; const double TWOPI_=6.283185307179586; double complex*Fd=NULL,*fd=NULL;
  if(USE_DBL){Fd=calloc(P,sizeof(double complex)); fd=malloc(P*sizeof(double complex)); for(slong i=0;i<P;i++)fd[i]=1;}
  acb_dft_pre_t pre; if(!USE_DBL) acb_dft_precomp_init(pre,M,prec);
  double complex*AF=NULL;
  if(USE_DBL && (double)Lt*M*16<=1.3e8){ AF=calloc((size_t)Lt*M,sizeof(double complex));
    for(slong n=0;n<N;n++){ double complex c=cd[n]; double e_=ed[n]; double complex*p=AF+md[n];
      for(int k=0;k<Lt;k++){ p[(size_t)k*M]+=c; c*=e_; } } }
  for(int k=0;k<Lt;k++){
    if(!USE_DBL){_acb_vec_zero(G,M); for(slong n=0;n<N;n++) acb_add(G+mi[n],G+mi[n],cur+n,prec);}
    if(USE_DBL){ /* rigorous: Higham Thm 24.2 bound + input rounding */
      static double complex*WC=NULL; static slong WM=0;
      if(!a&&WM==M){a=malloc(M*sizeof(double complex));W=WC;}
      if(!a){a=malloc(M*sizeof(double complex));W=malloc(M/2*sizeof(double complex)); if(WC)free(WC); WC=W; WM=M;
        arb_t ang,sn,cs; arb_init(ang);arb_init(sn);arb_init(cs);
        for(slong k=0;k<M/2;k++){arb_mul_si(ang,tp,-k,prec);arb_div_si(ang,ang,M,prec);arb_sin_cos(sn,cs,ang,prec);
          W[k]=arf_get_d(arb_midref(cs),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(sn),ARF_RND_NEAR);} arb_clear(ang);arb_clear(sn);arb_clear(cs);}
      double dG1=0,n2=0,Sk=0;
      if(AF) memcpy(a,AF+(size_t)k*M,M*sizeof(double complex)); else {memset(a,0,M*sizeof(double complex)); for(slong n=0;n<N;n++) a[md[n]]+=cd[n];}
      Sk=S0d*1.000001*ldexp(1.0,-k)*(1+rho)/(1-rho); /* |cd_n|<=w_n(1+rho0)..., |e_n|<=1/2 */ dG1=sqrt(2.0)*(rho+gc*(1+rho))*Sk*(1+1e-9);
      for(slong m=0;m<M;m++) n2+=creal(a[m])*creal(a[m])+cimag(a[m])*cimag(a[m]);
      n2=sqrt(n2)*(1+1e-9)+dG1;
      fft_d(a,M,W);
      double lg=log2((double)M), mu=2*uu, g4=4*uu/(1-4*uu), eta=mu+g4*(sqrt(2.0)+mu);
      double err=(lg*eta/(1-lg*eta))*sqrt((double)M)*n2*(1+1e-9)+dG1+ (sqrt((double)M)*n2)*4*uu; /* + rounding of outputs */
      { double Fk=1; for(int j=1;j<=k;j++) Fk*=(M_PI/2)/j; Fk*=1.000001; double ge=sqrt(2.0)*err;
        EF+=Fk*ge+Fk*(Sk+ge)*(4.0*(k+2))*uu*1.01; ESUM+=Fk*(Sk+ge)*sqrt(2.0);
        double x0=TWOPI_/(double)M;
        for(slong i=0;i<P;i++){ slong q=i-half; Fd[i]+=fd[i]*a[((q%M)+M)%M];
          double x=x0*(double)(-q); fd[i]=fd[i]*(_Complex_I*x)/(double)(k+1);} }
    } else { acb_dft_precomp(Gh,G,pre,prec);
    for(slong i=0;i<P;i++){ slong q=i-half; acb_addmul(F+i,fac+i,Gh+((q%M)+M)%M,prec);
       arb_mul_si(tmp,tp,-q,prec); arb_div_si(tmp,tmp,M*(k+1),prec); acb_mul_onei(z,fac+i); acb_mul_arb(fac+i,z,tmp,prec);} }
    if(USE_DBL){ if(!AF) for(slong n=0;n<N;n++) cd[n]*=ed[n]; rho=rho+uu*(1+rho)*1.0000001; }
    else for(slong n=0;n<N;n++) acb_mul_arb(cur+n,cur+n,e+n,prec);
  }
  if(!USE_DBL) acb_dft_precomp_clear(pre); free(AF);
  if(USE_DBL){ EF+=ESUM*(Lt+2)*uu*1.01; EPH+=EF;
    for(slong i=0;i<P;i++){arb_set_d(acb_realref(F+i),creal(Fd[i]));arb_set_d(acb_imagref(F+i),cimag(Fd[i]));}
    free(Fd);free(fd);}
  double t2=now();
  if(G_FOUT){_acb_vec_set(G_FOUT,F,P);} G_EPH=EPH;
   arb_t co[4],th3,v,th,c,cc,gg,zr; for(int k=0;k<4;k++)arb_init(co[k]); arb_init(th3);arb_init(v);arb_init(th);arb_init(c);arb_init(cc);arb_init(gg);arb_init(zr);
  double H=1.0; slong nb=(slong)ceil(P*delta/(2*H)); slong ncert=0;
  if(USE_DBL){arb_set_d(tmp,EPH); arb_add_error(trunc,tmp);}
  for(slong b=0;b<nb;b++){
    slong lo=(slong)ceil(b*2*H/delta), hi=(slong)ceil((b+1)*2*H/delta); if(hi>P)hi=P; if(lo>=hi)continue;
    double cu=(2*b+1)*H; arb_set_d(c,cu); arb_add(c,c,T,prec);
    theta_ser(co,c,0,prec); arb_t cb[4]; for(int k=0;k<4;k++)arb_init(cb[k]); theta_ser(cb,c,H,prec); arb_abs(th3,cb[3]);
    slong pp=USE_DBL?64:prec;
    if(USE_DBL){ arb_div(tmp,co[0],tp,prec); arb_floor(tmp,tmp,prec); arb_mul(tmp,tmp,tp,prec); arb_sub(co[0],co[0],tmp,prec);
      arb_set_d(tmp,cu-H); arb_add(tmp,tmp,T,prec); gabcke(gg,tmp,prec); }
    for(slong p=lo;p<hi;p++){ if(USE_DBL){
      arb_mul_si(v,d,p,pp); arb_sub_d(v,v,cu,pp);
      arb_mul(th,co[2],v,pp); arb_add(th,th,co[1],pp); arb_mul(th,th,v,pp); arb_add(th,th,co[0],pp);
      arb_abs(tmp,v); arb_pow_ui(tmp,tmp,3,pp); arb_mul(tmp,tmp,th3,pp); arb_add_error(th,tmp);
      acb_set(z,F+p); arb_add_error(acb_realref(z),trunc); arb_add_error(acb_imagref(z),trunc);
      acb_t e1; acb_init(e1); arb_sin_cos(acb_imagref(e1),acb_realref(e1),th,pp); acb_mul(z,z,e1,pp); acb_clear(e1);
      arb_mul_2exp_si(zr,acb_realref(z),1);
      arb_add(tmp,c,v,prec); int ok=c0term(cc,tmp,N,pp); arb_add(zr,zr,cc,pp); arb_add_error(zr,gg);
      sg[p]= (!ok||arb_contains_zero(zr))?0:(arb_is_positive(zr)?1:-1); if(sg[p])ncert++; continue; }
      arb_mul_si(v,d,p,prec); arb_sub_d(v,v,cu,prec);   /* v = p*delta - cu (exact delta) */
      arb_mul(th,co[2],v,prec); arb_add(th,th,co[1],prec); arb_mul(th,th,v,prec); arb_add(th,th,co[0],prec);
      arb_abs(tmp,v); arb_pow_ui(tmp,tmp,3,prec); arb_mul(tmp,tmp,th3,prec); arb_add_error(th,tmp);
      acb_set(z,F+p); arb_add_error(acb_realref(z),trunc); arb_add_error(acb_imagref(z),trunc); if(USE_DBL){arb_set_d(tmp,EPH); arb_add_error(acb_realref(z),tmp); arb_add_error(acb_imagref(z),tmp);}
      acb_t e1; acb_init(e1); arb_sin_cos(acb_imagref(e1),acb_realref(e1),th,prec); acb_mul(z,z,e1,prec); acb_clear(e1);
      arb_mul_2exp_si(zr,acb_realref(z),1);
      arb_add(tmp,c,v,prec); int ok=c0term(cc,tmp,N,prec); arb_add(zr,zr,cc,prec); gabcke(gg,tmp,prec); arb_add_error(zr,gg);
      sg[p]= (!ok||arb_contains_zero(zr))?0:(arb_is_positive(zr)?1:-1); if(sg[p])ncert++;
    }
    for(int k=0;k<4;k++)arb_clear(cb[k]);
  }
  double t3=now(); tm[0]=t1-t0; tm[1]=t2-t1; tm[2]=t3-t2; *Nout=N; *Lout=Lt; *ncout=ncert;
  if(a){free(a);} if(cd){free(cd);free(ed);free(md);} tm[3]=EPH;
  for(int k=0;k<4;k++)arb_clear(co[k]); arb_clear(th3);arb_clear(v);arb_clear(th);arb_clear(c);arb_clear(cc);arb_clear(gg);arb_clear(zr);
  arb_clear(pr);arb_clear(trunc);arb_clear(fl); acb_clear(z);
  _acb_vec_clear(cur,NC);_acb_vec_clear(G,M);_acb_vec_clear(Gh,M);_acb_vec_clear(F,P);_acb_vec_clear(fac,P);_arb_vec_clear(e,NE);free(mi);
  arb_clear(T);arb_clear(d);arb_clear(pi);arb_clear(tp);arb_clear(l);arb_clear(w);arb_clear(x);arb_clear(ph);arb_clear(tmp);arb_clear(S0);arb_clear(g);
}

static void zgrid(double Td,double delta,slong P,signed char*sg,slong*Nout,int*Lout,slong*ncout,double*tm){zgrid_c(NULL,Td,delta,P,sg,Nout,Lout,ncout,tm);}
typedef struct{double u; signed char s;} pt;
static int cmp_pt(const void*x,const void*y){double a=((pt*)x)->u,b=((pt*)y)->u; return a<b?-1:a>b;}
/* usage: fftzr T delta P a b Lw [signfile] [dbl]  -- coarse grid, then automatic in-process refinement of wide gaps */

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: dpcheck T delta P\n");return 2;}
  double Td=atof(argv[1]),delta=atof(argv[2]); slong P=atol(argv[3]);
  if(P<2||(P&(P-1))){fprintf(stderr,"error: P must be a power of two\n");return 2;}
  signed char*sg=malloc(P); slong N,ncert; int Lt; double tm[4];
  acb_ptr Fa=_acb_vec_init(P), Fd=_acb_vec_init(P);
  USE_DBL=0; G_FOUT=Fa; zgrid(Td,delta,P,sg,&N,&Lt,&ncert,tm); slong nca=ncert; double ta=tm[1];
  USE_DBL=1; G_FOUT=Fd; zgrid(Td,delta,P,sg,&N,&Lt,&ncert,tm); slong ncd=ncert;
  double maxdiff=0,maxrad=0; slong bad=0; arb_t t; arb_init(t);
  for(slong i=0;i<P;i++){ for(int c=0;c<2;c++){ arb_srcptr a=c?acb_imagref(Fa+i):acb_realref(Fa+i); arb_srcptr d=c?acb_imagref(Fd+i):acb_realref(Fd+i);
      arb_sub(t,d,a,128); double df=fabs(arf_get_d(arb_midref(t),ARF_RND_NEAR)); double rd=mag_get_d(arb_radref(a));
      if(df>maxdiff)maxdiff=df; if(rd>maxrad)maxrad=rd; if(df+rd>G_EPH) bad++; } }
  printf("{\"T\":%.0f,\"delta\":%g,\"P\":%ld,\"N\":%ld,\"L\":%d,\"cert_arb\":%ld,\"cert_dbl\":%ld,\"EPH\":%.3e,\"max_abs_diff\":%.3e,\"max_arb_rad\":%.3e,\"ratio\":%.3e,\"violations\":%ld}\n",
    Td,delta,P,N,Lt,nca,ncd,G_EPH,maxdiff,maxrad,maxdiff/G_EPH,bad);
  return 0;}
