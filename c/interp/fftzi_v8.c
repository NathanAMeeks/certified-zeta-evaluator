/* [NUM] batched certified Z(t) on uniform grid u_p=p*delta, p<P, t=T+u (Taylor-NUFFT + acb_dft) + Turing/Trudgian.
   usage: fftz T delta P a b Lw [signfile]   (T integer-valued double ok) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/acb_dft.h>
#include <flint/arb_hypgeom.h>
#include <complex.h>
#undef I
static double complex TB[1025]; static int tbi=0; static double DT_MAX=0; static long long DT_N=0,DT_BAD=0;
static int USE_DBL=0; static slong G_AN=0; static double *AH=NULL,*AL=NULL,ARAD=0;
/* speed3: on-disk cache of a_n=log(n)/2pi as double-double (AH+AL) with Arb-derived radius ARAD; FNV-1a checksum */
static unsigned long long fnv(const void*p,size_t n,unsigned long long h){const unsigned char*c=p;for(size_t i=0;i<n;i++){h^=c[i];h*=1099511628211ULL;}return h;}
static int ah_load(const char*fn,slong need){ FILE*f=fopen(fn,"rb"); if(!f)return 0; char mg[4]; long long nc=0; double ar=0; unsigned long long ck=0,h=1469598103934665603ULL;
  if(fread(mg,1,4,f)!=4||memcmp(mg,"AHC1",4)||fread(&nc,8,1,f)!=1||fread(&ar,8,1,f)!=1){fclose(f);return 0;}
  AH=realloc(AH,nc*8); AL=realloc(AL,nc*8); if(fread(AH,8,nc,f)!=(size_t)nc||fread(AL,8,nc,f)!=(size_t)nc||fread(&ck,8,1,f)!=1){fclose(f);return 0;} fclose(f);
  h=fnv(&nc,8,h); h=fnv(&ar,8,h); h=fnv(AH,nc*8,h); h=fnv(AL,nc*8,h); if(h!=ck){fprintf(stderr,"AH cache checksum mismatch\n");return 0;} ARAD=ar; G_AN=nc; return nc>=need?1:2; }
static void ah_save(const char*fn){ FILE*f=fopen(fn,"wb"); if(!f)return; long long nc=G_AN; unsigned long long h=1469598103934665603ULL; h=fnv(&nc,8,h); h=fnv(&ARAD,8,h); h=fnv(AH,nc*8,h); h=fnv(AL,nc*8,h);
  fwrite("AHC1",1,4,f); fwrite(&nc,8,1,f); fwrite(&ARAD,8,1,f); fwrite(AH,8,nc,f); fwrite(AL,8,nc,f); fwrite(&h,8,1,f); fclose(f); } static int MF=2;
/* speed2: saved coarse grid for certified band-limited interpolation (Lemma BL) */
static double complex*GF=NULL; static double GRF=0,GS0=0,GDEL=0,*GZ=NULL; static slong GP=0,GNN=0; static long long NINTERP=0; static double complex*GE=NULL; static double*GTP=NULL;
/* speed3: explicit textbook complex product (ac-bd, ad+bc): bit-identical to C's operator* on finite values, but no __muldc3 NaN fallback, so loops vectorize */
static inline double complex cmul(double complex a,double complex b){double ar=creal(a),ai=cimag(a),br=creal(b),bi=cimag(b); return CMPLX(ar*br-ai*bi, ar*bi+ai*br);}
/* [HEUR only, never in a certified path] 4-lane SoA complex dot product for targeting samples */
static inline double complex mac4(const double*gr,const double*gi,const double*wr,const double*wi,slong n){
  double ar[4]={0,0,0,0},ai[4]={0,0,0,0}; slong m=0;
  for(;m+4<=n;m+=4) for(int l=0;l<4;l++){ ar[l]+=gr[m+l]*wr[m+l]-gi[m+l]*wi[m+l]; ai[l]+=gr[m+l]*wi[m+l]+gi[m+l]*wr[m+l]; }
  for(;m<n;m++){ ar[0]+=gr[m]*wr[m]-gi[m]*wi[m]; ai[0]+=gr[m]*wi[m]+gi[m]*wr[m]; }
  return CMPLX((ar[0]+ar[1])+(ar[2]+ar[3]),(ai[0]+ai[1])+(ai[2]+ai[3])); }
/* iterative radix-2 DIT FFT, sign -1 */
static void fft_d(double complex*a,slong M,const double complex*W){
  for(slong i=1,j=0;i<M;i++){slong bit=M>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j){double complex t=a[i];a[i]=a[j];a[j]=t;}}
  for(slong len=2;len<=M;len<<=1){slong st=M/len;for(slong i=0;i<M;i+=len)for(slong k=0;k<len/2;k++){double complex u=a[i+k],v=cmul(a[i+k+len/2],W[k*st]);a[i+k]=u+v;a[i+k+len/2]=u-v;}}}

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
  arb_t x; arb_init(x); arb_set_d(g,-0.75); arb_pow(x,t,g,prec); arb_set_ui(g,127); arb_div_ui(g,g,1000,prec); arb_mul(g,g,x,prec); /* exact 0.127 ball (was double 0.127) */ arb_clear(x);}
/* C0 term; returns 0 if |cos 2pi p|<1e-3 */
static int c0term(arb_t res,const arb_t t,slong N,slong prec){
  arb_t a,p,pi2,num,den,x; arb_init(a);arb_init(p);arb_init(pi2);arb_init(num);arb_init(den);arb_init(x);
  arb_const_pi(pi2,prec); arb_mul_2exp_si(pi2,pi2,1); arb_div(a,t,pi2,prec); arb_sqrt(a,a,prec); arb_sub_si(p,a,N,prec);
  arb_mul(den,pi2,p,prec); arb_cos(den,den,prec); int ok=1;
  if(!(arb_is_nonnegative(p) && arb_lt(p,(arb_set_ui(x,1),x)))) ok=0; /* RS needs N=floor(sqrt(t/2pi)) at this t */
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
  slong N=(slong)floor(sqrt(Td/(2*M_PI))); slong half=P/2, M=1; int mf_=(!C&&getenv("GMF"))?atoi(getenv("GMF")):MF; while(M<mf_*P) M<<=1;
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
  if(USE_DBL){
    if(!tbi){arb_t an,sn,cs; arb_init(an);arb_init(sn);arb_init(cs);
      for(int j=0;j<=1024;j++){arb_mul_si(an,tp,-j,prec);arb_div_si(an,an,1024,prec);arb_sin_cos(sn,cs,an,prec);
        TB[j]=arf_get_d(arb_midref(cs),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(sn),ARF_RND_NEAR);} tbi=1; arb_clear(an);arb_clear(sn);arb_clear(cs);}
    cd=malloc(N*sizeof(double complex)); ed=malloc(N*sizeof(double)); md=malloc(N*sizeof(slong));
    double amax=log((double)N+1)/(2*M_PI)*(1+1e-12), dpsi=0, dX=0; const double TWOPI=6.283185307179586;
    slong AN=G_AN;
    if(AN<N){ arb_t an,r2; arb_init(an); arb_init(r2); AH=realloc(AH,N*sizeof(double)); AL=realloc(AL,N*sizeof(double)); ARAD=0;
      for(slong n=1;n<=N;n++){ arb_div(an,L_+n-1,tp,prec); double h=arf_get_d(arb_midref(an),ARF_RND_NEAR);
        arb_sub_d(r2,an,h,prec); double l=arf_get_d(arb_midref(r2),ARF_RND_NEAR); AH[n-1]=h; AL[n-1]=l;
        double rr=mag_get_d(arb_radref(r2))+fabs(l)*uu; if(rr>ARAD)ARAD=rr; } AN=N; G_AN=N; arb_clear(an); arb_clear(r2);}
    if(!(Td==floor(Td)&&fabs(Td)<9e15)) { fprintf(stderr,"T must be integer\n"); exit(1);}
    double hf=(double)half; dpsi=fabs(Td)*ARAD*1.000001+12*uu+((half&(half-1))?hf*uu*1.01:0); /* hf*X exact only if half=2^k */
    for(slong n=1;n<=N;n++){ double ah=AH[n-1]; double xr=delta*ah, X=xr-floor(xr);
      double xm=X*(double)M, m=nearbyint(xm); ed[n-1]=xm-m; md[n-1]=((((slong)m)%M)+M)%M;
      double p=Td*ah, pe=fma(Td,ah,-p); double f=(p-floor(p))+pe; f+=Td*AL[n-1]; {double q=hf*X; q-=floor(q); f+=q;} f-=floor(f);
      double jj=nearbyint(f*1024.0), r=f-jj/1024.0, x=r*TWOPI, x2=x*x;
      double c=1-x2*(0.5-x2*(1.0/24-x2/720)), sn=x*(1-x2*(1.0/6-x2*(1.0/120-x2/5040)));
      double w=1.0/sqrt((double)n); S0d+=w;
      cd[n-1]=w*cmul(TB[(int)jj],CMPLX(c,-sn)); }
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
  double LTF=getenv("LTF")?atof(getenv("LTF")):100.0; /* trunc target (x Gabcke); any value is rigorous: trunc is added to every ball */
  int Lt=1; for(;;Lt++){ arb_pow_ui(trunc,pr,Lt,prec); arb_fac_ui(fl,Lt,prec); arb_div(trunc,trunc,fl,prec); arb_exp(tmp,pr,prec);
     arb_mul(trunc,trunc,tmp,prec); arb_mul(trunc,trunc,S0,prec); if(arf_get_d(arb_midref(trunc),ARF_RND_UP)<=LTF*arf_get_d(arb_midref(g),ARF_RND_DOWN)) break;}
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
      n2=sqrt(n2)*(1+(double)(M+2)*uu)+dG1; /* gamma_{M+1} for the sum of squares; old 1e-9 only valid for M<=2^24 */
      fft_d(a,M,W);
      double lg=log2((double)M), mu=2*uu, g4=4*uu/(1-4*uu), eta=mu+g4*(sqrt(2.0)+mu);
      double err=(lg*eta/(1-lg*eta))*sqrt((double)M)*n2*(1+1e-9)+dG1+ (sqrt((double)M)*n2)*4*uu; /* + rounding of outputs */
      { double Fk=1; for(int j=1;j<=k;j++) Fk*=(M_PI/2)/j; Fk*=1.000001; double ge=sqrt(2.0)*err;
        EF+=Fk*ge+Fk*(Sk+ge)*(4.0*(k+2))*uu*1.01; ESUM+=Fk*(Sk+ge)*sqrt(2.0);
        double x0=TWOPI_/(double)M;
        for(slong i=0;i<P;i++){ slong q=i-half; Fd[i]+=cmul(fd[i],a[((q%M)+M)%M]);
          double x=x0*(double)(-q); double kk1=(double)(k+1); fd[i]=CMPLX((-cimag(fd[i])*x)/kk1,(creal(fd[i])*x)/kk1);} /* = fd*(i x)/(k+1) exactly as before (+-0 terms exact) */ }
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
  if(getenv("DUMPF")){FILE*fo=fopen(getenv("DUMPF"),"w"); for(slong i=0;i<P;i++){fprintf(fo,"%ld %.17g %.17g %.6g %.6g %.6g\n",i,arf_get_d(arb_midref(acb_realref(F+i)),ARF_RND_NEAR),arf_get_d(arb_midref(acb_imagref(F+i)),ARF_RND_NEAR),mag_get_d(arb_radref(acb_realref(F+i))),mag_get_d(arb_radref(acb_imagref(F+i))),USE_DBL?EPH:0.0);} fclose(fo); exit(0);}
  double t2=now();
   arb_t co[4],th3,v,th,c,cc,gg,zr; for(int k=0;k<4;k++)arb_init(co[k]); arb_init(th3);arb_init(v);arb_init(th);arb_init(c);arb_init(cc);arb_init(gg);arb_init(zr);
  double H=getenv("PH")?atof(getenv("PH")):(getenv("OLDPTS")?1.0:8.0); const int dbt=getenv("DBLTEST")!=NULL; slong nb=(slong)ceil(P*delta/(2*H)); slong ncert=0;
  if(USE_DBL){arb_set_d(tmp,EPH); arb_add_error(trunc,tmp);}
  if(USE_DBL&&!C){ arf_t ub; arf_init(ub); GP=P; GDEL=delta; GNN=N; GF=malloc(P*sizeof(double complex)); GZ=calloc(P,sizeof(double)); GE=calloc(P,sizeof(double complex)); GTP=calloc(P,sizeof(double));
    for(slong i=0;i<P;i++) GF[i]=arf_get_d(arb_midref(acb_realref(F+i)),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(acb_imagref(F+i)),ARF_RND_NEAR); /* exact: mids are doubles */
    arb_get_abs_ubound_arf(ub,trunc,prec); GRF=arf_get_d(ub,ARF_RND_UP)*(1+1e-12); arb_get_ubound_arf(ub,S0,prec); GS0=arf_get_d(ub,ARF_RND_UP)*(1+1e-12); arf_clear(ub); }
  for(slong b=0;b<nb;b++){
    slong lo=(slong)ceil(b*2*H/delta), hi=(slong)ceil((b+1)*2*H/delta); if(hi>P)hi=P; if(lo>=hi)continue;
    double cu=(2*b+1)*H; arb_set_d(c,cu); arb_add(c,c,T,prec);
    theta_ser(co,c,0,prec); arb_t cb[4]; for(int k=0;k<4;k++)arb_init(cb[k]); theta_ser(cb,c,H,prec); arb_abs(th3,cb[3]);
    slong pp=USE_DBL?64:prec;
    if(USE_DBL){ arb_div(tmp,co[0],tp,prec); arb_floor(tmp,tmp,prec); arb_mul(tmp,tmp,tp,prec); arb_sub(co[0],co[0],tmp,prec);
      arb_set_d(tmp,cu-H); arb_add(tmp,tmp,T,prec); gabcke(gg,tmp,prec); }
    /* ---- speed3 [LEM DP, candidate]: per-point step in binary64 with hand bounds (see SPEED2.md §8). Built with -ffp-contract=off. ---- */
    int dblok=0; double cod[3],cor[3],th3d=0,C0m=0,C0r=0,ggd=0,trd=0;
    if(USE_DBL&&!getenv("OLDPTS")){ arf_t ub_; arf_init(ub_); arb_t bb; arb_init(bb);
      for(int k=0;k<3;k++){ cod[k]=arf_get_d(arb_midref(co[k]),ARF_RND_NEAR); arb_sub_d(bb,co[k],cod[k],prec); arb_get_abs_ubound_arf(ub_,bb,prec); cor[k]=arf_get_d(ub_,ARF_RND_UP); }
      arb_get_abs_ubound_arf(ub_,th3,prec); th3d=arf_get_d(ub_,ARF_RND_UP);
      arb_set(bb,c); mag_set_d(arb_radref(bb),H*1.000001); /* t-ball covering the whole block */
      if(c0term(cc,bb,N,prec)){ C0m=arf_get_d(arb_midref(cc),ARF_RND_NEAR); arb_sub_d(bb,cc,C0m,prec); arb_get_abs_ubound_arf(ub_,bb,prec); C0r=arf_get_d(ub_,ARF_RND_UP);
        arb_get_abs_ubound_arf(ub_,gg,prec); ggd=arf_get_d(ub_,ARF_RND_UP); arb_get_abs_ubound_arf(ub_,trunc,prec); trd=arf_get_d(ub_,ARF_RND_UP);
        dblok=(C0r<1e-3&&isfinite(trd)&&trd<1e-3); }
      arf_clear(ub_); arb_clear(bb); }
    for(slong p=lo;p<hi;p++){ if(dblok){ const double u_=ldexp(1,-53), R2=1.4142135623730951, INV2PI=0.15915494309189535, TWOPI2=6.283185307179586;
      double pd=(double)p*delta, vd=pd-cu, ev=1.01*u_*(fabs(pd)+fabs(vd)), av=fabs(vd)+ev;
      double th_=cod[0]+vd*(cod[1]+vd*cod[2]);
      double eth=cor[0]+cor[1]*av+cor[2]*av*av+av*av*av*th3d + (fabs(cod[1])+cor[1])*ev+2*av*(fabs(cod[2])+cor[2])*ev*1.01 + (4*u_/(1-4*u_))*(fabs(cod[0])+av*fabs(cod[1])+av*av*fabs(cod[2])); eth*=1.25; /* gamma_4 Horner; 25% margin (Taylor term attained to 0.999 in referee tests) */
      double f=th_*INV2PI; f-=floor(f); double jj=nearbyint(f*1024.0), r=f-jj/1024.0, x=r*TWOPI2, x2=x*x;
      double cs=1-x2*(0.5-x2*(1.0/24-x2/720)), sn=x*(1-x2*(1.0/6-x2*(1.0/120-x2/5040)));
      double complex tb=TB[(int)jj]; double ec=creal(tb)*cs+cimag(tb)*sn, es=creal(tb)*sn-cimag(tb)*cs; /* conj(TB)*(cs+i sn) = e^{i th} */
      double ecs=(16+3.15)*u_+2.02*u_*fabs(th_)+1e-300; /* +pi*u: f-floor(f) inexact for -1<f<0 */
      double Fr=arf_get_d(arb_midref(acb_realref(F+p)),ARF_RND_NEAR), Fi=arf_get_d(arb_midref(acb_imagref(F+p)),ARF_RND_NEAR);
      double a_=Fr*ec, b_=Fi*es, d_=a_-b_, zr_=2*d_+C0m, aF=fabs(Fr)+fabs(Fi);
      double Ent=2*aF*(eth+ecs)*1.000001 + 2.02*u_*(fabs(a_)+fabs(b_)+fabs(d_))*1.01 + 1.01*u_*fabs(zr_) + C0r;
      double E=(Ent + 2*R2*trd + 2*R2*trd*(eth+ecs) + ggd)*(1+1e-12);
      sg[p]=(fabs(zr_)>E)?(zr_>0?1:-1):0; if(sg[p])ncert++;
      if(GZ&&!C){ GZ[p]=zr_; GE[p]=ec+_Complex_I*es; GTP[p]=cod[1]+2*vd*cod[2]; }
      if(dbt&&(p%13==0)){ /* adversarial: 256-bit Arb truth for 2Re(e^{i theta(t)}F^)+C0(t) at the same point */
        arb_t tt,th2[4],z2; acb_t w2; arb_init(tt);arb_init(z2);acb_init(w2); for(int k=0;k<4;k++)arb_init(th2[k]);
        arb_set_d(tt,delta); arb_mul_si(tt,tt,p,256); arb_add(tt,tt,T,256); theta_ser(th2,tt,0,256);
        acb_set_d_d(w2,Fr,Fi); acb_t e2; acb_init(e2); arb_sin_cos(acb_imagref(e2),acb_realref(e2),th2[0],256); acb_mul(w2,w2,e2,256); acb_clear(e2);
        arb_mul_2exp_si(z2,acb_realref(w2),1); int ok2=c0term(cc,tt,N,256); arb_add(z2,z2,cc,256);
        if(ok2){ double tr_=arf_get_d(arb_midref(z2),ARF_RND_NEAR); double rr=fabs(zr_-tr_)/Ent; if(rr>DT_MAX)DT_MAX=rr; if(rr>1)DT_BAD++; DT_N++; }
        arb_clear(tt);arb_clear(z2);acb_clear(w2); for(int k=0;k<4;k++)arb_clear(th2[k]); }
      continue; }
      if(USE_DBL){
      arb_mul_si(v,d,p,pp); arb_sub_d(v,v,cu,pp);
      arb_mul(th,co[2],v,pp); arb_add(th,th,co[1],pp); arb_mul(th,th,v,pp); arb_add(th,th,co[0],pp);
      arb_abs(tmp,v); arb_pow_ui(tmp,tmp,3,pp); arb_mul(tmp,tmp,th3,pp); arb_add_error(th,tmp);
      acb_set(z,F+p); arb_add_error(acb_realref(z),trunc); arb_add_error(acb_imagref(z),trunc);
      acb_t e1; acb_init(e1); arb_sin_cos(acb_imagref(e1),acb_realref(e1),th,pp); acb_mul(z,z,e1,pp);
      if(GE&&!C){ GE[p]=arf_get_d(arb_midref(acb_realref(e1)),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(acb_imagref(e1)),ARF_RND_NEAR);
        double vd=arf_get_d(arb_midref(v),ARF_RND_NEAR); GTP[p]=arf_get_d(arb_midref(co[1]),ARF_RND_NEAR)+2*vd*arf_get_d(arb_midref(co[2]),ARF_RND_NEAR); }
      acb_clear(e1);
      arb_mul_2exp_si(zr,acb_realref(z),1);
      arb_add(tmp,c,v,prec); int ok=c0term(cc,tmp,N,pp); arb_add(zr,zr,cc,pp); arb_add_error(zr,gg);
      sg[p]= (!ok||arb_contains_zero(zr))?0:(arb_is_positive(zr)?1:-1); if(sg[p])ncert++; if(GZ&&!C) GZ[p]=arf_get_d(arb_midref(zr),ARF_RND_NEAR); continue; }
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
/* [LEM BL, candidate] Certified band-limited (Gaussian-windowed sinc) interpolation of the RS main sum
   F(u)=sum_{n<=N} n^{-1/2} e^{-i(T+u)log n} from the coarse grid F(k*delta), k<GP (Platt 2017-style; own derivation, see SPEED2.md).
   g(u)=e^{i b0 u}F(u) has spectrum in [-b,b], b=b0=log(N)/2. K(x)=sinc(Ax)exp(-x^2/2h^2), A=pi/delta.
   |F(u)-sum_{|j|<=K} F(u_k) e^{-i b0 x_j} K(x_j)| <= S0*(Eal+Etr), Eal=erfc(z1)(1+1/(1-e^{-2 z1 D})), z1=(A-b)h/sqrt2, D=sqrt2*A*h,
   Etr=2 e^{-X^2/2h^2}/(A X (1-e^{-X delta/h^2})), X=(K-1/2)delta. Grid-ball and double-rounding errors added at runtime. */
static FILE*BLF=NULL; double BL_EAL=0,BL_ETR=0;
static void zinterp(double Td,double o,double dr,slong PRw,signed char*sr){
  static int init=0; static double Ad,Qd,hd,dd,Eal,Etr; static double complex wd; static slong K; static arb_t A,b0,i2h;
  const slong pr=128; const double uu=ldexp(1,-53);
  arb_t x,y,t1,t2; arb_init(x);arb_init(y);arb_init(t1);arb_init(t2);
  if(!init){ arb_init(A);arb_init(b0);arb_init(i2h); dd=GDEL; arb_const_pi(A,pr); arb_div_d(A,A,dd,pr); Ad=arf_get_d(arb_midref(A),ARF_RND_NEAR);
    arb_log_ui(b0,GNN,pr); arb_mul_2exp_si(b0,b0,-1); arb_sub(t1,A,b0,pr); double Amb=arf_get_d(arb_midref(t1),ARF_RND_DOWN);
    if(!(Amb>0)){fprintf(stderr,"delta too coarse for interpolation\n");exit(1);}
    double tgt=getenv("BLTGT")?atof(getenv("BLTGT")):92.0; /* X(A-b) target: error ~ e^{-tgt/2} */
    K=(slong)ceil(tgt/Amb/dd+0.5)+1; double X=(K-0.5)*dd; hd=sqrt(X/Amb);
    arb_set_d(x,hd); arb_mul(i2h,x,x,pr); arb_mul_2exp_si(i2h,i2h,1); arb_inv(i2h,i2h,pr); /* 1/(2h^2) */
    arb_mul_d(t2,t1,hd,pr); arb_sqrt_ui(y,2,pr); arb_div(t2,t2,y,pr); /* z1 */
    arb_mul_d(x,A,hd,pr); arb_mul(x,x,y,pr); /* D */ arb_mul(x,x,t2,pr); arb_mul_2exp_si(x,x,1); arb_neg(x,x); arb_exp(x,x,pr); arb_sub_ui(x,x,1,pr); arb_neg(x,x); arb_inv(x,x,pr); arb_add_ui(x,x,1,pr);
    arb_hypgeom_erfc(t2,t2,pr); arb_mul(t2,t2,x,pr); arf_t ub; arf_init(ub); arb_get_ubound_arf(ub,t2,pr); Eal=arf_get_d(ub,ARF_RND_UP);
    arb_set_d(x,X); arb_mul(t1,x,x,pr); arb_mul(t1,t1,i2h,pr); arb_neg(t1,t1); arb_exp(t1,t1,pr); /* e^{-X^2/2h^2} */
    arb_set_d(y,X*dd); arb_mul(y,y,i2h,pr); arb_mul_2exp_si(y,y,1); arb_neg(y,y); arb_exp(y,y,pr); arb_sub_ui(y,y,1,pr); arb_neg(y,y); /* 1-e^{-X d/h^2} (X*dd rounding: tiny, absorbed below) */
    arb_mul(y,y,A,pr); arb_mul(y,y,x,pr); arb_div(t1,t1,y,pr); arb_mul_2exp_si(t1,t1,1); arb_get_ubound_arf(ub,t1,pr); Etr=arf_get_d(ub,ARF_RND_UP)*1.01; arf_clear(ub);
    arb_set_d(x,dd); arb_mul(x,x,x,pr); arb_mul(x,x,i2h,pr); arb_mul_2exp_si(x,x,1); arb_neg(x,x); arb_exp(x,x,pr); Qd=arf_get_d(arb_midref(x),ARF_RND_NEAR);
    arb_mul_d(x,b0,-dd,pr); arb_sin_cos(t1,t2,x,pr); wd=arf_get_d(arb_midref(t2),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(t1),ARF_RND_NEAR);
    BL_EAL=Eal; BL_ETR=Etr; fprintf(stderr,"BL: delta=%g A-b=%.4g K=%ld h=%.4g Eal=%.3g Etr=%.3g S0=%.4g grid_rad=%.3g\n",dd,Amb,K,hd,Eal,Etr,GS0,GRF); init=1; }
  /* theta Taylor about window centre (as in zgrid_c) */
  static arb_t co[4],cb[4],c,T,gg; static int cinit=0; static double cuc=-1e300,cH=0,cTd=0; /* speed3: theta/Gabcke cache per 8-unit block */
  arb_t v,th,zr,cc,ub_; acb_t z,e1; arb_init(v);arb_init(th);arb_init(zr);arb_init(cc);arb_init(ub_); acb_init(z);acb_init(e1);
  if(!cinit){for(int k=0;k<4;k++){arb_init(co[k]);arb_init(cb[k]);} arb_init(c);arb_init(T);arb_init(gg); cinit=1;}
  double span=PRw*dr, uc, H;
  if(span<=4){ uc=floor(o/8.0)*8.0+4.0; H=4.0+span+dr; if(o+span>uc+H-dr) {uc=o+0.5*span; H=0.5*span+dr;} }
  else { uc=o+0.5*span; H=0.5*span+dr; }
  if(!(uc==cuc&&H<=cH&&Td==cTd)){ if(span<=4&&H<8.0+dr) H=8.0+dr; cuc=uc; cH=H; cTd=Td;
    arb_set_d(T,Td); arb_set_d(c,uc); arb_add(c,c,T,pr);
    theta_ser(co,c,0,pr); theta_ser(cb,c,H,pr); arb_abs(cb[3],cb[3]);
    arb_const_pi(t1,pr); arb_mul_2exp_si(t1,t1,1); arb_div(t2,co[0],t1,pr); arb_floor(t2,t2,pr); arb_mul(t2,t2,t1,pr); arb_sub(co[0],co[0],t2,pr);
    arb_set_d(t1,uc-H); arb_add(t1,t1,T,pr); gabcke(gg,t1,pr); /* Gabcke bound decreasing in t: take lower end of block */ }
  double eps_t=((double)(K+3)*(K+3)+20.0*(K+1))*uu, gam=(2*K+2)*uu/(1-(2*K+2)*uu);
  for(slong q=0;q<PRw;q++){ double ud=o+q*dr; slong k0=(slong)llround(ud/dd); sr[q]=0; if(k0-K<0||k0+K>=GP) continue;
    arb_set_d(x,dr); arb_mul_si(x,x,q,pr); arb_set_d(t1,o); arb_add(x,x,t1,pr); arb_set_d(t1,dd); arb_mul_si(t1,t1,k0,pr); arb_sub(x,x,t1,pr); /* x0=u-k0*delta */
    arb_mul(y,A,x,pr); arb_sin(t1,y,pr); double s=arf_get_d(arb_midref(t1),ARF_RND_NEAR); arb_sinc(t1,y,pr); double sc=arf_get_d(arb_midref(t1),ARF_RND_NEAR);
    arb_mul(t1,x,x,pr); arb_mul(t1,t1,i2h,pr); arb_neg(t1,t1); arb_exp(t1,t1,pr); double G0=arf_get_d(arb_midref(t1),ARF_RND_NEAR);
    arb_mul_2exp_si(t1,x,1); arb_add_d(t2,t1,dd,pr); arb_mul_d(t2,t2,dd,pr); arb_mul(t2,t2,i2h,pr); arb_neg(t2,t2); arb_exp(t2,t2,pr); double Rp=arf_get_d(arb_midref(t2),ARF_RND_NEAR);
    arb_neg(t1,t1); arb_add_d(t2,t1,dd,pr); arb_mul_d(t2,t2,dd,pr); arb_mul(t2,t2,i2h,pr); arb_neg(t2,t2); arb_exp(t2,t2,pr); double Rm=arf_get_d(arb_midref(t2),ARF_RND_NEAR);
    arb_mul(y,b0,x,pr); arb_neg(y,y); arb_sin_cos(t1,t2,y,pr); double complex E0=arf_get_d(arb_midref(t2),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(t1),ARF_RND_NEAR);
    double x0=arf_get_d(arb_midref(x),ARF_RND_NEAR);
    double complex S=GF[k0]*E0*(sc*G0); double SA=cabs(GF[k0])*fabs(sc*G0), L1=fabs(sc*G0);
    for(int side=0;side<2;side++){ double G=G0,R=side?Rm:Rp,sgn=1; double complex E=E0, w=side?conj(wd):wd;
      for(slong j=1;j<=K;j++){ G*=R; R*=Qd; E=cmul(E,w); sgn=-sgn; double xj=side?(x0-j*dd):(x0+j*dd); double Kj=sgn*s*G/(Ad*xj);
        double complex f=GF[side?k0+j:k0-j]; S+=cmul(f,E)*Kj; SA+=cabs(f)*fabs(Kj); L1+=fabs(Kj);} }
    double Eint=GS0*(Eal+Etr) + sqrt(2.0)*GRF*L1*1.01 + 2*(eps_t+gam*(1+eps_t))*SA*1.01;
    arb_set_d(acb_realref(z),creal(S)); arb_set_d(acb_imagref(z),cimag(S)); arb_set_d(t1,Eint); arb_add_error(acb_realref(z),t1); arb_add_error(acb_imagref(z),t1);
    arb_set_d(v,ud-uc); arb_set_d(t1,dr); arb_mul_si(t1,t1,q,pr); arb_set_d(t2,o); arb_sub_d(t2,t2,uc,pr); arb_add(v,t1,t2,pr); /* v=u-uc exact */
    arb_mul(th,co[2],v,pr); arb_add(th,th,co[1],pr); arb_mul(th,th,v,pr); arb_add(th,th,co[0],pr);
    arb_abs(t1,v); arb_pow_ui(t1,t1,3,pr); arb_mul(t1,t1,cb[3],pr); arb_add_error(th,t1);
    arb_sin_cos(acb_imagref(e1),acb_realref(e1),th,pr); acb_mul(z,z,e1,pr); arb_mul_2exp_si(zr,acb_realref(z),1);
    arb_add(t1,c,v,pr); int ok=c0term(cc,t1,GNN,pr); arb_add(zr,zr,cc,pr); arb_add_error(zr,gg);
    sr[q]=(!ok||arb_contains_zero(zr))?0:(arb_is_positive(zr)?1:-1); NINTERP++;
    if(BLF) fprintf(BLF,"%.10f %.17g %.6g\n",ud,arf_get_d(arb_midref(zr),ARF_RND_NEAR),mag_get_d(arb_radref(zr))); }
  arb_clear(v);arb_clear(th);arb_clear(zr);arb_clear(cc);arb_clear(ub_); acb_clear(z);acb_clear(e1);
  arb_clear(x);arb_clear(y);arb_clear(t1);arb_clear(t2);
}

static int cmp_pt(const void*x,const void*y){double a=((pt*)x)->u,b=((pt*)y)->u; return a<b?-1:a>b;}
/* usage: fftzr T delta P a b Lw [signfile] [dbl]  -- coarse grid, then automatic in-process refinement of wide gaps */
/* [LEM BL-DP, candidate] binary64 certified probe: Lemma BL sum + Lemma DP final step, no Arb per point.
   Evaluated as the BL interpolant at u'=k0*delta+x0d (x0d=fma(-k0,delta,u), |u'-u|<=ex=u|x0d|), plus |F(u)-F(u')|<=S0*logN*ex.
   Seeds: sinc and exp(-w) by Taylor (deg 24 / 12) with Horner bound gamma_{2n+3}*e^{|z|}; cis via the Arb-certified TB table (16u+2.02u|phi|). */
static void zinterp_d(double Td,double u,signed char*sr,double*zout,double*Eout,double*Entout){
  static int init=0; static double Ad,bd,Qd,hd,dd,Eal,Etr,i2hd,SC[13],EC[13],logN; static double complex wd; static slong K;
  const double uu=ldexp(1,-53); const slong pr=128; *sr=0; if(zout)*zout=0;
  if(!init){ signed char d_; zinterp(Td,u,1e-3,1,&d_); /* forces BL parameter init (Arb) */ NINTERP--;
    arb_t A,b0,t,f; arb_init(A);arb_init(b0);arb_init(t);arb_init(f); dd=GDEL; arb_const_pi(A,pr); arb_div_d(A,A,dd,pr); Ad=arf_get_d(arb_midref(A),ARF_RND_NEAR);
    arb_log_ui(b0,GNN,pr); logN=arf_get_d(arb_midref(b0),ARF_RND_UP)*(1+1e-15); arb_mul_2exp_si(b0,b0,-1); bd=arf_get_d(arb_midref(b0),ARF_RND_NEAR);
    arb_sub(t,A,b0,pr); double Amb=arf_get_d(arb_midref(t),ARF_RND_DOWN); double tgt=getenv("BLTGT")?atof(getenv("BLTGT")):92.0;
    K=(slong)ceil(tgt/Amb/dd+0.5)+1; double X=(K-0.5)*dd; hd=sqrt(X/Amb); /* identical to zinterp */
    arb_set_d(t,hd); arb_mul(t,t,t,pr); arb_mul_2exp_si(t,t,1); arb_inv(t,t,pr); i2hd=arf_get_d(arb_midref(t),ARF_RND_NEAR);
    arb_set_d(f,dd); arb_mul(f,f,f,pr); arb_mul(f,f,t,pr); arb_mul_2exp_si(f,f,1); arb_neg(f,f); arb_exp(f,f,pr); Qd=arf_get_d(arb_midref(f),ARF_RND_NEAR);
    arb_mul_d(f,b0,-dd,pr); arb_t s_,c_; arb_init(s_);arb_init(c_); arb_sin_cos(s_,c_,f,pr); wd=arf_get_d(arb_midref(c_),ARF_RND_NEAR)+_Complex_I*arf_get_d(arb_midref(s_),ARF_RND_NEAR);
    for(int k=0;k<=12;k++){ arb_fac_ui(f,2*k+1,pr); arb_inv(f,f,pr); if(k&1)arb_neg(f,f); SC[k]=arf_get_d(arb_midref(f),ARF_RND_NEAR);
      arb_fac_ui(f,k,pr); arb_inv(f,f,pr); if(k&1)arb_neg(f,f); EC[k]=arf_get_d(arb_midref(f),ARF_RND_NEAR); }
    /* Eal,Etr copied from zinterp's computation */ Eal=BL_EAL; Etr=BL_ETR;
    arb_clear(A);arb_clear(b0);arb_clear(t);arb_clear(f);arb_clear(s_);arb_clear(c_); init=1; }
  /* block data (Arb once per 8-unit block) */
  static double cuc=-1e300,cTd=0,cod[3],cor[3],th3d,C0m,C0r,ggd; static int cok=0;
  double uc=floor(u/8.0)*8.0+4.0, H=4.0+1e-6;
  if(!(uc==cuc&&Td==cTd)){ arb_t co[4],cb[4],c,t1,t2,cc; for(int k=0;k<4;k++){arb_init(co[k]);arb_init(cb[k]);} arb_init(c);arb_init(t1);arb_init(t2);arb_init(cc);
    arf_t ub; arf_init(ub); arb_set_d(c,Td); arb_add_d(c,c,uc,pr); theta_ser(co,c,0,pr); theta_ser(cb,c,H,pr);
    arb_const_pi(t1,pr); arb_mul_2exp_si(t1,t1,1); arb_div(t2,co[0],t1,pr); arb_floor(t2,t2,pr); arb_mul(t2,t2,t1,pr); arb_sub(co[0],co[0],t2,pr);
    for(int k=0;k<3;k++){ cod[k]=arf_get_d(arb_midref(co[k]),ARF_RND_NEAR); arb_sub_d(t1,co[k],cod[k],pr); arb_get_abs_ubound_arf(ub,t1,pr); cor[k]=arf_get_d(ub,ARF_RND_UP); }
    arb_get_abs_ubound_arf(ub,cb[3],pr); th3d=arf_get_d(ub,ARF_RND_UP);
    arb_set(t1,c); mag_set_d(arb_radref(t1),H); cok=c0term(cc,t1,GNN,pr);
    if(cok){ C0m=arf_get_d(arb_midref(cc),ARF_RND_NEAR); arb_sub_d(t2,cc,C0m,pr); arb_get_abs_ubound_arf(ub,t2,pr); C0r=arf_get_d(ub,ARF_RND_UP); }
    arb_set_d(t1,Td); arb_add_d(t1,t1,uc-H,pr); gabcke(t2,t1,pr); arb_get_abs_ubound_arf(ub,t2,pr); ggd=arf_get_d(ub,ARF_RND_UP);
    for(int k=0;k<4;k++){arb_clear(co[k]);arb_clear(cb[k]);} arb_clear(c);arb_clear(t1);arb_clear(t2);arb_clear(cc); arf_clear(ub); cuc=uc; cTd=Td; }
  if(!cok) return;
  slong k0=(slong)llround(u/dd); if(k0-K<0||k0+K>=GP) return;
  double x0=fma(-(double)k0,dd,u), ex=1.01*uu*fabs(x0);
  /* sinc, sin */ double y=Ad*x0, z=y*y, sc=SC[12]; for(int k=11;k>=0;k--) sc=sc*z+SC[k];
  double ey=2.02*uu*fabs(y); double esc=(29*uu/(1-29*uu))*exp(fabs(y))+0.5*ey+1e-22; double s=y*sc; double es=fabs(y)*esc+uu*fabs(s)*1.01+ey*1.0;
  if(fabs(y)>1.6) return;
  /* exp seeds */ double w0=x0*x0*i2hd, wp=(2*x0+dd)*dd*i2hd, wm=(dd-2*x0)*dd*i2hd; if(w0>0.5||wp>0.5||wm>0.5||wm<0||wp<0) return;
  double G0=EC[12],Rp=EC[12],Rm=EC[12]; for(int k=11;k>=0;k--){G0=G0*w0+EC[k]; Rp=Rp*wp+EC[k]; Rm=Rm*wm+EC[k]; /* EC signed: sum (-1)^k w^k/k! = e^{-w} */}
  /* cis(-b x0) via table */ double ph=-bd*x0, f=ph*0.15915494309189535; f-=floor(f); double jj=nearbyint(f*1024.0), r=f-jj/1024.0, xx=r*6.283185307179586, x2=xx*xx;
  double cs=1-x2*(0.5-x2*(1.0/24-x2/720)), sn=xx*(1-x2*(1.0/6-x2*(1.0/120-x2/5040))); double complex tb=TB[(int)jj]; double complex E0=CMPLX(creal(tb)*cs+cimag(tb)*sn, creal(tb)*sn-cimag(tb)*cs);
  double complex S=cmul(GF[k0],E0)*(sc*G0); double SA=cabs(GF[k0])*fabs(sc*G0), L1=fabs(sc*G0), SB=0;
  for(int side=0;side<2;side++){ double G=G0,R=side?Rm:Rp,sgn=1; double complex E=E0, w=side?conj(wd):wd;
    for(slong j=1;j<=K;j++){ G*=R; R*=Qd; E=cmul(E,w); sgn=-sgn; double xj=side?(x0-j*dd):(x0+j*dd); double g_=G/(Ad*xj); double Kj=sgn*s*g_;
      double complex fz=GF[side?k0+j:k0-j]; S+=cmul(fz,E)*Kj; double af=cabs(fz); SA+=af*fabs(Kj); SB+=af*fabs(g_); L1+=fabs(Kj);} }
  double eps_t=((double)(K+3)*(K+3)+60.0*(K+1)+28.0+4.1*fabs(ph))*uu, gam=(2*K+2)*uu/(1-(2*K+2)*uu);
  double Efp=2*(eps_t+gam*(1+eps_t))*SA*1.01 + es*SB*1.02 + esc*cabs(GF[k0])*G0*1.02 + GS0*logN*ex;
  double Eint=Efp + GS0*(Eal+Etr) + sqrt(2.0)*GRF*L1*1.01; /* |F(u) - S| <= Eint (complex modulus) */
  /* Lemma DP final step at u */
  double vd=u-uc, ev=1.01*uu*fabs(vd), av=fabs(vd)+ev;
  double th_=cod[0]+vd*(cod[1]+vd*cod[2]);
  double eth=cor[0]+cor[1]*av+cor[2]*av*av+av*av*av*th3d + (fabs(cod[1])+cor[1])*ev+2*av*(fabs(cod[2])+cor[2])*ev*1.01 + (4*uu/(1-4*uu))*(fabs(cod[0])+av*fabs(cod[1])+av*av*fabs(cod[2])); eth*=1.25;
  f=th_*0.15915494309189535; f-=floor(f); jj=nearbyint(f*1024.0); r=f-jj/1024.0; xx=r*6.283185307179586; x2=xx*xx;
  cs=1-x2*(0.5-x2*(1.0/24-x2/720)); sn=xx*(1-x2*(1.0/6-x2*(1.0/120-x2/5040))); tb=TB[(int)jj];
  double ec=creal(tb)*cs+cimag(tb)*sn, es_=creal(tb)*sn-cimag(tb)*cs, ecs=(16+3.15)*uu+2.02*uu*fabs(th_)+1e-300;
  double Fr=creal(S),Fi=cimag(S), a_=Fr*ec, b_=Fi*es_, d_=a_-b_, zr=2*d_+C0m, aF=fabs(Fr)+fabs(Fi);
  double Ent=2*aF*(eth+ecs)*1.000001 + 2.02*uu*(fabs(a_)+fabs(b_)+fabs(d_))*1.01 + 1.01*uu*fabs(zr) + C0r + 2*Efp*(1+eth+ecs);
  double E=(Ent + 2*(Eint-Efp)*(1+eth+ecs) + ggd)*(1+1e-12);
  *sr=(fabs(zr)>E)?(zr>0?1:-1):0; NINTERP++; if(zout)*zout=zr; if(Eout)*Eout=E; if(Entout)*Entout=Ent;
}
/* 256-bit Arb truth for the adversarial test: 2Re(e^{i theta(T+u)} * sum_k GF[k] e^{-i b x_k} K(x_k)) + C0(T+u), same K and support */
static double zinterp_truth(double Td,double u){
  const slong pr=256; static int init=0; static arb_t A,b0,i2h; static slong K;
  if(!init){ arb_init(A);arb_init(b0);arb_init(i2h); arb_const_pi(A,pr); arb_div_d(A,A,GDEL,pr); arb_log_ui(b0,GNN,pr); arb_mul_2exp_si(b0,b0,-1);
    arb_t t; arb_init(t); arb_sub(t,A,b0,pr); double Amb=arf_get_d(arb_midref(t),ARF_RND_DOWN); double tgt=getenv("BLTGT")?atof(getenv("BLTGT")):92.0;
    K=(slong)ceil(tgt/Amb/GDEL+0.5)+1; double X=(K-0.5)*GDEL, hd=sqrt(X/Amb); arb_set_d(i2h,hd); arb_mul(i2h,i2h,i2h,pr); arb_mul_2exp_si(i2h,i2h,1); arb_inv(i2h,i2h,pr); arb_clear(t); init=1; }
  slong k0=(slong)llround(u/GDEL); acb_t S,tm,e; arb_t x,y,g,th[4],tt,cc,z; acb_init(S);acb_init(tm);acb_init(e); arb_init(x);arb_init(y);arb_init(g);arb_init(tt);arb_init(cc);arb_init(z); for(int k=0;k<4;k++)arb_init(th[k]);
  for(slong k=k0-K;k<=k0+K;k++){ arb_set_d(x,u); arb_set_d(y,GDEL); arb_mul_si(y,y,k,pr); arb_sub(x,x,y,pr);
    arb_mul(y,A,x,pr); arb_sinc(g,y,pr); arb_mul(y,x,x,pr); arb_mul(y,y,i2h,pr); arb_neg(y,y); arb_exp(y,y,pr); arb_mul(g,g,y,pr);
    arb_mul(y,b0,x,pr); arb_neg(y,y); arb_sin_cos(acb_imagref(e),acb_realref(e),y,pr); acb_set_d_d(tm,creal(GF[k]),cimag(GF[k])); acb_mul(tm,tm,e,pr); acb_mul_arb(tm,tm,g,pr); acb_add(S,S,tm,pr); }
  arb_set_d(tt,Td); arb_add_d(tt,tt,u,pr); theta_ser(th,tt,0,pr); arb_sin_cos(acb_imagref(e),acb_realref(e),th[0],pr); acb_mul(S,S,e,pr);
  arb_mul_2exp_si(z,acb_realref(S),1); c0term(cc,tt,GNN,pr); arb_add(z,z,cc,pr); double r=arf_get_d(arb_midref(z),ARF_RND_NEAR);
  acb_clear(S);acb_clear(tm);acb_clear(e); arb_clear(x);arb_clear(y);arb_clear(g);arb_clear(tt);arb_clear(cc);arb_clear(z); for(int k=0;k<4;k++)arb_clear(th[k]); return r; }

static int USE_DP=1; static long long DPT_N=0,DPT_BAD=0,DPF_N=0; static double DPT_MAX=0;
static void probe(double Td,double u,double dr_,signed char*sr){ if(!USE_DP){zinterp(Td,u,dr_,1,sr);return;} double z=0,E=0,Ent=0; zinterp_d(Td,u,sr,&z,&E,&Ent);
  if(!*sr){ zinterp(Td,u,dr_,1,sr); DPF_N++; } /* fallback to Arb probe */
  if(getenv("DPTEST")&&Ent>0){ double tr=zinterp_truth(Td,u); double rr=fabs(z-tr)/Ent; if(rr>DPT_MAX)DPT_MAX=rr; if(rr>1)DPT_BAD++; DPT_N++; } }
/* v8 (Nathan A. Meeks' delta cache extender): extend AH/AL only for n=G_AN+1..Nt. Each entry is computed by the exact
   same Arb ops/precision as the fresh build in zgrid_c (arb_log_ui, /2pi, prec=PREC), so entries are bit-identical;
   ARAD = max over all stored n, which equals a fresh build's max over the same range (still a valid bound for every n<=N used). */
static double AHX_SEC=0; static slong AHX_N=0;
static int ah_extend(slong Nt){ if(G_AN>=Nt) return 0; double t0=now(); slong prec=PREC; arb_t tp,L,an,r2; arb_init(tp);arb_init(L);arb_init(an);arb_init(r2);
  arb_const_pi(tp,prec); arb_mul_2exp_si(tp,tp,1); double uu=ldexp(1,-53); AH=realloc(AH,Nt*sizeof(double)); AL=realloc(AL,Nt*sizeof(double));
  double ar=G_AN?ARAD:0; for(slong n=G_AN+1;n<=Nt;n++){ arb_log_ui(L,n,prec); arb_div(an,L,tp,prec); double h=arf_get_d(arb_midref(an),ARF_RND_NEAR);
    arb_sub_d(r2,an,h,prec); double l=arf_get_d(arb_midref(r2),ARF_RND_NEAR); AH[n-1]=h; AL[n-1]=l; double rr=mag_get_d(arb_radref(r2))+fabs(l)*uu; if(rr>ar)ar=rr; }
  AHX_N=Nt-G_AN; ARAD=ar; G_AN=Nt; arb_clear(tp);arb_clear(L);arb_clear(an);arb_clear(r2); AHX_SEC=now()-t0; return 1; }
int main(int argc,char**argv){
  double Td=atof(argv[1]),delta=atof(argv[2]); slong P=atol(argv[3]); double A=atof(argv[4]),B=atof(argv[5]),Lw=atof(argv[6]);
  slong prec=PREC; if(argc>8) USE_DBL=1; double t0=now(); if(getenv("NODP")) USE_DP=0;
  signed char*sg=malloc(P); slong N,ncert; int Lt; double tm[4];
  const char*ahc=getenv("AHCACHE"); int ahl=0; slong Nneed=(slong)floor(sqrt(Td/(2*M_PI)));
  if(ahc&&USE_DBL){ ahl=ah_load(ahc,Nneed); if(ahl!=1&&!getenv("NOAHX")){ /* [search-only] anticipate next sliding windows: slack AHXS (fraction of T) + chunk rounding AHXC */
      double sl=getenv("AHXS")?atof(getenv("AHXS")):0.0; slong ch=getenv("AHXC")?atol(getenv("AHXC")):65536; slong Nt=(slong)floor(sqrt(Td*(1+sl)/(2*M_PI))); if(ch>1) Nt=((Nt+ch-1)/ch)*ch; if(Nt<Nneed)Nt=Nneed;
      if(ah_extend(Nt)) { ah_save(ahc); ahl=1; fprintf(stderr,"AHX: extended +%ld to %ld in %.3fs\n",(long)AHX_N,(long)G_AN,AHX_SEC);} } }
  zgrid(Td,delta,P,sg,&N,&Lt,&ncert,tm);
  if(ahc&&USE_DBL&&ahl!=1) ah_save(ahc);
  if(getenv("GZDUMP")){FILE*f=fopen(getenv("GZDUMP"),"w"); for(slong p=0;p<P;p++) fprintf(f,"%.10f %.17g\n",p*delta,GZ[p]); fclose(f); return 0;}
  if(getenv("BLDUMP")){BLF=fopen(getenv("BLDUMP"),"w"); signed char*sr=malloc(20000); zinterp(Td,60.0,0.01,20000,sr); fclose(BLF); return 0;}
  slong M=1; while(M<2*P) M<<=1;
  arb_t T,pi,tmp,c; arb_init(T);arb_init(pi);arb_init(tmp);arb_init(c); arb_set_d(T,Td); arb_const_pi(pi,prec);
  arb_t co[4]; for(int k=0;k<4;k++)arb_init(co[k]);
  slong np=0,cap=P+(1<<16); pt*pts=malloc(cap*sizeof(pt));
  for(slong p=0;p<P;p++) if(sg[p]){pts[np].u=p*delta;pts[np].s=sg[p];np++;}
  /* ---- speed2 dip pass: targeting is heuristic (parabola through |Z| at 3 grid points); every added sign is certified via Lemma BL ---- */
  int useI=!getenv("NOINTERP"); double tdip0=now(); slong ndip=0; int ND=getenv("NDIV")?atoi(getenv("NDIV")):8; double tau=getenv("DIPTAU")?atof(getenv("DIPTAU")):3.0;
  int useE=useI&&!getenv("NOEXT"); double ktau=getenv("KTAU")?atof(getenv("KTAU")):1.0;
  if(useE){ /* [HEUR targeting] Z' at grid points from BL kernel derivative: F'(u_p)=-i b F_p + sum_{j!=0} F_{p-j} e^{-i b j d}(-1)^j G(jd)/(jd) */
    double b=0.5*log((double)GNN), Amb=M_PI/delta-b; slong KK=(slong)ceil(92.0/Amb/delta+0.5)+1; double X=(KK-0.5)*delta, h2=X/Amb;
    double complex*cw=malloc((2*KK+1)*sizeof(double complex)); for(slong j=-KK;j<=KK;j++) cw[j+KK]=j?cexp(-_Complex_I*b*j*delta)*((j&1)?-1.0:1.0)*exp(-(j*delta)*(j*delta)/(2*h2))/(j*delta):0;
    int scr5=getenv("SCR5")?atoi(getenv("SCR5")):(delta>0.15); double kt5=getenv("KT5")?atof(getenv("KT5")):0.5; double kflat=getenv("KFLAT")?atof(getenv("KFLAT")):0.1; double kgate=getenv("KGATE")?atof(getenv("KGATE")):0.25; slong ngate=0; double*Zpp=calloc(P,sizeof(double));
    double complex*cw2=malloc((2*KK+1)*sizeof(double complex)); /* [HEUR] F'' at grid points: d^2/du^2 of e^{-ibx}K(x) at x=j*delta */
    for(slong j=-KK;j<=KK;j++){ double x=j*delta,G=exp(-x*x/(2*h2)),sg_=(j&1)?-1.0:1.0;
      cw2[j+KK]= j? cexp(-_Complex_I*b*x)*(-2.0*_Complex_I*b*sg_*G/x + sg_*G*(-2.0/(x*x)-2.0/h2)) : (-b*b - M_PI*M_PI/(3*delta*delta) - 1.0/h2); }
    double*Zp=calloc(P,sizeof(double)); slong ncand=0; const char*hdbg=getenv("HDBG"); int scr=getenv("SCREEN")?atoi(getenv("SCREEN")):(delta>0.15); double nzt=getenv("NZT")?atof(getenv("NZT")):0.1;
    for(slong p=KK;p+KK<P;p++){ double u=p*delta; if(u<A-Lw-2*delta||u>B+Lw+2*delta) continue; double complex Fp=-_Complex_I*b*GF[p];
      for(slong j=-KK;j<=KK;j++) if(j) Fp+=cw[j+KK]*GF[p-j]; Zp[p]=2*creal(GE[p]*(_Complex_I*GTP[p]*GF[p]+Fp));
      if(scr5){ double complex Fpp=0; for(slong j=-KK;j<=KK;j++) Fpp+=cmul(cw2[j+KK],GF[p-j]); double tp=GTP[p];
        Zpp[p]=2*creal(GE[p]*(-tp*tp*GF[p]+2.0*_Complex_I*tp*Fp+Fpp)); } }
    free(cw); free(cw2); if(scr5&&!getenv("NZT")) nzt=0; /* v5: quintic(Z,Z',Z'') + gated 9-sample screen replaces the near-end rule */ signed char sr1[1]; slong cap2=0; double tz0=now(), tpr=0, tsm=0;
    slong nint_=0, nfast=0, nfasthit=0; for(slong p=KK;p+1+KK<P;p++){ double u=p*delta; if(u<A-Lw-delta||u>B+Lw+delta) continue; nint_++;
      double z0=GZ[p],z1=GZ[p+1],m0=delta*Zp[p],m1=delta*Zp[p+1]; double H[33]; int nsc=0; double tl[33]; int nl=0;
      double c0q=delta*delta*Zpp[p], c1q=delta*delta*Zpp[p+1];
      static double QB[33][6], HB[33][4]; static int qbi=0; /* v6: basis tables, same expressions => bit-identical H */
      if(!qbi){ for(int i=0;i<=32;i++){ double t=i/32.0,t2=t*t,t3=t2*t,t4=t3*t,t5=t4*t;
          QB[i][0]=(1-10*t3+15*t4-6*t5); QB[i][1]=(t-6*t3+8*t4-3*t5); QB[i][2]=0.5*(t2-3*t3+3*t4-t5); QB[i][3]=(10*t3-15*t4+6*t5); QB[i][4]=(-4*t3+7*t4-3*t5); QB[i][5]=0.5*(t3-2*t4+t5);
          HB[i][0]=(2*t3-3*t2+1); HB[i][1]=(t3-2*t2+t); HB[i][2]=(-2*t3+3*t2); HB[i][3]=(t3-t2); } qbi=1; }
      for(int i=0;i<=32;i++){
        H[i]= scr5? QB[i][0]*z0+QB[i][1]*m0+QB[i][2]*c0q+QB[i][3]*z1+QB[i][4]*m1+QB[i][5]*c1q
                  : HB[i][0]*z0+HB[i][1]*m0+HB[i][2]*z1+HB[i][3]*m1;
        if(i&&((H[i]<0)!=(H[i-1]<0))){ tl[nl++]=(i-0.5)/32.0; nsc++; } }
      int gridsc=(sg[p]&&sg[p+1]&&sg[p]!=sg[p+1])?1:0; double pts_t[16]; int npt=0;
      int cand=(sg[p]==0||sg[p+1]==0||nsc>gridsc);
      if(!cand&&scr&&(fabs(z0)<nzt*fabs(m0)||fabs(z1)<nzt*fabs(m1))) cand=1; /* [HEUR] a zero within ~nzt*delta of a grid end: a hidden partner may sit right next to it */
      if(!cand){ double mz=gridsc?fmax(fabs(z0),fabs(z1)):fmin(fabs(z0),fabs(z1)); for(int i=1;i<32;i++) if(fabs(H[i])<=fabs(H[i-1])&&fabs(H[i])<=fabs(H[i+1])&&((H[i-1]<0)==(H[i]<0))&&((H[i+1]<0)==(H[i]<0))&&fabs(H[i])<(scr5?kt5:ktau)*mz){cand=1;break;} }
      if(!cand&&scr5){ double mh=0; for(int i=0;i<=32;i++) if(fabs(H[i])>mh) mh=fabs(H[i]); double loc=fmax(fabs(GZ[p-1]),fabs(GZ[p+2])); if(mh<kflat*loc) cand=1; } /* [HEUR] |Z| flat near 0 across the whole interval */
      enum{na=64}; double za[na+1]; static double complex*WT=NULL; static double *WTr,*WTi,*GFr,*GFi; slong WM=2*KK+2;
        if(!WT){ double Ad=M_PI/delta; WT=malloc((na+1)*WM*sizeof(double complex));
          for(int i=0;i<=na;i++) for(slong m=-KK;m<=KK+1;m++){ double x=i*delta/na-m*delta; double Kx=(fabs(x)<1e-12)?1.0:sin(Ad*x)/(Ad*x); WT[i*WM+m+KK]=GF?cexp(-_Complex_I*b*x)*Kx*exp(-x*x/(2*h2)):0; }
          WTr=malloc((na+1)*WM*8); WTi=malloc((na+1)*WM*8); for(slong q=0;q<(na+1)*WM;q++){WTr[q]=creal(WT[q]);WTi[q]=cimag(WT[q]);}
          GFr=malloc(P*8); GFi=malloc(P*8); for(slong q=0;q<P;q++){GFr[q]=creal(GF[q]);GFi[q]=cimag(GF[q]);} }
      int gate=1; if(scr5&&scr){ double mn=1e300; for(int i=0;i<=32;i++) if(fabs(H[i])<mn) mn=fabs(H[i]); double loc=fmax(fmax(fabs(GZ[p-1]),fabs(GZ[p+2])),fmax(fabs(z0),fabs(z1))); if(gridsc){ mn=fmin(fabs(z0),fabs(z1)); for(int i=1;i<32;i++) if(fabs(H[i])<=fabs(H[i-1])&&fabs(H[i])<=fabs(H[i+1])&&((H[i-1]<0)==(H[i]<0))&&((H[i+1]<0)==(H[i]<0))&&fabs(H[i])<mn) mn=fabs(H[i]); }
        gate=(mn<kgate*loc); /* crossing itself excluded for gridsc=1 */ if(gate) ngate++; }
      double zraw8[9]; int have8=0; /* v6: raw 9-sample values reused by the fine sampler (identical formula, same p) */
      if(!cand&&scr&&gate){ /* [HEUR] coarse approximate screen: 9 BL samples per interval */ double zs8[9];
        for(int ii=0;ii<=8;ii++){ int i=ii*8; double complex S=mac4(GFr+p-KK,GFi+p-KK,WTr+i*WM,WTi+i*WM,WM);
          zs8[ii]=zraw8[ii]=2*creal(GE[p]*cexp(_Complex_I*GTP[p]*(i*delta/na))*S); } have8=1;
        { double c0_=z0-zs8[0], c1_=z1-zs8[8]; for(int ii=0;ii<=8;ii++) zs8[ii]+=c0_*(1-ii/8.0)+c1_*(ii/8.0); }
        int n8=0; for(int ii=1;ii<=8;ii++) if((zs8[ii]<0)!=(zs8[ii-1]<0)) n8++; if(n8>gridsc) cand=1;
        double mz=gridsc?fmax(fabs(z0),fabs(z1)):fmin(fabs(z0),fabs(z1));
        for(int ii=1;ii<8&&!cand;ii++) if(fabs(zs8[ii])<=fabs(zs8[ii-1])&&fabs(zs8[ii])<=fabs(zs8[ii+1])&&((zs8[ii-1]<0)==(zs8[ii]<0))&&((zs8[ii+1]<0)==(zs8[ii]<0))&&fabs(zs8[ii])<ktau*mz) cand=1; }
      double tq0=now(); static int fastH=-1; if(fastH<0) fastH=getenv("NOFASTH")?0:1;
      if(cand&&fastH&&gridsc==0&&nsc==2&&sg[p]&&sg[p+1]){ /* v7 [HEUR targeting]: quintic shows a clean pair -> certified probe at its midpoint first; accept only if the certified sign is opposite to both grid ends, else fall back to the fine sampler */
        double uu_=(p+0.5*(tl[0]+tl[1]))*delta; double tq1=now(); probe(Td,uu_,delta,sr1); ndip++; tpr+=now()-tq1; nfast++;
        if(sr1[0]&&sr1[0]!=sg[p]){ if(np+1>cap){cap=2*(np+1);pts=realloc(pts,cap*sizeof(pt));} pts[np].u=uu_;pts[np].s=sr1[0];np++; nfasthit++; cand=0; }
        else if(sr1[0]){ if(np+1>cap){cap=2*(np+1);pts=realloc(pts,cap*sizeof(pt));} pts[np].u=uu_;pts[np].s=sr1[0];np++; } }
      if(cand){ ncand++; /* [HEUR] approximate (uncertified) BL samples to locate roots; only certified evaluations count */
        /* two-level: 17 coarse samples, then full resolution only in segments near sign changes / |z| minima / small ends [HEUR] */
        char fl[na+1]; memset(fl,0,sizeof fl);
        #define ZA_AT(i) do{ double s_=(i)*delta/na; double complex S=mac4(GFr+p-KK,GFi+p-KK,WTr+(i)*WM,WTi+(i)*WM,WM); \
          za[i]=2*creal(GE[p]*cexp(_Complex_I*GTP[p]*s_)*S); fl[i]=1; }while(0)
        for(int i=0;i<=na;i+=4){ if(have8&&(i%8)==0){ za[i]=zraw8[i/8]; fl[i]=1; } else ZA_AT(i); }
        double c0_=GZ[p]-za[0], c1_=GZ[p+1]-za[na];
        for(int i=0;i<=na;i+=4) za[i]+=c0_*(1-(double)i/na)+c1_*((double)i/na);
        char rs[na/4]; memset(rs,0,sizeof rs);
        for(int g=0;g<na/4;g++){ double a0=za[4*g],a1=za[4*g+4]; if((a0<0)!=(a1<0)) rs[g]=1; }
        for(int g=0;g<=na/4;g++){ int i=4*g; double ai=fabs(za[i]); int lm=(g==0||ai<=fabs(za[i-4]))&&(g==na/4||ai<=fabs(za[i+4]));
          if(lm){ if(g>0)rs[g-1]=1; if(g<na/4)rs[g]=1; } }
        for(int g=0;g<na/4;g++) if(rs[g]) for(int i=4*g+1;i<4*g+4;i++){ ZA_AT(i); za[i]+=c0_*(1-(double)i/na)+c1_*((double)i/na); }
        #undef ZA_AT
        double rl[65]; int nr=0, lastI=0; for(int i=1;i<=na;i++){ if(!fl[i]) continue; if((za[i]<0)!=(za[lastI]<0)) rl[nr++]=0.5*(i+lastI)/na; lastI=i; }
        if(nr>gridsc) for(int r=0;r+1<nr&&npt<16;r++) pts_t[npt++]=0.5*(rl[r]+rl[r+1]);
        else if(nsc>gridsc) for(int r=0;r+1<nl&&npt<16;r++) pts_t[npt++]=0.5*(tl[r]+tl[r+1]); }
      if(hdbg){double c_=atof(hdbg); if(fabs(u-c_)<0.4){double mn=1e9; for(int i=0;i<=32;i++) if(fabs(H[i])<mn)mn=fabs(H[i]); fprintf(stderr,"u=%.3f z0=%.4g z1=%.4g m0=%.4g m1=%.4g nsc=%d g=%d npt=%d minH=%.3g\n",u,z0,z1,m0,m1,nsc,gridsc,npt,mn);}}
      double tq1=now(); tsm+=tq1-tq0; for(int q=0;q<npt;q++){ double uu_=(p+pts_t[q])*delta; probe(Td,uu_,delta,sr1); ndip++;
        if(sr1[0]){ if(np+1>cap){cap=2*(np+1);pts=realloc(pts,cap*sizeof(pt));} pts[np].u=uu_;pts[np].s=sr1[0];np++; } } tpr+=now()-tq1; }
    free(Zp); (void)cap2; free(Zpp); fprintf(stderr,"fastH: tried=%ld hit=%ld\n",nfast,nfasthit); fprintf(stderr,"ext: gated=%ld intervals=%ld candidates=%ld probes=%ld t_zp=%.3f t_loop=%.3f t_finesample=%.3f t_probe=%.3f\n",ngate,nint_,ncand,ndip,tz0-tdip0,now()-tz0,tsm,tpr); }
  if(useI&&!useE){ signed char*sr=malloc(2*ND+1);
    for(slong p=1;p+1<P;p++){ double u=p*delta; if(u<A-Lw-delta||u>B+Lw+delta) continue; int cand=0;
      if(sg[p]==0) cand=1;
      else if(sg[p-1]==sg[p]&&sg[p+1]==sg[p]){ double a_=fabs(GZ[p-1]),b_=fabs(GZ[p]),c_=fabs(GZ[p+1]);
        if(b_<a_&&b_<=c_){ double cu_=a_-2*b_+c_; double vmin=b_-(c_-a_)*(c_-a_)/(8*cu_); if(vmin<tau*b_) cand=1; } }
      if(!cand) continue; ndip++;
      zinterp(Td,(p-1)*delta,delta/ND,2*ND+1,sr);
      if(np+2*ND+1>cap){cap=2*(np+2*ND+1);pts=realloc(pts,cap*sizeof(pt));}
      for(int q=0;q<=2*ND;q++) if(sr[q]){pts[np].u=(p-1)*delta+q*(delta/ND);pts[np].s=sr[q];np++;} }
    free(sr); }
  qsort(pts,np,sizeof(pt),cmp_pt);
  /* speed3: narrow sign-change intervals inside the Turing zones by certified bisection (tightens the Turing integrals; targeting only) */
  if(useI){ double TW=getenv("TZW")?atof(getenv("TZW")):0.02; signed char sr1[1];
    { double bd[4]={A-Lw,A,B,B+Lw}; for(int q=0;q<4;q++) for(int tr_=0;tr_<8;tr_++){ double m=bd[q]+(tr_?((tr_&1)?1:-1)*1e-3*((tr_+1)/2):0); /* certified point at each zone boundary */
        probe(Td,m,1e-3,sr1); ndip++; if(sr1[0]){ if(np+1>cap){cap=2*(np+1);pts=realloc(pts,cap*sizeof(pt));} pts[np].u=m;pts[np].s=sr1[0];np++; break; } }
      qsort(pts,np,sizeof(pt),cmp_pt); }
    for(int rd=0;rd<6;rd++){ slong add=0, np0=np;
      for(slong i=1;i<np0;i++){ double l=pts[i-1].u,r=pts[i].u; if(pts[i].s==pts[i-1].s||r-l<=TW) continue;
        if(!((l>=B&&r<=B+Lw)||(l>=A-Lw&&r<=A))) continue;
        double m=0.5*(l+r); probe(Td,m,TW,sr1); ndip++; if(sr1[0]){ if(np+1>cap){cap=2*(np+1);pts=realloc(pts,cap*sizeof(pt));} pts[np].u=m;pts[np].s=sr1[0];np++;add++; } }
      if(!add) break; qsort(pts,np,sizeof(pt),cmp_pt); } }
  double tdip=now()-tdip0;
  slong nz=0,dis0=0; double *zs=malloc(np*sizeof(double)),*ze=malloc(np*sizeof(double)); { double lu=0; int ls=0;
  for(slong i=0;i<np;i++){ if(i>0&&pts[i].u-pts[i-1].u<1e-9){ if(pts[i].s!=pts[i-1].s) dis0++; continue; }
    if(ls&&pts[i].s!=ls){zs[nz]=lu;ze[nz]=pts[i].u;nz++;} lu=pts[i].u; ls=pts[i].s; } }
  /* Turing: N(b)<=1+(C(b+L)-In+Thint(b,b+L))/L ; N(a)>=1+(-C(a)+Im+Thint(a-L,a))/L */
  arb_t Ub,La,In,Im,I,R,h,Cb; arb_init(Ub);arb_init(La);arb_init(In);arb_init(Im);arb_init(I);arb_init(R);arb_init(h);arb_init(Cb);
  arb_zero(In);arb_zero(Im); slong found=0;
  for(slong j=0;j<nz;j++){ if(zs[j]>=B&&ze[j]<=B+Lw){arb_set_d(tmp,B+Lw-ze[j]);arb_add(In,In,tmp,prec);}
     if(zs[j]>=A-Lw&&ze[j]<=A){arb_set_d(tmp,zs[j]-(A-Lw));arb_add(Im,Im,tmp,prec);} if(zs[j]>=A&&ze[j]<=B)found++; }
  double ends[2][2]={{B,B+Lw},{A-Lw,A}};
  for(int s=0;s<2;s++){ double x0=ends[s][0],x1=ends[s][1],hh=(x1-x0)/2; arb_set_d(c,(x0+x1)/2); arb_add(c,c,T,prec);
     theta_ser(co,c,hh,prec); arb_t cp[4]; for(int k=0;k<4;k++)arb_init(cp[k]); theta_ser(cp,c,0,prec);
     arb_set_d(h,hh); arb_mul(I,cp[0],h,prec); arb_mul_2exp_si(I,I,1); arb_pow_ui(tmp,h,3,prec); arb_mul(tmp,tmp,cp[2],prec); arb_mul_ui(tmp,tmp,2,prec); arb_div_ui(tmp,tmp,3,prec); arb_add(I,I,tmp,prec);
     arb_abs(R,co[3]); arb_pow_ui(tmp,h,4,prec); arb_mul(R,R,tmp,prec); arb_mul_2exp_si(R,R,-1); arb_add_error(I,R); arb_div(I,I,pi,prec);
     arb_set_d(tmp,s==0?x1:x0+0); arb_add(tmp,tmp,T,prec); if(s==1){arb_set_d(tmp,A);arb_add(tmp,tmp,T,prec);} arb_log(Cb,tmp,prec); {arb_t q_; arb_init(q_); arb_set_ui(q_,59); arb_div_ui(q_,q_,1000,prec); arb_mul(Cb,Cb,q_,prec); arb_set_ui(q_,2067); arb_div_ui(q_,q_,1000,prec); arb_add(Cb,Cb,q_,prec); arb_clear(q_);} /* exact Trudgian constants (double 0.059 rounds DOWN) */
     if(s==0){arb_sub(Ub,Cb,In,prec);arb_add(Ub,Ub,I,prec);arb_div_d(Ub,Ub,Lw,prec);arb_add_ui(Ub,Ub,1,prec);}
     else{arb_neg(La,Cb);arb_add(La,La,Im,prec);arb_add(La,La,I,prec);arb_div_d(La,La,Lw,prec);arb_add_ui(La,La,1,prec);}
     for(int k=0;k<4;k++)arb_clear(cp[k]); }
  arf_t u_,l_; arf_init(u_);arf_init(l_); arb_get_ubound_arf(u_,Ub,prec); arb_get_lbound_arf(l_,La,prec);
  fmpz_t fu,flo; fmpz_init(fu);fmpz_init(flo); arf_get_fmpz(fu,u_,ARF_RND_FLOOR); arf_get_fmpz(flo,l_,ARF_RND_CEIL);
  fmpz_t mx; fmpz_init(mx); fmpz_sub(mx,fu,flo);
  double t4=now();
  slong mxs=fmpz_get_si(mx), found0=found;
  /* ---- in-process refinement. Heuristic only chooses WHERE to add points; every added sign is a certified Arb ball. ---- */
  double meansp=2*M_PI/log(Td/(2*M_PI));
  slong ng=0; double *gl=malloc(nz*sizeof(double)),*gr=malloc(nz*sizeof(double)),*gw=malloc(nz*sizeof(double)); char*done=calloc(nz,1);
  for(slong j=0;j+1<nz;j++) if(zs[j]>=A&&ze[j+1]<=B){gl[ng]=ze[j];gr[ng]=zs[j+1];gw[ng]=gr[ng]-gl[ng];ng++;}
  /* ---- deficit-targeted refinement (CP). Targeting is heuristic; soundness rests only on certified signs + Turing max. ---- */
  found=found0; slong nrefine=0, extra=0, disagree=0; double tref=0; int rounds=0; slong budget=(argc>9)?atol(argv[9]):4000000;
  double dr=delta/10; const slong PR=getenv("RPR")?atol(getenv("RPR")):256; zcache Cc={0}; int stopreason=0; /*0 proven,1 budget,2 dr floor,3 overcount*/
  double lr=log(Td/(2*M_PI))/(2*M_PI); /* local zero density, heuristic */
  slong maxrounds=useI?(getenv("MAXR")?atol(getenv("MAXR")):20):200; /* interp: fail fast, caller retries at delta/2 */ char*tried=calloc(ng?ng:1,1);
  while(found<mxs){
    if(found>mxs){stopreason=3;break;}
    if(extra>=budget||rounds>=maxrounds){stopreason=1;break;}
    rounds++;
    /* current sign-change list from merged points */
    slong nc=0; double *cl=malloc(np*sizeof(double)),*cr=malloc(np*sizeof(double));
    for(slong i=1;i<np;i++) if(pts[i-1].u>=A&&pts[i].u<=B&&pts[i].s!=pts[i-1].s){cl[nc]=pts[i-1].u;cr[nc]=pts[i].u;nc++;}
    /* deficit D_j = j - expected(A..mid_j); step score at gap j = mean D over K after - mean D over K before */
    slong K=40; double *D=malloc((nc+1)*sizeof(double));
    for(slong j=0;j<nc;j++){double m=0.5*(cl[j]+cr[j]); D[j]=j-(m-A)*lr*(1+0.0);} 
    ng=0; gl=realloc(gl,(nc+1)*sizeof(double)); gr=realloc(gr,(nc+1)*sizeof(double)); gw=realloc(gw,(nc+1)*sizeof(double));
    for(slong j=0;j+1<nc;j++){ double sb=0,sa=0; slong kb=0,ka=0;
      for(slong k=j;k>=0&&k>j-K;k--){sb+=D[k];kb++;} for(slong k=j+1;k<nc&&k<=j+K;k++){sa+=D[k];ka++;}
      double step=(kb&&ka)?(sb/kb-sa/ka):0; double Eg=(cl[j+1]-cr[j])*lr;
      gl[ng]=cr[j]; gr[ng]=cl[j+1]; gw[ng]=step+0.25*Eg; ng++; }
    /* also gaps at the window ends (before first / after last change) */
    free(D); free(cl); free(cr);
    /* pick top gaps not yet resolved at this dr: want = deficit + 4 */
    slong want=useI?8*(mxs-found)+64:(mxs-found)+4, picked=0; /* interp points are cheap: refine more gaps per round */ double *pl=malloc(want*sizeof(double)),*pr_=malloc(want*sizeof(double));
    for(slong w=0;w<want;w++){ slong bi=-1;
      for(slong g=0;g<ng;g++){ if(gw[g]<-1e9) continue; /* skip gaps already finer than dr: width<=2dr */
        if(gr[g]-gl[g]<=2*dr*1.0001) continue;
        /* skip gap already refined at current dr: detect via point density (has interior points spaced <=dr) */
        slong lo=0,hi=np; while(lo<hi){slong m=(lo+hi)/2; if(pts[m].u<=gl[g])lo=m+1; else hi=m;}
        if(lo<np && pts[lo].u<gr[g] && pts[lo].u-gl[g]<=dr*1.0001) continue;
        if(bi<0||gw[g]>gw[bi]) bi=g; }
      if(bi<0) break; pl[picked]=gl[bi]; pr_[picked]=gr[bi]; picked++; gw[bi]=-2e9; }
    if(picked==0){ if(dr/10<delta*1e-4){stopreason=2;free(pl);free(pr_);break;}
      dr/=10; if(Cc.init){_arb_vec_clear(Cc.l,Cc.N);_arb_vec_clear(Cc.w,Cc.N);_arb_vec_clear(Cc.x,Cc.N);_arb_vec_clear(Cc.e,Cc.N);free(Cc.mi);arb_clear(Cc.S0);Cc.init=0;}
      free(pl);free(pr_); continue; }
    for(slong w=0;w<picked;w++){
      double o=floor(pl[w]); slong PRw=PR; if(useI){o=pl[w]; PRw=(slong)ceil((pr_[w]-pl[w])/dr)+1;} else while(PRw*dr<ceil(pr_[w])-o) PRw<<=1; /* integer start: T+o exact in double */
      { signed char*sr=malloc(PRw); slong Nr,ncr; int Lr; double tr[4]; double s0=now();
        MF=getenv("RMF")?atoi(getenv("RMF")):16; if(useI){zinterp(Td,o,dr,PRw,sr);Lr=0;tr[0]=tr[1]=tr[2]=0;} else zgrid_c(&Cc,Td+o,dr,PRw,sr,&Nr,&Lr,&ncr,tr); MF=2; if(getenv("RDBG"))fprintf(stderr,"ref P=%ld L=%d %.3f %.3f %.3f\n",PRw,Lr,tr[0],tr[1],tr[2]); tref+=now()-s0; nrefine++;
        if(np+PRw>cap){cap=2*(np+PRw);pts=realloc(pts,cap*sizeof(pt));}
        for(slong q=0;q<PRw;q++) if(sr[q]){ double u=o+q*dr; if(u<A||u>B) continue;
          pts[np].u=u;pts[np].s=sr[q];np++;extra++; }
        free(sr); }
    }
    free(pl);free(pr_);
    qsort(pts,np,sizeof(pt),cmp_pt);
    found=0; disagree=0; for(slong i=1;i<np;i++){ if(pts[i].u-pts[i-1].u<1e-9){ if(pts[i].s!=pts[i-1].s) disagree++; continue; }
      if(pts[i-1].u>=A&&pts[i].u<=B&&pts[i].s!=pts[i-1].s) found++; }
    if(disagree){stopreason=4;break;}
    if(found>mxs){stopreason=3;break;}
  }
  free(tried); disagree+=dis0;
  double t5=now(); if(USE_DP) fprintf(stderr,"DP probes: arb_fallbacks=%lld\n",DPF_N); if(getenv("DPTEST")) fprintf(stderr,"DPTEST n=%lld max|err|/bound=%.3g violations=%lld\n",DPT_N,DPT_MAX,DPT_BAD); if(getenv("DBLTEST")) fprintf(stderr,"DBLTEST n=%lld max|err|/bound=%.3g violations=%lld\n",DT_N,DT_MAX,DT_BAD);
  if(getenv("ZDUMP")){FILE*f=fopen(getenv("ZDUMP"),"w"); for(slong i=1;i<np;i++) if(pts[i-1].u>=A&&pts[i].u<=B&&pts[i].s!=pts[i-1].s&&pts[i].u-pts[i-1].u>1e-9) fprintf(f,"%.9f %.9f\n",pts[i-1].u,pts[i].u); fclose(f);}
  printf("{\"T\":%.0f,\"P\":%ld,\"N\":%ld,\"L\":%d,\"M\":%ld,\"certified\":%ld,\"signchanges\":%ld,\"found_coarse\":%ld,\"Na_lower\":\"",Td,P,N,Lt,M,ncert,nz,found0);
  fmpz_print(flo); printf("\",\"Nb_upper\":\""); fmpz_print(fu); printf("\",\"max_possible\":"); fmpz_print(mx);
  printf(",\"refine_windows\":%ld,\"refine_points\":%ld,\"refine_disagreements\":%ld,\"found\":%ld,\"proven\":%d,\"status\":\"%s\",\"rounds\":%d,\"final_dr\":%g",nrefine,extra,disagree,found,(found==mxs&&disagree==0),(found==mxs&&disagree==0)?"PROVEN":(stopreason==1?"UNPROVEN_budget":stopreason==2?"UNPROVEN_drfloor":stopreason==3?"ERROR_overcount":"UNPROVEN_disagree"),rounds,dr);
  printf(",\"EPH\":%.3g,\"sec_setup\":%.3f,\"sec_nufft\":%.3f,\"sec_points\":%.3f,\"sec_refine\":%.3f,\"sec_total\":%.3f,\"interp\":%d,\"dips\":%ld,\"sec_dip\":%.3f,\"n_interp\":%lld,\"evals_per_zero\":%.3f,\"ms_per_zero\":%.4f}\n",tm[3],tm[0],tm[1],tm[2],tref,t5-t0,useI,ndip,tdip,NINTERP,(P+(double)(useI?NINTERP:extra))/(double)found,1e3*(t5-t0)/(double)found);
  if(argc>7){FILE*f=fopen(argv[7],"wb");fwrite(sg,1,P,f);fclose(f);}
  return 0;}
