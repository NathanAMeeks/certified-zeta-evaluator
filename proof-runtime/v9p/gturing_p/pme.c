/* G-Turing-P: Platt multieval signs for Turing mesh. Nathan A. Meeks & Grok. Independent of rsbatch.c (Arb/FLINT DFT path). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/acb_dirichlet.h>
#include <flint/thread_support.h>
int main(int argc,char**argv){
  /* usage: pme T W threads  -> prints "u sign" for u=n/A, |u|<=W ; sign 0 = uncertified */
  fmpz_t T,J; fmpz_init(T); fmpz_init(J); fmpz_set_str(T,argv[1],10);
  double W=atof(argv[2]); int th=atoi(argv[3]); flint_set_num_threads(th);
  slong prec=128, A=8,B=4096; if(argc>5){A=atol(argv[5]);B=atol(argv[6]);}
  double x=log(atof(argv[1])),x2=x*x,x3=x2*x,x4=x2*x2;
  double dlogJ=-0.4035+0.49086*x+0.00016299*x2-3.6139e-06*x3+2.9323e-08*x4;
  double dK=79.032-1.781*x+0.039243*x2-0.00094859*x3+7.3149e-06*x4;
  double dgrid=1186.9+130.17*x-7.4059*x2+0.17895*x3-0.001602*x4;
  double dh=178.66-15.127*x+0.93132*x2-0.02311*x3+0.00022146*x4;
  if(argc>4) prec=atol(argv[4]);
  arb_t h; arb_init(h); arb_set_d(h,dh); fmpz_set_si(J,(slong)exp(dlogJ));
  slong K=(slong)dK, sg=((slong)(dgrid/2))*2+1, N=A*B;
  arb_ptr out=_arb_vec_init(N);
  if(th>1) acb_dirichlet_platt_multieval_threaded(out,T,A,B,h,J,K,sg,prec);
  else acb_dirichlet_platt_multieval(out,T,A,B,h,J,K,sg,prec);
  fprintf(stderr,"J=%ld K=%ld sigma=%ld h=%g\n",(long)fmpz_get_si(J),(long)K,(long)sg,dh);
  for(slong i=0;i<N;i++){ double u=(double)(i-N/2)/A; if(fabs(u)>W) continue;
    int s= arb_is_positive(out+i)?1:(arb_is_negative(out+i)?-1:0);
    printf("%.4f %d %.6e\n",u,s,arf_get_d(arb_midref(out+i),ARF_RND_NEAR)); }
  return 0;}
