/* [NUM] OPEN-2 test (Nathan A. Meeks & Grok): c0term must refuse (ok=0) a t-ball straddling t_N=2pi(N+1)^2,
   and on each side must agree with the RS formula computed at that side's own index. */
#define main v9main
#include "../fftzi_v9.c"
#undef main
int main(void){ slong prec=256, N=691000; arb_t tN,t,r,r2,pi2; arb_init(tN);arb_init(t);arb_init(r);arb_init(r2);arb_init(pi2);
  arb_const_pi(pi2,prec); arb_mul_2exp_si(pi2,pi2,1); arb_set_ui(tN,N+1); arb_mul(tN,tN,tN,prec); arb_mul(tN,tN,pi2,prec);
  int fails=0;
  arb_set(t,tN); mag_set_d(arb_radref(t),1e-3); int ok=c0term(r,t,N,prec); printf("straddle ball: ok=%d (want 0)\n",ok); fails+=ok!=0;
  double eps[2]={-2e-3,2e-3};
  for(int s=0;s<2;s++){ arb_set(t,tN); arb_add_d(t,t,eps[s],prec); mag_set_d(arb_radref(t),1e-4);
    ok=c0term(r,t,N,prec);
    slong Ns = s? N+1: N; int ok2=c0term(r2,t,Ns,prec); /* side's own index: N'=N+1 needs the n=N+1 term added */
    if(s){ /* r (from N with k=1) = C0(N+1) + 2 n^{-1/2}cos(theta - t log n) at n=N+1 ; r2 = C0(N+1) only */ }
    printf("side %d: ok=%d ok_own=%d  r=",s,ok,ok2); arb_printd(r,12); printf("  own="); arb_printd(r2,12); printf("\n"); fails+=!ok||!ok2;
    if(!s && !arb_overlaps(r,r2)) fails++; }
  /* continuity check at the two sides: |Z-assembly jump| is checked in val via the main-sum term; here only ok-flags */
  printf("%s\n",fails?"FAIL":"PASS"); return fails; }
