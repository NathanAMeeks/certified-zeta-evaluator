#define main fftzi_main
#include "fftzi_v9.c"
#undef main
int main(int argc,char**argv){ slong prec=atol(argv[1]); double rad=atof(argv[2]); double x; arb_t p,r; arb_init(p);arb_init(r);
  while(scanf("%lf",&x)==1){ arb_set_d(p,x); if(rad>0) mag_set_d(arb_radref(p),rad); int ok=1; psi_c0(r,p,prec,&ok);
    printf("%.17g %d %s %.17g\n",x,ok,arf_get_str(arb_midref(r),90),mag_get_d(arb_radref(r))); } return 0; }
