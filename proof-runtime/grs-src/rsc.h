/* [NUM] rigorous sin/cos for |x|<=2^20 with a priori abs error <= RSC_E = 64*2^-53 each (Nathan A. Meeks & Grok audit fix).
   Reduction: k=rint(x*2/pi), r1=fma(-k,C1,x) (one rounding, <=u|r1|), r=r1-k*C2 (<=u|kC2|+u|r|), C2 err<=2^-107 => |r-r*|<=2.2u for |r|<=0.8.
   Poly: Taylor to deg 19/20 (trunc <=0.8^21/21!<1e-21); coeff rounding <=u*cosh(.8); r2 rounding <=1u; Horner gamma_22*sum|a_i||r|^i<=30u. Total <=36u<64u.
   Returns 0 (caller must refuse) if |x|>2^20 or |r|>0.8. Needs IEEE double, round-to-nearest, no -ffast-math. */
#include <math.h>
static inline int rsc(double x,double*s,double*c){
  if(!(fabs(x)<=1048576.0)) return 0;
  const double C1=0x1.921fb54442d18p0, C2=0x1.1a62633145c07p-54, I=0x1.45f306dc9c883p-1;
  double k=nearbyint(x*I); double r=fma(-k,C1,x); r=r-k*C2; if(!(fabs(r)<=0.8)) return 0;
  double r2=r*r;
  double ps=1.0/121645100408832000.0; /* 1/19! */
  ps=ps*r2-1.0/355687428096000.0; ps=ps*r2+1.0/1307674368000.0; ps=ps*r2-1.0/6227020800.0; ps=ps*r2+1.0/39916800.0;
  ps=ps*r2-1.0/362880.0; ps=ps*r2+1.0/5040.0; ps=ps*r2-1.0/120.0; ps=ps*r2+1.0/6.0; ps=ps*r2-1.0; /* = -(sin r)/r ... */
  double sn=-ps*r;
  double pc=1.0/2432902008176640000.0; /* 1/20! */
  pc=pc*r2-1.0/6402373705728000.0; pc=pc*r2+1.0/20922789888000.0; pc=pc*r2-1.0/87178291200.0; pc=pc*r2+1.0/479001600.0;
  pc=pc*r2-1.0/3628800.0; pc=pc*r2+1.0/40320.0; pc=pc*r2-1.0/720.0; pc=pc*r2+1.0/24.0; pc=pc*r2-0.5; pc=pc*r2+1.0;
  long q=((long)k)&3;
  if(q==0){*s=sn;*c=pc;} else if(q==1){*s=pc;*c=-sn;} else if(q==2){*s=-sn;*c=-pc;} else {*s=-pc;*c=sn;}
  return 1; }
#define RSC_E (64.0*0x1p-53)
