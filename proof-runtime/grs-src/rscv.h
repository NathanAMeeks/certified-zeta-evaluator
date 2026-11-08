/* Lane-friendly copy of rsc() (Nathan A. Meeks & Grok): identical operations and order, branches replaced by selects; returns ok flag. */
#pragma omp declare simd notinbranch
static inline int rscv(double x,double*s,double*c){
  const double C1=0x1.921fb54442d18p0, C2=0x1.1a62633145c07p-54, I=0x1.45f306dc9c883p-1;
  int ok=(fabs(x)<=1048576.0);
  double k=nearbyint(x*I); double r=fma(-k,C1,x); r=r-k*C2; ok&=(fabs(r)<=0.8);
  double r2=r*r;
  double ps=1.0/121645100408832000.0;
  ps=ps*r2-1.0/355687428096000.0; ps=ps*r2+1.0/1307674368000.0; ps=ps*r2-1.0/6227020800.0; ps=ps*r2+1.0/39916800.0;
  ps=ps*r2-1.0/362880.0; ps=ps*r2+1.0/5040.0; ps=ps*r2-1.0/120.0; ps=ps*r2+1.0/6.0; ps=ps*r2-1.0;
  double sn=-ps*r;
  double pc=1.0/2432902008176640000.0;
  pc=pc*r2-1.0/6402373705728000.0; pc=pc*r2+1.0/20922789888000.0; pc=pc*r2-1.0/87178291200.0; pc=pc*r2+1.0/479001600.0;
  pc=pc*r2-1.0/3628800.0; pc=pc*r2+1.0/40320.0; pc=pc*r2-1.0/720.0; pc=pc*r2+1.0/24.0; pc=pc*r2-0.5; pc=pc*r2+1.0;
  long q=((long)k)&3;
  double S=(q==0)?sn:((q==1)?pc:((q==2)?-sn:-pc));
  double C=(q==0)?pc:((q==1)?-sn:((q==2)?-pc:sn));
  *s=S; *c=C; return ok; }
