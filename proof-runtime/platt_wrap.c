/* [NUM] Thin C wrapper around FLINT 3 Platt isolation (Nathan A. Meeks & Grok). Candidate, unrefereed; not RH.
   Header-free: declares only the FLINT 3 ABI it uses (validated against FLINT 3.1.3 headers). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long slong; typedef long fmpz; typedef unsigned long ulong;
typedef struct { fmpz exp; slong size; ulong d[2]; } arf_struct;
typedef struct { arf_struct a, b; } arf_interval_struct;
typedef struct { fmpz exp; ulong man; } mag_struct; typedef struct { arf_struct mid; mag_struct rad; } arb_struct;
extern slong acb_dirichlet_platt_local_hardy_z_zeros(arb_struct*, const fmpz*, slong, slong);
extern void arb_get_lbound_arf(arf_struct*, const arb_struct*, slong); extern void arb_get_ubound_arf(arf_struct*, const arb_struct*, slong);
extern void fmpz_init(fmpz*); extern void fmpz_clear(fmpz*); extern int fmpz_set_str(fmpz*, const char*, int);
extern char* fmpz_get_str(char*, int, const fmpz*); extern void flint_free(void*);
extern void arf_get_fmpz_2exp(fmpz*, fmpz*, const arf_struct*);
extern slong acb_dirichlet_platt_isolate_local_hardy_z_zeros(arf_interval_struct*, const fmpz*, slong, slong);
extern void flint_set_num_threads(int);
static void pr(FILE* f, const arf_struct* x){ fmpz m, e; fmpz_init(&m); fmpz_init(&e); arf_get_fmpz_2exp(&m,&e,x);
  char*s=fmpz_get_str(NULL,10,&m),*t=fmpz_get_str(NULL,10,&e); fprintf(f," %s %s",s,t); flint_free(s); flint_free(t); fmpz_clear(&m); fmpz_clear(&e);}
/* usage: platt_iso n len [prec] [threads] [i|r]  (i=isolating intervals, r=refined balls) -> line 1: count; then per zero k (index n+k): ma ea mb eb  (a=ma*2^ea, b=mb*2^eb, exact) */
int main(int c, char** v){ if(c<3) return 2; fmpz n; fmpz_init(&n); if(fmpz_set_str(&n,v[1],10)) return 2;
  slong len=atol(v[2]), prec=c>3?atol(v[3]):256; if(c>4) flint_set_num_threads(atoi(v[4]));
  if(c>5 && v[5][0]=='r'){ /* refined mode: certified ball for each zero, printed as exact [lower,upper] */
    arb_struct* z=calloc(len,sizeof *z); slong k=acb_dirichlet_platt_local_hardy_z_zeros(z,&n,len,prec); printf("%ld\n",k);
    for(slong i=0;i<k;i++){ arf_struct l,u; memset(&l,0,sizeof l); memset(&u,0,sizeof u); arb_get_lbound_arf(&l,&z[i],prec); arb_get_ubound_arf(&u,&z[i],prec); pr(stdout,&l); pr(stdout,&u); printf("\n"); }
    return 0; }
  arf_interval_struct* r=calloc(len,sizeof *r); /* == arf_init (all-zero) */ slong k=acb_dirichlet_platt_isolate_local_hardy_z_zeros(r,&n,len,prec);
  printf("%ld\n",k); for(slong i=0;i<k;i++){ pr(stdout,&r[i].a); pr(stdout,&r[i].b); printf("\n"); }
  free(r); /* process exit; limbs reclaimed by OS */ fmpz_clear(&n); return 0; }
