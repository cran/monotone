
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void fake( int *n, double* x, double* w );
extern void fitm( int* n, double* x, double* w );
extern void wmrmnh( int* n, double* x, double* w );
extern void amalgm( int* n, double* x, double* w );
extern void pav( int* n, double* x, double* w );
extern void isoreg( int* n, double* x, double* w );
extern void iso_pava( int* n, double* x, double* w );
extern void isotonic( int* n, double* x, double* w );
extern void isomean( int* n, double* x, double* w );
extern void pooled_pava( int* n, double* x, double* w );
extern void linear_pava( int* n, double* x, double* w );
extern void inplace_pava( int* n, double* x, double* w );
extern void md_pava( int* n, double* x, double* w );
extern void reg_1d_l2( int* n, double* x, double* w );
extern void jbkpava( int* n, double* x, double* w );

extern void monotoneC( int* n, double* x, double* w );
extern void unimonotoneC( int* n, double* x, double* w );
extern void bimonotoneC( int* n, int* m, double* x, double* w, int* maxiter, double* eps );

static const R_CMethodDef CEntries[] = {
  {"fake",         ( DL_FUNC ) &fake,         3},
  {"fitm",         ( DL_FUNC ) &fitm,         3},
  {"wmrmnh",       ( DL_FUNC ) &wmrmnh,       3},
  {"amalgm",       ( DL_FUNC ) &amalgm,       3},
  {"pav",          ( DL_FUNC ) &pav,          3},
  {"isoreg",       ( DL_FUNC ) &isoreg,       3},
  {"iso_pava",     ( DL_FUNC ) &iso_pava,     3},
  {"isotonic",     ( DL_FUNC ) &isotonic,     3},
  {"isomean",      ( DL_FUNC ) &isomean,      3},
  {"pooled_pava",  ( DL_FUNC ) &pooled_pava,  3},
  {"linear_pava",  ( DL_FUNC ) &linear_pava,  3},
  {"inplace_pava", ( DL_FUNC ) &inplace_pava, 3},
  {"md_pava",      ( DL_FUNC ) &md_pava,      3},
  {"reg_1d_l2",    ( DL_FUNC ) &reg_1d_l2,    3},
  {"jbkpava",      ( DL_FUNC ) &jbkpava,      3},
  {"monotoneC",    ( DL_FUNC ) &monotoneC,    3},
  {"unimonotoneC", ( DL_FUNC ) &unimonotoneC, 3},
  {"bimonotoneC",  ( DL_FUNC ) &bimonotoneC,  6},
  {NULL, NULL, 0}
};

void R_init_monotone( DllInfo *dll )
{
  R_registerRoutines( dll, CEntries, NULL, NULL, NULL );
  R_useDynamicSymbols( dll, FALSE );
}
