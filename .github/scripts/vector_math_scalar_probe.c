/* Diagnostic only: interpose selected vector libm calls with scalar libm.
 * No solver source change; call counts verify whether the selected path ran.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
typedef double v4df __attribute__((vector_size(32)));
static unsigned long long counts[6];
#define UNARY(NAME, INDEX) v4df _ZGVdN4v_##NAME(v4df x) { \
  v4df y; counts[INDEX]++; \
  for (int i=0; i<4; ++i) y[i]=NAME(x[i]); return y; }
#ifdef SELECT_EXP
UNARY(exp,0)
#endif
#ifdef SELECT_LOG
UNARY(log,1)
#endif
#ifdef SELECT_LOG10
UNARY(log10,2)
#endif
#ifdef SELECT_POW
v4df _ZGVdN4vv_pow(v4df x, v4df y) {
 v4df z; counts[3]++; for (int i=0;i<4;++i) z[i]=pow(x[i],y[i]); return z;
}
#endif
#ifdef SELECT_SIN
UNARY(sin,4)
#endif
#ifdef SELECT_COS
UNARY(cos,5)
#endif
__attribute__((destructor)) static void report_counts(void) {
 const char *names[6]={"exp","log","log10","pow","sin","cos"};
 for(int i=0;i<6;++i) if(counts[i]) fprintf(stderr,"SCALAR_VECTOR_CALL %s %llu\n",names[i],counts[i]);
}
