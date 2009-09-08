extern void BFpolint(long double xa[],long double ya[],int n,long double x, long double *y,long double *dy);
extern void BFratint(long double xa[],long double ya[],int n, long double x,long double *y,long double *dy);
extern void BFinterpol(long double x,long double *Bx,long double *By,long double *Bz);
extern void BFderivs(long double x, long double *y, long double *dydx);
extern void BFrkck(long double y[], long double dydx[], int n, long double x, long double h, long double yout[],long double yerr[], void (*BFderivs)(long double, long double [], long double []));
extern void BFrkqs(long double y[], long double dydx[], int n, long double *x, long double htry, long double eps,long double yscal[], long double *hdid, long double *hnext, void (*BFderivs)(long double, long double [], long double []));
extern void BFodeintrk(long double ystart[], int nvar, long double x1, long double x2, long double eps, long double h1,long double hmin, int *nok, int *nbad,void (*BFderivs)(long double, long double [], long double []),void (*BFrkqs)(long double [], long double [], int, long double *, long double, long double, long double [],long double *, long double *, void (*)(long double, long double [], long double [])));
extern void BFderivs_new(long double x, long double *y, long double *dydx,long double *B);
extern void BFinterpol_new(long double x,long double *B);

void hunt(long double xx[], int n, long double x, int *jlo);
