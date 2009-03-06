extern void derivs(long double,long double *,long double *);

extern void rkck(long double y[], long double dydx[], int n, long double x, long double h, long double yout[], long double yerr[], void (*derivs)(long double, long double [], long double []));
	
extern void mmid(long double *y, long double *dydx, int nvar, long double xs, long double htot, int nstep, long double *yout,void (*derivs)(long double,long double*,long double*));

extern void rzextr(int iest,long double xest, long double *yest,long double *yz, long double *dy, int nv);

extern void pzextr(int iest, long double xest, long double *yest, long double *yz, long double *dy, int nv);

extern void bsstep(long double *y, long double *dydx, int nv, long double *xx, long double htry, long double eps,long double *yscal, long double *hdid, long double *hnext, void (*derivs)(long double, long double*, long double*));

extern void odeint(long double *ystart,int nvar, long double x1, long double x2, long double eps, long double h1, long double hmin,int *nok, int *nbad, void (*derivs)(long double, long double*, long double*), void (*bsstep)(long double *,long double *,int, long double*, long double, long double, long double*, long double*, long double*, void (*derivs)(long double,long double*,long double*)));

extern void rkqs(long double y[], long double dydx[], int n, long double *x, long double htry, long double eps, long double yscal[], long double *hdid, long double *hnext, void (*derivs)(long double, long double [], long double []));

extern void odeintrk(long double ystart[], int nvar, long double x1, long double x2, long double eps, long double h1, long double hmin, int *nok, int *nbad,void (*derivs)(long double, long double *, long double *),void (*rkqs)(long double [], long double [], int, long double *, long double, long double, long double [],long double *, long double *, void (*)(long double, long double [], long double [])));

extern void rosenbrock(long double *y,long double *dydx,int n, long double *x,long double htry,long double eps, long double *yscale,long double *hdid,long double *hnext,void (*derivs)(long double,long double *,long double *));
extern void ludcmp(long double **a,int n,int *indx,long double *d);
extern void lubksb(long double **a,int n,int *indx,long double *b);
extern void jacobn(long double x,long double *y,long double *dfdx,long double **dfdy,int n,long double *B);
