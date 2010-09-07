void bsstep(long double *y, long double *dydx, int nv, long double *xx, long double htry, long double eps,long double *yscal, long double *hdid, long double *hnext, void *params, void (*derivs)(long double, long double*, long double*, void*));

void rkqs(long double y[], long double dydx[], int n, long double *x, long double htry, long double eps, long double yscal[], long double *hdid, long double *hnext, void *params, void (*derivs)(long double, long double [], long double [], void*));
