#include "main.h"
//#include "vars.h"

#define MAXSTP 1.0e10
#define TINY 1.0e-30


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
//The value ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.
void rkqs(long double y[], long double dydx[], int n, long double *x, long double htry, long double eps,
long double yscal[], long double *hdid, long double *hnext,
void (*derivs)(long double, long double [], long double []))

/*
Fifth�order Rung�Kutta step with
monitoring of local truncation error to ensure accuracy and
adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
at the
starting value of the independent variable x. Also input are the stepsize to be attempted
htry, the required accurac eps, and the vector
yscal[1..n] against which the error is
scaled. On output, y and x are replaced b their new values, hdid is the stepsize that was
actuall accomplished, and hnext is the estimated next stepsize. derivs is the user�supplied
routine that computes the
rig ht�hand side derivatives.     */
{
	
	//cout << "RK";

int i;
long double errmax, h,htemp,xnew,*yerr,*ytemp;
yerr=dvector(1,n);
ytemp=dvector(1,n);
h=htry;                                                           //Set stepsize to the initial trial value.
for (;;) {
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);                            //Take a step .
        errmax=0.0;                                                     //Evaluate accuracy.
        for (i=1;i<=n-2;i++) errmax=fmaxl(errmax,fabsl(yerr[i]/yscal[i]));
        errmax /= eps;                                                 //Scale relative to required tolerance.
        if (errmax <= 1.0) break;                                //Step succeeded. Compute size of next step.

        htemp=SAFETY*h*powl(errmax,PSHRNK);                             //Truncation error too large, reduce stepsize.
        h=(h >= 0.0 ? fmaxl(htemp,0.1*h) : fminl(htemp,0.1*h));          //No more than a factor of 10.
        xnew=(*x)+h;
        if (xnew == *x)
               { h = 1.5 *h;
               printf("!");
               break; }     // avoid stepsize underflow
        if (xnew == *x)
                {nrerror("stepsize underflow in rkqs"); fprintf(LOGSCR,"stepsize underflow in rkqs\n");}
        }
if (errmax > ERRCON) *hnext=SAFETY*h*powl(errmax,PGROW);
else *hnext=5.0*h;                                             //No more than a factor of 5 increase.
*x += (*hdid=h);

for (i=1;i<=n;i++) y[i]=ytemp[i];
free_vector(ytemp,1,n);
free_vector(yerr,1,n);
}

//The routine rkqs calls the routine rkck to take a Cash�Karp Runge�Kutta step:
#include "nrutil.h"
void rkck(long double y[], long double dydx[], int n, long double x, long double h, long double yout[],
long double yerr[], void (*derivs)(long double, long double [], long double []))
/*
Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x,use
the fifth�order Cash�Karp Rung�Kutta method to advance the solution over an interval h
and return the incremented variables as yout[1..n]. Also return an estimate of the local
truncation error in yout
using the embedded fourth�order method. The user supplies the routine
derivs(x,y,dydx), which returns derivatives dydx at x.    */
{
int i;


static long double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
dc5 = -277.00/14336.0;
long double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
dc4=c4-13525.0/55296.0,dc6=c6-0.25;
long double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
ak2=dvector(1,n);
ak3=dvector(1,n);
ak4=dvector(1,n);
ak5=dvector(1,n);
ak6=dvector(1,n);
ytemp=dvector(1,n);
for (i=1;i<=n;i++)                                          //First step.
        ytemp[i]=y[i]+b21*h*dydx[i];
(*derivs)(x+a2*h,ytemp,ak2);                                //Second step.
for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
(*derivs)(x+a3*h,ytemp,ak3);                                //Third step.
for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
(*derivs)(x+a4*h,ytemp,ak4);                                //Fourth step.
for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
(*derivs)(x+a5*h,ytemp,ak5);                                //Fifth step.
for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
(*derivs)(x+a6*h,ytemp,ak6);                                //Sixth step.
for (i=1;i<=n;i++)                                          //Accumulate increments with proper weights.
        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
for (i=1;i<=n;i++) 
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
//Estimate error as difference between fourth and fifth order methods.
free_dvector(ytemp,1,n);
free_dvector(ak6,1,n);
free_dvector(ak5,1,n);
free_dvector(ak4,1,n);
free_dvector(ak3,1,n);
free_dvector(ak2,1,n);
}
/*
Noting that the above routines are all in single precision, don't be too greedy in
specifying eps. The punishment for excessive greediness is interesting and worthy of
Gilbert and Sullivan's Mikado: The routine can always achieve an apparent zero error 
by making the stepsize so small that quantities of order hy # add to quantities of order 
y as if they were zero. Then the routine chugs happily along taking infinitely many
infinitesimal steps and never changing the dependent variables one iota. (You guard 
against this catastrophic loss of your computer budget by signaling on abnormally 
small stepsizes or on the dependent variable vector remaining unchanged from step 
to step. On a personal workstation you guard against it by not taking too long a
lunch hour while your program is running.)

Here is a full�fledged ``driver'' for Runge�Kutta with adaptive stepsize control. 
We warmly recommend this routine, or one like it, for a variety of problems, notably
including garden�variety ODEs or sets of ODEs, and definite integrals (augmenting
the methods of Chapter 4). For storage of intermediate results (if you desire to
inspect them) we assume that the top�level pointer references *xp and **yp have been
validly initialized (e.g., by the utilities vector() and matrix()). Because steps
occur at unequal intervals results are only stored at intervals greater than dxsav. The
top�level variable kmax indicates the maximum number of steps that can be stored.
If kmax=0 there is no intermediate storage, and the pointers *xp and **yp need not
point to valid memory. Storage of steps stops if kmax is exceeded, except that the 
ending values are always stored. Again, these controls are merely indicative of what 
you might need. The routine odeint should be customized to the problem at hand.*/


// extern int kmax,kount;
// extern long double *xp,**yp,dxsav;

/*User storage for intermediate results. Preset kmax and dxsav in the
calling program. If kmax #=0 results are stored at approximate intervals
dxsav in the arra s xp[1..kount], yp[1..nvar] [1..kount],wherekount is
output b odeint.
Defining declarations for these variables, with
memor allocations xp[1..kmax] and yp[1..nvar][1..kmax] for the arra s, should be in 
the calling program.   */


void odeintrk(long double ystart[], int nvar, long double x1, long double x2, long double eps, long double h1,
long double hmin, int *nok, int *nbad,
void (*derivs)(long double, long double [], long double []),
void (*rkqs)(long double [], long double [], int, long double *, long double, long double, long double [],
long double *, long double *, void (*)(long double, long double [], long double [])))   /*
Rung�Kutta driver with adaptive stepsize control. Integate
starting values ystart[1..nvar]
from x1 to x2 with accurac eps,
storing intermediate results
ing lobal variables. h1 should
be set as a gessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On
output nok and nbad are the number
ofg ood and bad (but retried and fixed) steps taken, and
ystart is replaced b values at the end of the integation interval. derivs is the user�supplied
routine for
calculating the rigt�hand side derivative, while rkqs isthenameofthestepper
routine to be used.    */
{
	int nstp,i;
	long double xsav = 0,x = 0,hnext = 0,hdid = 0,h = 0;
	long double *yscal,*y,*dydx;
	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=(x2 > x1) ? fabsl(h1) : -fabsl(h1);
	//h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i]; 
	if (kmax > 0) xsav=x-dxsav*2.0;                                             //Assures storag of first step.
	for (nstp=1;nstp<=MAXSTP;nstp++) {                                          //Take at most MAXSTP steps.
		(*derivs)(x,y,dydx); 
		for (i=1;i<=nvar;i++) //Scaling used to monitor accuracy.
			yscal[i]=fabsl(y[i])+fabsl(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabsl(x-xsav) > fabsl(dxsav)) {
			xp[++kount]=x;                                                              //Store intermediate results.
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;                                        //If stepsize can overshoot, decrease.
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);                    // Successful Runge Kutte step
		if (hdid == h) ++(*nok); 
		else ++(*nbad);
		if ((x - x2)*(x2-x1) >= 0.0) {                                                //Are we done?
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;                                                              //Save final step.
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return;                                                                     //Normal exit.
		}
		if (fabsl(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	}
	nrerror("Too many steps in routine odeint");
}


//--------------------------------------------------------------------------
// Modified Midpoint Method

void mmid(long double *y, long double *dydx, int nvar, long double xs,
    long double htot, int nstep, long double *yout,
    void (*derivs)(long double,long double*,long double*))

{
	int n,i;
	long double x,swap,h2,h,*ym,*yn,*dvector(int,int);
	void free_dvector(long double*,int,int);

	ym=dvector(1,nvar);
	yn=dvector(1,nvar);
	h=htot/nstep;
	for (i=1;i<=nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	(*derivs)(x,yn,yout);
        if (stopall == 1) return;
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
		for (i=1;i<=nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		(*derivs)(x,yn,yout);
                if (stopall == 1) return;
	}
	for (i=1;i<=nvar;i++)
		yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	free_dvector(yn,1,nvar);
	free_dvector(ym,1,nvar);
}


//--------------------------------------------------------------------------
// Rational Extrapolation

#define IMAX 11
#define SHRINK 0.95
#define GROW 1.2

extern long double **d,*x;	/* defined in BSSTEP */

void rzextr(int iest,long double xest, long double *yest,long double *yz,long double *dy, int nv){
    int k,j;
    long double yy = 0,v = 0,ddy = 0,c = 0,b1 = 0,b = 0,*fx,*dvector(int,int);
    void free_dvector(long double*,int,int);

    fx= dvector(1, iest);
    x[iest]=xest;
    if (iest == 1)
    	for (j=1;j<=nv;j++) {
    	      yz[j]=yest[j];
	      d[j][1]=yest[j];
	      dy[j]=yest[j];
        }
    else {
  	for (k=1;k<iest;k++)
		fx[k+1]=x[iest-k]/xest;
	for (j=1;j<=nv;j++) {
        	yy=yest[j];
		v=d[j][1];
		c=yy;
		d[j][1]=yy;
		for (k=2;k<=iest;k++) {
                    b1=fx[k]*v;
		    b=b1-c;
                    if (b) {
		    	b=(c-v)/b;
		    	ddy=c*b;
		    	c=b1*b;
		    } else
		      	ddy=v;
                    if (k != iest) v=d[j][k];
		    d[j][k]=ddy;
		    yy += ddy;
                }
	       	dy[j]=ddy;
	      	yz[j]=yy;
        }
    }
    free_vector(fx,1,iest);
}


//--------------------------------------------------------------------------
// Polynomial Extrapolation

extern long double **d,*x;	/* defined in BSSTEP */

void pzextr(int iest, long double xest, long double *yest, long double *yz,
    long double *dy, int nv)

{
    int k1,j;
    long double q,f2,f1,delta,*c,*dvector(int,int);
    void free_dvector(long double*,int,int);

    c=dvector(1,nv);
    x[iest]=xest;
    for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
    if (iest == 1) {
    	for (j=1;j<=nv;j++) d[j][1]=yest[j];
    }
    else{
        for (j=1;j<=nv;j++) c[j]=yest[j];
        for (k1=1;k1<iest;k1++) {
            delta=1.0/(x[iest-k1]-xest);
	    f1=xest*delta;
	    f2=x[iest-k1]*delta;
	    for (j=1;j<=nv;j++) {
	    	q=d[j][k1];
	    	d[j][k1]=dy[j];
	    	delta=c[j]-q;
	    	dy[j]=f1*delta;
	    	c[j]=f2*delta;
	    	yz[j] += dy[j];
	    }
        }
	for (j=1;j<=nv;j++) d[j][iest]=dy[j];
    }
    free_vector(c,1,nv);
}

//--------------------------------------------------------------------------
// Bulirsch-Stoer Step Size Controler

#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0E-05
#define REDMIN 0.7
#define SCALMX 0.1

long double **d=0,*x=0;	/* defining declaration */

void bsstep(long double *y, long double *dydx, int nv, long double *xx,
    long double htry, long double eps, long double *yscal, long double *hdid,
    long double *hnext, void (*derivs)(long double, long double*, long double*))
{
    void mmid(long double*,long double*,int,long double,long double,int,long double*,
         void (*derivs) (long double, long double*, long double*));
    void rzextr(int,long double,long double*,long double*,long double*,int);
    int i = 0,iq = 0,k = 0,kk = 0,km = 0;
    static int first= 1,kmax,kopt;
    static long double epsold = -1.0,xnew = 0;
    long double epsl = 0,errmax = 0,fact = 0,h = 0,red = 0,scale = 0,work  = 0,wrkmin = 0,xest = 0;
    long double *err,*yerr,*ysav,*yseq;
    static long double a[IMAXX+1];
    static long double alf[KMAXX+1][KMAXX+1];
    static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
    int reduct,exitflag=0;

	 //cout << "BS";
	 
    d=dmatrix(1,nv,1,KMAXX);
    err=dvector(1,KMAXX);
    x=dvector(1,KMAXX);
    yerr=dvector(1,nv);
    ysav=dvector(1,nv);
    yseq=dvector(1,nv);
    if (eps != epsold) {
       *hnext = xnew = -1.0e29;
       epsl=SAFE1*eps;
       a[1]=nseq[1]+1;
       for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
       for (iq=2;iq<=KMAXX;iq++) {
           for (k=1;k<iq;k++)
               alf[k][iq] =powl(epsl,(a[k+1] -a[iq+1])/
                         ((a[iq+1]-a[1]+1.0)*(2*k+1)));
           }
       epsold=eps;
       for (kopt=2; kopt<KMAXX; kopt++)
           if(a[kopt+1]>a[kopt]*alf[kopt-1][kopt]) break;
       kmax=kopt;
    }
    h=htry;
    for (i=1;i<=nv;i++) ysav[i]=y[i];
    if (*xx != xnew || h != (*hnext)) {
       first=1;
       kopt=kmax;
    }
    reduct=0;
    for (;;) {
        for (k=1;k<=kmax;k++) {
            xnew= (*xx) +h;
            if (xnew == (*xx)) {nrerror("step size underflow in bsstep"); OutputState(y, 1);} 
            mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
            if (stopall == 1) return;
            xest = h/nseq[k];
            xest=xest*xest;
            pzextr(k,xest,yseq,y,yerr,nv);
            if (k != 1) {
               errmax=TINY;
               for (i=1;i<=nv;i++) errmax=fmaxl(errmax,fabsl(yerr[i]/yscal[i]));
               errmax /= eps;
               km=k-1;
               err[km]=powl(errmax/SAFE1,1.0/(2*km+1));
            }
            if (k != 1 && (k >= kopt-1 || first)) {
               if (errmax < 1.0) {
                  exitflag=1;
                  break;
               }
               if (k == kmax || k == kopt+1) {
                  red=SAFE2/err[km];
                  break;
               }
               else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                  red=1.0/err[km] ;
                  break;
               }
               else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                  red=alf[km][kmax-1]*SAFE2/err[km];
                  break;
               }
               else if (alf[km][kopt] < err[km]) {
                  red=alf[km][kopt-1]/err[km];
                  break;
               }
            }
        }
        if (exitflag) break;
        red=fminl(red,REDMIN);
        red=fmaxl(red,REDMAX);
        h *= red;
        reduct=1;
    }
    *xx=xnew;
    *hdid=h;
    first=0;
    wrkmin=1.0E35;
    for (kk=1;kk<=km;kk++) {
        fact=fmaxl(err[kk],SCALMX);
        work=fact*a[kk+1] ;
        if (work < wrkmin) {
           scale=fact;
           wrkmin=work;
           kopt=kk+1;
        }
    }
    *hnext=h/scale;
    if (kopt >= k && kopt != kmax && !reduct) {
       fact=fmaxl(scale/alf[kopt-1][kopt],SCALMX);
       if (a[kopt+1]*fact <= wrkmin) {
          *hnext=h/fact;
          kopt++;
       }
    }
    free_dvector(yseq,1,nv);
    free_dvector(ysav,1,nv);
    free_dvector(yerr,1,nv);
    free_dvector(x,1,KMAXX);
    free_dvector(err,1,KMAXX);
    free_dmatrix(d,1,nv,1,KMAXX);
}

#undef IMAX
#undef NUSE
#undef SHRINK
#undef GROW


//--------------------------------------------------------------------------
// Integration "Driver" for Integration with adaptive stepsize control

//int kmax=0,kount=0;  /* defining declaration */
//long double *xp=0,**yp=0,dxsav=0;  /* defining declaration */
//-----------------------------------------------------------------------------
//int kmax=0,kount=0;  /* defining declaration */
//long double *xp=0,**yp=0,dxsav=0;  /* defining declaration */
//-----------------------------------------------------------------------------
void odeint (long double *ystart,int nvar, long double x1, long double x2, long double eps,long double h1, long double hmin, int *nok, int *nbad,void (*derivs)(long double, long double*, long double*),void (*integrator)(long double *,long double *, int, long double*, long double,long double, long double*, long double*,long double*, void (*derivs)(long double,long double*,long double*))){
		int nstp,i;
		long double xsav = 0,x = 0,hnext = 0,hdid = 0,h = 0;
		long double *yscal = NULL,*y = NULL,*dydx = NULL;

		yscal=dvector(1,nvar);
		y=dvector(1,nvar);
		dydx=dvector(1,nvar);
		
		x=x1;
		h=(x2 > x1) ? fabsl(h1) : -fabsl(h1);
		*nok = (*nbad) = kount = 0;
		
		for (i=1;i<=nvar;i++) 
			y[i]=ystart[i];
		if (kmax > 0) 
			xsav=x-dxsav*2.0;                             //Assures storage of first step.
		for (nstp=1;nstp<=MAXSTP;nstp++){                          //Step loop: Take at most MAXSTP steps. 
			
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// only check for reflections of absorbtions when outside bulk of PENeLOPE storage volume
						if((y[3]>detz-wanddicke)||(y[1]<=StorVolrmin+wanddicke)||(y[1]>=StorVolrmax-wanddicke)||(y[3]<=StorVolzmin+wanddicke)||(y[3]>=StorVolzmax)  )
						{
							blankint = GeomCheck(y,y[1],y[2],y[3],y[4],y[5],y[6],x);
							if(stopall)
								return;
							//cout << "GEOMCHECK t: " << x <<  " r: " << y[1] << " z: " << y[3] << endl;
							//cin >> blankint;
						}
			// spin flip properties according to Vladimirsky and thumbrule
			if (spinflipcheck == 2){   
				// y[6]: phidot
				vlad = vladimirsky(y[1], Br, Bphi, Bz, dBrdr, dBrdphi, dBrdz, dBphidr, dBphidphi, dBphidz, dBzdr, dBzdphi, dBzdz, y[2], y[6], y[4]);
				frac = thumbrule(Br, Bphi, Bz, dBrdr, dBrdphi, dBrdz, dBphidr, dBphidphi, dBphidz, dBzdr, dBzdphi, dBzdz, y[2], y[6], y[4]);
				if (vlad > 1e-99){
					vladtotal = vladtotal * (1-vlad);
					if (vladtotal < 0.9999)
						printf(" VladShit at t= %.17LG\n",x);
				}
				if ((vlad > vladmax)&&(vlad > 1e-99))
					vladmax=log10l(vlad);
				if ((frac > thumbmax)&&(frac > 1e-99))
					thumbmax=log10l(frac);
			}
			// Ende Spinflipwahrscheinlichkeit
			
			// Trajectory length calculation
						if (TrajectoryLength)
						{
							trajlength =  sqrtl(powl(ytemp1-y[1],2)+powl(ytemp3-y[3],2)+powl(ytemp5*ytemp1-y[5]*y[1],2));
							trajlengthsum += trajlength;
							ytemp1 = y[1]; ytemp3 = y[3]; ytemp5 = y[5];
						}
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			(*derivs)(x,y,dydx);
			for (i=1;i<=nvar;i++)
				yscal[i]=fabsl(y[i])+fabsl(dydx[i]*h)+TINY;          //Scaling used to monitor accuracy.
			if (kmax > 0){
				if (fabsl(x-xsav) > fabsl(dxsav)){
					if (kount < kmax-1){
						
						// Vektoren und B-Feld zwischenspeichern
						if(SaveIntermediate)
						{
							xp[++kount]=x;                     //Store intermediate results.
							for (i=1;i<=nvar;i++) 
							yp[i][kount]=y[i];	
							BFeld(yp[1][kount],yp[5][kount],yp[3][kount],xp[kount]);
							Bp[1][kount]=Br;
							Bp[2][kount]=dBrdr;
							Bp[3][kount]=dBrdphi;
							Bp[4][kount]=dBrdz;
							Bp[5][kount]=Bphi;
							Bp[6][kount]=dBphidr;
							Bp[7][kount]=dBphidphi;
							Bp[8][kount]=dBphidz;
							Bp[9][kount]=Bz;
							Bp[10][kount]=dBzdr;
							Bp[11][kount]=dBzdphi;
							Bp[12][kount]=dBzdz;
							Bp[13][kount]=Bws;
							Ep[1][kount] = Er;
							Ep[2][kount] = Ez;
							//cout << "integrator: " << delx_n << " " << kount << " " << xp[kount] << " " << Br << " " << Bz << endl;
						}
						// ENDE Vektoren und B-Feld zwischenspeichern				
				
						xsav=x;
						
					}
				}
			}
		
			if ((x+h-x2)*(x+h-x1) > 0.0) 
				h=x2-x;                                          //If stepsize can overshoot, decrease.

			//(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);        // runge kutta
			 (*integrator)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);        // runge kutta or bulirsch stoer step
			
			
			if (hdid == h) 
				++(*nok); 
			else 
				++(*nbad);                                   //Are we done?
			if (stopall == 1) 
				return;
		
			if ((x-x2)*(x2-x1) >= 0.0){							
				for (i=1;i<=nvar;i++)   // save final step in ystart
					ystart[i]=y[i];
				if (kmax){  // war kmax, ist nun deaktiviert
					xp[++kount]=x;                                             //Save final step.
					for (i=1;i<=nvar;i++) 
						yp[i][kount]=y[i];					
					// B-Feld zwischenspeichern
					if(SaveIntermediate)
					{
					BFeld(yp[1][kount],yp[5][kount],yp[3][kount],xp[kount]);
					Bp[1][kount]=Br;
					Bp[2][kount]=dBrdr;
					Bp[3][kount]=dBrdphi;
					Bp[4][kount]=dBrdz;
					Bp[5][kount]=Bphi;
					Bp[6][kount]=dBphidr;
					Bp[7][kount]=dBphidphi;
					Bp[8][kount]=dBphidz;
					Bp[9][kount]=Bz;
					Bp[10][kount]=dBzdr;
					Bp[11][kount]=dBzdphi;
					Bp[12][kount]=dBzdz;
					Bp[13][kount]=Bws;
					Ep[1][kount] = Er;
					Ep[2][kount] = Ez;
					//cout << "Ja!";
					}
					// B-Feld zwischenspeichern Ende					
				}
				
				// Trajectory length calculation
					if (TrajectoryLength)
					{
						trajlength =  sqrtl(powl(ytemp1-y[1],2)+powl(ytemp3-y[3],2)+powl(ytemp5*ytemp1-y[5]*y[1],2));
						trajlengthsum += trajlength;
						ytemp1 = y[1]; ytemp3 = y[3]; ytemp5 = y[5];
					}

			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			return;
		}
		
		if (fabsl(hnext) <= hmin) 
			nrerror("Step size too small in ODEINT");
		h=hnext;
		

	}
	//nrerror("Too many steps in routine ODEINT");
}
#undef MAXSTP
#undef TINY

#define SAFETY 0.9
#define GROW 1.5
//#define PGROW -0.25
#define SHRNK 0.5
//#define PSHRNK (-1.0/3.0)
//#define ERRCON 0.1296
#define MAXTRY 40
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

/*
Fourth-order Rosenbrock step for integrating stiff ODE's with monitoring of local truncation error to adjust stepsize.
Input are the dependent variable vector y[1...n] and its derivative dydx[1...n] at the starting value of the independent vairable x.
Also input are the stepsize to be attempted htry, the required accuracy eps, and the vextor yscal[1...n] against which the error is
scaled. On output y and x are replaced by their new values, hdid is the stepsize that was actually accomplished, and hnext is the
estimated next stepsize. derivs is a user-supplied routine that computes the derivatives of the right-hand sied with respect to x,
while jacobn is a user-supplied routine that computes the Jacobi matrix of derivatives of the right-hand side with respect to the
components of y.
*//*
void rosenbrock(long double *y,long double *dydx,int n, long double *x,long double htry,long double eps, long double *yscale,long double *hdid,long double *hnext,void (*derivs)(long double,long double *,long double *))
{
	int i,j=0,jtry,*indx;
	long double d,errmax,h,xsav;
	long double *g1,*g2,*g3,*g4,*ysav,*B,*dfdx,*dysav,*err;
	long double **a,**dfdy;
	
	indx = ivector(1,n);
	a = save_matrix(1,n,1,n);
	dfdx = dvector(1,n);
	dfdy = save_matrix(1,n,1,n);
	dysav = vector(1,n);
	err = vector(1,n);
	g1 = vector(1,n);
	g2 = vector(1,n);
	g3 = vector(1,n);
	g4 = vector(1,n);
	ysav = vector(1,n);
	B = vector(1,3);
	xsav = *x;
	h = htry;
	
	BFinterpol_new(*x,B);
	jacobn(*x,y,dfdx,dfdy,n,B);

	for(i=1;i<=n;i++)
	{
		ysav[i] = y[i];
		dysav[i] = dydx[i];
	}
	
	for(jtry=1;jtry<=MAXTRY;jtry++)
	{
		for(i=1;i<=n;i++)
		{
			for(j=1;j<=n;j++)
			{
				a[i][j] = (-1)*dfdy[i][j];
			}
			a[i][i] += 1.0/(GAM*h);
		}
		ludcmp(a,n,indx,&d);
		for(i=1;i<=n;i++)
		{
			//g1[i] = dysav[i]+h*C1X*dfdx[i];
			g1[i] = dysav[i]+h;
		}
		lubksb(a,n,indx,g1);
		for(i=1;i<=n;i++)
		{
			y[i] = ysav[i]+A21*g1[i];
		}
		*x = xsav+A2X*h;
		(*BFderivs_new)(*x,y,dydx,B);
		for(i=1;i<=n;i++)
		{
			//g2[i] = dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
			g2[i] = dydx[i]+C21*g1[i]/h;
		}
		lubksb(a,n,indx,g2);
		for(i=1;i<=n;i++)
		{
			y[i] = ysav[i]+A31*g1[i]+A32*g2[i];
		}
		*x = xsav+A3X*h;
		(*BFderivs_new)(*x,y,dydx,B);
		for(i=1;i<=n;i++)
		{
			//g3[i] = dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
			g3[i] = dydx[i]+(C31*g1[i]+C32*g2[i])/h;
		}
		lubksb(a,n,indx,g3);
		for(i=1;i<=n;i++)
		{
			//g4[i] = dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
			g4[i] = dydx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
		}
		lubksb(a,n,indx,g4);
		for(i=1;i<=n;i++)
		{
			y[i] = ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i] = E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
		}
		*x = xsav+h;
		if(*x == xsav) 
			nrerror("stepsize not significant in rosenbrock!");
		errmax = 0.0;
		for(i=1;i<=n;i++)
		{
			errmax = fmaxl(errmax,fabsl(err[i]/yscale[i]));
		}
		errmax /= eps;
		if(errmax <= 1.0)
		{
			*hdid = h;
			*hnext = (errmax > ERRCON ? SAFETY*h*powl(errmax,PGROW) : GROW+h);
			free_vector(ysav,1,n);
			//printf("ysav free\n");
			free_vector(g4,1,n);
			//printf("g4 free\n");
			free_vector(g3,1,n);
			//printf("g3 free\n");
			free_vector(g2,1,n);
			//printf("g2 free\n");
			free_vector(g1,1,n);
			//printf("g1 free\n");
			free_vector(err,1,n);
			//printf("err free\n");
			free_vector(dysav,1,n);
			//printf("dysav free\n");
			//free_vector(dfdx,1,n); // we don't need this one
			//printf("dfdx free\n");
			free_matrix(dfdy,1,n,1,n);
			//printf("dfdy free\n");
			free_matrix(a,1,n,1,n);
			//printf("a free\n");
			free_ivector(indx,1,n);
			//printf("indx free\n");
			free_vector(B,1,3);
			//printf("B free\n");
			return;
		}
		else
		{
			*hnext = SAFETY*h*powl(errmax,PSHRNK);
			h = (h >=0.0 ? fmaxl(*hnext,SHRNK*h) : fminl(*hnext,SHRNK*h));
		}
	}
	nrerror("exeeded MAXTRY in rosenbrock!");
}

void jacobn(long double x,long double *y,long double *dfdx,long double **dfdy,int n,long double *B)
{
	int i;
	for(i=0;i<=n;i++)
		dfdx[i] = 0.0;
	dfdx = NULL;
	dfdy[1][1] =  0.0;
	dfdy[1][2] = -1*B[3]*gamma_n;
	dfdy[1][3] =  B[2]*gamma_n;
	dfdy[2][1] =  B[3]*gamma_n;
	dfdy[2][2] =  0.0;
	dfdy[2][3] = -1*B[1]*gamma_n;
	dfdy[3][1] = -1*B[2]*gamma_n;
	dfdy[3][2] =  B[1]*gamma_n;
	dfdy[3][3] =  0.0;
	
}

void jacobn_bf(long double x,long double *y,long double *dfdx,long double **dfdy,int n,long double *B)
{
	int i;
	//for(i=0;i<=n;i++)
	//	dfdx[i] = 0.0;
	dfdx = NULL;
	dfdy[1][1] =  0.0;
	dfdy[1][2] = -1*B[3]*gamma_n;
	dfdy[1][3] =  B[2]*gamma_n;
	dfdy[2][1] =  B[3]*gamma_n;
	dfdy[2][2] =  0.0;
	dfdy[2][3] = -1*B[1]*gamma_n;
	dfdy[3][1] = -1*B[2]*gamma_n;
	dfdy[3][2] =  B[1]*gamma_n;
	dfdy[3][3] =  0.0;
	
}

void lubksb(long double **a,int n,int *indx,long double *b)
{
	int i,ii=0,ip,j;
	long double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

#undef TINY
#define TINY 1.0e-20;

void ludcmp(long double **a,int n,int *indx,long double *d)
{
	int i,imax,j,k;
	long double big,dum,sum,temp;
	long double *vv,*vector();
	void nrerror(),free_vector();

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabsl(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabsl(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}

#undef TINY
*/

