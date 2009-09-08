//---------------------------------------------------------------------------
// user supplied equations for BruteForce integration of Bloch Equation

#include "main.h"
//#include "vars.h"

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 1.0e10
#define TINY 1.0e-30


void BFpolint(long double xa[],long double ya[],int n,long double x, long double *y,long double *dy){
	int i,m,ns=1;
	long double den,dif,dift,ho,hp,w;
	long double *c,*d; //*dvector();

	dif=fabsl(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabsl(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");    // only occurs if 2 input values xa have same value
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}

#define FREERETURN {free_dvector(d,1,n);free_dvector(c,1,n);return;}

void BFratint(long double xa[],long double ya[],int n, long double x,long double *y,long double *dy){
	
	int m,i,ns=1;
	long double w,t,hh,h,dd,*c,*d; //*dvector();
	
	c=dvector(1,n);
	d=dvector(1,n);
	hh=fabsl(x-xa[1]);
	for (i=1;i<=n;i++) {
		h=fabsl(x-xa[i]);
		if (h == 0.0) {
			*y=ya[i];
			*dy=0.0;
			FREERETURN
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) nrerror("Error in routine RATINT");
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	FREERETURN
}

#undef FREERETURN

void BFinterpol(long double x,long double *Bx,long double *By,long double *Bz){
	int ordnung = 5;              // ordnung means that polynomial or rational function is of order (ordnung-1) 
	long double dy;
	hunt(BFtime, offset, x, &BFindex);       // x lies between BFindex and (BFindex+1)

	if (BFindex < 1)
        BFindex = 1;            // start interpolation on index 1
	if (BFindex > (offset-ordnung+1))
        BFindex = (offset-ordnung+1);   //do not interpolate out of boundary of array values	
	
	/*
	for(int hp=0;hp<=3000;hp++)
	{
		fprintf(LOGSCR,"%d t %LG Bx %LG By %LG Bz %LG \n",hp,BFtime[hp],BFField[2][hp],BFField[1][hp],BFField[3][hp]);
	}*/
	
	BFpolint(&BFtime[BFindex-1],&BFField[1][BFindex-1],ordnung,x, Bx,&dy);   // interpolation of Bx on BFtime[Bfindex to BFindex+ordnung-1]
	BFpolint(&BFtime[BFindex-1],&BFField[2][BFindex-1],ordnung,x, By,&dy);
	BFpolint(&BFtime[BFindex-1],&BFField[3][BFindex-1],ordnung,x, Bz,&dy);
	
	//*Bx=10787.1388395727981340375148505 * x * x * x-379485.201099151545103425695947 * x * x + 4425265.48802544343656564311135 * x - 17089608.2309508308489746319934;
	//*By=16511.3272584589837067791527764 * x * x * x-594543.289553430286823139813885 * x * x + 7118405.94609221980769652287442 * x - 28331085.0061183671391854839742;
	//*Bz=1944281.06055634049312257407394 * x * x * x-69404618.8278242196709266829602 * x * x + 822964791.430415938909499801820 * x - 3239972252.28819208686600696119;
		
	//*Bx=B1 * cosl(x*gamma_n*1.0e-3);
	//*By=B1 * sinl(x*gamma_n*1.0e-3);
	//*Bz=1.0e-3;
/*
	if ((*Bx!=0) && (*By!=0) && (*Bz!=0))
	{
		Bxdev = fabsl(((B1 * cosl(x*gamma_n*1e-3)) - *Bx)/ *Bx);
		Bydev = fabsl(((B1 * sinl(x*gamma_n*1e-3)) - *By)); /// *By);
		Bzdev = fabsl((1.0e-3 - *Bz)/ *Bz);

		maxBxdev = fmaxl(maxBxdev,Bxdev);
		maxBydev = fmaxl(maxBydev,Bydev);
		maxBzdev = fmaxl(maxBzdev,Bzdev);
	}  
*/
return;
}


void BFderivs(long double x, long double *y, long double *dydx){

	long double Bx, By, Bz;
//	long double tmptime;
	//tmptime = clock();
	BFinterpol(x, &Bx, &By, &Bz);
	//timer3 += (long double)(((long double)clock() - (long double)tmptime)  / CLOCKS_PER_SEC);
	
	// dIdt = gamma_n I_n x B
	/*dydx[1]= gamma_n * (By * y[3] - Bz * y[2]);
	dydx[2]= gamma_n * (Bz * y[1] - Bx * y[3]);
	dydx[3]= gamma_n * (Bx * y[2] - By * y[1]);*/
	
	dydx[1]= gamma_n * (y[2] * Bz - y[3] * By);
	dydx[2]= gamma_n * (y[3] * Bx - y[1] * Bz);
	dydx[3]= gamma_n * (y[1] * By - y[2] * Bx);
	
	return;
}



/*
Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x,use
the fifth�order Cash�Karp Rung�Kutta method to advance the solution over an interval h
and return the incremented variables as yout[1..n]. Also return an estimate of the local
truncation error in yout
using the embedded fourth�order method. The user supplies the routine
derivs(x,y,dydx), which returns derivatives dydx at x.    */

void BFrkck(long double y[], long double dydx[], int n, long double x, long double h, long double yout[],long double yerr[], void (*BFderivs)(long double, long double [], long double [])){

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
	(*BFderivs)(x+a2*h,ytemp,ak2);                                //Second step.
	for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*BFderivs)(x+a3*h,ytemp,ak3);                                //Third step.
	for(i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*BFderivs)(x+a4*h,ytemp,ak4);                                //Fourth step.
	for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*BFderivs)(x+a5*h,ytemp,ak5);                                //Fifth step.
	for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*BFderivs)(x+a6*h,ytemp,ak6);                                //Sixth step.
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

void BFrkqs(long double y[], long double dydx[], int n, long double *x, long double htry, long double eps,long double yscal[], long double *hdid, long double *hnext, void (*BFderivs)(long double, long double [], long double [])){
	int i;
	long double errmax, h,htemp,xnew,*yerr,*ytemp;
	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;                                                           //Set stepsize to the initial trial value.
	for (;;) {
        BFrkck(y,dydx,n,*x,h,ytemp,yerr,BFderivs);                            //Take a step .
        errmax=0.0;                                                     //Evaluate accuracy.
        for (i=1;i<=n;i++) 
			errmax=fmaxl(errmax,fabsl(yerr[i]/yscal[i]));
		errmax /= eps;                                                 //Scale relative to required tolerance.
        if (errmax <= 1.0) 
			break;                                //Step succeeded. Compute size of next step.
        htemp=SAFETY*h*powl(errmax,PSHRNK);                             //Truncation error too large, reduce stepsize.
        h=(h >= 0.0 ? fmaxl(htemp,0.1*h) : fminl(htemp,0.1*h));          //No more than a factor of 10.
        xnew=(*x)+h;
        if (xnew == *x){ 
			h = 1.5 *h; 
			break; 
		}     // avoid stepsize underflow
        if (xnew == *x){
			nrerror("stepsize underflow in rkqs"); 
			fprintf(LOGSCR,"stepsize underflow in rkqs\n");
			}
	}
	if (errmax > ERRCON) 
		*hnext=SAFETY*h*powl(errmax,PGROW);
	else 
		*hnext=5.0*h;                                             //No more than a factor of 5 increase.
	*x += (*hdid=h);
	for (i=1;i<=n;i++) 
		y[i]=ytemp[i];
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
}


/*
Rung�Kutta driver with adaptive stepsize control. Integrate
starting values ystart[1..nvar] from x1 to x2 with accurac eps, storing intermediate results in global variables. h1 should
be set as a gessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On output nok and nbad are the number of good 
and bad (but retried and fixed) steps taken, and ystart is replaced b values at the end of the integation interval. derivs is the user�supplied
routine for calculating the rigt�hand side derivative, while rkqs isthenameofthestepper routine to be used.
*/
void BFodeintrk(long double ystart[], int nvar, long double x1, long double x2, long double eps, long double h1,long double hmin, int *nok, int *nbad,void (*BFderivs)(long double, long double [], long double []),void (*integrator)(long double [], long double [], int, long double *, long double, long double, long double [],long double *, long double *, void (*)(long double, long double [], long double [])))
{
	int nstp,i;
	long double xsav = 0,x = 0,hnext = 0,hdid = 0,h = 0;
	long double *yscal,*y,*dydx;
//	long double tmptime;
	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=(x2 > x1) ? fabsl(h1) : -fabsl(h1);
	//h=SIGN(h1,x2-x1);
	*nok = (*nbad) = BFkount = 0;
	for (i=1;i<=nvar;i++) 
		y[i]=ystart[i]; 
	if (BFkmax > 0) 
		xsav=x-BFdxsav*2.0;                                             //Assures storag of first step.
	for (nstp=1;nstp<=MAXSTP;nstp++) {                                          //Take at most MAXSTP steps.
		(*BFderivs)(x,y,dydx); 
		for (i=1;i<=nvar;i++)//Scaling used to monitor accuracy.
			yscal[i]=fabsl(y[i])+fabsl(dydx[i]*h)+TINY;
		if ((BFkmax > 0) && (BFkount < BFkmax-1) && (fabsl(x-xsav) > fabsl(BFdxsav))) {
			BFxp[++BFkount]=x;                                                              //Store intermediate results.
			for (i=1;i<=nvar;i++) 
				BFyp[i][BFkount]=y[i];
				//cout << "|";
			//RP
			//tmptime = clock();
			BFinterpol(BFxp[BFkount], &BFypFields[1][BFkount], &BFypFields[2][BFkount], &BFypFields[3][BFkount]);
			//timer3 += (long double)(((long double)clock() - (long double)tmptime)  / CLOCKS_PER_SEC);
			//RP End
			xsav=x;
		}
		if (BFkount==BFkmax)
			cout << "Too many points in intermediate array...!" << endl;
		if ((x+h-x2)*(x+h-x1) > 0.0) 
			h=x2-x;                                        //If stepsize can overshoot, decrease.
		(*integrator)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,BFderivs);                    // Successful Runge Kutte step
		if (hdid == h) 
			++(*nok); 
		else 
			++(*nbad);
		if ((x - x2)*(x2-x1) >= 0.0) {                                                //Are we done?
			for (i=1;i<=nvar;i++) 
				ystart[i]=y[i];
			if (BFkmax) {
				BFxp[++BFkount]=x;                                                              //Save final step.
				for (i=1;i<=nvar;i++) 
					BFyp[i][BFkount]=y[i];
				//RP
				//tmptime = clock();
				BFinterpol(BFxp[BFkount], &BFypFields[1][BFkount], &BFypFields[2][BFkount], &BFypFields[3][BFkount]);
				//timer3 += (long double)(((long double)clock() - (long double)tmptime)  / CLOCKS_PER_SEC);
				//RP End
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


// Routine um in einer 1D Array xx den Index zu finden, so dass der Wert x in
// der Array zwischen den Indizes n und n+1 liegt  (jlo ist ein first guess des indexes und der r�ckgabewert)
//Given an array
// xx[1..n], and given a value x,returns a value jlo such that x is between
//xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or decreasing.
// jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as the
//initialguess for jlo on output.
 void hunt(long double xx[], int n, long double x, int *jlo){
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

