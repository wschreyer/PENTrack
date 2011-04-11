//---------------------------------------------------------------------------
// user supplied equations for BruteForce integration of Bloch Equation

#include <cmath>
#include <sstream>
#include <vector>
#include <iomanip>

#include "nr/nr3.h"
#include "nr/stepper.h"
#include "nr/stepperdopr853.h"
#include "nr/odeint.h"

#include "bruteforce.h"
#include "globals.h"
#include "fields.h"

using namespace std;

FILE *BFLOG = NULL; 
//parameters for different experiment phases
int flipspin = 0;
vector<long double> BFtimes;

const long double BFdxsav=5e-7;  // timestep for intermediate output
long double BFBminmem=10.0;
long BFZeilencount, intsteps = 0;
int BFFilecount = 0; // to control output filename of BF
long double BFTargetB=0.1; // Babs < BFTargetB => integrate,
VecDoub *I_n = NULL; // spin vector

// struct used to calculate derivatives
struct TBFderivs{
	Doub xx1, xx2, c[4], cx[4], cy[4], cz[4]; // interpolation start time, end time, coefficients for cubic spline
	Doub b[3], R11, R12, R22;
	// constructor, pass integration step and field values
	TBFderivs(long double x1, VecDoub_I &y1, long double B1[4][4], long double x2, VecDoub_I &y2, long double B2[4][4]){
		xx1 = x1;
		xx2 = x2;

		long double dBxdt1 = B1[0][1]*y1[3] + B1[0][2]*y1[4] + B1[0][3]*y1[5];
		long double dBydt1 = B1[1][1]*y1[3] + B1[1][2]*y1[4] + B1[1][3]*y1[5];
		long double dBzdt1 = B1[2][1]*y1[3] + B1[2][2]*y1[4] + B1[2][3]*y1[5];

		long double dBxdt2 = B2[0][1]*y2[3] + B2[0][2]*y2[4] + B2[0][3]*y2[5];
		long double dBydt2 = B2[1][1]*y2[3] + B2[1][2]*y2[4] + B2[1][3]*y2[5];
		long double dBzdt2 = B2[2][1]*y2[3] + B2[2][2]*y2[4] + B2[2][3]*y2[5];
		// calculate coefficients of cubic splines Bi(x) = ci[0]*x^3 + ci[1]*x^2 + ci[2]*x + ci[3]   (i = x,y,z)
		// with boundary conditions Bi(x1) = Bi1, Bi(x2) = Bi2, dBidt(x1) = dBidt1, dBidt(x2) = dBidt2
		long double h = x2-x1;
		cx[0] = 2*B1[0][0] - 2*B2[0][0] + dBxdt1*h + dBxdt2*h;
		cx[1] = 3*B2[0][0] - 3*B1[0][0] - 2*dBxdt1*h - dBxdt2*h;
		cx[2] = dBxdt1*h;
		cx[3] = B1[0][0]; 
		cy[0] = 2*B1[1][0] - 2*B2[1][0] + dBydt1*h + dBydt2*h;
		cy[1] = 3*B2[1][0] - 3*B1[1][0] - 2*dBydt1*h - dBydt2*h;
		cy[2] = dBydt1*h;
		cy[3] = B1[1][0]; 
		cz[0] = 2*B1[2][0] - 2*B2[2][0] + dBzdt1*h + dBzdt2*h;
		cz[1] = 3*B2[2][0] - 3*B1[2][0] - 2*dBzdt1*h - dBzdt2*h;
		cz[2] = dBzdt1*h;
		cz[3] = B1[2][0]; 
	};
	// interpolate B
	void Binterp(long double x, long double B[3]){ 
		x = (x-xx1)/(xx2-xx1);
		B[0] = ((cx[0]*x + cx[1])*x + cx[2])*x + cx[3]; 
		B[1] = ((cy[0]*x + cy[1])*x + cy[2])*x + cy[3];
		B[2] = ((cz[0]*x + cz[1])*x + cz[2])*x + cz[3];
	};
	// integrator calls BFderivs(x,y,dydx) to get derivatives
	void operator()(Doub x, VecDoub_I &y, VecDoub_O &dydx){
		long double B[3];
		Binterp(x,B);
		dydx[0] = -gamma_n * (y[1] * B[2] - y[2] * B[1]);
		dydx[1] = -gamma_n * (y[2] * B[0] - y[0] * B[2]);
		dydx[2] = -gamma_n * (y[0] * B[1] - y[1] * B[0]);
	};
};

// print integration steps to file
void PrintBFStep(Output &out, TBFderivs &BFderivs){
	for (int i = 0; i < out.count; i++){
		if (BFZeilencount>50000){
			fclose(BFLOG);
			BFLOG = NULL;
		}
	
		if (!BFLOG){
			BFFilecount++;
			ostringstream BFoutfile1;
			BFoutfile1 << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "BF" << setw(3) << BFFilecount << setw(0) << ".out";
			if(!(BFLOG = fopen(BFoutfile1.str().c_str(),"w"))) 
			{
				perror("fopen");
				exit(1);
			}
			fprintf(BFLOG,"t Babs Polar logPolar Ix Iy Iz Bx By Bz\n");
			printf(" ##");
			printf("%s", BFoutfile1.str().c_str());
			printf("## ");
			BFZeilencount=1;
		}
	
		long double B[3];
		BFderivs.Binterp(out.xsave[i],B);
		long double BFBws =  sqrtl(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
		long double BFpol = (out.ysave[0][i]*B[0] + out.ysave[1][i]*B[1] + out.ysave[2][i]*B[2])/BFBws;
		long double BFlogpol = 0;
		if (BFpol<0.5) 
			BFlogpol = log10l(0.5-BFpol);
		else if (BFpol==0.5) 
			BFlogpol = 0.0;
		fprintf(BFLOG,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG\n",
						out.xsave[i],BFBws,BFpol,BFlogpol,2*out.ysave[0][i],2*out.ysave[1][i],2*out.ysave[2][i],B[0]/BFBws,B[1]/BFBws,B[2]/BFBws);
		
		BFZeilencount++;
	}
}

// calculate spin flip probability by integrating Bloch equations
long double BruteForceIntegration(long double x1, VecDoub_I &y1, long double B1[4][4], long double x2, VecDoub_I &y2, long double B2[4][4], long double xmax){
	bool BruteForce1 = false, BruteForce2 = false;;
	for (unsigned int i = 0; i < BFtimes.size(); i += 2){
		BruteForce1 |= (x1 >= BFtimes[i] && x1 < BFtimes[i+1]);
		BruteForce2 |= (x2 >= BFtimes[i] && x2 < BFtimes[i+1]);
	}
	if (BruteForce1 || BruteForce2){
		BFBminmem = min(BFBminmem,min(B1[3][0],B2[3][0])); // save lowest value for info
		
		// check if this value is worth for Bloch integration
		if (B1[3][0] < BFTargetB || B2[3][0] < BFTargetB){
			if (I_n == NULL){ 
				I_n = new VecDoub(3,0.0); // fill spin vector with values
				if (B1[3][0] > 0){
					(*I_n)[0] = B1[0][0]/B1[3][0]*0.5;
					(*I_n)[1] = B1[1][0]/B1[3][0]*0.5;
					(*I_n)[2] = B1[2][0]/B1[3][0]*0.5;
				}
				else
					(*I_n)[2] = 0.5;
				printf(" BF starttime %.6LG ",x1);
				fflush(stdout);
			}
	
			TBFderivs BFderivs(x1,y1,B1,x2,y2,B2); // create derivs-struct
			int nsave = 0;
			if (ausgabewunsch==OUTPUT_EVERYTHINGandSPIN || ausgabewunsch==OUTPUT_ENDPOINTSandSPIN)
				nsave = (int)((x2-x1)/BFdxsav);
			Output out(nsave); // create output object (save integrations steps in BFdxsav intervals)
			// create integrator<stepper<derivobject> >, params: spin vector, start time, end time, atol, rtol, first step size, min step size, output object, derivobject
			Odeint<StepperDopr853<TBFderivs> > ode(*I_n, x1, x2, 1e-13, 0, 1e-9, 0, out, BFderivs);
			ode.integrate(); // integrate
			intsteps += ode.nok + ode.nbad; // add up integration steps
			if ((ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)||(ausgabewunsch==OUTPUT_ENDPOINTSandSPIN))
				PrintBFStep(out,BFderivs); // print integrations steps
			
			if (B2[3][0] > BFTargetB || !BruteForce2){
				// output of polarisation after BF int completed
				// calculate polarisation at end of step BFpol = (I_n*B/|B|) in [-1/2,1/2]
				long double BFpol = ((*I_n)[0]*B2[0][0] + (*I_n)[1]*B2[1][0] + (*I_n)[2]*B2[2][0])/B2[3][0];
		
				printf(" BF endtime %LG, BFflipprop %.17LG, BFpol %.17LG, intsteps taken %ld, Bmin %LG\n",
						x2, 1-(BFpol+0.5), BFpol, intsteps, BFBminmem);
		
				BFBminmem=10; // reset values when done
				intsteps = 0;
				delete I_n;
				I_n = NULL;
	
				return BFpol + 0.5;
			}
		}
	}	
	
	return 1;	
}



