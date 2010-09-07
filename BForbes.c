#include <iostream>
#include <cmath>
#include "BForbes.h"
#include "globals.h"
#include "nrutil.h"

using namespace std;

TF1 *ForbesBrInt1;
TF1 *ForbesBzInt1;
TF1 *ForbesBzInt2;
TF1 *ForbesBzInt3;
TF1 *ForbesBzInt4;
TF1 *ForbesBrR;
TF1 *ForbesBrZ;
TF1 *ForbesBzR;
TF1 *ForbesBzZ;

// coil data for Forbes method
long double rFo, phiFo, zFo, aFo,  bFo,  R_0Fo,  J_0Fo, zoffsetFo;
int sign1 = 1, sign2 = 1;
long double C1a = 0.02, C1b=0.015, C1R_0 = 0.515, C1J_0=3e8, C1zoffset = 0;
long double aF[100], bF[100], R_0[100], zoffset[100], J_0[100];
int CoilNr=0;   // number of coils read in

// Forbes: "Rapid computation of static fields produced by thick circular solenoids"


// read the conductor data from a .cond file written by Vectorfield Opera
int ReadMagnets(const char *path)
{
	long double CoilAngle1=-1, CoilAngle2=-1, CoilAngle3=-1;
	long double dummy, CoilZOffset, Symmetry;
	long double CoilPoint1R, CoilPoint1Z,  CoilPoint2R,  CoilPoint2Z,  CoilPoint3R,  CoilPoint3Z,  CoilPoint4R,  CoilPoint4Z;
	
	
	FILE *cfg = NULL;
	
	char line[1024];
	cfg = fopen(path, "r");
	if(cfg == NULL){  
		Log("File coils.cond could not be opened!!!");
		return 1;
	}
	
	do // read the file
	{
		//for (i=1;i<=1024;i++)			
		//	line[i]='0'; // reset string
		
		fgets(line,1024,cfg); // get a line
		
		//if((line[0]=='C')||(line[0]=='Q')) continue; // if its a comment goto next line
			
		if((line[0]=='D')&&(line[7]=='G')&&(line[8]=='R')&&(line[9]=='A'))
		{cout << "RACETRACk COIL DETECTED" << endl; continue; }// discard racetrack coils
			
		if((line[0]=='D')&&(line[7]=='G')&&(line[8]=='S')&&(line[9]=='O'))
		{
			//printf(line);
			fgets(line,1024,cfg); // 1st line of coil definition: local coordinate system 1, only used for racetrack so far
		
			fgets(line,1024,cfg); // 2nd line of coil definition: local coordinat system 2 offsets
			sscanf(line,"%25Lf %25Lf %25Lf", &dummy, &dummy, &CoilZOffset);
		
			fgets(line,1024,cfg); // 3rd line of coil definition: local coordinat system 2 angles
			sscanf(line,"%25Lf %25Lf %25Lf", &CoilAngle1, &CoilAngle2, &CoilAngle3);
			if((CoilAngle1!=90)||(CoilAngle2!=0)||(CoilAngle3!=90))  // check if coils have right orientation: z as rotational axis
			{
				cout << "CoilAngle1 "	<< CoilAngle1 << "CoilAngle2 "<< CoilAngle2 << "CoilAngle3 " << CoilAngle3 << endl;
				cout << "Reading of coil number "	<< CoilNr << " failed!!! Wrong orientation!!! Aborting!!" << endl;
				return 1;
			}
		
			fgets(line,1024,cfg); // 4th line of coil definition: cross section points 1 and 2 
			sscanf(line,"%25Lf %25Lf %25Lf %25Lf", &CoilPoint1R, &CoilPoint1Z, &CoilPoint2R, &CoilPoint2Z);
			
			fgets(line,1024,cfg); // 5th line of coil definition: cross section points 3 and 4 
			sscanf(line,"%25Lf %25Lf %25Lf %25Lf", &CoilPoint3R, &CoilPoint3Z, &CoilPoint4R, &CoilPoint4Z);
			if (CoilPoint3Z!=CoilPoint1Z)
				aF[CoilNr]=fabs((CoilPoint3Z-CoilPoint1Z)/2.0);
			else aF[CoilNr]=fabs((CoilPoint2Z-CoilPoint1Z)/2.0);
				
			if (CoilPoint2R!=CoilPoint1R)
				bF[CoilNr]=fabs((CoilPoint2R-CoilPoint1R)/2.0);
			else bF[CoilNr]=fabs((CoilPoint3R-CoilPoint1R)/2.0);
				
			if (CoilPoint1R<CoilPoint2R)
				R_0[CoilNr]= CoilPoint1R+(CoilPoint2R-CoilPoint1R)/2.0;
			
			zoffset[CoilNr]=CoilZOffset+CoilPoint1Z+(CoilPoint3Z-CoilPoint1Z)/2.0;
			
			fgets(line,1024,cfg); // 6th line of coil definition: Curvatures -> not used
			
			fgets(line,1024,cfg); // 7th line of coil definition: current density, symmetry, label
			sscanf(line,"%25Lf %25Lf ", &J_0[CoilNr], &Symmetry);
			if ((Symmetry!=0)&&(Symmetry!=1))
			{
				cout << "Reading of coil number "	<< CoilNr << " failed!!! Wrong symmetry!!! Aborting!!" << endl;
				return 1;
			}
			
			fgets(line,1024,cfg); // 8th line of coil definition: reflection codes -> not used
			
			fgets(line,1024,cfg); // 9th line of coil definition: accuracy of field calculation in Opera -> not used		
			
			cout << "Coil " << CoilNr << " a " << aF[CoilNr] << " b " << bF[CoilNr] << " R_0 " << R_0[CoilNr] << " zoffset " << zoffset[CoilNr] << " J_0 " << J_0[CoilNr] << endl ;
			CoilNr++;
		}
			
			// printf(line);
		
	}while(!feof(cfg));
	
	Log("%i solenoids were read in!!! \n",CoilNr);
	
	return 0;
}

// parameter array for Coils:
//[R_0(radius of center),a (height),b (radial dimension),r (point to calculate field),z (point to calculate),sign1 (sign1 in integregral of Br),sign2 (sign2 in integregral of Br),J_0 (current density),mu0 (clear))]
long double BForbes(long double r, long double phi, long double z, long double t, long double Bi[3], long double dBidrj[3][3])
{
	if (r==0)
			r=1e-8;   // solve singularity
	
	
	ForbesBrInt1 = new TF1("ForbesBrInt1",BrFoIntRoot,0,4,9);
	ForbesBrR = new TF1("ForbesBrR",BrRootR,0,3,9);	
	ForbesBrZ = new TF1("ForbesBrZ",BrRootZ,-2,3,9);	
	ForbesBzInt1 = new TF1("ForbesBzInt1",BzFoInt1Root,0,4,9);
	ForbesBzInt2 = new TF1("ForbesBzInt2",BzFoInt2Root,0,4,9);
	ForbesBzInt3 = new TF1("ForbesBzInt3",BzFoInt3Root,0,4,9);
	ForbesBzInt4 = new TF1("ForbesBzInt4",BzFoInt4Root,0,4,9);
	ForbesBzR = new TF1("ForbesBzR",BzRootR,0,3,9);	
	ForbesBzZ = new TF1("ForbesBzZ",BzRootZ,-2,3,9);	
	//ForbesBrInt1->SetParameters(R_0Fo,aFo,bFo,rFo,zFo,sign1,sign2,J_0Fo,mu0);
	
						
	//long double Integral = ForbesBrInt1->Integral(0,1);
	//cout << "Integral = " << Integral << endl;
	
	//long double a = height/2.0, b = rthickness/2.0, R_0 = innerradius + height/2.0;
	//Bphi = RodFieldMultiplicator * BFeldSkal * Bphi;
	//dBphidr = RodFieldMultiplicator * BFeldSkal * dBphidr;
	//dBphidz = RodFieldMultiplicator * BFeldSkal * dBphidz;
	
	int i;
	for (i=0;i<CoilNr;i++){
	//for (i=7;i<8;i++){
		OneCoilRoot(r, phi, z, aF[i], bF[i], R_0[i], J_0[i], zoffset[i], Bi, dBidrj);	
	//OneCoil(r, phi, z, aF[i], bF[i], R_0[i], J_0[i], zoffset[i]);		
	}	
	
	//OneCoilRoot(r, phi, z, aF[7], bF[7], R_0[7], J_0[7], zoffset[7]);	
	//OneCoilRoot(r, phi, z, aF[25], bF[25], R_0[25], J_0[25], zoffset[25]);
	
	Bi[1] = 0.0;
	dBidrj[0][1] = 0.0;
	dBidrj[1][0] = 0.0;
	dBidrj[1][1] = 0.0;
	dBidrj[1][2] = 0.0;
	dBidrj[2][1] = 0.0;
			
	return Bi[0];
	
}


// ROOT VERSION


// Br of one coil in root nomenclature
Double_t BrRootR(Double_t *x, Double_t *par)
{
	
	// compute Br with Forbes method	
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBrInt1->SetParameter(5,1);  // set sign1 to 1
	ForbesBrInt1->SetParameter(6,1); // set sign2 to 1
	long double t3 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBrInt1->SetParameter(5,-1);
	ForbesBrInt1->SetParameter(6,1);
	long double t4 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBrInt1->SetParameter(5,1);
	ForbesBrInt1->SetParameter(6,-1);
	long double t5 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBrInt1->SetParameter(5,-1);
	ForbesBrInt1->SetParameter(6,-1);
	long double t6 = ForbesBrInt1->Integral(0,0.3141592654e1);
		
	return -par[8] * par[7] / 0.3141592654e1 * (t3 - t4 - t5 + t6) / 0.2e1;
}

// Br of one coil in root nomenclature
Double_t BrRootZ(Double_t *x, Double_t *par)
{
	
	
	// compute Br with Forbes method	
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBrInt1->SetParameter(5,1);  // set sign1 to 1
	ForbesBrInt1->SetParameter(6,1); // set sign2 to 1
	long double t3 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBrInt1->SetParameter(5,-1);
	ForbesBrInt1->SetParameter(6,1);
	long double t4 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBrInt1->SetParameter(5,1);
	ForbesBrInt1->SetParameter(6,-1);
	long double t5 = ForbesBrInt1->Integral(0,0.3141592654e1);
	ForbesBrInt1->SetParameters(par);
	ForbesBrInt1->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBrInt1->SetParameter(5,-1);
	ForbesBrInt1->SetParameter(6,-1);
	long double t6 = ForbesBrInt1->Integral(0,0.3141592654e1);
		
	return -par[8] * par[7] / 0.3141592654e1 * (t3 - t4 - t5 + t6) / 0.2e1;
}


// Bz of one coil in root nomenclature
Double_t BzRootR(Double_t *x, Double_t *par)
{
	
	// compute Bz with Forbes method	
	ForbesBzInt1->SetParameters(par);
	ForbesBzInt1->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBzInt2->SetParameters(par);
	ForbesBzInt2->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBzInt3->SetParameters(par);
	ForbesBzInt3->SetParameter(3,x[0]);  // set the parameter of the integral to r
	ForbesBzInt4->SetParameters(par);
	ForbesBzInt4->SetParameter(3,x[0]);  // set the parameter of the integral to r
		
	
	long double t5 = ForbesBzInt1->Integral(0,0.3141592654e1);
	long double t6 = ForbesBzInt2->Integral(0,0.3141592654e1);
	long double t7 = ForbesBzInt3->Integral(0,0.3141592654e1);
	long double t8 = ForbesBzInt4->Integral(0,0.3141592654e1);
	
	return par[8] * par[7] / 0.3141592654e1 / x[0] * (t5 - t6 - t7 + t8) / 0.2e1;	
}

// Bz of one coil in root nomenclature
Double_t BzRootZ(Double_t *x, Double_t *par)
{
	
	// compute Bz with Forbes method	
	ForbesBzInt1->SetParameters(par);
	ForbesBzInt1->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBzInt2->SetParameters(par);
	ForbesBzInt2->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBzInt3->SetParameters(par);
	ForbesBzInt3->SetParameter(4,x[0]);  // set the parameter of the integral to z
	ForbesBzInt4->SetParameters(par);
	ForbesBzInt4->SetParameter(4,x[0]);  // set the parameter of the integral to z
		
	
	long double t5 = ForbesBzInt1->Integral(0,0.3141592654e1);
	long double t6 = ForbesBzInt2->Integral(0,0.3141592654e1);
	long double t7 = ForbesBzInt3->Integral(0,0.3141592654e1);
	long double t8 = ForbesBzInt4->Integral(0,0.3141592654e1);
	
	return mu0 * J_0Fo / 0.3141592654e1 / par[3] * (t5 - t6 - t7 + t8) / 0.2e1;	
}

// integral in Br fuction
Double_t BrFoIntRoot(Double_t *beta, Double_t *par)
{
	long double R_0Fo = par[0], aFo = par[1], bFo = par[2], rFo = par[3], zFo = par[4], sign1= par[5], sign2 = par[6];
	long double cb = cos(beta[0]), sb = sin(beta[0]);	
	long double wurzel = sqrtl(powl(R_0Fo + sign1 * bFo - rFo * cb, 0.2e1) + rFo * rFo * powl(sb, 0.2e1) + powl(zFo+ sign2 * aFo, 0.2e1));
	return cb * (wurzel + rFo * cb * logl(R_0Fo + sign1 * bFo - rFo* cb + wurzel));
	//return sin(beta[0]);
}

// integral 1 in Bz fuction
Double_t BzFoInt1Root(Double_t *beta, Double_t *par)
{
	long double x = beta[0], R_0Fo = par[0], aFo = par[1], bFo = par[2], rFo = par[3], zFo = par[4];
	long double cb = cos(x), sb = sin(x);	
	return cb  * (0.2e1 * rFo *  cb  * (aFo - zFo) * log(R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo + bFo, 0.2e1)) * log((sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) - aFo + zFo) / (sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) + aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan(((R_0Fo + bFo - rFo *  cb ) * (R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) / rFo /  sb   / (aFo - zFo)) + atan((aFo - zFo + R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) / rFo /  sb  ) - atan((-aFo + zFo + R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) / rFo /  sb  )) + 0.5e0 * (aFo - zFo) * sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)));
}

Double_t BzFoInt2Root(Double_t *beta, Double_t *par)
{
	long double x = beta[0], R_0Fo = par[0], aFo = par[1], bFo = par[2], rFo = par[3], zFo = par[4];
	long double cb = cos(x), sb = sin(x);	
	return  cb  * (0.2e1 * rFo *  cb  * (aFo - zFo) * log(R_0Fo - bFo - rFo *  cb  + sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo - bFo, 0.2e1)) * log((sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) - aFo + zFo) / (sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) + aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan(((R_0Fo - bFo - rFo *  cb ) * (R_0Fo - bFo - rFo *  cb  + sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) / rFo /  sb   / (aFo - zFo)) + atan((aFo - zFo + R_0Fo - bFo - rFo *  cb  + sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) / rFo /  sb  ) - atan((-aFo + zFo + R_0Fo - bFo - rFo *  cb  + sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) / rFo /  sb  )) + 0.5e0 * (aFo - zFo) * sqrt(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)));
}

Double_t BzFoInt3Root(Double_t *beta, Double_t *par)
{
	long double x = beta[0], R_0Fo = par[0], aFo = par[1], bFo = par[2], rFo = par[3], zFo = par[4];
	long double cb = cos(x), sb = sin(x);	
	return cb  * (0.2e1 * rFo *  cb  * (-aFo - zFo) * log(R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo + bFo, 0.2e1)) * log((sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) + aFo + zFo) / (sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) - aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan(((R_0Fo + bFo - rFo *  cb ) * (R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) / rFo /  sb   / (-aFo - zFo)) + atan((-aFo - zFo + R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) / rFo /  sb  ) - atan((aFo + zFo + R_0Fo + bFo - rFo *  cb  + sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) / rFo /  sb  )) + 0.5e0 * (-aFo - zFo) * sqrt(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)));
}

Double_t BzFoInt4Root(Double_t *beta, Double_t *par)
{
	long double x = beta[0], R_0Fo = par[0], aFo = par[1], bFo = par[2], rFo = par[3], zFo = par[4];
	long double cb = cos(x), sb = sin(x);	
	return cb * (0.2e1 * rFo * cb * (-aFo - zFo) * log(R_0Fo - bFo - rFo * cb + sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl(cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo - bFo, 0.2e1)) * log((sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) + aFo + zFo) / (sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) - aFo - zFo)) + rFo * rFo * cb *  sb   * (-0.2e1 * atan(((R_0Fo - bFo - rFo * cb ) * (R_0Fo - bFo - rFo * cb + sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) / rFo /  sb   / (-aFo - zFo)) + atan((-aFo - zFo + R_0Fo - bFo - rFo * cb + sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) / rFo /  sb  ) - atan((aFo + zFo + R_0Fo - bFo - rFo * cb + sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) / rFo /  sb  )) + 0.5e0 * (-aFo - zFo) * sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)));
}


long double OneCoilRoot(long double r, long double phi, long double z, long double a, long double b, long double R_0, long double J_0, long double zoffset,
						long double Bi[3], long double dBidrj[3][3])
{
	
	ForbesBrR->SetParameters(R_0,a,b,r,z-zoffset,1,1,J_0,mu0);
	ForbesBrZ->SetParameters(R_0,a,b,r,z-zoffset,1,1,J_0,mu0);
	/*for (int i=0;i<9;i++)
	{
	cout << endl << "Parameter " << i << " = " << ForbesBrR->GetParameter(i) << endl;
	}*/
	Bi[0] += ForbesBrR->Eval(r,0,0,0);
	ForbesBzR->SetParameters(R_0,a,b,r,z-zoffset,1,1,J_0,mu0);
	ForbesBzZ->SetParameters(R_0,a,b,r,z-zoffset,1,1,J_0,mu0);
	Bi[2] += ForbesBzR->Eval(r,0,0,0);	
	
	dBidrj[0][0] += ForbesBrR->Derivative(r);
	dBidrj[0][2] += ForbesBrZ->Derivative(z-zoffset);
	dBidrj[2][0] += ForbesBzR->Derivative(r);
	dBidrj[2][2] += ForbesBzZ->Derivative(z-zoffset);
	
	return Bi[0];

}

// END ROOT VERSION
































// NUMERICAL RECIPES VERSION

long double BrFoInt(long double beta)
{
	long double cb = cosl(beta), sb = sinl(beta);	
	long double wurzel = sqrt(powl(R_0Fo + sign1 * bFo - rFo * cb, 0.2e1) + rFo * rFo * powl(sb, 0.2e1) + powl(zFo+ sign2 * aFo, 0.2e1));
	return cb * (wurzel + rFo * cb * log(R_0Fo + sign1 * bFo - rFo* cb + wurzel));
}

long double BzFoInt1(long double beta)
{
	long double cb = cosl(beta), sb = sinl(beta);
	return cb  * (0.2e1 * rFo *  cb  * (aFo - zFo) * logl(R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo + bFo, 0.2e1)) * logl((sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) - aFo + zFo) / (sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) + aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan2l(((R_0Fo + bFo - rFo *  cb ) * (R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) , rFo *  sb   * (aFo - zFo)) + atan2l((aFo - zFo + R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) , rFo *  sb  ) - atan2l((-aFo + zFo + R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) ,rFo *  sb  )) + 0.5e0 * (aFo - zFo) * sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)));
}

long double BzFoInt2(long double beta)
{
	long double cb = cosl(beta), sb = sinl(beta);
	return  cb  * (0.2e1 * rFo *  cb  * (aFo - zFo) * logl(R_0Fo - bFo - rFo *  cb  + sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo - bFo, 0.2e1)) * logl((sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) - aFo + zFo) / (sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)) + aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan2l(((R_0Fo - bFo - rFo *  cb ) * (R_0Fo - bFo - rFo *  cb  + sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) , rFo *  sb   * (aFo - zFo)) + atan2l((aFo - zFo + R_0Fo - bFo - rFo *  cb  + sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) , rFo *  sb  ) - atan2l((-aFo + zFo + R_0Fo - bFo - rFo *  cb  + sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1))) , rFo *  sb  )) + 0.5e0 * (aFo - zFo) * sqrtl(powl(R_0Fo - bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo - aFo, 0.2e1)));
}

long double BzFoInt3(long double beta)
{
	long double cb = cosl(beta), sb = sinl(beta);
	return cb  * (0.2e1 * rFo *  cb  * (-aFo - zFo) * logl(R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl( cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo + bFo, 0.2e1)) * logl((sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) + aFo + zFo) / (sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) - aFo - zFo)) + rFo * rFo *  cb  *  sb   * (-0.2e1 * atan2l(((R_0Fo + bFo - rFo *  cb ) * (R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) , rFo /  sb   / (-aFo - zFo)) + atan2l((-aFo - zFo + R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) , rFo *  sb  ) - atan2l((aFo + zFo + R_0Fo + bFo - rFo *  cb  + sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) , rFo *  sb  )) + 0.5e0 * (-aFo - zFo) * sqrtl(powl(R_0Fo + bFo - rFo *  cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)));
}

long double BzFoInt4(long double beta)
{
	long double cb = cosl(beta), sb = sinl(beta);
	return cb * (0.2e1 * rFo * cb * (-aFo - zFo) * logl(R_0Fo - bFo - rFo * cb + sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + 0.25e0 * (0.3e1 * rFo * rFo * powl(cb , 0.2e1) - 0.3e1 * rFo * rFo * powl( sb  , 0.2e1) - powl(R_0Fo - bFo, 0.2e1)) * logl((sqrt(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) + aFo + zFo) / (sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)) - aFo - zFo)) + rFo * rFo * cb *  sb   * (-0.2e1 * atan2l(((R_0Fo - bFo - rFo * cb ) * (R_0Fo - bFo - rFo * cb + sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) + rFo * rFo * powl( sb  , 0.2e1)) , rFo *  sb   * (-aFo - zFo)) + atan2l((-aFo - zFo + R_0Fo - bFo - rFo * cb + sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) , rFo *  sb  ) - atan2l((aFo + zFo + R_0Fo - bFo - rFo * cb + sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1))) , rFo *  sb  )) + 0.5e0 * (-aFo - zFo) * sqrtl(powl(R_0Fo - bFo - rFo * cb , 0.2e1) + rFo * rFo * powl( sb  , 0.2e1) + powl(zFo + aFo, 0.2e1)));
}

long double OneCoilBrR(long double rT)
{
	
	rFo =	rT;
	// compute Br with Forbes method	
	sign1 = +1;  sign2 =+1;  
	long double t3 = qromb(BrFoInt,0,pi);	
	sign1 = -1;  sign2 =+1;  
	long double t4 = qromb(BrFoInt,0,pi);	
	sign1 = +1;  sign2 =-1;  
	long double t5 = qromb(BrFoInt,0,pi);	
	sign1 = -1;  sign2 =-1;  
	long double t6 = qromb(BrFoInt,0,pi);	
	
	return -mu0 * J_0Fo / 0.3141592654e1 * (t3 - t4 - t5 + t6) / 0.2e1;	
	
}

long double OneCoilBzR(long double rT)
{
	rFo =	rT;
	// compute Bz with Forbes method
	
	long double t5 = qromb(BzFoInt1,0,pi);
	long double t6 = qromb(BzFoInt2,0,pi);
	long double t7 = qromb(BzFoInt3,0,pi);
	long double t8 = qromb(BzFoInt4,0,pi);
	
	return mu0 * J_0Fo / 0.3141592654e1 / rFo * (t5 - t6 - t7 + t8) / 0.2e1;
	
}

long double OneCoilBrZ(long double zT)
{
	zFo =zT;
	// compute Br with Forbes method	
	sign1 = +1;  sign2 =+1;  
	long double t3 = qromb(BrFoInt,0,pi);	
	sign1 = -1;  sign2 =+1;  
	long double t4 = qromb(BrFoInt,0,pi);	
	sign1 = +1;  sign2 =-1;  
	long double t5 = qromb(BrFoInt,0,pi);	
	sign1 = -1;  sign2 =-1;  
	long double t6 = qromb(BrFoInt,0,pi);	
	return -mu0 * J_0Fo / 0.3141592654e1 * (t3 - t4 - t5 + t6) / 0.2e1;	
	
}

long double OneCoilBzZ(long double zT)
{
	zFo =zT;
	// compute Bz with Forbes method
	
	long double t5 = qromb(BzFoInt1,0,pi);
	long double t6 = qromb(BzFoInt2,0,pi);
	long double t7 = qromb(BzFoInt3,0,pi);
	long double t8 = qromb(BzFoInt4,0,pi);
	
	return mu0 * J_0Fo / 0.3141592654e1 / rFo * (t5 - t6 - t7 + t8) / 0.2e1;
}

long double OneCoil(long double r, long double phi, long double z, long double a, long double b, long double R_0, long double J_0, long double zoffset,
					long double Bi[3], long double dBidrj[3][3])
{
	rFo = r, phiFo = phi, zFo = z-zoffset, aFo = a; bFo = b; R_0Fo = R_0; J_0Fo = J_0, zoffsetFo = zoffset;
	//				^ compensate coil position in z
		long double err;
	
	Bi[0] += OneCoilBrR(rFo);
	Bi[2] += OneCoilBzR(rFo);
	
	
	dBidrj[0][0] += dfridr(OneCoilBrR, rFo, 0.01, &err);
	dBidrj[0][2] += dfridr(OneCoilBrZ, zFo, 0.01, &err);
	dBidrj[2][0] += dfridr(OneCoilBzR, rFo, 0.01, &err);
	dBidrj[2][2] += dfridr(OneCoilBzZ, zFo, 0.01, &err);
	return Bi[2];

}

#define FUNC(x) ((*func)(x))
long double trapzd(long double (*func)(long double), long double a, long double b, int n)
//This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
//as a pointer to the function to be integrated between limits a and b, also input. When called with
//n=1, the routine returns the crudest estimate of b
//a f(x)dx. Subsequent calls with n=2,3,...
//(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
{
	long double x,tnm,sum,del;
	static long double s;
	int it,j;
	if (n == 1) 
	{
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} 
	else 
	{
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm; //This is the spacing of the points to be added.
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm); //This replaces s by its refined value.
		return s;
	}
}


#include <math.h>
#define EPS 1.0e-8
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

//Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
//JMAX limits the total number of steps; K is the number of points used in the extrapolation.
long double qromb(long double (*func)(long double), long double a, long double b)
//Returns the integral of the function func from a to b. Integration is performed by Romberg�s
//method of order 2K, where, e.g., K=2 is Simpson�s rule.
{
	void BFpolint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy);
	long double trapzd(long double(*func)(long double), long double a, long double b, int n);
	void nrerror(char error_text[]);
	long double ss,dss;
	long double s[JMAXP],h[JMAXP+1];  //These store the successive trapezoidal approximations and their relative stepsizes 
	int j; 
	h[1]=1.0;
	
	for (j=1;j<=JMAX;j++) 
	{
		s[j]=trapzd(func,a,b,j);
		if (j >= K) 
		{
			BFpolint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabsl(dss) <= EPS*fabsl(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
		//This is a key step: The factor is 0.25 even though the stepsize is decreased by only
		//0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1),
		// not just a polynomial in h.
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;  //Never get here.
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define CON 1.4 //Stepsize is decreased by CON at each iteration.
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10 //Sets maximum size of tableau.
#define SAFE 2.0 //Return when error is SAFE worse than the best so far.
long double dfridr(long double (*func)(long double), long double x, long double h, long double *err)
//Returns the derivative of a function func at a point x by Ridders� method of polynomial
//extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but
//rather should be an increment in x over which func changes substantially. An estimate of the
//error in the derivative is returned as err.
{
int i,j;
long double errt,fac,hh,**a,ans;
if (h == 0.0) nrerror("h must be nonzero in dfridr.");
a=matrix(1,NTAB,1,NTAB);
hh=h;
a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
*err=BIG;
for (i=2;i<=NTAB;i++) {
//Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of
//extrapolation.
hh /= CON;
a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh); //Try new, smaller stepsize
fac=CON2; 
for (j=2;j<=i;j++) { //Compute extrapolations of various orders, requiring
				//no new function evaluations.
a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
fac=CON2*fac;
errt=fmaxl(fabsl(a[j][i]-a[j-1][i]),fabsl(a[j][i]-a[j-1][i-1]));
	
//The error strategy is to compare each new extrapolation to one order lower, both
//at the present stepsize and the previous one.
if (errt <= *err) { //If error is decreased, save the improved answer.
*err=errt;
ans=a[j][i];
}
}
if (fabsl(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
//If higher order is worse by a significant factor SAFE, then quit early.
}
free_matrix(a,1,NTAB,1,NTAB);
return ans;
}


// END NUMERICAL RECIPES VERSION
