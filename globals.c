#include <cstdio>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "globals.h"

//set the duration of the experiment
long double FillingTime = 0;				// filling time, entrance open
long double CleaningTime = 0;				// cleaning without field
long double RampUpTime = 0;					// ramping up coils
long double FullFieldTime = 1000;			// storing in full field
long double RampDownTime = 5;				// ramping down coils
long double EmptyingTime = 0;				// emptying without field
long double StorageTime = 1500.0;			// time when ramping down shall start, if xend > storage time, let neutron decay

int jobnumber = 0;
string inpath = ".", outpath = ".";
int ausgabewunsch=5; // Ausgabewunsch

// print percent of (x-x1)/len
void percent(long double x, long double x1, long double len, int &perc){
	// write status to console
	// one point per 2 percent of endtime
	if((x - x1)*100/len > perc)
	{
		perc++;
		if(perc%10==0)
		{
			printf("%i%%",perc);
			fflush(stdout);
		}
		if((perc%10!=0)&&(perc%2==0)) 
		{
			printf(".");
			fflush(stdout);
		}
	}
}

// return experiment stage (filling, cleaning, rampup, storage, fullfield, rampdown, ...)
int ExpPhase(long double t){
	if ((t < FillingTime)&&(FillingTime>0))
	{      // filling in neutrons
		return 1;
	}	
	else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
	{      // spectrum cleaning
		return 2;
	}
	else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
	{   // ramping up field
		return 3;
	}
	else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
	{  // fullfield time
		return 4;
	}
	else if ((t >= (RampUpTime+CleaningTime+FillingTime+FullFieldTime)) && (t < (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime)) && (RampDownTime != 0))
	{   // ramping down field
		return 5;
	}
	else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
	{      // emptying of neutrons into detector
		return 6;
	}
	return 0;
}

// rotate vector into new coordinate sys whose z-axis lies on NORMALIZED vector n (active transformation)
void RotateVector(long double v[3], long double n[3])
{
	long double cosalpha = n[2], sinalpha = sqrt(1 - cosalpha*cosalpha);	// rotation angle (angle between z and n)
	if (sinalpha > 1e-30){ // when normal not parallel to z-axis rotate new velocity into the coordinate system where the normal is the z-axis
		long double a[2] = {-n[1]/sinalpha, n[0]/sinalpha};	// rotation axis (z cross n), a[2] = 0
		long double vtemp[3] = {v[0],v[1],v[2]};
		// rotate velocity vector
		v[0] = (cosalpha + a[0]*a[0]*(1 - cosalpha))*	vtemp[0] +  a[0]*a[1]*(1 - cosalpha)*				vtemp[1] + a[1]*sinalpha*	vtemp[2];
		v[1] =  a[1]*a[0]*(1 - cosalpha)*				vtemp[0] + (cosalpha + a[1]*a[1]*(1 - cosalpha))*	vtemp[1] - a[0]*sinalpha*	vtemp[2];
		v[2] = -a[1]*sinalpha*							vtemp[0] +  a[0]*sinalpha*							vtemp[1] + cosalpha*		vtemp[2];
	}
	else if (cosalpha < 0){
		v[0] = -v[0];
		v[1] = -v[1];
		v[2] = -v[2];
	}
}

//======== Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta ======================================================
void BOOST(long double beta[3], long double p[4]){
   //Boost this Lorentz vector (copy&paste from ROOT)
   long double b2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
   long double gamma = 1.0 / sqrt(1.0 - b2);
   long double bp = beta[0]*p[1] + beta[1]*p[2] + beta[2]*p[3];
   long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

   p[1] = (p[1] + gamma2*bp*beta[0] + gamma*beta[0]*p[0]);
   p[2] = (p[2] + gamma2*bp*beta[1] + gamma*beta[1]*p[0]);
   p[3] = (p[3] + gamma2*bp*beta[2] + gamma*beta[2]*p[0]);
   p[0] = (gamma*(p[0] + bp));
}
//======== end of BOOST ====================================================================================================

// energy distribution of protons (0 < 'Energie' < 750 eV)
// proton recoil spectrum from "Diplomarbeit M. Simson"
// result always < 1!
long double ProtonSpectrum(long double E){
	long double DeltaM = m_n - m_p;
	long double Xi = m_e / DeltaM;
	long double Sigma = 1 - 2*E*m_n / pow(DeltaM, 2) / pow(c_0, 2);
	long double g1 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 				4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	long double g2 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 	4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	long double littlea = -0.1017;
	return 0.5*(g1 + littlea * g2);
}

// energy distribution of electrons (0 < 'Energie' < 782 keV)
// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
// result always < 1!
long double ElectronSpectrum(long double E){
	long double Qvalue = 0.782; //[MeV]
	return 8.2*sqrt(E*1e-6*E*1e-6 + 2*E*1e-6*m_e*c_0*c_0*1e-6) * pow(Qvalue - E*1e-6, 2) * (E*1e-6 + m_e*c_0*c_0*1e-6);
}
