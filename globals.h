#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <cstdio>

#define ID_UNKNOWN 0 ///< standard kennz flag for particles
#define ID_NOT_FINISH -1 ///< kennz flag for particles which reached ::StorageTime
#define ID_HIT_BOUNDARIES -2 ///< kennz flag for particles which left bounding box of TParticle::geom
#define ID_NRERROR -3 ///< kennz flag for particles which produces a serious numerical error (step size underflow, missed reflection, ...)
#define ID_DECAYED -4 ///< kennz flag for particles which reached TParticle::xend
#define ID_INITIAL_NOT_FOUND -5 ///< kennz flag for particles which had a too low total energy to find a initial spot in the source volume

#define NEUTRON 1 ///< TParticle::protneut of neutrons
#define PROTON 2 ///< TParticle::protneut of protons
#define BF_ONLY 3 ///< set protneut in configuration to this value to print out a ramp heating analysis
#define BF_CUT 4 ///< set protneut in configuration to this value to print out a planar slice through electric/magnetic fields
#define ELECTRON 6 ///< TParticle::protneut of electrons
#define GEOMETRY 7 ///< set protneut in configuration to this value to print out a sampling of the geometry

#define POLARISATION_GOOD 1 ///< Printed to endlog for low field seeking neutron
#define POLARISATION_BAD 2 ///< Printed to endlog for high field seeking neutron
#define POLARISATION_NONE 3 ///< Printed to endlog for neutron without magnetic moment

#define OUTPUT_EVERYTHING 1 ///< configuration value to print endlog+tracklog
#define OUTPUT_ENDPOINTS 2 ///< configuration value to print endlog only
#define OUTPUT_NOTHING 5
#define OUTPUT_EVERYTHINGandSPIN 3 ///< configuration value to print endlog+tracklog+spintrack
#define OUTPUT_ENDPOINTSandSPIN 4 ///< configuration value to print endlog+spintrack

// physical constants
const long double pi = 3.1415926535897932384626; ///< Pi
const long double ele_e = 1.602176487E-19; ///< elementary charge [C]
const long double gravconst = 9.80665; ///< g [m/s]
const long double conv = pi/180.; ///< deg to rad conversion factor
const long double mu0 = 4*pi*1e-7; ///< magnetic permeability [Vs/Am]
const long double m_n = 1.674927211E-27/ele_e; ///< neutron mass [eV/c^2]
const long double m_p = 1.672621637E-27/ele_e; ///< proton mass [eV/c^2]
const long double m_e = 9.10938215e-31/ele_e; ///< electron mass [eV/c^2]
const long double c_0 = 299792458; ///< light speed [m/s]
const long double hbar = 1.05457266e-34; ///< planck constant [Js]
const long double mu_nSI = -0.96623641e-26;	///< Neutron magnetic moment [J/T]
const long double gamma_n = -1.83247185e8; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]

const long double lengthconv = 0.01; ///< length conversion factor cgs -> SI [cm -> m]
const long double Bconv = 1e-4; ///< magnetic field conversion factor cgs -> SI [G -> T]
const long double Econv = 1e2; ///< electric field conversion factor cgs -> SI [V/cm -> V/m]

long double StorageTime = 1500.; ///< max. simulation time

int jobnumber = 0; ///< job number, read from command line paramters, used for parallel calculations
string inpath = "."; ///< path to configuration files, read from command line paramters
string outpath = "."; ///< path where the log file should be saved to, read from command line parameters
int outputopt = 5; ///< output options chosen by user

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


#endif /*GLOBALS_H_*/
