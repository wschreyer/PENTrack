#ifndef GLOBALS_H_
#define GLOBALS_H_

#define KENNZAHL_UNKNOWN 0 ///< standard kennz flag for particles
#define KENNZAHL_NOT_FINISH 1 ///< kennz flag for particles which reached ::StorageTime
#define KENNZAHL_HIT_BOUNDARIES 2 ///< kennz flag for particles which left bounding box of TParticle::geom
#define KENNZAHL_NRERROR 3 ///< kennz flag for particles which produces a serious numerical error (step size underflow, missed reflection, ...)
#define KENNZAHL_DECAYED 4 ///< kennz flag for particles which reached TParticle::xend
#define KENNZAHL_INITIAL_NOT_FOUND 5 ///< kennz flag for particles which had a too low total energy to find a initial spot in the source volume

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

#include <string>

using namespace std;

void percent(long double x, long double x1, long double len, int &perc); ///< print percent of (x-x1)/len
void RotateVector(long double v[3], long double n[3]);///< rotate vector into new coordinate system whose z-axis lies on NORMALIZED vector n (active transformation)
void BOOST(long double beta[3], long double p[4]); ///< Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta
long double ProtonSpectrum(long double E); ///< proton recoil spectrum from neutron beta decay
long double ElectronSpectrum(long double E); ///< electron recoil spectrum from neutron beta decay

void CylToCart(long double v_r, long double v_phi, long double phi, long double &v_x, long double &v_y); ///< Translate cylindrical VECTOR (not point) components into cartesian components


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
const long double hquer = 1.05457266e-34; ///< planck constant [Js]
const long double mu_nSI = -0.96623641e-26;	///< Neutron magnetic moment [J/T]
const long double gamma_n = -1.83247185e8; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]

const long double lengthconv = 0.01; ///< length conversion factor cgs -> SI [cm -> m]
const long double Bconv = 1e-4; ///< magnetic field conversion factor cgs -> SI [G -> T]
const long double Econv = 1e2; ///< electric field conversion factor cgs -> SI [V/cm -> V/m]

extern long double StorageTime; ///< max. simulation time

extern int jobnumber; ///< job number, read from command line paramters, used for parallel calculations
extern string inpath; ///< path to configuration files, read from command line paramters
extern string outpath; ///< path where the log file should be saved to, read from command line parameters
extern int ausgabewunsch; ///< user choice for ausgabewunsch


#endif /*GLOBALS_H_*/
