#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <map>

#define ID_UNKNOWN 0 ///< standard kennz flag for particles
#define ID_NOT_FINISH -1 ///< kennz flag for particles which reached ::StorageTime
#define ID_HIT_BOUNDARIES -2 ///< kennz flag for particles which left bounding box of TParticle::geom
#define ID_NRERROR -3 ///< kennz flag for particles which produces a serious numerical error (step size underflow, missed reflection, ...)
#define ID_DECAYED -4 ///< kennz flag for particles which reached TParticle::tau
#define ID_INITIAL_NOT_FOUND -5 ///< kennz flag for particles which had a too low total energy to find a initial spot in the source volume

#define PARTICLE 1 ///< set particletype in configuration to this value to simulate particles
#define BF_ONLY 3 ///< set particletype in configuration to this value to print out a ramp heating analysis
#define BF_CUT 4 ///< set particletype in configuration to this value to print out a planar slice through electric/magnetic fields
#define GEOMETRY 7 ///< set particletype in configuration to this value to print out a sampling of the geometry

// physical constants
static const long double pi = 3.1415926535897932384626L; ///< Pi
static const long double ele_e = 1.602176487E-19L; ///< elementary charge [C]
static const long double gravconst = 9.80665L; ///< g [m/s]
static const long double conv = pi/180.L; ///< deg to rad conversion factor
static const long double mu0 = 4*pi*1e-7L; ///< magnetic permeability [Vs/Am]
static const long double m_n = 1.674927211E-27L/ele_e; ///< neutron mass [eV/c^2]
static const long double m_p = 1.672621637E-27L/ele_e; ///< proton mass [eV/c^2]
static const long double m_e = 9.10938215e-31L/ele_e; ///< electron mass [eV/c^2]
static const long double c_0 = 299792458.L; ///< light speed [m/s]
static const long double hbar = 1.05457266e-34L; ///< planck constant [Js]
static const long double mu_nSI = -0.96623641e-26L;	///< Neutron magnetic moment [J/T]
static const long double gamma_n = -1.83247185e8L; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]

static const long double lengthconv = 0.01; ///< length conversion factor cgs -> SI [cm -> m]
static const long double Bconv = 1e-4; ///< magnetic field conversion factor cgs -> SI [G -> T]
static const long double Econv = 1e2; ///< electric field conversion factor cgs -> SI [V/cm -> V/m]

extern long long int jobnumber; ///< job number, read from command line paramters, used for parallel calculations
extern std::string inpath; ///< path to configuration files, read from command line paramters
extern std::string outpath; ///< path where the log file should be saved to, read from command line parameters

/**
 * Print progress bar.
 *
 * Prints a point every 2% and a number every 10%
 *
 * @param percentage Progress of current action (0..1)
 * @param lastprint Put in and return last printed percentage
 */
void PrintPercent(double percentage, int &lastprint);

/**
 * Rotate a vector.
 *
 * Rotate vector into new coordinate system whose z-axis lies on NORMALIZED vector n (active transformation)
 */
void RotateVector(long double v[3], const long double n[3]);

/**
 * Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta.
 *
 * Copy & paste from ROOT
 *
 * @param beta Vector v/c of moving reference frame
 * @param p Four-vector
 */
void BOOST(long double beta[3], long double p[4]);

/**
 * Energy distribution of protons from free neutron beta decay (0 < E < 750 eV)
 *
 * From diploma thesis M. Simson.
 *
 * @param E energy
 *
 * @return Returns a probability between 0 and 1 that a proton with energy E is created.
 */
double ProtonBetaSpectrum(double E);


/**
 * Energy distribution of electrons from free neutron decay (0 < E < 782 keV)
 *
 * From "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
 *
 * @param E energy
 *
 * @return Returns probability between 0 and 1 that an electron with energy E is created.
 */
double ElectronBetaSpectrum(double E);


typedef std::map<std::string, std::map<std::string, std::string> > TConfig; ///< map of sections containing a map of key-value pairs

/**
 * Read variables from *.in file into map.
 *
 * in file is assumed to have structure
 * [section]
 * key arbitrary string
 * key value
 * [section]
 * ...
 *
 * @param inpath Path to in file.
 * @param vars Return TConfig map
 */
void ReadInFile(const char *inpath, TConfig &vars);

#endif /*GLOBALS_H_*/
