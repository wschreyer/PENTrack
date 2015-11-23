#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <map>

#define ID_UNKNOWN 0 ///< standard flag for particles
#define ID_NOT_FINISH -1 ///< flag for particles which reached ::StorageTime
#define ID_HIT_BOUNDARIES -2 ///< flag for particles which left bounding box of TParticle::geom
#define ID_ODEINT_ERROR -3 ///< flag for particles which produced a numerical error during ODe integration
#define ID_DECAYED -4 ///< flag for particles which reached TParticle::tau
#define ID_INITIAL_NOT_FOUND -5 ///< flag for particles which had a too low total energy to find a initial spot in the source volume
#define ID_CGAL_ERROR -6 ///< flag for particles which produced an error during geometry collision checks
#define ID_GEOMETRY_ERROR -7 ///< flag for particles which produced an error while tracking material boundaries along the trajectory
#define ID_ABSORBED_IN_MATERIAL 1 ///< flag for particles that were absorbed inside a material
#define ID_ABSORBED_ON_SURFACE 2 ///< flag for particles that were absorbed on a material surface

#define PARTICLE 1 ///< set simtype in configuration to this value to simulate particles
#define BF_ONLY 3 ///< set simtype in configuration to this value to print out a ramp heating analysis
#define BF_CUT 4 ///< set simtype in configuration to this value to print out a planar slice through electric/magnetic fields
#define GEOMETRY 7 ///< set simtype in configuration to this value to print out a sampling of the geometry
#define MR_THETA_OUT_ANGLE 8 ///< set simtype in configuration to this value to output a 2d histogram of the MR model's diffuse reflection 
			     ///< probability for every solid angle
#define MR_THETA_I_ENERGY 9 ///< set simtype in configuration to this value to output a 2d histogram of the MR models' diffuse reflection
			    ///< probability for theta_i vs neutron energy

// physical constants
static const long double pi = 3.1415926535897932384626L; ///< Pi
static const long double ele_e = 1.602176487E-19L; ///< elementary charge [C]
static const long double gravconst = 9.80665L; ///< g [m/s]
static const long double boltzconst = 1.38064852E-23L/ele_e; /// Boltzmann's constant [ev/K] from http://physics.nist.gov/cgi-bin/cuu/Value?sigma
static const long double avogadroconst = 6.022140857E23L; /// Avogadro's constant [ mol^-1 ] from http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avagadro%27s+number
static const long double conv = pi/180.L; ///< deg to rad conversion factor
static const long double mu0 = 4*pi*1e-7L; ///< magnetic permeability [Vs/Am]
static const long double m_n = 1.674927211E-27L/ele_e; ///< neutron mass [eV/c^2]
static const long double m_p = 1.672621637E-27L/ele_e; ///< proton mass [eV/c^2]
static const long double m_e = 9.10938215e-31L/ele_e; ///< electron mass [ eV/c^2 ]
static const long double m_hg = 198.96828064/(1000*avogadroconst)/ele_e; ///< mercury-199 mass [ev/c^2] (http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Hg)
static const long double m_xe = 128.9047808611/(1000*avogadroconst)/ele_e; ///< xenon-129 mass [ev/c^2] (http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Xe)
static const long double c_0 = 299792458.L; ///< light speed [m/s]
static const long double hbar = 1.05457266e-34L; ///< planck constant [Js]
static const long double mu_nSI = -0.96623641e-26L;	///< Neutron magnetic moment [J/T]
static const long double mu_hgSI = 2.555118e-27L; ///< mercury-199 magnetic moment [J/T]
static const long double mu_xeSI = -3.392939e-27L; ///< Xenon-129 avogadroconstmagnetic moment [J/T]
static const long double gamma_n = -1.83247185e8L; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]
static const long double gamma_hg = 4.76901003e7L; ///< from: http://www.sciencedirect.com/science/article/pii/S0370269314007692 [ 1/Ts ]
static const long double gamma_xe = -7.399707336e7L; ///< from: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio [ 1/Ts ]
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
 * Rotate vector into new coordinate with basis vectors x and z (active transformation)
 */
void RotateVector(double v[3], const double z[3], const double x[3] = NULL);

/**
 * Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta.
 *
 * Copy & paste from ROOT
 *
 * @param beta Vector v/c of moving reference frame
 * @param p Four-vector
 */
void BOOST(double beta[3], double p[4]);

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

/**
 * Energy distribution of comagnetometer gasses from Maxwell-Boltzmann distribution. 
 *
 * From "http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html"
 *
 * @param T temperature of the gas in Kelvin
 * @param E energy in eV
 *
 * @return Returns probability between 0 and 1 that a gas molecule with energy E being created. 
 */
double MaxwellBoltzSpectrum (double T, double E); 

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
