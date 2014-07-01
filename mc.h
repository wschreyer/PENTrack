/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <cstdlib>
#include <string>
#include <map>

#include <boost/random.hpp>

#include "muParser.h"

/**
 * For each section in particle.in such a struct is created containing all user options
 */
struct TParticleConfig{
	double tau; ///< lifetime
	double tmax; ///< max. simulation time
	double lmax; ///< max. trajectory length
	int polarization; ///< initial polarization
	double Emin; ///< min. initial energy
	double Emax; ///< max. initial energy
	mu::Parser spectrum; ///< Parsed energy spectrum given by user
	double phi_v_min; ///< Parsed minimum for initial azimuthal angle of velocity given by user
	double phi_v_max; ///< Parsed maximum for initial azimuthal angle of velocity given by user
	mu::Parser phi_v; ///< Parsed initial azimuthal angle distribution of velocity given by user
	double theta_v_min; ///< Parsed minimum for initial polarl angle of velocity given by user
	double theta_v_max; ///< Parsed maximum for initial polar angle of velocity given by user
	mu::Parser theta_v; ///< Parsed initial polar angle distribution of velocity given by user
};

/**
 * Class to generate random numbers in several distributions.
 */
class TMCGenerator{
private:
	boost::mt19937_64 rangen; ///< random number generator
	double xvar; ///< x variable for formula parser
	std::map<std::string, TParticleConfig> pconfigs;

public:
	uint64_t seed; ///< initial random seed

	/**
	 * Constructor.
	 *
	 * Create random seed and read infile.
	 *
	 * @param infile Path to configuration file
	 */
	TMCGenerator(const char *infile);

	/**
	 * Destructor.
	 *
	 * Delete random number generator.
	 */
	~TMCGenerator();

	/// return uniformly distributed random number in [min..max]
	double UniformDist(double min, double max);
	

	/// return sine distributed random number in [min..max] (in rad!)
	double SinDist(double min, double max);
	

	/// return sin(x)*cos(x) distributed random number in [min..max] (0 <= min/max < pi/2!)
	double SinCosDist(double min, double max);
	

	/// return x^2 distributed random number in [min..max]
	double SquareDist(double min, double max);
	

	/// return linearly distributed random number in [min..max]
	double LinearDist(double min, double max);


	/// return sqrt(x) distributed random number in [min..max]
	double SqrtDist(double min, double max);
	

	/**
	 * Create isotropically distributed 3D angles.
	 *
	 * @param phi Azimuth
	 * @param theta Polar angle
	 */
	void IsotropicDist(long double &phi, long double &theta);
	

	/// energy distribution of UCNs
	double NeutronSpectrum();

	/**
	 * Energy distribution for each particle type
	 */
	double Spectrum(const std::string &particlename);


	/**
	 * Angular velocity distribution for each particle type
	 * 
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 * @param phi_v Returns velocity azimuth
	 * @param theta_v Returns velocity polar angle
	 */
	void AngularDist(const std::string &particlename, long double &phi_v, long double &theta_v);

	/**
	 * Lifetime of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns lifetimes using an exponentially decaying or flat distribution, depending on user choice in particle.in
	 */
	double LifeTime(const std::string &particlename);
	

	/**
	 * Max. trajectory length of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns max. trajectory length, depending on user choice in particle.in
	 */
	double MaxTrajLength(const std::string &particlename);
	

	/**
	 * Initial polarisation of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns deiced of fixed polarisation (-1, 0, 1), depending on user choice in particle.in
	 */
	int DicePolarisation(const std::string &particlename);


	/**
	 * Simulate neutron beta decay.
	 *
	 * Calculates velocities of neutron decay products proton and electron.
	 *
	 * Reaction(s):
	 *
	 *     n0  ->  p+  +  R-
	 *
	 *     R-  ->  e-  +  nue
	 * 
	 * Procedure:
	 * (1) Dice the energy of the decay proton according to #ProtonBetaSpectrum in the rest frame of the neutron and
	 *     calculate the proton momentum via the energy momentum relation.
	 * (2) Dice isotropic orientation of the decay proton.
	 * (3) Calculate 4-momentum of rest R- via 4-momentum conservation.
	 * (4) Get fixed electron energy from two body decay of R-.
	 * (5) Dice isotropic electron orientation in the rest frame of R-.
	 * (6) Lorentz boost electron 4-momentum into moving frame of R-.
	 * (7) Calculate neutrino 4-momentum via 4-momentum conservation.
	 * (8) Boost all 4-momentums into moving neutron frame.
	 * 
	 * Cross-check:
	 * (9) Print neutrino 4-momentum invariant mass (4-momentum square, should be zero).
	 * 
	 * @param v_n Velocity of decayed neutron
	 * @param E_p Returns proton kinetic energy
	 * @param E_e Returns electron kinetic energy
	 * @param phi_p Returns azimuth of proton velocity vector
	 * @param phi_p Returns azimuth of electron velocity vector
	 * @param theta_p Returns polar angle of proton velocity vector
	 * @param theta_e Returns polar angle of electron velocity vector
	 */
	void NeutronDecay(long double v_n[3], double &E_p, double &E_e, double &phi_p, double &phi_e, double &theta_p, double &theta_e);
};

#endif /*MC_H_*/
