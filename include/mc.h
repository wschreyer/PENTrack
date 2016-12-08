/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <cstdlib>
#include <string>
#include <map>
#include <random>

/**
 * Class to generate random numbers in several distributions.
 */
class TMCGenerator{
private:
	mutable std::mt19937_64 rangen;
	uint64_t seed; ///< initial random seed
	mutable std::piecewise_linear_distribution<double> ProtonBetaSpectrumDist; ///< predefined energy distribution of protons from free-neutron decay
	mutable std::piecewise_linear_distribution<double> ElectronBetaSpectrumDist; ///< predefined energy distribution of electron from free-neutron decay

public:
	/**
	 * Constructor.
	 *
	 * Create random seed and read infile.
	 *
	 * @param infile Path to configuration file
	 * @param aseed Set a custom seed (optional), if set to zero seed is determined from high-resolution clock
	 */
	TMCGenerator(const uint64_t aseed = 0);

	/// return uniformly distributed random number in [min..max]
	double UniformDist(const double min, const double max) const;
	

	/// return sine distributed random number in [min..max] (in rad!)
	double SinDist(const double min, const double max) const;
	

	/// return sin(x)*cos(x) distributed random number in [min..max] (0 <= min/max < pi/2!)
	double SinCosDist(const double min, const double max) const;
	

	/// return x^2 distributed random number in [min..max]
	double SquareDist(const double min, const double max) const;
	

	/// return linearly distributed random number in [min..max]
	double LinearDist(const double min, const double max) const;


	/// return sqrt(x) distributed random number in [min..max]
	double SqrtDist(const double min, const double max) const;
	

	/// return normally distributed random number
	double NormalDist(const double mean, const double sigma) const;


	/// return exp(-exponent*x) distributed random number
	double ExpDist(const double exponent) const;


	/**
	 * Create isotropically distributed 3D angles.
	 *
	 * @param phi Azimuth
	 * @param theta Polar angle
	 */
	void IsotropicDist(double &phi, double &theta) const;
	

	/**
	 * Initial polarisation of different particles
	 *
	 * @param polarisation Projection of spin vector onto magnetic field vector
	 *
	 * @return Returns random polarisation (-1,1), weighted by spin projection onto magnetic field vector
	 */
	double DicePolarisation(const double polarisation) const;


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
	 * (1) Dice the energy of the decay proton according to globals.h#ProtonBetaSpectrum in the rest frame of the neutron and
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
	 * @param phi_e Returns azimuth of electron velocity vector
	 * @param theta_p Returns polar angle of proton velocity vector
	 * @param theta_e Returns polar angle of electron velocity vector
	 * @param pol_p Returns polarisation of proton spin
	 * @param pol_e Returns polarisation of electron spin
	 */
	void NeutronDecay(const double v_n[3], double &E_p, double &E_e, double &phi_p, double &phi_e, double &theta_p, double &theta_e, double &pol_p, double &pol_e) const;


	/**
	 * Outputs velocity spectrum for TOF expt given a chosen energy spectrum
	 * (guide is along y -axis)
	 *
	 * @param Ekin  - kinetic energy (does not get changed)
	 * @param phi - azimuthal associated with velocity
	 * @param theta - polar associated with velocity
	 */
	void tofDist(double &Ekin, double &phi, double &theta) const;

	/**
	 * Create a piecewise linear distribution from a muParser function
	 *
	 * @param func Formula string to parse
	 * @param range_min Lower range limit of distribution
	 * @param range_max Upper range limit of distribution
	 *
	 * @return Return piecewise linear distribution
	 */
	std::piecewise_linear_distribution<double> ParseDist(const std::string &func, const double range_min, const double range_max) const;

	/**
	 * Create a piecewise linear distribution from a function with single parameter
	 *
	 * @param func Function with single parameter
	 * @param range_min Lower range limit of distribution
	 * @param range_max Upper range limit of distribution
	 *
	 * @return Return piecewise linear distribution
	 */
	std::piecewise_linear_distribution<double> ParseDist(double (*func)(double), const double range_min, const double range_max) const;

	double SampleDist(std::piecewise_linear_distribution<double> &dist) const{
		return dist(rangen);
	}
};

#endif /*MC_H_*/
