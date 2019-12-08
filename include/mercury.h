/**
 * \file
 * Mercury class definition. Mercury-199 is used as a comagnetometer in the EDM experiment at TRIUMF. 
 */

#ifndef MERCURY_H_
#define MERCURY_H_

#include "particle.h"

extern const char* NAME_MERCURY; ///< name of TMERCURY class

/**
 * Mercury-199 particle class.
 *
 * Simulates a mercury atom including gravitation and Lorentz-force
 */
struct TMercury: TParticle{
public:
	/**
	 * Create mercury-199 atom.
	 *
	 * Wraps basic TParticle constructor
	 *
	 * @param number Particle number
	 * @param t Starting time
	 * @param x Initial x coordinate
	 * @param y Initial y coordinate
	 * @param z Initial z coordinate
	 * @param E Initial kinetic energy
	 * @param phi Initial azimuth of velocity vector
	 * @param theta Initial polar angle of velocity vector
	 * @param polarisation polarisation (+/-1)
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TMercury(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
			TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);

protected:
	/**
	 * This method is executed, when a particle encounters a material boundary.
	 * The atom is reflected specularly or diffusely depending on the material type. 
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(TStep &stepper, const double normal[3],
			const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;

	/**
	 * This method is executed on each step.
	 *
	 * Do nothing. 
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(TStep &stepper,
			const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Mercury decay (not used)
	 *
	 * For parameter doc see TParticle::Decay
	 */
	void Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;
};

#endif // MERCURY_H_
