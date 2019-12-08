/**
 * \file
 * Proton class definition.
 */

#ifndef PROTON_H_
#define PROTON_H_

extern const char* NAME_PROTON; ///< name of TProton class

#include "globals.h"
#include "particle.h"

/**
 * Proton particle class.
 *
 * Simulates a proton including gravitation and Lorentz-force
 */
struct TProton: TParticle{
public:
	/**
	 * Create proton
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
	 * @param polarisation polarisation of particle (+/-1)
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TProton(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
			TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);

protected:
	/**
	 * This method is executed, when a particle crosses a material boundary.
	 *
	 * Nothing happens to protons.
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(TStep &stepper, const double normal[3],
			const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * This method is executed on each step.
	 *
	 * Protons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(TStep &stepper,
			const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Proton decay (not used)
	 *
	 * For parameter doc see TParticle::Decay
	 */
	void Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;
};

#endif // PROTON_H_
