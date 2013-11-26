/**
 * \file
 * Proton class definition.
 */

#ifndef PROTON_H_
#define PROTON_H_

static const char* NAME_PROTON = "proton";

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
	 * Proton constructor, set start values randomly according to *.in files.
	 *
	 * Sets start time TParticle::tstart according to [SOURCE] in geometry.in.
	 * Sets start energy according to TMCGenerator::Spectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 *
	 * @param number Particle number
	 * @param ageometry TGeometry in which the particle will be simulated
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TFieldManager used to calculate energies (can be NULL)
	 */
	TProton(int number, TGeometry &ageometry, TSource &src,
				TMCGenerator &mcgen, TFieldManager *afield): TParticle(NAME_PROTON, ele_e, m_p, 0){
		Init(number, ageometry, src, mcgen, afield);
	};

	/**
	 * Generic proton constructor
	 *
	 * @param number Particle number
	 * @param t Start time
	 * @param atau Life time
	 * @param x Start x coordinate
	 * @param y Start y coordinate
	 * @param z Start z coordinate
	 * @param vx Velocity x component
	 * @param vy Velocity y component
	 * @param vz Velocity z component
	 * @param pol Start polarisation
	 * @param trajl Max. simulated trajectory length
	 * @param ageometry TGeometry in which the particle will be simulated
	 * @param afield Electric and magnetic fields (needed to determine total energy)
	 */
	TProton(int number, long double t, long double atau, long double x, long double y, long double z,
			long double vx, long double vy, long double vz, int pol, long double trajl,
			TGeometry &ageometry, TFieldManager *afield): TParticle(NAME_PROTON, ele_e, m_p, 0){
		InitV(number, t, atau, x, y, z, vx, vy, vz, pol, trajl, ageometry, afield);
	}

protected:
	/**
	 * This method is executed, when a particle crosses a material boundary.
	 *
	 * Nothing happens to protons.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the proton is leaving
	 * @param entering Solid that the proton is entering (can be modified by method)
	 * @return Returns true if particle was reflected
	 */
	bool OnHit(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], const solid *leaving, const solid *&entering){
		return false;
	};


	/**
	 * This method is executed on each step.
	 *
	 * Protons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid in which the proton is at the moment
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, const solid *currentsolid){
		if (currentsolid->ID != geom->defaultsolid.ID){
			x2 = x1;
			y2 = y1;
			ID = currentsolid->ID;
			return true;
		}
		return false;
	};


	/**
	 * Proton decay (not used)
	 */
	void Decay(){

	};

};


#endif // PROTON_H_
