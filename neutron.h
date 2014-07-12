/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

#include "particle.h"

extern const char* NAME_NEUTRON;

/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Create neutron
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
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 */
	TNeutron(int number, double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Each reflection has a certain chance of being diffuse.
	 * Diffuse reflection angles are cosine-distributed around the normal vector of the surface.
	 * Logs surface hits if reflectlog is set in config.in.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering (can be modified by method)
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);


	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * Caclulates adiabacity condition (TNeutron::frac).
	 * Estimates spin flip probability according to Vladimirskii (TNeutron::vlad)
	 * and by doing a bruteforce integration of the Bloch equation describing spin precession.
	 * Additionally it can print out a spatial neutron distribution matrix.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid);


	/**
	 * Neutron decay.
	 *
	 * Create decay proton and electron and add them to secondaries list
	 */
	void Decay();


	/**
	 * Potential energy, other than electromagnetic or gravitational, which is added to total energy of neutron
	 *
	 * @param t Time
	 * @param y Position vector
	 * @param polarisation Particle polarisation
	 * @param field TFieldManager for electromagnetic potential
	 * @param solids List of solids for optional material potential
	 *
	 * @return Returns potential energy plus Fermi-Potential of solid
	 */
	value_type Epot(value_type t, state_type y, int polarisation, TFieldManager *field, map<solid, bool> &solids);

};


#endif // NEUTRON_H_
