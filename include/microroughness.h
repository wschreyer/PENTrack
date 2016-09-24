/*
 * microroughness.h
 *
 *  Created on: Dec 10, 2015
 *      Author: wolfgang
 */

#ifndef INCLUDE_MICROROUGHNESS_H_
#define INCLUDE_MICROROUGHNESS_H_

#include "geometry.h"

namespace MR{
	/**
	 * Check if the MicroRoughness is model is applicable to the current interaction
	 *
	 * @param v veclocity right before hitting the material
	 * @param normal normal vector of the material boundary
	 * @param leaving Material that the particle is leaving through the boundary
	 * @param entering Material, that the particle entering through the boundary
	 *
	 * @return Returns true, if the MicroRoughness model can be used.
	 */
	bool MRValid(const double v[3], const double normal[3], const solid &leaving, const solid &entering);

	/**
	 * Return MicroRoughness model probability distribution for scattering angles theta, phi
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param integral Compute phi-integral of diffuse scattering probability distribution
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid, which the neutron would leave if it were transmitted
	 * @param entering Solid, which the neutron would enter if it were transmitted
	 * @param theta Polar angle of scattered velocity vector (0 < theta < pi/2)
	 * @param phi Azimuthal angle of scattered velocity vector (0 < phi < 2*pi), ignored if integral == true
	 *
	 * @return Returns probability of reflection/transmission, in direction (theta, phi) or its integration over phi=0..2pi if integral == true
	 */
	double MRDist(const bool transmit, const bool integral, const double v[3], const double normal[3], const solid &leaving, const solid &entering, const double theta, const double phi);

	/**
	 * Calculate total diffuse scattering probability according to MicroRoughness model by doing numerical theta-integration of MRDist with integral = true
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid, which the neutron would leave if it were transmitted
	 * @param entering Solid, which the neutron would enter if it were transmitted
	 *
	 * @return Returns probability of reflection/transmission
	 */
	double MRProb(const bool transmit, const double v[3], const double normal[3], const solid &leaving, const solid &entering);

	/*
	 * Calculate maximum of MicroRoughness model distribution by doing numerical minimization of TNeutron::NegMRDist
	 *
	 * @param transmit True, if the particle is transmitted through the material boundary
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of material boundary
	 * @param leaving Material, that the particle is leaving at this material boundary
	 * @param entering Material, that the particle is entering at this material boundary
	 *
	 * @return Returns maximal value of MicroRoughness model distribution in range (theta = 0..pi/2, phi = 0..2pi)
	 */
	double MRDistMax(const bool transmit, const double v[3], const double normal[3], const solid &leaving, const solid &entering);
};


#endif /* INCLUDE_MICROROUGHNESS_H_ */
