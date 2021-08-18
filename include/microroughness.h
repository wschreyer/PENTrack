/**
 * \file
 * Routines to calculate the microroughness scattering distribution, its integral, and its maximum.
 * See A. Steyerl, Effect of surface roughness on the total reflexion and transmission of slow neutrons, Z. Phys. A 254, 169–188 (1972)
 * https://doi.org/10.1007/BF01380066
 */

#ifndef INCLUDE_MICROROUGHNESS_H_
#define INCLUDE_MICROROUGHNESS_H_

namespace MR{
	/**
	 * Check if the MicroRoughness is model is applicable to the current interaction, see equations (17) and (18) in Steyerl's publication
	 *
	 * @param v veclocity right before hitting the material
	 * @param normal normal vector of the material boundary
	 * @param Estep potential step at the material boundary
	 * @param RMSroughness root-mean-square roughness of the surface at the material boundary
	 * @param correlationLength correlation length of the surface at the boundary
	 *
	 * @return Returns true, if the MicroRoughness model can be used.
	 */
	bool MRValid(const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength);

	/**
	 * Return MicroRoughness model probability distribution for scattering angles theta, phi. See equations (20) and (21) in Steyerl's publication
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param integral Compute phi-integral of diffuse scattering probability distribution
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of hit surface
	 * @param Estep potential step at the material boundary
	 * @param RMSroughness root-mean-square roughness of the surface at the material boundary
	 * @param correlationLength correlation length of the surface at the boundary
	 * @param theta Polar angle of scattered velocity vector (0 < theta < pi/2)
	 * @param phi Azimuthal angle of scattered velocity vector (0 < phi < 2*pi), ignored if integral == true
	 *
	 * @return Returns probability of reflection/transmission, in direction (theta, phi) or its integration over phi=0..2pi if integral == true
	 */
	double MRDist(const bool transmit, const bool integral, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength, const double theta, const double phi);

	/**
	 * Calculate total diffuse scattering probability according to MicroRoughness model by doing numerical theta-integration of MRDist with integral = true
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of hit surface
	 * @param Estep potential step at the material boundary
	 * @param RMSroughness root-mean-square roughness of the surface at the material boundary
	 * @param correlationLength correlation length of the surface at the boundary
	 *
	 * @return Returns probability of reflection/transmission
	 */
	double MRProb(const bool transmit, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength);

	/*
	 * Calculate maximum of MicroRoughness model distribution by doing numerical minimization of TNeutron::NegMRDist
	 *
	 * @param transmit True, if the particle is transmitted through the material boundary
	 * @param v velocity right before surface hit
	 * @param normal Normal vector of material boundary
	 * @param Estep potential step at the material boundary
	 * @param RMSroughness root-mean-square roughness of the surface at the material boundary
	 * @param correlationLength correlation length of the surface at the boundary
	 *
	 * @return Returns maximal value of MicroRoughness model distribution in range (theta = 0..pi/2, phi = 0..2pi)
	 */
	double MRDistMax(const bool transmit, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength);
};


#endif /* INCLUDE_MICROROUGHNESS_H_ */
