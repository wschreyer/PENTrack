/**
 * \file
 * The EDM fields used in the Ramsey cycle. 
 */

#ifndef STATICEDMFIELD_H_
#define STATICEDMFIELD_H_

#include "field.h"

/**
 * Magnetic field definition for the static B_0 field experienced by neutrons throughout the Ramsey cycle.
 */

struct TEDMStaticB0GradZField: public TField{
	double edmB0z0; ///< z component of B0 field at the origin in T
	double edmdB0z0dz; ///< z derivative of z component of the B0 field: dB0zdz at the origin in T

	/**
	 * Constructor, requires B0z component at centre of the EDM cell and the derivative of the B0z component along the z-axis.
	 */
	TEDMStaticB0GradZField(double abz, double adB0zdz);

	void BField(double x, double y, double z, double t, double B[4][4]);
	
	/**
	 * Adds no electric field. 
	 * Required because we are inheriting from abstract class TField which has virtual void EField function.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Electric potential
	 * @param Ei Electric field components
	 */
	void EField(double x, double y, double z, double t, double &V, double Ei[3]){};
};



#endif /*STATICEDMFIELD_H_*/
