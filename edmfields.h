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
	
	/**
	 * The calculation for the B0 components are derived from a first order approximation of
	 * an axially symmetric static magnetic field. See for example: http://physics.princeton.edu/~mcdonald/examples/axial.pdf
	 *
	 * @param x, the x-coordinate in the field's coordinate system
	 * @param y, the x-coordinate in the field's coordinate system
	 * @param z, the x-coordinate in the field's coordinate system
	 * @param t, the time (not used since the field is presumed to be static)
	 * @param B, Magnetic field component matrix to which the values are added
	**/
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
	void EField(double x, double y, double z, double t, double &V, double Ei[3]) {};
};

/**
 * The static homogenous field that is parallel to the B0 field.
 * 
 */
struct TEDMStaticEField: public TField {
	double magEField; 	

	/**
 	* Constructor requires just the magnitude of the field. 
 	*/
	TEDMStaticEField(double amagEField);

	/**
	 * Adds no magnetic field.
	 *
	 * @param x, the x-coordinate in the field's coordinate system
	 * @param y, the x-coordinate in the field's coordinate system
	 * @param z, the x-coordinate in the field's coordinate system
	 * @param t, the time (not used since the field is presumed to be static)
	 * @param B, Magnetic field component matrix to which the values are added
	**/
	void BField(double x, double y, double z, double t, double B[4][4]) {};
	
	/**
	 * The static E field should be parallel to the B0 field in Ramsey cycle.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Electric potential
	 * @param Ei Electric field components
	 */
	void EField(double x, double y, double z, double t, double &V, double Ei[3]); 


};


#endif /*STATICEDMFIELD_H_*/
