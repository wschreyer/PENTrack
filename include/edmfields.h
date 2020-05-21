/**
 * \file
 * The EDM fields used in the Ramsey cycle. 
 */

#ifndef STATICEDMFIELD_H_
#define STATICEDMFIELD_H_

#include "field.h"

/**
 * Nearly homogeneous magnetic field with small gradient.
 *
 * Typically used in nEDM experiments
 */
class TEDMStaticB0GradZField: public TField{
private:
	double edmB0xoff; ///< the x-coordinate offset from the origin
	double edmB0yoff; ///< the y-coordinate offset from the origin
	double edmB0zoff; ///< the z-coordinate offset from the origin
	double pol_ang1; ///< the polar angle theta
	double azm_ang2; ///< the azimuthal angle phi
	double edmB0z0; ///< the z component of the B0 field at the origin in Tesla
	double edmdB0z0dz; ///< the z derivative of z component of the B0 field: dB0zdz at the origin in Tesla
	double Rot1[3][3]; ///< rotation matrix
	double Rot2[3][3]; ///< rotation matrix
	double Rot3[3][3]; ///< rotation matrix
	double dB[9]; ///< Theory values for derivatives in the BField frame
	double Bd[6][3]; ///< change in B field in this coordinate system to be transformed into original coordinate system
	
public:
	/**
	 * Magnetic field definition for the static B_0 field experienced by neutrons throughout the Ramsey cycle.
	 *
	 * Constructor requires:
	 *
	 * @param xoff the x-coordinate offset from the origin
	 * @param yoff the y-coordinate offset from the origin
	 * @param zoff the z-coordinate offset from the origin
	 * @param ang1 the polar angle theta
	 * @param ang2 the azimuthal angle phi
	 * @param abz the z component of the B0 field at the origin in Tesla
	 * @param adB0zdz the z derivative of z component of the B0 field: dB0zdz at the origin in Tesla
	 */
	TEDMStaticB0GradZField(const double xoff, const double yoff, const double zoff, const double ang1, const double ang2, const double abz, const double adB0zdz);
	
	/**
	 * The calculation for the B0 components are derived from a first order approximation of
	 * an axially symmetric static magnetic field. See for example: http://physics.princeton.edu/~mcdonald/examples/axial.pdf
	 *
	 * @param x the x-coordinate in the field's coordinate system
	 * @param y the x-coordinate in the field's coordinate system
	 * @param z the x-coordinate in the field's coordinate system
	 * @param t the time (not used since the field is presumed to be static)
	 * @param B Returns magnetic-field components
	 * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
	**/
	void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const override;
	
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
	 **/
	void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override {};
};

/**
 * The static homogenous field that is parallel to the B0 field.
 * 
 */
class TEDMStaticEField: public TField {
private:
	double exMag; ///< X component of electric field
	double eyMag; ///< Y component of electric field
	double ezMag; ///< Z component of electric field

public:
	/**
 	* Constructor requires just the magnitude of the field.
 	*
 	* @param aexMag X component of electric field
 	* @param aeyMag Y component of electric field
 	* @param aezMag Z component of electric field
 	*/
	TEDMStaticEField(const double aexMag, const double aeyMag, const double aezMag);

	/**
	 * Adds no magnetic field.
	 *
	 * @param x the x-coordinate in the field's coordinate system
	 * @param y the x-coordinate in the field's coordinate system
	 * @param z the x-coordinate in the field's coordinate system
	 * @param t the time (not used since the field is presumed to be static)
	 * @param B Returns magnetic-field components
	 * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
	**/
	void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = nullptr) const override{};
	
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
	void EField(const double x, const double y, const double z, const double t,
			double &V, double Ei[3]) const override;

};


#endif /*STATICEDMFIELD_H_*/
