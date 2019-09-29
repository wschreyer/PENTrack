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
	double BoundaryWidth; ///< Distance from the edges where folding begins
	double xmax; ///< maximum x value for the field to permiate
	double xmin; ///< minimum x value for the field to permiate
	double ymax; ///< maximum y value for the field to permiate
	double ymin; ///< minimum y value for the field to permiate
	double zmax; ///< maximum z value for the field to permiate
	double zmin; ///< minimum z value for the field to permiate
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
	 * @param bW Distance from the edges where folding begins
	 * @param _xmax maximum x value for the field to permiate
	 * @param _xmin minimum x value for the field to permiate
	 * @param _ymax maximum y value for the field to permiate
	 * @param _ymin minimum y value for the field to permiate
	 * @param _zmax maximum z value for the field to permiate
	 * @param _zmin minimum z value for the field to permiate
	 * @param Bscale Formula to scale magnetic field
	 */
	TEDMStaticB0GradZField(const double xoff, const double yoff, const double zoff, const double ang1, const double ang2, const double abz, const double adB0zdz, const double bW,
			const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale);
	
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
	
private:
	/**
	 * Folds the BField component (Bx, By, or Bz) and its spatial derivative components 
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param Bxi The Bfield component (Bx,By, or Bz) to be folded
	 * @param dBScaled The three derivatives of each Bfield component; Only the ones pertaining to Bxi are updated in a call
	 * @param i Integer to specify which Bxi component you passed; x = 0, y = 1, z = 2
	 **/
	void FieldSmthr(const double x, const double y, const double z, double Bxi[3], double dBScaled[9], const int i) const;
	
	/**
	 * Computes how much scale the field and its derivatives in the x y and z direction based off of how close they are to their respective boundaries
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param compFactors First three elements to compress in the x,y,z directions respectively. The next 3 are to indicate upper boundary (a value of 1) or a lower boundary (a value of 0)
	 **/
	void CompressionFactor(const double x, const double y, const double z, double *compFactors) const;
	
	/**
	 * Smooth function used to scale the field at the edges
	 *
	 * @param x function parameter (x = 0..1)
	 *
	 * @return Returns number between 0 and 1, smoothly rising with x
	 **/
	double SmthrStp(const double x) const;
	
	/**
	 * Derivative of SmthStp
	 *
	 * @param x function parameter (x = 0..1)
	 *
	 * @return Returns derivative of SmthrStp at parameter x
	 **/

	double SmthrStpDer(const double x) const;
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
 	* @param Escale Formula used to scale electric field
 	*/
	TEDMStaticEField(const double aexMag, const double aeyMag, const double aezMag, const std::string &Escale);

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
	void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = nullptr) const{};
	
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
			double &V, double Ei[3]) const;

};


#endif /*STATICEDMFIELD_H_*/
