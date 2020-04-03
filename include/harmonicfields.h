/**
 * \file
 * Header file for the implementation of a magnetic field determined by 
 * coefficients, provided as inputs by the user, of an expansion in
 * terms of harmonic polynomials. This implementation is based on the
 * work described in "Computation_summary.pdf" (!!! unfinished pls
 * update!)
 */

#ifndef HARMONICFIELD_H_
#define HARMONICFIELD_H_

#include "field.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/roots.hpp>

/**
 * Nearly homogeneous magnetic field with small gradient.
 *
 * Typically used in nEDM experiments
 */
class HarmonicExpandedBField: public TField{
private:

	double xoff; ///< the x-coordinate offset from the origin
	double yoff; ///< the y-coordinate offset from the origin
	double zoff; ///< the z-coordinate offset from the origin
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

	// the rotational axis does not have to be normalized for input, there is 
	// a normalization to make it a unit vector within the program
	double axis_x; // x-component of rotational axis
	double axis_y; // y-component of rotational axis
	double axis_z; // z-component of rotational axis
	double angle;  // angle through which to rotate
	double G[24]; ///< G parameters; coefficients of the harmonic expansion

public:
	/**
	 * Magnetic field definition for the static B_0 field experienced by neutrons throughout the Ramsey cycle.
	 *
	 * Constructor requires:
	 *
	 * @param _xoff the x-coordinate offset from the origin
	 * @param _yoff the y-coordinate offset from the origin
	 * @param _zoff the z-coordinate offset from the origin
	 * @param bW Distance from the edges where folding begins
	 * @param _xmax maximum x value for the field to permiate
	 * @param _xmin minimum x value for the field to permiate
	 * @param _ymax maximum y value for the field to permiate
	 * @param _ymin minimum y value for the field to permiate
	 * @param _zmax maximum z value for the field to permiate
	 * @param _zmin minimum z value for the field to permiate
	 * @param Bscale Formula to scale magnetic field
	 * @param _axis_x the x-component of rotational axis
	 * @param _axis_y the y-component of rotational axis
	 * @param _axis_z the z-component of rotational axis
	 * @param _angle the angle through which to rotate
	 * @param G parameters; coefficients of the harmonic expansion
	 */

	HarmonicExpandedBField(const double _xoff, const double _yoff, const double _zoff, const double bW,
			const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale,
			const double _axis_x, const double _axis_y, const double _axis_z, const double _angle, 
			const double G0, const double G1, const double G2, const double G3, const double G4, const double G5, const double G6, const double G7, const double G8, 
			const double G9, const double G10, const double G11, const double G12, const double G13, const double G14, const double G15, const double G16, const double G17, 
			const double G18, const double G19, const double G20, const double G21, const double G22, const double G23);

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

#endif /*HARMONICFIELD_H_*/