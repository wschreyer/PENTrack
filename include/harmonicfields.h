/**
 * \file
 * Header file for the implementation of a magnetic field in terms of harmonic polynomials.
 * Polynomial coefficients are provided as inputs by the user.
 * For more information about the harmonic polynomial decomposition see https://doi.org/10.1103/PhysRevA.99.042112
 */

#ifndef HARMONICFIELD_H_
#define HARMONICFIELD_H_

#include "field.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/roots.hpp>

/**
 * Magetic field described by harmonic polynomials and their coefficients.
 *
 * Typically used in nEDM experiments
 * See https://doi.org/10.1103/PhysRevA.99.042112
 */
class HarmonicExpandedBField: public TField{
private:

	double xoff; ///< the x-coordinate offset from the origin
	double yoff; ///< the y-coordinate offset from the origin
	double zoff; ///< the z-coordinate offset from the origin

	// the rotational axis does not have to be normalized for input, there is 
	// a normalization to make it a unit vector within the program
	double axis_x; ///< x-component of rotational axis (NOT IMPLEMENTED YET)
	double axis_y; ///< y-component of rotational axis (NOT IMPLEMENTED YET)
	double axis_z; ///< z-component of rotational axis (NOT IMPLEMENTED YET)
	double angle;  ///< angle through which to rotate (NOT IMPLEMENTED YET)
	std::vector<std::tuple<int, int, double> > Glm; ///< Coefficients (l, m, Glm) of the harmonic expansion

public:
	/**
	 * Construct field class.
	 *
	 * @param _xoff the x-coordinate offset from the origin
	 * @param _yoff the y-coordinate offset from the origin
	 * @param _zoff the z-coordinate offset from the origin
	 * @param _axis_x the x-component of rotational axis (NOT IMPLEMENTED YET)
	 * @param _axis_y the y-component of rotational axis (NOT IMPLEMENTED YET)
	 * @param _axis_z the z-component of rotational axis (NOT IMPLEMENTED YET)
	 * @param _angle the angle through which to rotate (NOT IMPLEMENTED YET)
	 * @param Glm parameters; coefficients of the harmonic expansion
	 */

	HarmonicExpandedBField(const double _xoff, const double _yoff, const double _zoff,
			const double _axis_x, const double _axis_y, const double _axis_z, const double _angle,
			const std::vector<std::tuple<int, int, double> > _Glm);

	/**
	 * Calculate magnetic field vector and spatial derivatives from harmonic polynomials and their coefficients up to fifth order
	 * See https://doi.org/10.1103/PhysRevA.99.042112
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
	 * Old implementation with non-lazy evaluation up to third order
	*/
	void BField_old(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const;
	/**
	 * New implementation using autodifferentiation to calculate any order. Very slow at the moment.
	*/
	void BField_new(const double _x, const double _y, const double _z, const double t, double B[3], double dBidxj[3][3]) const;
	
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

#endif /*HARMONICFIELD_H_*/
