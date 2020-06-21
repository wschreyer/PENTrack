/*
 * field.h
 *
 *  Contains virtual base class for field calculation methods
 */

#ifndef FIELD_H_
#define FIELD_H_

#include <array>
#include <memory>

#include "exprtk.hpp"


/**
 * Virtual base class for all field calculation methods
 */
class TField{
public:
	/**
	 * Calculate magnetic field at a given position and time.
	 *
	 * Has to be implemented by all derived field calculation classes.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param B Returns magnetic field vector
	 * @param dBidxj Return spatial derivatives of magnetic field components (optional)
	 */
	virtual void BField (const double x, const double y, const double z, const double t,
            double B[3], double dBidxj[3][3]) const = 0;

	/**
	 * Calculate electric field and potential at a given position.
	 *
	 * Has to be implemented by all derived field calculation classes.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Return electric potential (!=0 only if a map with potential was loaded)
	 * @param Ei Returns electric field vector
	 */
	virtual void EField (const double x, const double y, const double z, const double t,
            double &V, double Ei[3]) const = 0;

};


/**
 * Class to calculate a time-dependent field-scaling factor based on a formula string
 */
class TFieldScaler{
private:
	exprtk::expression<double> scaler; ///< formula interpreter for field scaling
	std::unique_ptr<double> tvar; ///< time variable for use in scaling-formula parsers. It needs to be a pointer to make sure the reference in the exprtk expression will not be invalidated when copying
public:
	/**
	 * Calculate time-dependent scaling factor from parsed formula
	 *
	 * @param t Time
	 *
	 * @return Return scaling factor
	 */
	double scalingFactor(const double t) const;

	/**
	 * Scale scalar field F with gradient dFdxi by calculated scaling factor
	 * 
	 * @param t Time used to calculate scaling factor
	 * @param F Value of scalar field
	 * @param dFdxi Gradient of scalar field
	 */
	void scaleScalarField(const double t, double &F, double dFdxi[3]) const;

	/**
	 * Scale vector field F with Jacobian dFidxj by calculated scaling factor
	 * 
	 * @param t Time used to calculate scaling factor
	 * @param F Value of vector field
	 * @param dFidxj Jacobian matrix of vector field
	 */
	void scaleVectorField(const double t, double F[3], double dFidxj[3][3]) const;

	/**
	 * Default constructor. Sets scaling factor to zero.
	 */
	TFieldScaler(): TFieldScaler("0") {};

	/**
	 * Constructor, set time-dependent scaling formula. Only variable allowed in the formula is "t" for time.
	 *
	 * @param scalingFormula String containing formula describing time-dependence of field
	 */
	TFieldScaler(const std::string &scalingFormula);
};


/**
 * Virtual base class defining a spatial boundary for field calculations. Can smoothly transition fields to zero at the boundary.
 * 
 */
class TFieldBoundary{
public:
	/**
	 * Check if valid boundaries are set
	 * 
	 * @return Returns true if boundaries are valid, false otherwise
	 */
	virtual bool hasBounds() const = 0;

	/**
	 * Check if coordinates are within the boundaries
	 * 
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * 
	 * @return Returns true if coordinates are within the boundaries
	 */
	virtual bool inBounds(const double x, const double y, const double z) const = 0;

	/**
	 * Smoothly scale field at the edges of the boundary region
	 *
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * @param F Field value at (x,y,z)
	 * @param dFdxi Field derivatives at (x,y,z)
	 */
	virtual void scaleScalarFieldAtBounds(const double x, const double y, const double z, double &F, double dFdxi[3]) const = 0;

	/**
	 * Smoothly reduce a vector field at the edges of the boundary region
	 *
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * @param F Field vector at (x,y,z)
	 * @param dFidxj Field derivatives at (x,y,z)
	 */
	virtual void scaleVectorFieldAtBounds(const double x, const double y, const double z, double F[3], double dFidxj[3][3]) const = 0;
};


/**
 * Class defining a spatial boundary box for field calculations. Can smoothly transition fields to zero at the boundary.
 * 
 */
class TFieldBoundaryBox: public TFieldBoundary{
private:
	double xmin, xmax, ymin, ymax, zmin, zmax; ///< variables for boundary box
	double boundaryWidth; ///< width of boundary along the bounding box that the field should be scaled in

	/**
	 * Smooth function used to scale the field at the edges
	 *
	 * @param x function parameter (x = 0..1)
	 *
	 * @return Returns number between 0 and 1, smoothly rising with x
	 */
	double smthrStp(const double x) const;

	/**
	 * Derivative of SmthStpDer
	 *
	 * @param x function parameter (x = 0..1)
	 *
	 * @return Returns derivative of smthrStp at parameter x
	 */
	double smthrStpDer(const double x) const;
public:
	/**
	 * Check if valid boundaries are set
	 * 
	 * @return Returns true if boundaries are valid, false otherwise
	 */
	bool hasBounds() const override;

	/**
	 * Check if coordinates are within the boundaries
	 * 
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * 
	 * @return Returns true if coordinates are within the boundaries xmin-xmax, ymin-ymax, zmin-zmax
	 */
	bool inBounds(const double x, const double y, const double z) const override;

	/**
	 * Smoothly scale field at the edges of the boundary region
	 *
	 * If coordinates are within BoundaryWidth of the edges of the boundary box,
	 * the field and its derivatives are scaled by the SmthrStp and SmthrStpDer functions.
	 *
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * @param F Field value at (x,y,z)
	 * @param dFdxi Field derivatives at (x,y,z)
	 */
	void scaleScalarFieldAtBounds(const double x, const double y, const double z, double &F, double dFdxi[3]) const override;

	/**
	 * Smoothly reduce a vector field at the edges of the boundary region
	 *
	 * If coordinates are within BoundaryWidth of the edges of the boundary box,
	 * the field and its derivatives are scaled by the SmthrStp and SmthrStpDer functions.
	 *
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param z z coordinate
	 * @param F Field vector at (x,y,z)
	 * @param dFidxj Field derivatives at (x,y,z)
	 */
	void scaleVectorFieldAtBounds(const double x, const double y, const double z, double F[3], double dFidxj[3][3]) const override;

	/**
	 * Default constructor. No scaling is applied
	 */
	TFieldBoundaryBox(): xmin(0.), xmax(0.), ymin(0.), ymax(0.), zmin(0.), zmax(0.), boundaryWidth(0.){}

	/**
	 * Constructor setting boundary parameters
	 * 
	 * @param _xmax Maximum x coordinate of bounding box
	 * @param _xmin Minimum x coordinate of bounding box
	 * @param _ymax Maximum x coordinate of bounding box
	 * @param _ymin Minimum y coordinate of bounding box
	 * @param _zmax Maximum x coordinate of bounding box
	 * @param _zmin Minimum z coordinate of bounding box
	 * @param _boundaryWidth If coordinates fall within this distance from the boundary, the field will be scaled to smoothly transition to no field outside the boundary
	 */
	TFieldBoundaryBox(const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const double _boundaryWidth);
};


/**
 * Class combining TField with TFieldBoundary and TFieldScaler for a fully featured field calculation
 */
class TFieldContainer{
private:
	std::unique_ptr<TField> field; ///< Class derived from TField
	TFieldScaler BScaler; ///< Scaler class for magnetic field
	TFieldScaler EScaler; ///< Scaler class for electric field
	std::unique_ptr<TFieldBoundary> boundary; ///< Class derived from TFieldBoundary
public:
	/**
	 * Constructor for a field using a box-shaped boundary
	 * 
	 * @param _field Class derived from TField
	 * @param BScalingFormula String containing the formula to calculate a time-dependent magnetic-field scaling factor
	 * @param EScalingFormula String containing the formula to calculate a time-dependent electric-field scaling factor
	 * @param xmax Maximum x coordinate of bounding box
	 * @param xmin Minimum x coordinate of bounding box
	 * @param ymax Maximum x coordinate of bounding box
	 * @param ymin Minimum y coordinate of bounding box
	 * @param zmax Maximum x coordinate of bounding box
	 * @param zmin Minimum z coordinate of bounding box
	 * @param boundaryWidth If coordinates fall within this distance from the boundary, the field will be scaled to smoothly transition to no field outside the boundary
	 */
	TFieldContainer(std::unique_ptr<TField> &&_field, const std::string &BScalingFormula, const std::string &EScalingFormula,
					const double xmax, const double xmin, const double ymax, const double ymin, const double zmax, const double zmin, const double boundaryWidth):
						field(std::move(_field)), BScaler(TFieldScaler(BScalingFormula)), EScaler(TFieldScaler(EScalingFormula)), 
						boundary(std::unique_ptr<TFieldBoundary>(new TFieldBoundaryBox(xmax, xmin, ymax, ymin, zmax, zmin, boundaryWidth))) {}

	/**
	 * Constructor for field without boundary
	 * 
	 * @param _field Class derived from TField
	 * @param BScalingFormula String containing the formula to calculate a time-dependent magnetic-field scaling factor (optional)
	 * @param EScalingFormula String containing the formula to calculate a time-dependent electric-field scaling factor (optional)
	 */
	TFieldContainer(std::unique_ptr<TField> &&_field, const std::string &BScalingFormula = "1", const std::string &EScalingFormula = "1"):
		TFieldContainer(std::move(_field), BScalingFormula, EScalingFormula, 0., 0., 0., 0., 0., 0., 0.) {}


	/**
	 * Calculate magnetic field at coordinates x,y,z taking into account time-dependent scaling and boundary
	 * 
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param B Returns magnetic field vector
	 * @param dBidxj Return spatial derivatives of magnetic field components (optional)
	 */
	void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const;


	/**
	 * Calculate electric potential and field at coordinates x,y,z taking into account time-dependent scaling and boundary
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Return electric potential (!=0 only if a map with potential was loaded)
	 * @param Ei Returns electric field vector
	 */
	void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const;
};





#endif /* FIELD_H_ */
