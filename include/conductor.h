/**
 * \file
 * Class to calculate magnetic fields produced by finite, straight conductors.
 */

#ifndef RACETRACK_H_
#define RACETRACK_H_

#include "field.h"

/**
 * Virtual base class for all straight wire classes.
 */
class TConductorField: public TField{
private:
	double I; ///< current through wire
	double SW1x; ///< x coordinate of start point
	double SW1y; ///< y coordinate of start point
	double SW1z; ///< z coordinate of start point
	double SW2x; ///< x coordinate of end point
	double SW2y; ///< y coordinate of end point
	double SW2z; ///< z coordinate of end point
public:
	/**
	 * Constructor
	 *
	 * Sets current through wire TConductorField::I and field scaling formula
	 */
	TConductorField(const double SW1xx, const double SW1yy, const double SW1zz, const double SW2xx, const double SW2yy, const double SW2zz, const double aI);

	/**
	 * Compute magnetic field of a straight, finite conductor.
	 *
	 * For parameter doc see TField::BField.
	 */
	void BField(const double x, const double y, const double z, const double t,
			double B[3], double dBidxj[3][3]) const override;

	/**
	 * Conductors produce no electric field.
	 *
	 * For parameter doc see TField::EField.
	 */
	void EField(const double x, const double y, const double z, const double t,
			double &V, double Ei[3]) const override {};
};


#endif /*RACETRACK_H_*/
