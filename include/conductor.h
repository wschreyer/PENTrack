/**
 * \file
 * Several optimized routines to calculate magnetic fields produced by (in)finite straight wires.
 */

#ifndef RACETRACK_H_
#define RACETRACK_H_

#include "field.h"

/**
 * Virtual base class for all straight wire classes.
 */
struct TConductorField: public TField{
	double I; ///< current through wire

	/**
	 * Constructor
	 *
	 * Sets current through wire TConductorField::I and field scaling formula
	 */
	TConductorField(double aI, std::string Bscale);

	/**
	 * Compute magnetic field.
	 *
	 * Virtual: has to be implented for all derived classes.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param B Magnetic field components
	 */
	virtual void BField(double x, double y, double z, double t, double B[4][4]) = 0;

	/**
	 * Adds no electric field.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Electric potential
	 * @param Ei Electric field components
	 * @param dEidxj Spatial derivatives of electric field components
	 */
	void EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3] = NULL){};
};

/**
 * Calculates magnetic field by a straight wire SW1->SW2
 */

struct TFiniteWire: public TConductorField{
	double SW1x; ///< x coordinate of start point
	double SW1y; ///< y coordinate of start point
	double SW1z; ///< z coordinate of start point
	double SW2x; ///< x coordinate of end point
	double SW2y; ///< y coordinate of end point
	double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire, current, and field scaling formula.
	 */
	TFiniteWire(double SW1xx, double SW1yy, double SW1zz, double SW2xx, double SW2yy, double SW2zz, double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to x-axis on x-z-plane (y = 0)
 */

struct TFiniteWireX: public TConductorField{
	double SW1x; ///< x coordinate of start point
	double SW2x; ///< x coordinate of end point
	double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire, current, and field scaling formula.
	 */
	TFiniteWireX(double SW1xx, double SW2xx, double SWzz, double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to y-axis on y-z-plane (x = 0)
 */

struct TFiniteWireY: public TConductorField{
	double SW1y; ///< y coordinate of start point
	double SW2y; ///< y coordinate of end point
	double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire, current, and field scaling formula.
	 */
	TFiniteWireY(double SW1yy, double SW2yy, double SWzz, double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};


/**
 * Calculates magnetic field by a vertical finite straight wire
 */
struct TFiniteWireZ: public TConductorField{
	double SWx; ///< x coordinate of wire
	double SWy; ///< y coordintate of wire
	double SW1z; ///< z coordinate of start point
	double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire, current, and field scaling formula.
	 */
	TFiniteWireZ(double SWxx, double SWyy, double SW1zz, double SW2zz, double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};


/**
 * Calculates magnetic field by a finite straight wire on z-axis
 */
struct TFiniteWireZCenter: public TConductorField{
	double SW1z; ///< z coordinate of start point
	double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire, current, and field scaling formula.
	 */
	TFiniteWireZCenter(double SW1zz, double SW2zz, double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};


/**
 * Calculates magnetic field by four racetrack coils.
 *
 * Calculates magnetic field by four racetrack coils on x-z, (-x)-z, y-z and (-y)-z planes,
 * meeting at z-axis, with bottom SW1z and top SW2z and width SWr
 */
struct TFullRacetrack: public TConductorField{
	double SW1z; ///< z coordinate of coil bottoms
	double SW2z; ///< z coordinate of coil tops
	double SWr; ///< width of coils

	/**
	 * Constructor, requires coordinates, size, current, and field scaling formula.
	 */
	TFullRacetrack(double SW1zz, double SW2zz, double SWrr, double aI, std::string Bscale);

	void BField(double x, double y, double z, double t, double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire parallel to z-axis
 */
struct TInfiniteWireZ: public TConductorField{
	double lx; ///< x coordinate of wire
	double ly; ///< y coordinate of wire

	/**
	 * Constructor, requires coordinates, current, and field scaling formula.
	 */
	TInfiniteWireZ(double lxx, double lyy, double aI, std::string Bscale);

	void BField(double x,double y,double z, double t, double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire on z-axis
 */
struct TInfiniteWireZCenter: public TConductorField{
	/**
	 * Constructor, requires current and field scaling formula.
	 */
	TInfiniteWireZCenter(double aI, std::string Bscale);

	void BField(const double x, const double y, const double z, double t, double B[4][4]);
};

#endif /*RACETRACK_H_*/
