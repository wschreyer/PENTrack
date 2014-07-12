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
	long double I; ///< current through wire
	/**
	 * Constructor
	 *
	 * Sets current through wire TConductorField::I
	 */
	TConductorField(long double aI);

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
	virtual void BField(long double x, long double y, long double z, long double t, long double B[4][4]) = 0;

	/**
	 * Adds no electric field.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Electric potential
	 * @param Ei Electric field components
	 */
	void EField(long double x, long double y, long double z, long double t, long double &V, long double Ei[3]){};
};

/**
 * Calculates magnetic field by a straight wire SW1->SW2
 */

struct TFiniteWire: public TConductorField{
	long double SW1x; ///< x coordinate of start point
	long double SW1y; ///< y coordinate of start point
	long double SW1z; ///< z coordinate of start point
	long double SW2x; ///< x coordinate of end point
	long double SW2y; ///< y coordinate of end point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWire(long double SW1xx, long double SW1yy, long double SW1zz, long double SW2xx, long double SW2yy, long double SW2zz, long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to x-axis on x-z-plane (y = 0)
 */

struct TFiniteWireX: public TConductorField{
	long double SW1x; ///< x coordinate of start point
	long double SW2x; ///< x coordinate of end point
	long double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireX(long double SW1xx, long double SW2xx, long double SWzz, long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to y-axis on y-z-plane (x = 0)
 */

struct TFiniteWireY: public TConductorField{
	long double SW1y; ///< y coordinate of start point
	long double SW2y; ///< y coordinate of end point
	long double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireY(long double SW1yy, long double SW2yy, long double SWzz, long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};


/**
 * Calculates magnetic field by a vertical finite straight wire
 */
struct TFiniteWireZ: public TConductorField{
	long double SWx; ///< x coordinate of wire
	long double SWy; ///< y coordintate of wire
	long double SW1z; ///< z coordinate of start point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireZ(long double SWxx, long double SWyy, long double SW1zz, long double SW2zz, long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};


/**
 * Calculates magnetic field by a finite straight wire on z-axis
 */
struct TFiniteWireZCenter: public TConductorField{
	long double SW1z; ///< z coordinate of start point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireZCenter(long double SW1zz, long double SW2zz, long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};


/**
 * Calculates magnetic field by four racetrack coils.
 *
 * Calculates magnetic field by four racetrack coils on x-z, (-x)-z, y-z and (-y)-z planes,
 * meeting at z-axis, with bottom SW1z and top SW2z and width SWr
 */
struct TFullRacetrack: public TConductorField{
	long double SW1z; ///< z coordinate of coil bottoms
	long double SW2z; ///< z coordinate of coil tops
	long double SWr; ///< width of coils

	/**
	 * Constructor, requires coordinates, size and current.
	 */
	TFullRacetrack(long double SW1zz, long double SW2zz, long double SWrr, long double aI);

	void BField(long double x, long double y, long double z, long double t, long double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire parallel to z-axis
 */
struct TInfiniteWireZ: public TConductorField{
	long double lx; ///< x coordinate of wire
	long double ly; ///< y coordinate of wire

	/**
	 * Constructor, requires coordinates and current.
	 */
	TInfiniteWireZ(long double lxx, long double lyy, long double aI);

	void BField(long double x,long double y,long double z, long double t, long double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire on z-axis
 */
struct TInfiniteWireZCenter: public TConductorField{
	/**
	 * Constructor, requires current.
	 */
	TInfiniteWireZCenter(long double aI);

	void BField(const long double x, const long double y, const long double z, long double t, long double B[4][4]);
};

#endif /*RACETRACK_H_*/
