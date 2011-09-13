/**
 * \file
 * Several optimized routines to calculate magnetic fields produced by (in)finite straight wires.
 */

#ifndef RACETRACK_H_
#define RACETRACK_H_

/**
 * Virtual base class for all straight wire classes.
 */
struct TRacetrack{
	long double I_rt; ///< current through wire
	TRacetrack(long double I): I_rt(I){}; ///< constructor
	virtual ~TRacetrack(){}; ///< destructor

	/**
	 * Returns magnetic field produced by wire.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param B Magnetic field components matrix
	 */
	virtual void BFeld(long double x, long double y, long double z, long double B[4][4]) = 0;
};

/**
 * Calculates magnetic field by a straight wire SW1->SW2
 */

struct TFiniteWire: public TRacetrack{
	long double SW1x; ///< x coordinate of start point
	long double SW1y; ///< y coordinate of start point
	long double SW1z; ///< z coordinate of start point
	long double SW2x; ///< x coordinate of end point
	long double SW2y; ///< y coordinate of end point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWire(long double SW1xx, long double SW1yy, long double SW1zz, long double SW2xx, long double SW2yy, long double SW2zz, long double I): TRacetrack(I), SW1x(SW1xx), SW1y(SW1yy), SW1z(SW1zz), SW2x(SW2xx), SW2y(SW2yy), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to x-axis on x-z-plane (y = 0)
 */

struct TFiniteWireX: public TRacetrack{
	long double SW1x; ///< x coordinate of start point
	long double SW2x; ///< x coordinate of end point
	long double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireX(long double SW1xx, long double SW2xx, long double SWzz, long double I): TRacetrack(I), SW1x(SW1xx), SW2x(SW2xx), SWz(SWzz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

/**
 * Calculates magnetic field by a finite straight wire parallel to y-axis on y-z-plane (x = 0)
 */

struct TFiniteWireY: public TRacetrack{
	long double SW1y; ///< y coordinate of start point
	long double SW2y; ///< y coordinate of end point
	long double SWz; ///< z coordinate of wire

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireY(long double SW1yy, long double SW2yy, long double SWzz, long double I): TRacetrack(I), SW1y(SW1yy), SW2y(SW2yy), SWz(SWzz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};


/**
 * Calculates magnetic field by a vertical finite straight wire
 */
struct TFiniteWireZ: public TRacetrack{
	long double SWx; ///< x coordinate of wire
	long double SWy; ///< y coordintate of wire
	long double SW1z; ///< z coordinate of start point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireZ(long double SWxx, long double SWyy, long double SW1zz, long double SW2zz, long double I): TRacetrack(I), SWx(SWxx), SWy(SWyy), SW1z(SW1zz), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};


/**
 * Calculates magnetic field by a finite straight wire on z-axis
 */
struct TFiniteWireZCenter: public TRacetrack{
	long double SW1z; ///< z coordinate of start point
	long double SW2z; ///< z coordinate of end point

	/**
	 * Constructor, requires coordinates of start and end point of wire and current.
	 */
	TFiniteWireZCenter(long double SW1zz, long double SW2zz, long double I): TRacetrack(I), SW1z(SW1zz), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};


/**
 * Calculates magnetic field by four racetrack coils.
 *
 * Calculates magnetic field by four racetrack coils on x-z, (-x)-z, y-z and (-y)-z planes,
 * meeting at z-axis, with bottom SW1z and top SW2z and width SWr
 */
struct TFullRacetrack: public TRacetrack{
	long double SW1z; ///< z coordinate of coil bottoms
	long double SW2z; ///< z coordinate of coil tops
	long double SWr; ///< width of coils

	/**
	 * Constructor, requires coordinates, size and current.
	 */
	TFullRacetrack(long double SW1zz, long double SW2zz, long double SWrr, long double I): TRacetrack(I), SW1z(SW1zz), SW2z(SW2zz), SWr(SWrr){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire parallel to z-axis
 */
struct TInfiniteWireZ: public TRacetrack{
	long double lx; ///< x coordinate of wire
	long double ly; ///< y coordinate of wire

	/**
	 * Constructor, requires coordinates and current.
	 */
	TInfiniteWireZ(long double lxx, long double lyy, long double I): TRacetrack(I), lx(lxx), ly(lyy){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};


/**
 * Calculates magnetic field by an infinite straight wire on z-axis
 */
struct TInfiniteWireZCenter: public TRacetrack{
	/**
	 * Constructor, requires current.
	 */
	TInfiniteWireZCenter(long double I): TRacetrack(I){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

#endif /*RACETRACK_H_*/
