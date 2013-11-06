/*
 * field.h
 *
 *  Contains virtual base class for field calculation methods
 */

#ifndef FIELD_H_
#define FIELD_H_

/**
 * Virtual base class for all field calculation methods
 */
class TField{
public:
	/**
	 * Add magnetic field at a given position and time.
	 *
	 * Add field components to the given field matrix B:
	 *	Bx,		dBxdx,	dBxdy,	dBxdz;
	 *	By,		dBydx,	dBydy,	dBydz;
	 *	Bz,		dBzdx,	dBzdy,	dBzdz;
	 *	Babs,	dBdx,	dBdy,	dBdz;
	 * Has to be implemented by all derived field calculation classes.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param B Magnetic field component matrix to which the values are added
	 */
	virtual void BField (long double x, long double y, long double z, long double t, long double B[4][4]) = 0;

	/**
	 * Add electric field and potential at a given position.
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
	virtual void EField (long double x, long double y, long double z, long double t, long double &V, long double Ei[3]) = 0;

	/**
	 * Virtual destructor
	 */
	virtual ~TField(){ };
};




#endif /* FIELD_H_ */
