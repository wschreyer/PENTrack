/*
 * field.h
 *
 *  Contains virtual base class for field calculation methods
 */

#ifndef FIELD_H_
#define FIELD_H_

#include <cstddef>

#include "muParser.h"

/**
 * Virtual base class for all field calculation methods
 */
class TField{
private:
	mu::Parser Bscaler; ///< formula parser for time-dependent magnetic field scaling formula
	mu::Parser Escaler; ///< formula parser for time-dependent electric field scaling formula
	double tvar; ///< time variable for use in scaling formular parsers

protected:
	/**
	 * Calculate magnetic field scaling from parsed formula
	 *
	 * @param t Time
	 *
	 * @return Return scaling factor
	 */
	double BScaling(double t){
		tvar = t;
		try{
			return Bscaler.Eval();
		}
		catch (mu::Parser::exception_type &e)
		{
			std::cout << e.GetMsg() << std::endl;
		}
		return 0;
	};

	/**
	 * Calculate electric field scaling from parsed formula
	 *
	 * @param t Time
	 *
	 * @return Return scaling factor
	 */
	double EScaling(double t){
		tvar = t;
		try{
			return Escaler.Eval();
		}
		catch (mu::Parser::exception_type &e)
		{
			std::cout << e.GetMsg() << std::endl;
		}
		return 0;
	};

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
	 * @param B Returns magnetic field component matrix
	 */
	virtual void BField (double x, double y, double z, double t, double B[4][4]) = 0;

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
	 * @param dEidxj Returns spatial derivatives of electric field components (optional)
	 */
	virtual void EField (double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3] = NULL) = 0;

	/**
	 * Generic constructor, should be called by every derived class.
	 */
	TField(std::string Bscale, std::string Escale){
		Bscaler.DefineVar("t", &tvar);
		Bscaler.SetExpr(Bscale);
		Escaler.DefineVar("t", &tvar);
		Escaler.SetExpr(Escale);
	}

	/**
	 * Virtual destructor
	 */
	virtual ~TField(){ };
};




#endif /* FIELD_H_ */
