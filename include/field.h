/*
 * field.h
 *
 *  Contains virtual base class for field calculation methods
 */

#ifndef FIELD_H_
#define FIELD_H_

#include <cstddef>
#include <iostream>

#include "exprtk.hpp"

/**
 * Virtual base class for all field calculation methods
 */
class TField{
private:
	exprtk::expression<double> Bscaler;
	exprtk::expression<double> Escaler;
	mutable double tvar; ///< time variable for use in scaling-formula parsers

protected:
	/**
	 * Calculate magnetic field scaling from parsed formula
	 *
	 * @param t Time
	 *
	 * @return Return scaling factor
	 */
	double BScaling(const double t) const{
		tvar = t;
		return Bscaler.value();
	};

	/**
	 * Calculate electric field scaling from parsed formula
	 *
	 * @param t Time
	 *
	 * @return Return scaling factor
	 */
	double EScaling(const double t) const{
		tvar = t;
		return Escaler.value();
	};

public:
	/**
	 * Add magnetic field at a given position and time.
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
			double B[3], double dBidxj[3][3] = NULL) const = 0;

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
	virtual void EField (const double x, const double y, const double z, const double t,
			double &V, double Ei[3], double dEidxj[3][3] = NULL) const = 0;

	/**
	 * Generic constructor, should be called by every derived class.
	 *
	 * @param Bscale String containing formula describing time-dependence of magnetic field
	 * @param Escale String containing formula describing time-dependence of electric field
	 */
	TField(const std::string &Bscale, const std::string &Escale){
		exprtk::symbol_table<double> symbol_table;
		symbol_table.add_variable("t",tvar);
		symbol_table.add_constants();
		Bscaler.register_symbol_table(symbol_table);
		exprtk::parser<double> parser;
		parser.compile(Bscale, Bscaler);

		Escaler.register_symbol_table(symbol_table);
		parser.compile(Escale, Escaler);
	}

	/**
	 * Virtual destructor
	 */
	virtual ~TField(){ };
};




#endif /* FIELD_H_ */
