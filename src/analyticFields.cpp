/**
 * \file
 * Static, analytical B fields
 *
 * Time-dependent scaling not implemented
 */

#include "analyticFields.h"

#include <iostream>

//TExponentialBFieldX constructor
TExponentialFieldX::TExponentialFieldX(const double _a1, const double _a2, const double _a3,
									const double _c1, const double _c2, const double _xmax, const double _xmin,
									const double _ymax, const double _ymin, const double _zmax, const double _zmin)
									: TField("1", "0") {
	a1 = _a1; a2 = _a2; a3 = _a3;
	c1 = _c1; c2 = _c2;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TExponentialFieldX::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	if (withinBounds(x, y, z)) {
		B[0] = a1 * exp(- a2* x + a3) + c1; //Bx contribution
		B[1] = y * a1 * a2 / 2 * exp(- a2* x + a3) + c2; //By contribution
		B[2] = z* a1 * a2 / 2 * exp(- a2* x + a3) + c2; //Bz contribution
		if (dBidxj != nullptr){
			dBidxj[0][0] = -a2 * a1 * exp(- a2* x + a3); // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = -a2 * x * y * a1 * a2 / 2 * exp(- a2* x + a3) + c2; //dBydx
			dBidxj[1][1] = a1 * a2 / 2 * exp(- a2* x + a3); //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = -a2 * x * z * a1 * a2 / 2 * exp(- a2* x + a3) + c2; //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = a1 * a2 / 2 * exp(- a2* x + a3); //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != nullptr){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TExponentialFieldX::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}

//TLinearFieldZ constructor
TLinearFieldZ::TLinearFieldZ(const double _a1, const double _a2, const double _xmax, const double _xmin,
									const double _ymax, const double _ymin, const double _zmax, const double _zmin)
									: TField("1", "0") {
	a1 = _a1; a2 = _a2;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TLinearFieldZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	if (withinBounds(x, y, z)) {
		B[0] = 0; //Bx contribution
		B[1] = 0; //By contribution
		B[2] = a1*x + a2; //Bz contribution
		if (dBidxj != nullptr){
			dBidxj[0][0] = 0; // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = 0; //dBydx
			dBidxj[1][1] = 0; //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = a1; //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = 0; //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != nullptr){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TLinearFieldZ::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}
