/**
 * \file
 * Static, analytical B fields
 *
 */

#include "analyticFields.h"

#include <iostream>

//TExponentialBFieldX constructor
TExponentialFieldX::TExponentialFieldX(const double _a1, const double _a2, const double _a3, const double _c1, const double _c2,
									const double _xmax, const double _xmin, const double _ymax, const double _ymin,
									const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; a2 = _a2; a3 = _a3;
	c1 = _c1; c2 = _c2;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TExponentialFieldX::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0) ) {
		B[0] = Bscale * (a1 * exp(- a2* x + a3) + c1); //Bx contribution
		B[1] = Bscale * (y * a1 * a2 / 2 * exp(- a2* x + a3) + c2); //By contribution
		B[2] = Bscale * (z* a1 * a2 / 2 * exp(- a2* x + a3) + c2); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = Bscale * (-a2 * a1 * exp(- a2* x + a3)); // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = Bscale * (-a2 * y * a1 * a2 / 2 * exp(- a2* x + a3)); //dBydx
			dBidxj[1][1] = Bscale * (a1 * a2 / 2 * exp(- a2* x + a3)); //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = Bscale * (-a2 * z * a1 * a2 / 2 * exp(- a2* x + a3)); //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = Bscale * (a1 * a2 / 2 * exp(- a2* x + a3)); //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
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
TLinearFieldZ::TLinearFieldZ(const double _a1, const double _a2, const double _xmax, const double _xmin, const double _ymax,
									const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; a2 = _a2;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TLinearFieldZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0)) {
		B[0] = 0; //Bx contribution
		B[1] = 0; //By contribution
		B[2] = Bscale * (a1*x + a2); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = 0; // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = 0; //dBydx
			dBidxj[1][1] = 0; //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = Bscale * a1; //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = 0; //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
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

//TB0GradZ constructor
TB0GradZ::TB0GradZ(const double _a1, const double _a2, const double _z0, const double _xmax, const double _xmin, const double _ymax,
									const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; a2 = _a2; z0 = _z0;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TB0GradZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0)) {
		B[0] = Bscale * (-(a1*z + a2)/2 * x); //Bx contribution
		B[1] = Bscale * (-(a1*z + a2)/2 * y); //By contribution
		B[2] = Bscale * (a1*z*z/2 + a2*z + z0); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = Bscale * (-(a1*z + a2) / 2); // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = Bscale * (- a1*x / 2); //dBxdz
			dBidxj[1][0] = 0; //dBydx
			dBidxj[1][1] = Bscale * (-(a1*z + a2) / 2); //dBydy
			dBidxj[1][2] = Bscale * (- a1*y / 2); //dBydz
			dBidxj[2][0] = 0; //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = Bscale * (a1*z + a2); //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TB0GradZ::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}

//TB0GradX2 constructor
TB0GradX2::TB0GradX2(const double _a1, const double _a2, const double _a3, const double _z0, const double _xmax, const double _xmin, const double _ymax,
									const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; a2 = _a2; a3 = _a3; z0 = _z0;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;

}

void TB0GradX2::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0)) {
		B[0] = Bscale * ( -a1/6*x*x*x - a2/4*x*x -a3/2*x ); //Bx contribution
		B[1] = Bscale * ( -( a1*x*x + a2*x + a3 ) / 2 * y ); //By contribution
		B[2] = Bscale * ((a1*x*x + a2*x + a3) * z + z0); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = Bscale * (-( a1*x*x + a2*x + a3 ) / 2); // dBxdx
			dBidxj[0][1] = 0; //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = Bscale * (-a1*x -a2/2) * y; //dBydx
			dBidxj[1][1] = Bscale * (-( a1*x*x + a2*x + a3 ) / 2); //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = Bscale * (2*a1*x + a2)*z; //dBzdx
			dBidxj[2][1] = 0; //dBzdy
			dBidxj[2][2] = Bscale * ( a1*x*x + a2*x + a3 ); //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TB0GradX2::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}

//TB0GradXY constructor
TB0GradXY::TB0GradXY(const double _a1, const double _a2, const double _z0, const double _xmax, const double _xmin, const double _ymax,
									const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; a2 = _a2; z0 = _z0;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;
}

void TB0GradXY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0)) {
		B[0] = Bscale * (-a1*x*x*y/4 - a2*x/2); //Bx contribution
		B[1] = Bscale * (-a1*x*y*y/4 - a2*y/2); //By contribution
		B[2] = Bscale * (a1*x*y*z + a2*z +z0); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = Bscale * (-(a1*x*y +a2) / 2); // dBxdx
			dBidxj[0][1] = Bscale * (-a1*x*x/4); //dBxdy
			dBidxj[0][2] = 0; //dBxdz
			dBidxj[1][0] = Bscale * (-a1*y*y/4); //dBydx
			dBidxj[1][1] = Bscale * (-(a1*x*y +a2) / 2); //dBydy
			dBidxj[1][2] = 0; //dBydz
			dBidxj[2][0] = Bscale * a1*y*z; //dBzdx
			dBidxj[2][1] = Bscale * a1*x*z; //dBzdy
			dBidxj[2][2] = Bscale * (a1*x*y + a2); //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TB0GradXY::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}


//TB0_XY constructor
TB0_XY::TB0_XY(const double _a1, const double _z0, const double _xmax, const double _xmin, const double _ymax,
									const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
									: TField(Bscale, "0") {
	a1 = _a1; z0 = _z0;
	xmax = _xmax; xmin = _xmin;
	ymax = _ymax; ymin = _ymin;
	zmax = _zmax; zmin = _zmin;
}

void TB0_XY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	double Bscale = BScaling(t);	// Time dependent scaling

	if (withinBounds(x, y, z) && (Bscale != 0)) {
		B[0] = Bscale * (a1 * y * z); //Bx contribution
		B[1] = Bscale * (a1 * x * z); //By contribution
		B[2] = Bscale * (a1 * x * y + z0); //Bz contribution
		if (dBidxj != NULL){
			dBidxj[0][0] = 0; // dBxdx
			dBidxj[0][1] = Bscale * (a1 * z); //dBxdy
			dBidxj[0][2] = Bscale * (a1 * y); //dBxdz
			dBidxj[1][0] = Bscale * (a1 * z); //dBydx
			dBidxj[1][1] = 0; //dBydy
			dBidxj[1][2] = Bscale * (a1 * x); //dBydz
			dBidxj[2][0] = Bscale * (a1 * y); //dBzdx
			dBidxj[2][1] = Bscale * (a1 * x); //dBzdy
			dBidxj[2][2] = 0; //dBzdz
		}
	} else {
		//Field is 0 outside of bounds
		for (int i = 0; i < 3; i++){
			B[i] = 0;
			if (dBidxj != NULL){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] = 0;
			}
		}
	}
}


int TB0_XY::withinBounds(const double x, const double y, const double z) const{
	return ((x >= this->xmin) and (x <= this->xmax) and (y >= this->ymin)
					and (y <= this->ymax)	and (z >= this->zmin) and (z <= this->zmax));
}
