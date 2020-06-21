/**
 * \file
 * Static, analytical B fields
 *
 */

#include "analyticFields.h"

using namespace std;

//TExponentialBFieldX constructor
TExponentialFieldX::TExponentialFieldX(const double _a1, const double _a2, const double _a3, const double _c1, const double _c2){
	a1 = _a1; a2 = _a2; a3 = _a3;
	c1 = _c1; c2 = _c2;
}

void TExponentialFieldX::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = a1 * exp(- a2* x + a3) + c1; //Bx contribution
	B[1] = y * a1 * a2 / 2 * exp(- a2* x + a3) + c2; //By contribution
	B[2] = z* a1 * a2 / 2 * exp(- a2* x + a3) + c2; //Bz contribution
	if (dBidxj != NULL){
		dBidxj[0][0] = -a2 * a1 * exp(- a2* x + a3); // dBxdx
		dBidxj[0][1] = 0; //dBxdy
		dBidxj[0][2] = 0; //dBxdz
		dBidxj[1][0] = -a2 * y * a1 * a2 / 2 * exp(- a2* x + a3); //dBydx
		dBidxj[1][1] = a1 * a2 / 2 * exp(- a2* x + a3); //dBydy
		dBidxj[1][2] = 0; //dBydz
		dBidxj[2][0] = -a2 * z * a1 * a2 / 2 * exp(- a2* x + a3); //dBzdx
		dBidxj[2][1] = 0; //dBzdy
		dBidxj[2][2] = a1 * a2 / 2 * exp(- a2* x + a3); //dBzdz
	}
}


//TLinearFieldZ constructor
TLinearFieldZ::TLinearFieldZ(const double _a1, const double _a2){
	a1 = _a1; a2 = _a2;
}

void TLinearFieldZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = 0; //Bx contribution
	B[1] = 0; //By contribution
	B[2] = (a1*x + a2); //Bz contribution
	if (dBidxj != NULL){
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
}


//TB0GradZ constructor
TB0GradZ::TB0GradZ(const double _a1, const double _a2, const double _z0){
	a1 = _a1; a2 = _a2; z0 = _z0;
}

void TB0GradZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = (-(a1*z + a2)/2 * x); //Bx contribution
	B[1] = (-(a1*z + a2)/2 * y); //By contribution
	B[2] = (a1*z*z/2 + a2*z + z0); //Bz contribution
	if (dBidxj != NULL){
		dBidxj[0][0] = (-(a1*z + a2) / 2); // dBxdx
		dBidxj[0][1] = 0; //dBxdy
		dBidxj[0][2] = (- a1*x / 2); //dBxdz
		dBidxj[1][0] = 0; //dBydx
		dBidxj[1][1] = (-(a1*z + a2) / 2); //dBydy
		dBidxj[1][2] = (- a1*y / 2); //dBydz
		dBidxj[2][0] = 0; //dBzdx
		dBidxj[2][1] = 0; //dBzdy
		dBidxj[2][2] = (a1*z + a2); //dBzdz
	}
}


//TB0GradX2 constructor
TB0GradX2::TB0GradX2(const double _a1, const double _a2, const double _a3, const double _z0){
	a1 = _a1; a2 = _a2; a3 = _a3; z0 = _z0;
}

void TB0GradX2::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = ( -a1/6*x*x*x - a2/4*x*x -a3/2*x ); //Bx contribution
	B[1] = ( -( a1*x*x + a2*x + a3 ) / 2 * y ); //By contribution
	B[2] = ((a1*x*x + a2*x + a3) * z + z0); //Bz contribution
	if (dBidxj != NULL){
		dBidxj[0][0] = (-( a1*x*x + a2*x + a3 ) / 2); // dBxdx
		dBidxj[0][1] = 0; //dBxdy
		dBidxj[0][2] = 0; //dBxdz
		dBidxj[1][0] = (-a1*x -a2/2) * y; //dBydx
		dBidxj[1][1] = (-( a1*x*x + a2*x + a3 ) / 2); //dBydy
		dBidxj[1][2] = 0; //dBydz
		dBidxj[2][0] = (2*a1*x + a2)*z; //dBzdx
		dBidxj[2][1] = 0; //dBzdy
		dBidxj[2][2] = ( a1*x*x + a2*x + a3 ); //dBzdz
	}
}


//TB0GradXY constructor
TB0GradXY::TB0GradXY(const double _a1, const double _a2, const double _z0){
	a1 = _a1; a2 = _a2; z0 = _z0;
}

void TB0GradXY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = (-a1*x*x*y/4 - a2*x/2); //Bx contribution
	B[1] = (-a1*x*y*y/4 - a2*y/2); //By contribution
	B[2] = (a1*x*y*z + a2*z +z0); //Bz contribution
	if (dBidxj != NULL){
		dBidxj[0][0] = (-(a1*x*y +a2) / 2); // dBxdx
		dBidxj[0][1] = (-a1*x*x/4); //dBxdy
		dBidxj[0][2] = 0; //dBxdz
		dBidxj[1][0] = (-a1*y*y/4); //dBydx
		dBidxj[1][1] = (-(a1*x*y +a2) / 2); //dBydy
		dBidxj[1][2] = 0; //dBydz
		dBidxj[2][0] = a1*y*z; //dBzdx
		dBidxj[2][1] = a1*x*z; //dBzdy
		dBidxj[2][2] = (a1*x*y + a2); //dBzdz
	}
}


//TB0_XY constructor
TB0_XY::TB0_XY(const double _a1, const double _z0) {
	a1 = _a1; z0 = _z0;
}

void TB0_XY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	B[0] = (a1 * y * z); //Bx contribution
	B[1] = (a1 * x * z); //By contribution
	B[2] = (a1 * x * y + z0); //Bz contribution
	if (dBidxj != NULL){
		dBidxj[0][0] = 0; // dBxdx
		dBidxj[0][1] = (a1 * z); //dBxdy
		dBidxj[0][2] = (a1 * y); //dBxdz
		dBidxj[1][0] = (a1 * z); //dBydx
		dBidxj[1][1] = 0; //dBydy
		dBidxj[1][2] = (a1 * x); //dBydz
		dBidxj[2][0] = (a1 * y); //dBzdx
		dBidxj[2][1] = (a1 * x); //dBzdy
		dBidxj[2][2] = 0; //dBzdz
	}
}


TCustomBField::TCustomBField(const std::string &_Bx, const std::string &_By, const std::string &_Bz){
	tvar = unique_ptr<double>(new double(0.0));
	xvar = unique_ptr<double>(new double(0.0));
	yvar = unique_ptr<double>(new double(0.0));
	zvar = unique_ptr<double>(new double(0.0));
	exprtk::symbol_table<double> symbol_table;
	symbol_table.add_variable("t",*tvar);
	symbol_table.add_variable("x",*xvar);
	symbol_table.add_variable("y",*yvar);
	symbol_table.add_variable("z",*zvar);
	symbol_table.add_constants();
	exprtk::parser<double> parser;

	std::array<std::string, 3> expr{_Bx, _By, _Bz};
	for (int i = 0; i < 3; ++i){
		Bexpr[i].register_symbol_table(symbol_table);
		if (not parser.compile(expr[i], Bexpr[i])){
			throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing CustomBField formula '" + expr[i] + "': " + parser.get_error(0).diagnostic);
		}
	}
}

void TCustomBField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	*xvar = x;
	*yvar = y;
	*zvar = z;
	*tvar = t;
	B[0] = Bexpr[0].value();
	B[1] = Bexpr[1].value();
	B[2] = Bexpr[2].value();
//	std::cout << B[0] << " " << B[1] << " " << B[2] << " ";
	
	if (dBidxj != nullptr){
		dBidxj[0][0] = exprtk::derivative(Bexpr[0], *xvar);
		dBidxj[0][1] = exprtk::derivative(Bexpr[0], *yvar);
		dBidxj[0][2] = exprtk::derivative(Bexpr[0], *zvar);
		dBidxj[1][0] = exprtk::derivative(Bexpr[1], *xvar);
		dBidxj[1][1] = exprtk::derivative(Bexpr[1], *yvar);
		dBidxj[1][2] = exprtk::derivative(Bexpr[1], *zvar);
		dBidxj[2][0] = exprtk::derivative(Bexpr[2], *xvar);
		dBidxj[2][1] = exprtk::derivative(Bexpr[2], *yvar);
		dBidxj[2][2] = exprtk::derivative(Bexpr[2], *zvar);
//		std::cout << dBidxj[0][0] << " " << dBidxj[0][1] << " " << dBidxj[0][2] << std::endl;
	}
//	std::cout << std::endl;

}