
/**
 * \file
 * The EDM fields used in the Ramsey Cycle
*/

#include <cmath>
#include "edmfields.h"

TEDMStaticB0GradZField::TEDMStaticB0GradZField(double abz, double adB0zdz): edmB0z0(abz), edmdB0z0dz(adB0zdz) {}; //end constructor

void TEDMStaticB0GradZField::BField(double x, double y, double z, double t, double B[4][4]){
	
	//Calculations for the first row of the matrix i.e. all the x parts
	B[0][0] += -x/2*edmdB0z0dz; // Bx
	B[0][1] += -edmdB0z0dz/2; //dBxdx
	B[0][2] += 0; B[0][3] += 0; //dBxdy, dBxdz

    	//Calculations for the second row of the matrix i.e. all the y parts
    	B[1][0] += -y/2*edmdB0z0dz; // By
    	B[1][1] += 0; B[1][3] += 0; //dBydx, dBy/dz
    	B[1][2] += -edmdB0z0dz/2; //dBydy

    	//Calculations for the third row of the matrix i.e. all the z parts
   	B[2][0] += edmB0z0 + edmdB0z0dz*z; // Bz
   	B[2][1] += 0; B[2][2] += 0; //dBzdx, dBzdy
    	B[2][3] += edmdB0z0dz; //dBzdz

    	//Calculations for the fourth row of the matrix i.e. all |B| parts
    	B[3][0] += sqrt((x*x + y*y + 4*z*z)*edmdB0z0dz*edmdB0z0dz/4 + 2*edmB0z0*edmdB0z0dz*z + edmB0z0*edmB0z0); //(should give the same results as calculation below)
//    	B[3][0] += sqrt(pow(B[0][0], 2) + pow(B[1][0], 2) + pow(B[2][0], 2)); // |B| (not correct when there will be other fields)
    	B[3][1] += (x*pow(edmdB0z0dz,2))/sqrt((x*x + y*y + 4*z*z)*edmdB0z0dz*edmdB0z0dz/4 + 2*edmB0z0*edmdB0z0dz*z + edmB0z0*edmB0z0); // d|B|/dx
    	B[3][2] += (y*pow(edmdB0z0dz,2))/sqrt((x*x + y*y + 4*z*z)*edmdB0z0dz*edmdB0z0dz/4 + 2*edmB0z0*edmdB0z0dz*z + edmB0z0*edmB0z0); //d|B|/dy
    	B[3][3] += (pow(edmdB0z0dz, 2)*z + edmB0z0*edmdB0z0dz)/sqrt((x*x + y*y + 4*z*z)*edmdB0z0dz*edmdB0z0dz/4 + 2*edmB0z0*edmdB0z0dz*z + edmB0z0*edmB0z0); //d|B|/dz

} // end TEDMStaticB0GradZField::BField

TEDMStaticEField::TEDMStaticEField ( double amagEField): magEField(amagEField) {}; //end TEDMStaticEField constructor

void TEDMStaticEField::EField (double x, double y, double z, double t, double &V, double Ei[3]) {
	Ei[0] += 0; Ei[1] += 0;
	Ei[2] += magEField;

	//V = integrate(E*r) => V = E_z * z
	V += Ei[2]*z;
} //end TEDMStaticEField::EField 
