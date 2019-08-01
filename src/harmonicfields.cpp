
/**
 * \file
 * Implementation of a magnetic field determined by 
 * coefficients, provided as inputs by the user, of an expansion in
 * terms of harmonic polynomials. This implementation is based on the
 * work described in "Computation_summary.pdf" (!!! unfinished pls
 * update!)
 * 
 * !!! Note that this implementation does not include a calculation of the
 * B field gradients/derivatives for arbitrary fields. 
*/

#include "harmonicfields.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <cmath>

using namespace std;

// HarmonicExpandedBField constructor
HarmonicExpandedBField::HarmonicExpandedBField(const double _xoff, const double _yoff, const double _zoff, const double ang1, const double ang2,
		const double abz, const double adB0zdz, const bool AC, const double frq, const double tstart1, const double tend1, const double pshift, const double bW,
		const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale, 
		const double G0, const double G1, const double G2, const double G3, const double G4, const double G5, const double G6, const double G7, const double G8, 
		const double G9, const double G10, const double G11, const double G12, const double G13, const double G14, const double G15, const double G16, const double G17, 
		const double G18, const double G19, const double G20, const double G21, const double G22, const double G23)
			: TField(Bscale, "0") {

	xoff = _xoff;
	yoff = _yoff;
	zoff = _zoff;
	pol_ang1 = ang1;
	azm_ang2 = ang2;
	edmB0z0 = abz;
	edmdB0z0dz = adB0zdz;
	ac = AC;
	f = frq;
	on1 = tstart1;
	off1 = tend1;
	phase = pshift;
	BoundaryWidth = bW;
	xmax = _xmax;
	xmin = _xmin;
	ymax = _ymax;
	ymin = _ymin;
	zmax = _zmax;
	zmin = _zmin;

    G[0]  = G0;
	G[1]  = G1;
	G[2]  = G2;
	G[3]  = G3;
	G[4]  = G4;
	G[5]  = G5;
	G[6]  = G6;
	G[7]  = G7;
	G[8]  = G8;
	G[9]  = G9;
	G[10] = G10;
	G[11] = G11;
	G[12] = G12;
	G[13] = G13;
	G[14] = G14;
	G[15] = G15;
	G[16] = G16;
	G[17] = G17;
	G[18] = G18;
	G[19] = G19;
	G[20] = G20;
	G[21] = G21;
	G[22] = G22;
	G[23] = G23;

	// !!! Does this rotation of the field itself actually work generally?
	// It is untested.

	// //Rotation Matrix for Bfield
	// Rot1[0][0] = cos(ang1) * cos(ang2);
	// Rot1[1][0] = -sin(ang2);
	// Rot1[2][0] = sin(ang1) * cos(ang2);
	// Rot1[0][1] = cos(ang1) * sin(ang2);
	// Rot1[1][1] = cos(ang2);
	// Rot1[2][1] = sin(ang1) * sin(ang2);
	// Rot1[0][2] = -sin(ang1);
	// Rot1[1][2] = 0;
	// Rot1[2][2] = cos(ang1);
	
	// // generate inverted indicies rotation matrix
	// for(int i=0;i<3;i++){
	// 	for(int j=0;j<3;j++){
	// 		Rot2[j][i]=Rot1[i][j];
	// 	}
	// }
	
	//Rotation matrix for dB
	// Rot3[0][0]=cos(ang1)*cos(ang2);
	// Rot3[0][1]=-sin(ang2);
	// Rot3[0][2]=sin(ang1)*cos(ang2);
	// Rot3[1][0]=cos(ang1)*sin(ang2);
	// Rot3[1][1]=cos(ang2);
	// Rot3[1][2]=sin(ang1)*sin(ang2);
	// Rot3[2][0]=-sin(ang1);
	// Rot3[2][1]=0;
	// Rot3[2][2]=cos(ang1);

	// // Theory values for derivatives in the BField frame
	// dB[0] = -edmdB0z0dz/2; 			//dBxdx
	// dB[1] = 0;				//dBydx
	// dB[2] = 0; 				//dBzdx
	// dB[3] = 0;  				//dBxdy
	// dB[4] = -edmdB0z0dz/2;			//dBydy
	// dB[5] = 0;				//dBzdy
	// dB[6] = 0;				//dBxdz
	// dB[7] = 0;				//dBydz
	// dB[8] = edmdB0z0dz;			//dBzdz
	
	// for (int i = 0; i < 6; i++){
	// 	for (int j = 0; j < 3; j++) {
	// 		Bd[i][j] = 0;
	// 	}
	// }

	// // !!! The dB implementation has been commented everywhere.	
	// // change in B field in this coordinate system to be transformed into 
	// // original coordinate system. 
	// //Bd[0]= Bdx/d(x,y,z) 
	// //Bd[1]= Bdy/d(x,y,z) 
	// //Bd[2]= Bdz/d(x,y,z)
	// Bd[0][0] = dB[0];
	// Bd[0][1] = dB[3];
	// Bd[0][2] = dB[6];
	// Bd[1][0] = dB[1];
	// Bd[1][1] = dB[4];
	// Bd[1][2] = dB[7];
	// Bd[2][0] = dB[2];
	// Bd[2][1] = dB[5];
	// Bd[2][2] = dB[8];
	
	// //apply rotation matrix to Change in Bfield parameters
	// for(int k=0;k<3;k++){
	// 	for(int i=0;i<3;i++){
	// 	   for(int j=0;j<3;j++){
	// 		   Bd[k+3][i]+=Rot3[i][j]*Bd[k][j];
	// 		}
	// 	}
	// }

	// //rearrangements of vectors for 2nd transform of Bfield. 
	// // Bdx = Bd(x',y',z')/dx, where z' is in the direction of the B field
	// Bd[0][0]=Bd[3][0];
	// Bd[0][1]=Bd[4][0];
	// Bd[0][2]=Bd[5][0];
	// Bd[1][0]=Bd[3][1];	//Bdy = Bd(x',y',z')/dy
	// Bd[1][1]=Bd[4][1];
	// Bd[1][2]=Bd[5][1];
	// Bd[2][0]=Bd[3][2];	//Bdz = Bd(x',y',z')/dz
	// Bd[2][1]=Bd[4][2];
	// Bd[2][2]=Bd[5][2];

	// //apply rotation matrix again
	// for(int k=0;k<3;k++){
	// 	for(int i=0;i<3;i++){
	// 		Bd[k+3][i]=0;
	// 		for(int j=0;j<3;j++){
	// 			Bd[k+3][i]+=Rot3[i][j]*Bd[k][j];
	// 		}
	// 	}
	// }

	// // !!! The dB implementation has been commented everywhere.
	// // Update dB to rotated values
	// dB[0]=Bd[3][0];
	// dB[1]=Bd[3][1];
	// dB[2]=Bd[3][2];
	// dB[3]=Bd[4][0];
	// dB[4]=Bd[4][1];
	// dB[5]=Bd[4][2];
	// dB[6]=Bd[5][0];
	// dB[7]=Bd[5][1];
	// dB[8]=Bd[5][2];
	}

void HarmonicExpandedBField::BField(const double x, const double y, 
const double z, const double t, double B[3], double dBidxj[3][3]) const{

	B[0] = 	  G[2] 	\
			+ G[3]  * (y) \
			+ G[5]  * (-x / 2) \
			+ G[6]  * (z) \
			+ G[7]  * (x) \
			+ G[8]  * (2 * x * y) \
			+ G[9]  * (2 * y * z) \
			+ G[10] * (-x * y / 2) \
			+ G[11] * (-x * z) \
			+ G[12] * (-(3 * pow(x,2) + pow(y,2) - 4 * pow(z,2)) / 4) \
			+ G[13] * (2 * x * z) \
			+ G[14] * ((pow(x,2) - pow(y,2))) \
			+ G[15] * ((3 * pow(x,2) * y - pow(y,3))) \
			+ G[16] * (6 * x * y * z) \
			+ G[17] * (-(3 * pow(x,2) * y + pow(y,3) - 6 * y * pow(z,2)) / 2) \
			+ G[18] * (-(3 / 2) * (x * y * z)) \
			+ G[19] * ((3 / 8) * (pow(x,3) + x * pow(y,2) - 4 * x * pow(z,2))) \
			+ G[20] * ((-1 / 4) * (9 * pow(x,2) * z + 3 * pow(y,2) * z - 4 * pow(z,3))) \
			+ G[21] * (-pow(x,3) + 3 * x * pow(z,2)) \
			+ G[22] * (3 * (pow(x,2) * z - pow(y,2) * z)) \
			+ G[23] * (pow(x,3) - 3 * x * pow(y,2));
			
	B[1] =	  G[0] 	\
			+ G[3]  * (x) \
			+ G[4]  * (z) \
			+ G[5]  * (-y / 2) \
			+ G[7]  * (-y) \
			+ G[8]  * (pow(x,2) - pow(y,2)) \
			+ G[9]  * (2 * x * z) \
			+ G[10] * ((pow(x,2) + 3 * pow(y,2) - 4 * pow(z,2)) / 4) \
			+ G[11] * (-y * z) \
			+ G[12] * (-x * y / 2) \
			+ G[13] * (-2 * y * z) \
			+ G[14] * (-2 * x * y) \
			+ G[15] * (pow(x,3) - 3 * x * pow(y,2)) \
			+ G[16] * (3 * (pow(x,2) * z - pow(y,2) * z)) \
			+ G[17] * (-(pow(x,3) + 3 * x * pow(y,2) - 6 * x * pow(z,2)) / 2) \
			+ G[18] * ( (-1 / 4) * (3 * pow(x,2) * z + 9 * pow(y,2) * z - 4 * pow(z,3))) \
			+ G[19] * ((3 / 8) * (pow(x,2) * y + pow(y,3) - 4 * y * pow(z,2))) \
			+ G[20] * ((-3 / 2) * x * y * z) \
			+ G[21] * (-3 * y * pow(z,2) + pow(y,3)) \
			+ G[22] * (-6 * x * y * z) \
			+ G[23] * (-3 * pow(x,2) * y + pow(y,3));

	B[2] =	  G[1] 	\
			+ G[4]  * (y) \
			+ G[5]  * (z) \
			+ G[6]  * (x) \
			+ G[9]  * (2 * x * y) \
			+ G[10] * (2 * y * z) \
			+ G[11] * (pow(z,2) - (1 / 2) * (pow(x,2) + pow(y,2))) \
			+ G[12] * (2 * x * z) \
			+ G[13] * (pow(x,2) - pow(y,2)) \
			+ G[16] * (3 * pow(x,2) * y - pow(y,3))\
			+ G[17] * (6 * x * y * z)  \
			+ G[18] * (3 * y * pow(z,2) - (3/4) * (pow(x,2) * y + pow(y,3))) \
			+ G[19] * (pow(z,3) - (3 / 2) * z * (pow(x,2) + pow(y,2))) \
			+ G[20] * (3 * x * pow(z,2) - (3/4) * (pow(x,3) + x * pow(y,2))) \
			+ G[21] * (3 * (pow(x,2) * z - pow(y,2) * z)) \
			+ G[22] * (pow(x,3) - 3 * x * pow(y,2));

		// !!! What are we doing? This loop doesn't sum the
		// polynomials!

		// !!! Print statements used during testing.
		// std::cout << "B[0]: " << B[0] << "\n\n";
		// std::cout << "B[1]: " << B[1] << "\n\n";
		// std::cout << "B[2]: " << B[2] << "\n\n";

		// for(int i=0;i<3;i++){
		// 	if(std::isnan(B[i])){
				
		// 		printf("Nan found!\n");
		// 		std::cout << "x: " << x << "\n";
		// 		std::cout << "y: " << y << "\n";
		// 		std::cout << "z: " << z << "\n";
		// 		std::cout << "Bx: " << B[0] << "\n";
		// 		std::cout << "By: " << B[1] << "\n";
		// 		std::cout << "Bz: " << B[2] << "\n\n";

		// 	}

		// 	else{
		// 		// // printf("Not a nan!\n\n");
		// 		// std::cout << "Bx: " << B[0] << "\n";
		// 		// std::cout << "By: " << B[1] << "\n";
		// 		// std::cout << "Bz: " << B[2] << "\n\n";
		// 	}
		
		// }

		// // !!! What is the utility of all this other stuff? How would one even 
		// // go about testing it?

		// double Bscale = BScaling(t);
		// if((ac && (t<on1 || t>off1)) || Bscale == 0){
		// 	return;
		// }
		// else{
		
		// 
		// 	//point to be rotated into Bfield reference frame
		// 	double t1[3]={x-edmB0xoff,y-edmB0yoff,z-edmB0zoff};
		// 	double t2[3]={0};

		// 	//rotate specified point to BField coordinate system
		// 	for(int i=0;i<3;i++){
		// 		for(int j=0;j<3;j++){
		// 			t2[i]+=Rot1[i][j]*t1[j];
		// 		}
		// 	}

		// 	//compute Bfield compoenents
		// 	double BF1[3];
		// 	double BF2[3]={0};

		// This just calculates the B field values based on a 
		// simple, 0th order gradient in the z-direction.
		// This would be sufficient for analyzing the magneto-gravitational
		// dephasing effect.
		// 	BF1[0] = -t2[0]/2*edmdB0z0dz;		// Bx
		// 	BF1[1] = -t2[1]/2*edmdB0z0dz;		// By
		// 	BF1[2] = edmB0z0 + edmdB0z0dz*t2[2];	// Bz

		
		// 	//rotate Bfield components back to global coordinates
		// 	for(int i=0;i<3;i++){
		// 		for(int j=0;j<3;j++){
		// 			BF2[i]+=Rot2[i][j]*BF1[j];
		// 		}
		// 	}

			// double dBScaled[9];
		// !!! It's unclear what the point of folding is?
		// 	//Initialize a local instance of dB to be folded
		// 	for(int i = 0; i < 9; i++)
		// 		dBScaled[i] = dB[i];

		// 	//apply AC scaling if necessary
		// 	if(ac){
		// 		double scalar;
		// 		scalar=sin((f*t+phase)*2*M_PI);
		// 		for (int i = 0; i < 3; i++)
		// 			BF2[i] *= scalar;
		// 		if (dBidxj != NULL){
		// 			for (int i = 0; i < 9; i++)
		// 				dBScaled[i]*=scalar;
		// 		}
		// 	}
			
			// //Fold the field near its boundaries
			// for (int i = 0; i < 3; i++){
			// 	FieldSmthr(x, y, z, BF2, dBScaled, i);
			// }

			// !!! Copied from the above, BF2 -> B
			//Fold the field near its boundaries
			// for (int i = 0; i < 3; i++){
			// 	FieldSmthr(x, y, z, B, dBScaled, i);
			// }

			// //set BField to EDM component
			// B[0] = BF2[0]*Bscale;
			// B[1] = BF2[1]*Bscale;
			// B[2] = BF2[2]*Bscale;

			// !!! Here I've repeated the above lines but I don't have the BF2
			// variable in use.
			// B[0] = B[0]*Bscale;
			// B[1] = B[1]*Bscale;
			// B[2] = B[2]*Bscale;

			// if (dBidxj != NULL){
			// 	//set BField gradient to EDM component
			// 	dBidxj[0][0] = dBScaled[0]*Bscale;
			// 	dBidxj[1][0] = dBScaled[1]*Bscale;
			// 	dBidxj[2][0] = dBScaled[2]*Bscale;
			// 	dBidxj[0][1] = dBScaled[3]*Bscale;
			// 	dBidxj[1][1] = dBScaled[4]*Bscale;
			// 	dBidxj[2][1] = dBScaled[5]*Bscale;
			// 	dBidxj[0][2] = dBScaled[6]*Bscale;
			// 	dBidxj[1][2] = dBScaled[7]*Bscale;
			// 	dBidxj[2][2] = dBScaled[8]*Bscale;
			// }
		// }
}



void HarmonicExpandedBField::FieldSmthr(const double x, const double y, const double z, double Bxi[3], double dBScaled[9], const int xi) const{
	
		// Fscale = P(x')*P(y')*P(z') in the boundary where P(xi') = SmthrStp(xi-xi_min / BoundaryWidth) for the lower boundary and SmthrStp(xi-xi_max / BoundaryWidth) for the upper boundary
		double compressionArray[6] = {1,1,1,0,0,0};
		double Fscale = 1; 
		double dBadd[3] = {0, 0, 0};
		
		// Throw an error if two edges to be scaled overlap, otherwise compute how far into each boundary (x,y,z) is
		try{
			CompressionFactor(x,y,z,compressionArray);
		}
		catch (const std::invalid_argument& e){
			std::cout << "max-min distance has to be at least twice the BoundaryWidth!" << "\n";
			exit(-1);
		}
			
		// compute the dBadd (d(Fscale)/dxi) term for all three directions
		for (int i = 0; i < 3; i++){
			if (compressionArray[i] != 1){
				Fscale *=SmthrStp(compressionArray[i]);
				if(compressionArray[i] == 0)
					dBadd[i] = 0;
				else if(compressionArray[i+3] == 1)
					dBadd[i] = -Bxi[xi]*SmthrStpDer(compressionArray[i])/(SmthrStp(compressionArray[i])*BoundaryWidth);
				else if(compressionArray[i+3] == 0)
					dBadd[i] = Bxi[xi]*SmthrStpDer(compressionArray[i])/(SmthrStp(compressionArray[i])*BoundaryWidth);

			}
		}

		// If one or more of the components was scaled, scale the Bfield its relevant directional components
		// Bxi' (Scaled) = Bxi (unscaled) * P(x')*P(y')*P(z')
		// dBxidxj' (Scaled) = P(x')*P(y')*P(z') * dBxixj (Unscaled) + Bxi * (d/dxj P(x')*P(y')*P(z'))
		// Note that for the supper boundary d/dxj P(xj') = -P'(xj')/BoundaryWidth and for the lower boundary d/dxj P(xj') = P'(xj')/BoundaryWidth
		if (Fscale != 1){
			Bxi[xi] *= Fscale; // scale field value
			int j = 0;
				for (int i = 0; i < 7; i = i + 3){
					dBScaled[xi + i] = dBScaled[xi + i]*Fscale + dBadd[j]*Fscale; // scale derivatives according to product rule
					j++;
			}
		}
	}
	
// Compute the percentage into the boundary widthy the coordinate is in the x,y,z directions
void HarmonicExpandedBField::CompressionFactor(const double x, const double y, const double z, double *compFactors) const{
	
	// Do nothing if BoundaryWidth is set to zero
	if (BoundaryWidth != 0){
	
		// throw error if the min/max distances are insufficient (min > max or not twice BoundaryWidth apart)
		if ((xmax - xmin) < 2*BoundaryWidth || (ymax - ymin) < 2*BoundaryWidth || (zmax - zmin) < 2*BoundaryWidth){
			throw std::invalid_argument( "max - min has to be at least twice the BoundaryWidth!" );
		}

		// Scale x,y,z if they are within BoundaryWidth of max or min
		if ((x < xmax) && (x >= (xmax - BoundaryWidth))){
			compFactors[0] = (xmax - x)/BoundaryWidth;
			compFactors[3] = 1;
		}

		if ((x > xmin) && (x <= (xmin + BoundaryWidth)))
			compFactors[0] = (x - xmin)/BoundaryWidth;


		if ((y < ymax) && (y >= (ymax - BoundaryWidth))){
			compFactors[1] = (ymax - y)/BoundaryWidth;
			compFactors[4] = 1;
			}

		if ((y > ymin) && (y <= (ymin + BoundaryWidth)))
			compFactors[1] = (y - ymin)/BoundaryWidth;


		if ((z < zmax) && (z >= (zmax - BoundaryWidth))){
			compFactors[2] = (zmax - z)/BoundaryWidth;
			compFactors[5] = 1;
			}

		if ((z > zmin) && (z <= (zmin + BoundaryWidth)))
			compFactors[2] = (z - zmin)/BoundaryWidth;


		// Nullify the field outside of the min/max
		if (x >= xmax || x <= xmin)
			compFactors[0] = 0;

		if (y >= ymax || y <= ymin)
			compFactors[1] = 0;

		if (z >= zmax || z <= zmin)
			compFactors[2] = 0;
	}
}

double HarmonicExpandedBField::SmthrStp(const double x) const{
	return 6*pow(x, 5) - 15*pow(x, 4) + 10*pow(x, 3);
}

double HarmonicExpandedBField::SmthrStpDer(const double x) const{
	return 30*pow(x, 4) - 60*pow(x, 3) + 30*pow(x,2);
}