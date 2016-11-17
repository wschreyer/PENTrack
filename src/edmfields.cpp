
/**
 * \file
 * The EDM fields used in the Ramsey Cycle
*/

#include <cmath>
#include "edmfields.h"
#include <iostream>
#include <stdlib.h>


//TEDMStaticB0GradZField constructor
TEDMStaticB0GradZField::TEDMStaticB0GradZField(const double xoff, const double yoff, const double zoff, const double ang1, const double ang2,
		const double abz, const double adB0zdz, const bool AC, const double frq, const double tstart1, const double tend1, const double pshift, const double bW,
		const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const std::string &Bscale)
			: TField(Bscale, "0") {
	edmB0xoff = xoff;
	edmB0yoff = yoff;
	edmB0zoff = zoff;
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

	//Rotation Matrix for Bfield
	Rot1[0][0]=cos(ang1)*cos(ang2);
	Rot1[1][0]=-sin(ang2);
	Rot1[2][0]=sin(ang1)*cos(ang2);
	Rot1[0][1]=cos(ang1)*sin(ang2);
	Rot1[1][1]=cos(ang2);
	Rot1[2][1]=sin(ang1)*sin(ang2);
	Rot1[0][2]=-sin(ang1);
	Rot1[1][2]=0;
	Rot1[2][2]=cos(ang1);
	
	//generate inverted indicies rotation matrix
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			Rot2[j][i]=Rot1[i][j];
		}
	}
	
	//Rotation matrix for dB
	Rot3[0][0]=cos(ang1)*cos(ang2);
	Rot3[0][1]=-sin(ang2);
	Rot3[0][2]=sin(ang1)*cos(ang2);
	Rot3[1][0]=cos(ang1)*sin(ang2);
	Rot3[1][1]=cos(ang2);
	Rot3[1][2]=sin(ang1)*sin(ang2);
	Rot3[2][0]=-sin(ang1);
	Rot3[2][1]=0;
	Rot3[2][2]=cos(ang1);


	//Theory values for derivatives in the BField frame
	dB[0] = -edmdB0z0dz/2; 			//dBxdx
	dB[1] = 0;				//dBydx
	dB[2] = 0; 				//dBzdx
	dB[3] = 0;  				//dBxdy
	dB[4] = -edmdB0z0dz/2;			//dBydy
	dB[5] = 0;				//dBzdy
	dB[6] = 0;				//dBxdz
	dB[7] = 0;				//dBydz
	dB[8] = edmdB0z0dz;			//dBzdz
	
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 3; j++) {
			Bd[i][j] = 0;
		}
	}
	
	// change in B field in this coordinate system to be transformed into original coordinate system. //Bd[0]= Bdx/d(x,y,z) //Bd[1]= Bdy/d(x,y,z) //Bd[2]= Bdz/d(x,y,z)
	Bd[0][0] = dB[0];
	Bd[0][1] = dB[3];
	Bd[0][2] = dB[6];
	Bd[1][0] = dB[1];
	Bd[1][1] = dB[4];
	Bd[1][2] = dB[7];
	Bd[2][0] = dB[2];
	Bd[2][1] = dB[5];
	Bd[2][2] = dB[8];
	

	
	//apply rotation matrix to Change in Bfield parameters
	for(int k=0;k<3;k++){
		for(int i=0;i<3;i++){
		   for(int j=0;j<3;j++){
			   Bd[k+3][i]+=Rot3[i][j]*Bd[k][j];
			}
		}
	}


	//rearrangements of vectors for 2nd transform of Bfield. Bdx = Bd(x',y',z')/dx, where z' is in the direction of the B field
	Bd[0][0]=Bd[3][0];
	Bd[0][1]=Bd[4][0];
	Bd[0][2]=Bd[5][0];
	Bd[1][0]=Bd[3][1];	//Bdy = Bd(x',y',z')/dy
	Bd[1][1]=Bd[4][1];
	Bd[1][2]=Bd[5][1];
	Bd[2][0]=Bd[3][2];	//Bdz = Bd(x',y',z')/dz
	Bd[2][1]=Bd[4][2];
	Bd[2][2]=Bd[5][2];

	//apply rotation matrix again
	for(int k=0;k<3;k++){
		for(int i=0;i<3;i++){
			Bd[k+3][i]=0;
			for(int j=0;j<3;j++){
				Bd[k+3][i]+=Rot3[i][j]*Bd[k][j];
			}
		}
	}

	// Update dB to rotated values
	dB[0]=Bd[3][0];
	dB[1]=Bd[3][1];
	dB[2]=Bd[3][2];
	dB[3]=Bd[4][0];
	dB[4]=Bd[4][1];
	dB[5]=Bd[4][2];
	dB[6]=Bd[5][0];
	dB[7]=Bd[5][1];
	dB[8]=Bd[5][2];

	}

void TEDMStaticB0GradZField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

	if(ac==true && (t<on1 || t>off1)){
		return;
	}
	else{
		
		//point to be rotated into Bfield reference frame
		double t1[3]={x-edmB0xoff,y-edmB0yoff,z-edmB0zoff};
		double t2[3]={0};

		//rotate specified point to BField coordinate system
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				t2[i]+=Rot1[i][j]*t1[j];
			}
		}

		//compute Bfield compoenents
		double BF1[3];
		double BF2[3]={0};

		BF1[0] = -t2[0]/2*edmdB0z0dz;		// Bx
		BF1[1] = -t2[1]/2*edmdB0z0dz;		// By
		BF1[2] = edmB0z0 + edmdB0z0dz*t2[2];	// Bz

	
		//rotate Bfield components back to global coordinates
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				BF2[i]+=Rot2[i][j]*BF1[j];
			}
		}

		double dBScaled[9];
		//Initialize a local instance of dB to be folded
		for(int i = 0; i < 9; i++)
			dBScaled[i] = dB[i];

		//apply AC scaling if necessary
		if(ac){
			double scalar;
			scalar=sin((f*t+phase)*2*M_PI);
			for (int i = 0; i < 3; i++)
				BF2[i] *= scalar;
			if (dBidxj != NULL){
				for (int i = 0; i < 9; i++)
					dBScaled[i]*=scalar;
			}
		}
		
		//Fold the field near its boundaries
		for (int i = 0; i < 3; i++){
			FieldSmthr(x, y, z, BF2, dBScaled, i);
		}

		//set BField to EDM component
		double Bscale = BScaling(t);
		B[0] = BF2[0]*Bscale;
		B[1] = BF2[1]*Bscale;
		B[2] = BF2[2]*Bscale;

		if (dBidxj != NULL){
			//set BField gradient to EDM component
			dBidxj[0][0] = dBScaled[0]*Bscale;
			dBidxj[1][0] = dBScaled[1]*Bscale;
			dBidxj[2][0] = dBScaled[2]*Bscale;
			dBidxj[0][1] = dBScaled[3]*Bscale;
			dBidxj[1][1] = dBScaled[4]*Bscale;
			dBidxj[2][1] = dBScaled[5]*Bscale;
			dBidxj[0][2] = dBScaled[6]*Bscale;
			dBidxj[1][2] = dBScaled[7]*Bscale;
			dBidxj[2][2] = dBScaled[8]*Bscale;
		}
	}
}

void TEDMStaticB0GradZField::FieldSmthr(const double x, const double y, const double z, double Bxi[3], double dBScaled[9], const int xi) const{
	
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
void TEDMStaticB0GradZField::CompressionFactor(const double x, const double y, const double z, double *compFactors) const{
	
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

double TEDMStaticB0GradZField::SmthrStp(const double x) const{
	return 6*pow(x, 5) - 15*pow(x, 4) + 10*pow(x, 3);
}

double TEDMStaticB0GradZField::SmthrStpDer(const double x) const{
	return 30*pow(x, 4) - 60*pow(x, 3) + 30*pow(x,2);
}

//TEDMStaticEField constructor
TEDMStaticEField::TEDMStaticEField(const double aexMag, const double aeyMag, const double aezMag, const std::string &Escale): TField("1", Escale){
	exMag = aexMag;
	eyMag = aeyMag;
	ezMag = aezMag;
}

void TEDMStaticEField::EField (const double x, const double y, const double z, const double t, double &V, double Ei[3], double dEidxj[3][3]) const{
	double Escale = EScaling(t);
	Ei[0] = exMag*Escale;
	Ei[1] = eyMag*Escale;
	Ei[2] = ezMag*Escale;

	V = Ei[2]*z*Escale;
	
}
