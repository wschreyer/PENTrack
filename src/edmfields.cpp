
/**
 * \file
 * The EDM fields used in the Ramsey Cycle
*/

#include "edmfields.h"

#include <iostream>

//TEDMStaticB0GradZField constructor
TEDMStaticB0GradZField::TEDMStaticB0GradZField(const double xoff, const double yoff, const double zoff, const double ang1, const double ang2,
		const double abz, const double adB0zdz) {
	edmB0xoff = xoff;
	edmB0yoff = yoff;
	edmB0zoff = zoff;
	pol_ang1 = ang1;
	azm_ang2 = ang2;
	edmB0z0 = abz;
	edmdB0z0dz = adB0zdz;

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

	//set BField to EDM component
	B[0] = BF2[0];
	B[1] = BF2[1];
	B[2] = BF2[2];

	if (dBidxj != nullptr){
		//set BField gradient to EDM component
		dBidxj[0][0] = dB[0];
		dBidxj[1][0] = dB[1];
		dBidxj[2][0] = dB[2];
		dBidxj[0][1] = dB[3];
		dBidxj[1][1] = dB[4];
		dBidxj[2][1] = dB[5];
		dBidxj[0][2] = dB[6];
		dBidxj[1][2] = dB[7];
		dBidxj[2][2] = dB[8];
	}
}


//TEDMStaticEField constructor
TEDMStaticEField::TEDMStaticEField(const double aexMag, const double aeyMag, const double aezMag){
	exMag = aexMag;
	eyMag = aeyMag;
	ezMag = aezMag;
}

void TEDMStaticEField::EField (const double x, const double y, const double z, const double t, double &V, double Ei[3]) const{
	Ei[0] = exMag;
	Ei[1] = eyMag;
	Ei[2] = ezMag;

	V = Ei[2]*z;
	
}
