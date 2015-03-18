/**
 * \file
 * Bicubic interpolation of axisymmetric field tables.
 */

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include "field_2d.h"
#include "globals.h"

using namespace std;


/**
 * Translate VECTOR (not point) components from cylindrical to cartesian coordinates.
 *
 * @param v_r Radial component of vector
 * @param v_phi Azimuthal component of vector
 * @param phi Azimuth of vector origin
 * @param v_x Returns x component of vector
 * @param v_y Returns y component of vector
 */
void CylToCart(double v_r, double v_phi, double phi, double &v_x, double &v_y){
	v_x = v_r*cos(phi) - v_phi*sin(phi);
	v_y = v_r*sin(phi) + v_phi*cos(phi);
}


void TabField::ReadTabFile(const char *tabfile, double Bscale, double Escale){
	ifstream FIN(tabfile, ifstream::in);
	if (!FIN.is_open()){
		cout << "\nCould not open " << tabfile << "!\n";
		exit(-1);
	}
	cout << "\nReading " << tabfile << "!\n";
	int intval;
	string line;
	FIN >> m >> intval >> n;
	rind.setlength(m);
	zind.setlength(n);

	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	bool skipy = true;
	if (line.substr(0,12) == " 2 Y [LENGU]")  skipy = false;
	getline(FIN,line);
	if (!skipy) getline(FIN,line);

	if (line.find("RBX") != string::npos){
		BrTab.setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("RBY") != string::npos){
		BphiTab.setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("RBZ") != string::npos){
		BzTab.setlength(m*n);
		getline(FIN,line);
	}

	if (line.find("EX") != string::npos){
		ErTab.setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("EY") != string::npos){
		EphiTab.setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("EZ") != string::npos){
		EzTab.setlength(m*n);
		getline(FIN,line);
	}

	if (line.find("RV") != string::npos){	// file contains potential?
		VTab.setlength(m*n);
		getline(FIN,line);
	}

	if (!FIN || line.substr(0,2) != " 0"){
		printf("%s not found or corrupt! Exiting...\n",tabfile);
		exit(-1);
	}

	int ri = 0,zi = -1, perc = 0;
	double r, z, val;
	while (FIN.good()){
		FIN >> r;
		if (!skipy) FIN >> val;
		FIN >> z;
		if (!FIN) break;
		r *= lengthconv;
		z *= lengthconv;
		if (zi >= 0 && z < zind[zi]){
			ri++;
			zi = 0;
		}
		else zi++;

		// status if read is displayed
		PrintPercent((float)zi*ri/(n*m), perc);


		rind[ri] = r;
		zind[zi] = z;
		int i2 = zi * m + ri;
		if (BrTab.length() > 0){
			FIN >> val;
			BrTab[i2] = val*Bconv*Bscale;
		}
		if (BphiTab.length() > 0){
			FIN >> val;
			BphiTab[i2] = val*Bconv*Bscale;
		}
		if (BzTab.length() > 0){
			FIN >> val;
			BzTab[i2] = val*Bconv*Bscale;
		}
		if (ErTab.length() > 0){
			FIN >> val;
			ErTab[i2] = val*Econv*Escale;
		}
		if (EphiTab.length() > 0){
			FIN >> val;
			EphiTab[i2] = val*Econv*Escale;
		}
		if (EzTab.length() > 0){
			FIN >> val;
			EzTab[i2] = val*Econv*Escale;
		}
		if (VTab.length() > 0){
			FIN >> val;
			VTab[i2] = val*Escale;
		}
		FIN >> ws;
	}

	if (Bscale == 0){
		BrTab.setlength(0);
		BphiTab.setlength(0);
		BzTab.setlength(0);
	}
	if (Escale == 0){
		ErTab.setlength(0);
		EphiTab.setlength(0);
		EzTab.setlength(0);
		VTab.setlength(0);
	}

	printf("\n");
	if (ri+1 != (int)m || zi+1 != (int)n){
		printf("The header says the size is %u by %u, actually it is %i by %i! Exiting...\n", m, n, ri+1, zi+1);
		exit(-1);
	}
	FIN.close();
}


void TabField::CheckTab(){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
	printf("The arrays are %u by %u.\n",m,n);
	printf("The r values go from %G to %G\n", rind[0], rind[m-1]);
	printf("The z values go from %G to %G.\n", zind[0], zind[n-1]);
	printf("rdist = %G zdist = %G\n", rind[1] - rind[0], zind[1] - zind[0]);

	double Babsmax = 0, Babsmin = 9e99, Babs;
	double Vmax = 0, Vmin = 9e99;
	for (int j=0; j < m; j++)
	{
		for (int k=0; k < n; k++)
		{
			Babs = 0;
			int i2 = k*m + j;
			if (BrTab.length() > 0) Babs += BrTab[i2]*BrTab[i2];
			if (BphiTab.length() > 0) Babs += BphiTab[i2]*BphiTab[i2];
			if (BzTab.length() > 0) Babs += BzTab[i2]*BzTab[i2];
			Babsmax = max(sqrt(Babs),Babsmax);
			Babsmin = min(sqrt(Babs),Babsmin);
			if (VTab.length() > 0)
			{
				Vmax = max(VTab[i2],Vmax);
				Vmin = min(VTab[i2],Vmin);
			}

		}
	}

	printf("The input table file has values of |B| from %G T to %G T and values of V from %G V to %G V\n",Babsmin,Babsmax,Vmin,Vmax);
}


double TabField::BFieldScale(double t){
	if (t < NullFieldTime || t >= NullFieldTime + RampUpTime + FullFieldTime + RampDownTime)
		return 0;
	else if (t >= NullFieldTime && t < NullFieldTime + RampUpTime){
		// ramping up field smoothly with cosine
		//return 0.5 - 0.5*cos(pi*(t - NullFieldTime)/RampUpTime);

		// linear ramp
		return (t - NullFieldTime)/RampUpTime;
	}
	else if (t >= NullFieldTime + RampUpTime && t < NullFieldTime + RampUpTime + FullFieldTime)
		return 1;
	else if (t >= NullFieldTime + RampUpTime + FullFieldTime && t < NullFieldTime + RampUpTime + FullFieldTime + RampDownTime){
		// ramping down field smoothly with cosine
		//return 0.5 + 0.5*cos(pi*(t - (RampUpTime + NullFieldTime + FullFieldTime)) / RampDownTime);

		// linear ramp
		return (1 - (t - RampUpTime - NullFieldTime - FullFieldTime)/RampDownTime);
	}
	else return 0;
}


TabField::TabField(const char *tabfile, double Bscale, double Escale,
		double aNullFieldTime, double aRampUpTime, double aFullFieldTime, double aRampDownTime){
	NullFieldTime = aNullFieldTime;
	RampUpTime = aRampUpTime;
	FullFieldTime = aFullFieldTime;
	RampDownTime = aRampDownTime;

	ReadTabFile(tabfile,Bscale,Escale); // open tabfile and read values into arrays

	CheckTab(); // print some info

	printf("Starting Preinterpolation ... ");
	if (BrTab.length() > 0){
		cout << "Br ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BrTab, 1, Brc);
	}
	if (BphiTab.length() > 0){
		cout << "Bhi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BphiTab, 1, Bphic);
	}
	if (BzTab.length() > 0){
		cout << "Bz ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BzTab, 1, Bzc);
	}
	if (VTab.length() > 0){
		cout << "V ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, VTab, 1, Vc);
		ErTab.setlength(0); // if there is a potential, we don't need the E-vector
		EphiTab.setlength(0);
		EzTab.setlength(0);
	}
	if (ErTab.length() > 0){
		cout << "Er ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ErTab, 1, Erc);
	}
	if (EphiTab.length() > 0){
		cout << "Ephi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, EphiTab, 1, Ephic);
	}
	if (EzTab.length() > 0){
		cout << "Ez ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, EzTab, 1, Ezc);
	}
}

TabField::~TabField(){
}


void TabField::BField(double x, double y, double z, double t, double B[4][4]){
	double r = sqrt(x*x+y*y);
	double Bscale = BFieldScale(t);
	if (Bscale != 0 && r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
		// bicubic interpolation
		double Br = 0, dBrdr = 0, dBrdz = 0, Bphi = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0;
		double Bx = 0, By = 0, Bz = 0, dBxdz = 0, dBydz = 0, dBzdz = 0;
		double dummy;
		double phi = atan2(y,x);
		if (BrTab.length() > 0){
			alglib::spline2ddiff(Brc, r, z, Br, dBrdr, dBrdz, dummy);
		}
		if (BphiTab.length() > 0){
			alglib::spline2ddiff(Bphic, r, z, Bphi, dBphidr, dBphidz, dummy);
		}
		CylToCart(Br,Bphi,phi,Bx,By);
		B[0][0] += Bx*Bscale;
		B[1][0] += By*Bscale;
		if (r > 0){
			B[0][1] += Bscale*(dBrdr*cos(phi)*cos(phi) - dBphidr*cos(phi)*sin(phi) + (Br*sin(phi)*sin(phi) + Bphi*cos(phi)*sin(phi))/r);
			B[0][2] += Bscale*(dBrdr*cos(phi)*sin(phi) - dBphidr*sin(phi)*sin(phi) - (Br*cos(phi)*sin(phi) + Bphi*cos(phi)*cos(phi))/r);
			B[1][1] += Bscale*(dBrdr*cos(phi)*sin(phi) + dBphidr*cos(phi)*cos(phi) - (Br*cos(phi)*sin(phi) - Bphi*sin(phi)*sin(phi))/r);
			B[1][2] += Bscale*(dBrdr*sin(phi)*sin(phi) + dBphidr*cos(phi)*sin(phi) + (Br*cos(phi)*cos(phi) - Bphi*cos(phi)*sin(phi))/r);
		}
		CylToCart(dBrdz,dBphidz,phi,dBxdz,dBydz);
		B[0][3] += dBxdz*Bscale;
		B[1][3] += dBydz*Bscale;
		if (BzTab.length() > 0){
			alglib::spline2ddiff(Bzc, r, z, Bz, dBzdr, dBzdz, dummy);
		}
		B[2][0] += Bz*Bscale;
		B[2][1] += dBzdr*cos(phi)*Bscale;
		B[2][2] += dBzdr*sin(phi)*Bscale;
		B[2][3] += dBzdz*Bscale;
	}
}


void TabField::EField(double x, double y, double z, double t, double &V, double Ei[3]){
	double Vloc, dVdrj[3], dummy;
	if (VTab.length() > 0){ // prefer E-field from potential over pure E-field interpolation
		double r = sqrt(x*x+y*y);
		if (r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
			// bicubic interpolation
			alglib::spline2ddiff(Vc, r, z, Vloc, dVdrj[0], dVdrj[2], dummy);
			double phi = atan2(y,x);
			V += Vloc;
			Ei[0] += -dVdrj[0]*cos(phi);
			Ei[1] += -dVdrj[0]*sin(phi);
			Ei[2] += -dVdrj[2];
		}
	}
	else if (ErTab.length() > 0 || EphiTab.length() > 0 || EzTab.length() > 0){
		double r = sqrt(x*x+y*y);
		if (r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
			// bicubic interpolation
			double Er = 0, Ephi = 0, Ex = 0, Ey = 0, Ez = 0;
			double phi = atan2(y,x);
			if (ErTab.length() > 0)
				Er = alglib::spline2dcalc(Erc, r, z);
			if (EphiTab.length() > 0)
				Ephi = alglib::spline2dcalc(Ephic, r, z);
			if (EzTab.length() > 0)
				Ez = alglib::spline2dcalc(Ezc, r, z);
			CylToCart(Er,Ephi,phi,Ex,Ey);
			Ei[0] += Ex;
			Ei[1] += Ey;
			Ei[2] += Ez;
		}
	}
}
