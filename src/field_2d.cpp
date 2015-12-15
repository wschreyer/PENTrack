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


void TabField::ReadTabFile(const char *tabfile, double Bscale, double Escale, alglib::real_1d_array &rind, alglib::real_1d_array &zind,
		alglib::real_1d_array BTabs[3], alglib::real_1d_array ETabs[3], alglib::real_1d_array &VTab){
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
		BTabs[0].setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("RBY") != string::npos){
		BTabs[1].setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("RBZ") != string::npos){
		BTabs[2].setlength(m*n);
		getline(FIN,line);
	}

	if (line.find("EX") != string::npos){
		ETabs[0].setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("EY") != string::npos){
		ETabs[1].setlength(m*n);
		getline(FIN,line);
	}
	if (line.find("EZ") != string::npos){
		ETabs[2].setlength(m*n);
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
		if (BTabs[0].length() > 0){
			FIN >> val;
			BTabs[0][i2] = val*Bconv*Bscale;
		}
		if (BTabs[1].length() > 0){
			FIN >> val;
			BTabs[1][i2] = val*Bconv*Bscale;
		}
		if (BTabs[2].length() > 0){
			FIN >> val;
			BTabs[2][i2] = val*Bconv*Bscale;
		}
		if (ETabs[0].length() > 0){
			FIN >> val;
			ETabs[0][i2] = val*Econv*Escale;
		}
		if (ETabs[1].length() > 0){
			FIN >> val;
			ETabs[1][i2] = val*Econv*Escale;
		}
		if (ETabs[2].length() > 0){
			FIN >> val;
			ETabs[2][i2] = val*Econv*Escale;
		}
		if (VTab.length() > 0){
			FIN >> val;
			VTab[i2] = val*Escale;
		}
		FIN >> ws;
	}

	if (Bscale == 0){
		BTabs[0].setlength(0);
		BTabs[1].setlength(0);
		BTabs[2].setlength(0);
	}
	if (Escale == 0){
		ETabs[0].setlength(0);
		ETabs[0].setlength(0);
		ETabs[0].setlength(0);
		VTab.setlength(0);
	}

	printf("\n");
	if (ri+1 != (int)m || zi+1 != (int)n){
		printf("The header says the size is %u by %u, actually it is %i by %i! Exiting...\n", m, n, ri+1, zi+1);
		exit(-1);
	}
	FIN.close();
}


void TabField::CheckTab(alglib::real_1d_array &rind, alglib::real_1d_array &zind,
		alglib::real_1d_array BTabs[3], alglib::real_1d_array ETabs[3], alglib::real_1d_array &VTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
	r_mi = rind[0];
	z_mi = zind[0];
	rdist = rind[1] - rind[0];
	zdist = zind[1] - zind[0];
	printf("The arrays are %u by %u.\n",m,n);
	printf("The r values go from %G to %G\n", rind[0], rind[m-1]);
	printf("The z values go from %G to %G.\n", zind[0], zind[n-1]);
	printf("rdist = %G zdist = %G\n", rdist, zdist);

	double Babsmax = 0, Babsmin = 9e99, Babs;
	double Vmax = 0, Vmin = 9e99;
	for (int j=0; j < m; j++)
	{
		for (int k=0; k < n; k++)
		{
			Babs = 0;
			int i2 = k*m + j;
			if (BTabs[0].length() > 0) Babs += BTabs[0][i2]*BTabs[0][i2];
			if (BTabs[1].length() > 0) Babs += BTabs[1][i2]*BTabs[1][i2];
			if (BTabs[2].length() > 0) Babs += BTabs[2][i2]*BTabs[2][i2];
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
		double aNullFieldTime, double aRampUpTime, double aFullFieldTime, double aRampDownTime,
		double alengthconv, double aBconv, double aEconv){
	NullFieldTime = aNullFieldTime;
	RampUpTime = aRampUpTime;
	FullFieldTime = aFullFieldTime;
	RampDownTime = aRampDownTime;
	lengthconv = alengthconv;
	Bconv = aBconv;
	Econv = aEconv;
	alglib::real_1d_array rind, zind, BTabs[3], ETabs[3], VTab;

	ReadTabFile(tabfile, Bscale, Escale, rind, zind, BTabs, ETabs, VTab); // open tabfile and read values into arrays

	CheckTab(rind, zind, BTabs, ETabs, VTab); // print some info

	printf("Starting Preinterpolation ... ");
	if (BTabs[0].length() > 0){
		cout << "Br ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[0], 1, Brc);
	}
	if (BTabs[1].length() > 0){
		cout << "Bhi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[1], 1, Bphic);
	}
	if (BTabs[2].length() > 0){
		cout << "Bz ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[2], 1, Bzc);
	}
	if (ETabs[0].length() > 0){
		cout << "Er ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[0], 1, Erc);
		VTab.setlength(0); // ignore potential if electric field map found
	}
	if (ETabs[1].length() > 0){
		cout << "Ephi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[1], 1, Ephic);
		VTab.setlength(0); // ignore potential if electric field map found
	}
	if (ETabs[2].length() > 0){
		cout << "Ez ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[2], 1, Ezc);
		VTab.setlength(0); // ignore potential if electric field map found
	}
	if (VTab.length() > 0){
		cout << "V ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, VTab, 1, Vc);
	}
}

TabField::~TabField(){
}


void TabField::BField(double x, double y, double z, double t, double B[4][4]){
	double r = sqrt(x*x+y*y);
	double Bscale = BFieldScale(t);
	if (Bscale != 0 && r >= r_mi && r <= r_mi + rdist*(m - 1) && z >= z_mi && z <= z_mi + zdist*(n - 1)){
		// bicubic interpolation
		double Br = 0, dBrdr = 0, dBrdz = 0, Bphi = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0;
		double Bx = 0, By = 0, Bz = 0, dBxdz = 0, dBydz = 0, dBzdz = 0;
		double dummy;
		double phi = atan2(y,x);
		if (Brc.c_ptr()->m*Brc.c_ptr()->n > 0){
			alglib::spline2ddiff(Brc, r, z, Br, dBrdr, dBrdz, dummy);
		}
		if (Bphic.c_ptr()->m*Bphic.c_ptr()->n > 0){
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
		if (Bzc.c_ptr()->m*Bzc.c_ptr()->n > 0){
			alglib::spline2ddiff(Bzc, r, z, Bz, dBzdr, dBzdz, dummy);
		}
		B[2][0] += Bz*Bscale;
		B[2][1] += dBzdr*cos(phi)*Bscale;
		B[2][2] += dBzdr*sin(phi)*Bscale;
		B[2][3] += dBzdz*Bscale;
	}
}


void TabField::EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3]){
	double r = sqrt(x*x+y*y);
	if (r >= r_mi && r <= r_mi + rdist*(m - 1) && z >= z_mi && z <= z_mi + zdist*(n - 1)){
		if (Erc.c_ptr()->m*Erc.c_ptr()->n > 0 || Ephic.c_ptr()->m*Ephic.c_ptr()->n > 0 || Ezc.c_ptr()->m*Ezc.c_ptr()->n > 0){ // prefer E-field interpolation over potential interpolation
			double phi = atan2(y,x);
			if (dEidxj == NULL){
				alglib::real_1d_array Er, Ephi, Ez;
				double Ex, Ey;
				if (Erc.c_ptr()->m*Erc.c_ptr()->n > 0)
					alglib::spline2dcalcv(Erc, r, z, Er); // use cheaper interpolation if derivatives not required
				if (Ephic.c_ptr()->m*Ephic.c_ptr()->n > 0)
					alglib::spline2dcalcv(Ephic, r, z, Ephi);
				if (Ezc.c_ptr()->m*Ezc.c_ptr()->n > 0)
					alglib::spline2dcalcv(Ezc, r, z, Ez);
				CylToCart(Er[0], Ephi[0], phi, Ex, Ey); // convert r,phi components to x,y components
				Ei[0] += Ex; // add electric field
				Ei[1] += Ey;
				Ei[2] += Ez[0];
			}
			else{
				double Er = 0, Ephi = 0, Ex = 0, Ey = 0, Ez = 0, dummy;
				double dErdr, dErdz, dEphidr, dEphidz, dEzdr, dEzdz, dExdz, dEydz;
				//bicubic interpolation
				if (Erc.c_ptr()->m*Erc.c_ptr()->n > 0)
					alglib::spline2ddiff(Erc, r, z, Er, dErdr, dErdz, dummy);
				if (Ephic.c_ptr()->m*Ephic.c_ptr()->n > 0)
					alglib::spline2ddiff(Ephic, r, z, Ephi, dEphidr, dEphidz, dummy);
				if (Ezc.c_ptr()->m*Ezc.c_ptr()->n > 0)
					alglib::spline2ddiff(Ezc, r, z, Ez, dEzdr, dEzdz, dummy);
				CylToCart(Er,Ephi,phi,Ex,Ey); // convert r,phi components to x,y components
				Ei[0] += Ex; // add electric field
				Ei[1] += Ey;
				Ei[2] += Ez;
				if (dEidxj != NULL){
					if (r > 0){ // convert r,phi derivatives to x,y derivatives
						dEidxj[0][0] += (dErdr*cos(phi)*cos(phi) - dEphidr*cos(phi)*sin(phi) + (Er*sin(phi)*sin(phi) + Ephi*cos(phi)*sin(phi))/r);
						dEidxj[0][1] += (dErdr*cos(phi)*sin(phi) - dEphidr*sin(phi)*sin(phi) - (Er*cos(phi)*sin(phi) + Ephi*cos(phi)*cos(phi))/r);
						dEidxj[1][0] += (dErdr*cos(phi)*sin(phi) + dEphidr*cos(phi)*cos(phi) - (Er*cos(phi)*sin(phi) - Ephi*sin(phi)*sin(phi))/r);
						dEidxj[1][1] += (dErdr*sin(phi)*sin(phi) + dEphidr*cos(phi)*sin(phi) + (Er*cos(phi)*cos(phi) - Ephi*cos(phi)*sin(phi))/r);
					}
					CylToCart(dErdz,dEphidz,phi,dExdz,dEydz); // convert z-derivatives of r,phi components into z-derivatives of x,y components
					dEidxj[0][2] += dExdz;
					dEidxj[1][2] += dEydz;
					dEidxj[2][0] += dEzdr*cos(phi);
					dEidxj[2][1] += dEzdr*sin(phi);
					dEidxj[2][2] += dEzdz;
				}
			}
		}
		else{
			double Vloc, dVdrj[3], dummy;
			// bicubic interpolation
			alglib::spline2ddiff(Vc, r, z, Vloc, dVdrj[0], dVdrj[2], dummy);
			double phi = atan2(y,x);
			V += Vloc;
			Ei[0] += -dVdrj[0]*cos(phi);
			Ei[1] += -dVdrj[0]*sin(phi);
			Ei[2] += -dVdrj[2];
			// !!!!!!!!!!! calculation of electric field derivatives not yet implemented for potential interpolation !!!!!!!!!!!!!!!!
		}

	}
}
