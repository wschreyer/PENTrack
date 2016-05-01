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


void TabField::ReadTabFile(const char *tabfile, alglib::real_1d_array &rind, alglib::real_1d_array &zind,
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
			BTabs[0][i2] = val*Bconv;
		}
		if (BTabs[1].length() > 0){
			FIN >> val;
			BTabs[1][i2] = val*Bconv;
		}
		if (BTabs[2].length() > 0){
			FIN >> val;
			BTabs[2][i2] = val*Bconv;
		}
		if (ETabs[0].length() > 0){
			FIN >> val;
			ETabs[0][i2] = val*Econv;
		}
		if (ETabs[1].length() > 0){
			FIN >> val;
			ETabs[1][i2] = val*Econv;
		}
		if (ETabs[2].length() > 0){
			FIN >> val;
			ETabs[2][i2] = val*Econv;
		}
		if (VTab.length() > 0){
			FIN >> val;
			VTab[i2] = val;
		}
		FIN >> ws;
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

TabField::TabField(const char *tabfile, std::string Bscale, std::string Escale, double alengthconv, double aBconv, double aEconv)
	: TField(Bscale, Escale){
	lengthconv = alengthconv;
	Bconv = aBconv;
	Econv = aEconv;
	alglib::real_1d_array rind, zind, BTabs[3], ETabs[3], VTab;

	ReadTabFile(tabfile, rind, zind, BTabs, ETabs, VTab); // open tabfile and read values into arrays

	CheckTab(rind, zind, BTabs, ETabs, VTab); // print some info

	printf("Starting Preinterpolation ... ");
	fBrc = fBphic = fBzc = fErc = fEphic = fEzc = fVc = false;
	if (BTabs[0].length() > 0){
		cout << "Br ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[0], 1, Brc);
		fBrc = true;
	}
	if (BTabs[1].length() > 0){
		cout << "Bhi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[1], 1, Bphic);
		fBphic = true;
	}
	if (BTabs[2].length() > 0){
		cout << "Bz ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, BTabs[2], 1, Bzc);
		fBzc = true;
	}
	if (ETabs[0].length() > 0){
		cout << "Er ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[0], 1, Erc);
		VTab.setlength(0); // ignore potential if electric field map found
		fErc = true;
	}
	if (ETabs[1].length() > 0){
		cout << "Ephi ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[1], 1, Ephic);
		VTab.setlength(0); // ignore potential if electric field map found
		fEphic = true;
	}
	if (ETabs[2].length() > 0){
		cout << "Ez ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, ETabs[2], 1, Ezc);
		VTab.setlength(0); // ignore potential if electric field map found
		fEzc = true;
	}
	if (VTab.length() > 0){
		cout << "V ... ";
		cout.flush();
		alglib::spline2dbuildbicubicv(rind, m, zind, n, VTab, 1, Vc);
		fVc = true;
	}
}

TabField::~TabField(){
}


void TabField::BField(double x, double y, double z, double t, double B[4][4]){
	double r = sqrt(x*x+y*y);
	double Bscale = BScaling(t);
	if (Bscale != 0 && r >= r_mi && r <= r_mi + rdist*(m - 1) && z >= z_mi && z <= z_mi + zdist*(n - 1)){
		// bicubic interpolation
		double Br = 0, dBrdr = 0, dBrdz = 0, Bphi = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0;
		double Bx = 0, By = 0, Bz = 0, dBxdz = 0, dBydz = 0, dBzdz = 0;
		double dummy;
		double phi = atan2(y,x);
		if (fBrc){
			alglib::spline2ddiff(Brc, r, z, Br, dBrdr, dBrdz, dummy);
		}
		if (fBphic){
			alglib::spline2ddiff(Bphic, r, z, Bphi, dBphidr, dBphidz, dummy);
		}
		CylToCart(Br,Bphi,phi,Bx,By);
		B[0][0] = Bx*Bscale;
		B[1][0] = By*Bscale;
		if (r > 0){
			B[0][1] = Bscale*(dBrdr*cos(phi)*cos(phi) - dBphidr*cos(phi)*sin(phi) + (Br*sin(phi)*sin(phi) + Bphi*cos(phi)*sin(phi))/r);
			B[0][2] = Bscale*(dBrdr*cos(phi)*sin(phi) - dBphidr*sin(phi)*sin(phi) - (Br*cos(phi)*sin(phi) + Bphi*cos(phi)*cos(phi))/r);
			B[1][1] = Bscale*(dBrdr*cos(phi)*sin(phi) + dBphidr*cos(phi)*cos(phi) - (Br*cos(phi)*sin(phi) - Bphi*sin(phi)*sin(phi))/r);
			B[1][2] = Bscale*(dBrdr*sin(phi)*sin(phi) + dBphidr*cos(phi)*sin(phi) + (Br*cos(phi)*cos(phi) - Bphi*cos(phi)*sin(phi))/r);
		}
		CylToCart(dBrdz,dBphidz,phi,dBxdz,dBydz);
		B[0][3] = dBxdz*Bscale;
		B[1][3] = dBydz*Bscale;
		if (fBzc){
			alglib::spline2ddiff(Bzc, r, z, Bz, dBzdr, dBzdz, dummy);
		}
		B[2][0] = Bz*Bscale;
		B[2][1] = dBzdr*cos(phi)*Bscale;
		B[2][2] = dBzdr*sin(phi)*Bscale;
		B[2][3] = dBzdz*Bscale;
	}
}


void TabField::EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3]){
	double r = sqrt(x*x+y*y);
	double Escale = EScaling(t);
	if (Escale != 0 && r >= r_mi && r <= r_mi + rdist*(m - 1) && z >= z_mi && z <= z_mi + zdist*(n - 1)){
		if (fErc || fEphic || fEzc){ // prefer E-field interpolation over potential interpolation
			double phi = atan2(y,x);
			if (dEidxj == NULL){
				alglib::real_1d_array Er("[0]"), Ephi("[0]"), Ez("[0]");
				double Ex, Ey;
				if (fErc)
					alglib::spline2dcalcvbuf(Erc, r, z, Er); // use cheaper interpolation if derivatives not required
				if (fEphic)
					alglib::spline2dcalcvbuf(Ephic, r, z, Ephi);
				if (fEzc)
					alglib::spline2dcalcvbuf(Ezc, r, z, Ez);
				CylToCart(Er[0], Ephi[0], phi, Ex, Ey); // convert r,phi components to x,y components
				Ei[0] = Ex*Escale; // set electric field
				Ei[1] = Ey*Escale;
				Ei[2] = Ez[0]*Escale;
			}
			else{
				double Er = 0, Ephi = 0, Ex = 0, Ey = 0, Ez = 0, dummy;
				double dErdr, dErdz, dEphidr, dEphidz, dEzdr, dEzdz, dExdz, dEydz;
				//bicubic interpolation
				if (fErc)
					alglib::spline2ddiff(Erc, r, z, Er, dErdr, dErdz, dummy);
				if (fEphic)
					alglib::spline2ddiff(Ephic, r, z, Ephi, dEphidr, dEphidz, dummy);
				if (fEzc)
					alglib::spline2ddiff(Ezc, r, z, Ez, dEzdr, dEzdz, dummy);
				CylToCart(Er,Ephi,phi,Ex,Ey); // convert r,phi components to x,y components
				Ei[0] = Ex*Escale; // set electric field
				Ei[1] = Ey*Escale;
				Ei[2] = Ez*Escale;
				if (r > 0){ // convert r,phi derivatives to x,y derivatives
					dEidxj[0][0] = Escale*(dErdr*cos(phi)*cos(phi) - dEphidr*cos(phi)*sin(phi) + (Er*sin(phi)*sin(phi) + Ephi*cos(phi)*sin(phi))/r);
					dEidxj[0][1] = Escale*(dErdr*cos(phi)*sin(phi) - dEphidr*sin(phi)*sin(phi) - (Er*cos(phi)*sin(phi) + Ephi*cos(phi)*cos(phi))/r);
					dEidxj[1][0] = Escale*(dErdr*cos(phi)*sin(phi) + dEphidr*cos(phi)*cos(phi) - (Er*cos(phi)*sin(phi) - Ephi*sin(phi)*sin(phi))/r);
					dEidxj[1][1] = Escale*(dErdr*sin(phi)*sin(phi) + dEphidr*cos(phi)*sin(phi) + (Er*cos(phi)*cos(phi) - Ephi*cos(phi)*sin(phi))/r);
				}
				CylToCart(dErdz,dEphidz,phi,dExdz,dEydz); // convert z-derivatives of r,phi components into z-derivatives of x,y components
				dEidxj[0][2] = dExdz*Escale;
				dEidxj[1][2] = dEydz*Escale;
				dEidxj[2][0] = dEzdr*cos(phi)*Escale;
				dEidxj[2][1] = dEzdr*sin(phi)*Escale;
				dEidxj[2][2] = dEzdz*Escale;
			}
		}
		else if (fVc){
			double Vloc, dVdrj[3], dummy;
			// bicubic interpolation
			alglib::spline2ddiff(Vc, r, z, Vloc, dVdrj[0], dVdrj[2], dummy);
			double phi = atan2(y,x);
			V = Vloc*Escale;
			Ei[0] = -dVdrj[0]*cos(phi)*Escale;
			Ei[1] = -dVdrj[0]*sin(phi)*Escale;
			Ei[2] = -dVdrj[2]*Escale;
			// !!!!!!!!!!! calculation of electric field derivatives not yet implemented for potential interpolation !!!!!!!!!!!!!!!!
		}

	}
}
