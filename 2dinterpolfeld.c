// this file contains functions to read in E and B fields, interpolate
// them and calculate derivatives

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "nr/nr3.h"
#include "nr/interp_1d.h"

#include "2dinterpolfeld.h"
#include "globals.h"

using namespace std;

// coefficients for bicubic interpolation
void bcucof(VecDoub_I &y, VecDoub_I &y1, VecDoub_I &y2, VecDoub_I &y12,
	const Doub d1, const Doub d2, Doub c[4][4]) {
	static Int wt_d[16*16]=
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
		2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
		0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
		-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
		9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
		-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
		2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
		-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
		4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
	Int l,k,j,i;
	Doub xx,d1d2=d1*d2;
	VecDoub cl(16),x(16);
	static MatInt wt(16,16,wt_d);
	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++) {
		xx=0.0;
		for (k=0;k<16;k++) xx += wt[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) c[i][j]=cl[l++];
}

// bicubic interpolation with predetermined coefficients
void bcuint_new(Doub c[4][4],
	const Doub x1l, const Doub x1u, const Doub x2l, const Doub x2u,
	const Doub x1, const Doub x2, Doub &ansy, Doub &ansy1, Doub &ansy2) {
	Int i;
	Doub t,u,d1=x1u-x1l,d2=x2u-x2l;
	if (x1u == x1l || x2u == x2l)
		throw("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}
	ansy1 /= d1;
	ansy2 /= d2;
}


void CylToCart(long double v_r, long double v_phi, long double phi, long double &v_x, long double &v_y){
	v_x = v_r*cos(phi) - v_phi*sin(phi);
	v_y = v_r*sin(phi) + v_phi*cos(phi);	
}

long double lengthconv = 0.01 , Bconv = 1e-4, Econv = 1e2;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 

// constructor
TabField::TabField(const char *tabfile, long double Bscale, long double Escale,
		long double aNullFieldTime, long double aRampUpTime, long double aFullFieldTime, long double aRampDownTime){
	NullFieldTime = aNullFieldTime;
	RampUpTime = aRampUpTime;
	FullFieldTime = aFullFieldTime;
	RampDownTime = aRampDownTime;

	MatDoub BrTab, BphiTab, BzTab;	// Br/Bz values
	MatDoub ErTab , EphiTab, EzTab;	// Er/Ez values
	MatDoub VTab; // potential values
	ReadTabFile(tabfile,Bscale,Escale,BrTab,BphiTab,BzTab,ErTab,EphiTab,EzTab,VTab); // open tabfile and read values into arrays
	
	CheckTab(BrTab,BphiTab,BzTab,ErTab,EphiTab,EzTab,VTab); // print some info

	printf("Starting Preinterpolation ... ");
	if (BrTab.nrows() > 0){
		printf("Br ... ");
		fflush(stdout);
		PreInterpol(Brc,BrTab); // precalculate interpolation coefficients for B field
	}
	if (BphiTab.nrows() > 0){
		printf("Bphi ... ");
		fflush(stdout);
		PreInterpol(Bphic,BphiTab);
	}
	if (BzTab.nrows() > 0){
		printf("Bz ... ");
		fflush(stdout);
		PreInterpol(Bzc,BzTab);
	}
	if (VTab.nrows() > 0){
		printf("V ... ");
		fflush(stdout);
		PreInterpol(Vc,VTab);
		ErTab.resize(0,0); // if there is a potential, we don't need the E-vector
		EphiTab.resize(0,0);
		EzTab.resize(0,0);
	}
	if (ErTab.nrows() > 0){
		printf("Er ... ");
		fflush(stdout);
		PreInterpol(Erc,ErTab); // precalculate interpolation coefficients for E field
	}
	if (EphiTab.nrows() > 0){
		printf("Ephi ... ");
		fflush(stdout);
		PreInterpol(Ephic,EphiTab);
	}
	if (EzTab.nrows() > 0){
		printf("Ez ... ");
		fflush(stdout);
		PreInterpol(Ezc,EzTab);
	}
	printf("Done (%.f MB)\n",float(Brc.nrows()*Brc.ncols() + Bphic.nrows()*Bphic.ncols() + Bzc.nrows()*Bzc.ncols()
								 + Erc.nrows()*Erc.ncols() + Ezc.nrows()*Ezc.ncols() + Vc.nrows()*Vc.ncols())
								*16*sizeof(Doub)/1024/1024);
}


// read table file given in the constructor parameter
void TabField::ReadTabFile(const char *tabfile, long double Bscale, long double Escale,
		MatDoub_O &BrTab, MatDoub_O &BphiTab, MatDoub_O &BzTab, MatDoub_O &ErTab, MatDoub_O &EphiTab, MatDoub_O &EzTab, MatDoub_O &VTab){
	ifstream FIN(tabfile, ifstream::in);
	if (!FIN.is_open()){
		printf("\nCould not open %s!\n",tabfile);
		exit(-1);
	}
	printf("\nReading %s ",tabfile);
	int intval;
	string line;
	FIN >> m >> intval >> n;
	
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	bool skipy = true;
	if (line.substr(0,12) == " 2 Y [LENGU]")  skipy = false;	
	getline(FIN,line);
	if (!skipy) getline(FIN,line);

	if (line.find("RBX") != string::npos){
		BrTab.resize(m,n);
		getline(FIN,line);
	}
	if (line.find("RBY") != string::npos){
		BphiTab.resize(m,n);
		getline(FIN,line);
	}
	if (line.find("RBZ") != string::npos){
		BzTab.resize(m,n);
		getline(FIN,line);
	}

	if (line.find("EX") != string::npos){
		ErTab.resize(m,n);
		getline(FIN,line);
	}
	if (line.find("EY") != string::npos){
		EphiTab.resize(m,n);
		getline(FIN,line);
	}
	if (line.find("EZ") != string::npos){
		EzTab.resize(m,n);
		getline(FIN,line);
	}
	
	if (line.find("RV") != string::npos){	// file contains potential?
		VTab.resize(m,n);
		getline(FIN,line);
	}

	if (!FIN || line.substr(0,2) != " 0"){
		printf("%s not found or corrupt! Exiting...\n",tabfile);	
		exit(-1);
	}
	
	VecDoub rind(m), zind(n);
	int ri = 0,zi = -1, perc = 0;
	long double r, z, val;
	while (FIN.good()){
		FIN >> r;
		if (!skipy) FIN >> val;
		FIN >> z;
		if (!FIN) break;
		r *= lengthconv;
		z *= lengthconv;
		if (zi > 0 && z < zind[zi-1]){
			ri++;
			zi = 0;
		}
		else zi++;
		
		// status if read is displayed
		percent(zi*ri ,0, n*m, perc);
		
		rind[ri] = r;
		zind[zi] = z;
		if (BrTab.nrows() > 0){
			FIN >> val;
			BrTab[ri][zi] = val*Bconv*Bscale;
		}
		if (BphiTab.nrows() > 0){
			FIN >> val;
			BphiTab[ri][zi] = val*Bconv*Bscale;
		}
		if (BzTab.nrows() > 0){
			FIN >> val;
			BzTab[ri][zi] = val*Bconv*Bscale;
		}
		if (ErTab.nrows() > 0){
			FIN >> val;
			ErTab[ri][zi] = val*Econv*Escale;
		}
		if (EphiTab.nrows() > 0){
			FIN >> val;
			EphiTab[ri][zi] = val*Econv*Escale;
		}
		if (EzTab.nrows() > 0){
			FIN >> val;
			EzTab[ri][zi] = val*Econv*Escale;
		}
		if (VTab.nrows() > 0){
			FIN >> val;
			VTab[ri][zi] = val*Escale;
		}
		FIN >> ws;
	}
	
	if (Bscale == 0){
		BrTab.resize(0,0);
		BphiTab.resize(0,0);
		BzTab.resize(0,0);
	}
	if (Escale == 0){
		ErTab.resize(0,0);
		EphiTab.resize(0,0);
		EzTab.resize(0,0);
		VTab.resize(0,0);
	}

	r_mi = rind[0];
	z_mi = zind[0];
	rdist = rind[1] - rind[0];
	zdist = rind[1] - rind[0];
	
	printf("\n");
	if (ri+1 != m || zi+1 != n){
		printf("The header says the size is %i by %i, actually it is %i by %i! Exiting...\n", m, n, ri+1, zi+1);
		exit(-1);
	}
	FIN.close();
}


// print info about tablefile, calculate min/max fields
void TabField::CheckTab(MatDoub_I &BrTab, MatDoub_I &BphiTab, MatDoub_I &BzTab, MatDoub_I &ErTab, MatDoub_I &EphiTab, MatDoub_I &EzTab, MatDoub_I &VTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB	
	printf("The arrays are %d by %d.\n",m,n);
	printf("The r values go from %LG to %LG\n",r_mi,r_mi + rdist*(m-1));
	printf("The z values go from %LG to %LG.\n",z_mi,z_mi + rdist*(n-1));
	printf("rdist = %LG zdist = %LG\n",rdist,zdist);
		
	long double Babsmax = 0, Babsmin = 9e99, Babs;
	long double Vmax = 0, Vmin = 9e99;
	for (int j=0; j < m; j++)
	{
		for (int k=0; k < n; k++)
		{
			Babs = 0;
			if (BrTab.nrows() > 0) Babs += BrTab[j][k]*BrTab[j][k];
			if (BphiTab.nrows() > 0) Babs += BphiTab[j][k]*BphiTab[j][k];
			if (BzTab.nrows() > 0) Babs += BzTab[j][k]*BzTab[j][k];
			Babsmax = max(sqrt(Babs),Babsmax);
			Babsmin = min(sqrt(Babs),Babsmin);
			if (VTab.nrows() > 0)
			{
				Vmax = max(VTab[j][k],Vmax);
				Vmin = min(VTab[j][k],Vmin);
			}
							
		}
	}
	
	printf("The input table file has values of |B| from %LG T to %LG T and values of V from %LG V to %LG V\n",Babsmin,Babsmax,Vmin,Vmax);
}


// calculate derivatives of Tab and write them to Tab1/2/12
void TabField::CalcDerivs(MatDoub_I &Tab, MatDoub_O &Tab1, MatDoub_O &Tab2, MatDoub_O &Tab12)
{
	VecDoub x(m), y(m);
	for (int zi = 0; zi < n; zi++){
		for (int ri = 0; ri < m; ri++){
			x[ri] = r_mi + ri*rdist;
			y[ri] = Tab[ri][zi];
		}
		Spline_interp spli(x,y); // splineinterpolate field values F in r-direction
		for (int ri = 0; ri < m-1; ri++)
			Tab1[ri][zi] = (y[ri+1] - y[ri])/rdist - (spli.y2[ri+1] - spli.y2[ri])*rdist/6 - spli.y2[ri]*rdist/2; // get dF/dr from spline coefficients spli.y2
		Tab1[m-1][zi] = (y[m-1] - y[m-2])/rdist - (spli.y2[m-1] - spli.y2[m-2])*rdist/6 + spli.y2[m-1]*rdist/2;
	}
	x.resize(n);
	y.resize(n);
	for (int ri = 0; ri < m; ri++){
		for (int zi = 0; zi < n; zi++){
			x[zi] = z_mi + zi*zdist;
			y[zi] = Tab[ri][zi];
		}
		Spline_interp spli(x,y); // splineinterpolate field values F in z-diretion
		for (int zi = 0; zi < n-1; zi++)
			Tab2[ri][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get dF/dz from spline coefficients spli.y2
		Tab2[ri][n-1] = (y[n-1] - y[n-2])/zdist - (spli.y2[n-1] - spli.y2[n-2])*zdist/6 + spli.y2[n-1]*zdist/2;
	}
	for (int ri = 0; ri < m; ri++){
		for (int zi = 0; zi < n; zi++){
			x[zi] = z_mi + zi*zdist;
			y[zi] = Tab1[ri][zi];
		}
		Spline_interp spli(x,y); // splineinterpolate dF/dr in z direction
		for (int zi = 0; zi < n-1; zi++)
			Tab12[ri][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get cross derivatives d2F/drdz from spline coefficients
		Tab12[ri][n-1] = (y[n-1] - y[n-2])/zdist - (spli.y2[n-1] - spli.y2[n-2])*zdist/6 + spli.y2[n-1]*zdist/2;
	}
}


// calculate interpolation coefficients to speed up interpolation
void TabField::PreInterpol(NRmatrix<Doub[4][4]> &coeff, MatDoub_I &Tab){
	MatDoub Tab1(m,n), Tab2(m,n), Tab12(m,n);	// dBi/dr, dBi/dz, d2Bi/drdz
	CalcDerivs(Tab,Tab1,Tab2,Tab12);
	
	// allocating space for the preinterpolation
	// The B*c are 4D arrays with m x n x 4 x 4 fields
	coeff.resize(m,n);

	int indr, indz;
	VecDoub yyy(4), yyy1(4), yyy2(4), yyy12(4);        //rectangle with values at the 4 corners for interpolation
	
	for (indr=0; indr < m-1; indr++)
	{
		for (indz=0; indz < n-1; indz++)
		{			
			// fill rectancle with values and derivatives
			yyy[0] = Tab[indr][indz];
			yyy[1] = Tab[indr+1][indz];
			yyy[2] = Tab[indr+1][indz+1];
			yyy[3] = Tab[indr][indz+1];
			
			yyy1[0] = Tab1[indr][indz];
			yyy1[1] = Tab1[indr+1][indz];
			yyy1[2] = Tab1[indr+1][indz+1];
			yyy1[3] = Tab1[indr][indz+1];
			
			yyy2[0] = Tab2[indr][indz];
			yyy2[1] = Tab2[indr+1][indz];
			yyy2[2] = Tab2[indr+1][indz+1];
			yyy2[3] = Tab2[indr][indz+1];
			
			yyy12[0] = Tab12[indr][indz];
			yyy12[1] = Tab12[indr+1][indz];
			yyy12[2] = Tab12[indr+1][indz+1];
			yyy12[3] = Tab12[indr][indz+1];
			
															
			// determine koefficients of interpolation
			bcucof(yyy,yyy1,yyy2,yyy12,rdist,zdist,coeff[indr][indz]);
									
		}
	}
}


// interpolate B-field (return true when calculation was successful)
bool TabField::BInterpol(long double t, long double x, long double y, long double z, long double B[4][4]){
	long double r = sqrt(x*x+y*y);
	long double Bscale = BFieldScale(t);
	if (Bscale != 0 && (r - r_mi)/rdist > 0 && (r - r_mi - (m-1)*rdist)/rdist < 0 && (z - z_mi)/zdist > 0 && (z - z_mi - (n-1)*zdist)/zdist < 0){
		// bicubic interpolation
		int indr = (int)((r - r_mi)/rdist);
		int indz = (int)((z - z_mi)/zdist);
		long double rl = r_mi + indr*rdist;
		long double zl = z_mi + indz*zdist;
		long double Br = 0, dBrdr = 0, dBrdz = 0, Bphi = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0;
		long double phi = atan2(y,x);
		if (Brc.nrows() > 0)   bcuint_new(Brc[indr][indz], 	 rl, rl+rdist, zl, zl+zdist, r, z, Br, dBrdr, dBrdz);
		if (Bphic.nrows() > 0) bcuint_new(Bphic[indr][indz], rl, rl+rdist, zl, zl+zdist, r, z, Bphi, dBphidr, dBphidz);
		CylToCart(Br,Bphi,phi,B[0][0],B[1][0]);
		if (r > 0){
			B[0][1] = dBrdr*cos(phi)*cos(phi) - dBphidr*cos(phi)*sin(phi) + (Br*sin(phi)*sin(phi) - Bphi*cos(phi)*sin(phi))/r;
			B[1][1] = dBrdr*cos(phi)*sin(phi) - dBphidr*cos(phi)*cos(phi) - (Br*cos(phi)*sin(phi) - Bphi*sin(phi)*sin(phi))/r;
			B[0][2] = dBrdr*cos(phi)*sin(phi) + dBphidr*sin(phi)*sin(phi) - (Br*cos(phi)*sin(phi) + Bphi*cos(phi)*cos(phi))/r;
			B[1][2] = dBrdr*sin(phi)*sin(phi) + dBphidr*cos(phi)*sin(phi) + (Br*cos(phi)*cos(phi) - Bphi*cos(phi)*sin(phi))/r;
		}
		CylToCart(dBrdz,dBphidz,phi,B[0][3],B[1][3]);
		if (Bzc.nrows() > 0)   bcuint_new(Bzc[indr][indz],   rl, rl+rdist, zl, zl+zdist, r, z, B[2][0], dBzdr, B[2][3]);
		B[2][1] = dBzdr*cos(phi);
		B[2][2] = dBzdr*sin(phi);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				B[i][j] *= Bscale;
		return (Brc.nrows() > 0 || Bphic.nrows() > 0 || Bzc.nrows() > 0);
	}
	return false;
}

//B-Feld nach Wunsch skalieren, BFeldSkal wird an alle B Komp und deren Ableitungen multipliziert
long double TabField::BFieldScale(long double t){
	if (t < NullFieldTime || t >= NullFieldTime + RampUpTime + FullFieldTime + RampDownTime)
		return 0;
	else if (t >= NullFieldTime && t < NullFieldTime + RampUpTime){
		// ramping up field smoothly with cosine
		//result = (0.5 - 0.5*cos(pi*(t - CleaningTime - FillingTime)/RampUpTime)) * BFeldSkalGlobal;

		// linear ramp
		return (t - NullFieldTime)/RampUpTime;
	}
	else if (t >= NullFieldTime + RampUpTime && t < NullFieldTime + RampUpTime + FullFieldTime)
		return 1;
	else if (t >= NullFieldTime + RampUpTime + FullFieldTime && t < NullFieldTime + RampUpTime + FullFieldTime + RampDownTime){
		// ramping down field smoothly with cosine
		//result = (0.5 + 0.5*cos(pi*(t - (RampUpTime + CleaningTime + FillingTime + FullFieldTime)) / RampDownTime)) * BFeldSkalGlobal;

		// linear ramp
		return (1 - (t - RampUpTime - NullFieldTime - FullFieldTime)/RampDownTime);
	}
	else return 0;
};


// interpolate E-field (return true when calculation was successful)
bool TabField::EInterpol(long double x, long double y, long double z, long double &V, long double Ei[3]){
	long double dummy, dVdrj[3];
	if (Vc.nrows() > 0){ // prefer E-field from potential over pure E-field interpolation
		long double r = sqrt(x*x+y*y);
		if ((r - r_mi)/rdist > 0 && (r - r_mi - (m-1)*rdist)/rdist < 0 && (z - z_mi)/zdist > 0 && (z - z_mi - (n-1)*zdist)/zdist < 0){
			// bicubic interpolation
			int indr = (int)((r - r_mi)/rdist);
			int indz = (int)((z - z_mi)/zdist);
			long double rl = r_mi + indr*rdist;
			long double zl = z_mi + indz*zdist;
			long double phi = atan2(y,x);
			bcuint_new(Vc[indr][indz], rl, rl+rdist, zl, zl+zdist, r, z, V, dVdrj[0], dVdrj[2]);
			CylToCart(-dVdrj[0],0,phi,Ei[0],Ei[1]);
			Ei[2] = -dVdrj[2];
			return true;
		}
	}
	else if (Erc.nrows() > 0 || Ephic.nrows() > 0 || Ezc.nrows() > 0){
		long double r = sqrt(x*x+y*y);
		if ((r - r_mi)/rdist > 0 && (r - r_mi - (m-1)*rdist)/rdist < 0 && (z - z_mi)/zdist > 0 && (z - z_mi - (n-1)*zdist)/zdist < 0){
			// bicubic interpolation
			int indr = (int)((r - r_mi)/rdist);
			int indz = (int)((z - z_mi)/zdist);
			long double rl = r_mi + indr*rdist;
			long double zl = z_mi + indz*zdist;
			long double Er = 0, Ephi = 0;
			long double phi = atan2(y,x);
			if (Erc.nrows() > 0)   bcuint_new(Erc[indr][indz],   rl, rl+rdist, zl, zl+zdist, r, z, Er, dummy, dummy);
			if (Ephic.nrows() > 0) bcuint_new(Ephic[indr][indz], rl, rl+rdist, zl, zl+zdist, r, z, Ephi, dummy, dummy);
			if (Ezc.nrows() > 0)   bcuint_new(Ezc[indr][indz],   rl, rl+rdist, zl, zl+zdist, r, z, Ei[2], dummy, dummy);
			CylToCart(Er,Ephi,phi,Ei[0],Ei[1]);
			return true;
		}
	}
	return false;
}



