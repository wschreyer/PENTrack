/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */


#include "field_3d.h"

#include <fstream>
#include <iostream>

#include "interpolation.h"

#include "tricubic.h"
#include "globals.h"


/**
 * calculate consecutive index of entry with x-index xi, y-index yi and z-index zi
 *
 * @param xi x index
 * @param yi y index
 * @param zi z index
 * @param xsize size of array in x dimension
 * @param ysize size of array in y dimension
 * @param zsize size of array in z dimension
 */
int INDEX_3D(int xi, int yi, int zi, int xsize, int ysize, int zsize){
	return zi*xsize*ysize + yi*xsize + xi;
}


void TabField3::ReadTabFile(const std::string &tabfile,
		std::vector<double> &BxTab, std::vector<double> &ByTab, std::vector<double> &BzTab, std::vector<double> &VTab){
	std::ifstream FIN(tabfile, std::ifstream::in);
	if (!FIN.is_open()){
		std::cout << "\nCould not open " << tabfile << "!\n";
		exit(-1);
	}
	std::cout << "\nReading " << tabfile << " ";
	std::string line;
	FIN >> xl >> yl >> zl;

	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);

	if (line.find("BX") != std::string::npos){
		BxTab.resize(xl*yl*zl);
		getline(FIN,line);
	}
	if (line.find("BY") != std::string::npos){
		ByTab.resize(xl*yl*zl);
		getline(FIN,line);
	}
	if (line.find("BZ") != std::string::npos){
		BzTab.resize(xl*yl*zl);
		getline(FIN,line);
	}

	if (line.find("V") != std::string::npos){	// file contains potential?
		VTab.resize(xl*yl*zl);
		getline(FIN,line);
	}

	if (!FIN || line.substr(0,2) != " 0"){
		std::cout << tabfile << " not found or corrupt! Exiting...\n";
		exit(-1);
	}

	std::vector<double> xind(xl), yind(yl), zind(zl);
	int xi = 0, yi = 0, zi = -1, perc = 0;
	double x, y, z, val;
	while (FIN.good()){
		FIN >> x;
		FIN >> y;
		FIN >> z;
		if (!FIN) break;
		x *= lengthconv;
		y *= lengthconv;
		z *= lengthconv;

		if (zi >= 0 && z < zind[zi]){
			if (yi >= 0 && y < yind[yi]){
				xi++;
				yi = 0;
			}
			else yi++;
			zi = 0;
		}
		else zi++;

		int i3 = INDEX_3D(xi, yi, zi, xl, yl, zl);

		// status if read is displayed
		PrintPercent((float)i3/(xl*yl*zl), perc);

		xind[xi] = x;
		yind[yi] = y;
		zind[zi] = z;
		if (BxTab.size() > 0){
			FIN >> val;
			BxTab[i3] = val*Bconv;
		}
		if (ByTab.size() > 0){
			FIN >> val;
			ByTab[i3] = val*Bconv;
		}
		if (BzTab.size() > 0){
			FIN >> val;
			BzTab[i3] = val*Bconv;
		}
		if (VTab.size() > 0){
			FIN >> val;
			VTab[i3] = val;
		}
		FIN >> std::ws;

	}

	std::cout << "\n";
	if (xi+1 != xl || yi + 1 != yl || zi+1 != zl){
		std::cout << "The header says the size is " << xl << " by " << yl << " by " << zl << ", actually it is " << xi+1 << " by " << yi+1 << " by " << zi+1 << "! Exiting...\n";
		exit(-1);
	}
	FIN.close();

	xdist = xind[1] - xind[0];
	ydist = yind[1] - yind[0];
	zdist = zind[1] - zind[0];
	x_mi = xind[0];
	y_mi = yind[0];
	z_mi = zind[0];

};


void TabField3::CheckTab(const std::vector<double> &BxTab, const std::vector<double> &ByTab,
		const std::vector<double> &BzTab, const std::vector<double> &VTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
	std::cout << "The arrays are " << xl << " by " << yl << " by " << zl << ".\n";
	std::cout << "The x values go from " << x_mi << " to " << x_mi + xdist*(xl-1) << "\n";
	std::cout << "The y values go from " << y_mi << " to " << y_mi + ydist*(yl-1) << "\n";
	std::cout << "The z values go from " << z_mi << " to " << z_mi + zdist*(zl-1) << ".\n";
	std::cout << "xdist = " << xdist << " ydist = " << ydist << " zdist = " << zdist << "\n";

	double Babsmax = 0, Babsmin = 9e99, Babs;
	double Vmax = 0, Vmin = 9e99;
	for (int j=0; j < xl; j++)
	{
		for (int k=0; k < yl; k++)
		{
			for (int l=0; l < zl; l++)
			{
				Babs = 0;
				int i3 = INDEX_3D(j, k, l, xl, yl, zl);
				if (BxTab.size() > 0)
					Babs += BxTab[i3] * BxTab[i3];
				if (ByTab.size() > 0)
					Babs += ByTab[i3] * ByTab[i3];
				if (BzTab.size() > 0)
					Babs += BzTab[i3] * BzTab[i3];
				Babsmax = std::max(sqrt(Babs),Babsmax);
				Babsmin = std::min(sqrt(Babs),Babsmin);
				if (VTab.size() > 0)
				{
					Vmax = std::max(VTab[i3],Vmax);
					Vmin = std::min(VTab[i3],Vmin);
				}

			}
		}
	}

	std::cout << "The input table file has values of |B| from " << Babsmin << " T to " << Babsmax << " T and values of V from " << Vmin << " V to " << Vmax << " V\n";
};


void TabField3::CalcDerivs(const std::vector<double> &Tab, std::vector<double> &Tab1, std::vector<double> &Tab2, std::vector<double> &Tab3,
		std::vector<double> &Tab12, std::vector<double> &Tab13, std::vector<double> &Tab23, std::vector<double> &Tab123) const
{
	alglib::real_1d_array x, y, diff;

	x.setlength(xl);
	y.setlength(xl);
	diff.setlength(xl);
	for (int zi = 0; zi < zl; zi++){
		for (int yi = 0; yi < yl; yi++){
			for (int xi = 0; xi < xl; xi++){
				x[xi] = x_mi + xi*xdist;
				y[xi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)]; // get derivatives dF/dx from spline interpolation
			}
			alglib::spline1dgriddiffcubic(x, y, diff);
			for (int xi = 0; xi < xl; xi++)
				Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[xi];
		}
	}

	x.setlength(yl);
	y.setlength(yl);
	diff.setlength(yl);
	for (int xi = 0; xi < xl; xi++){
		for (int zi = 0; zi < zl; zi++){
			for (int yi = 0; yi < yl; yi++){
				x[yi] = y_mi + yi*ydist;
				y[yi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)]; // get derivatives dF/dy from spline interpolation
			}
			alglib::spline1dgriddiffcubic(x, y, diff);
			for (int yi = 0; yi < yl; yi++)
				Tab2[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[yi];
		}
	}

	x.setlength(zl);
	y.setlength(zl);
	diff.setlength(zl);
	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get derivatives dF/dz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab3[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	x.setlength(yl);
	y.setlength(yl);
	diff.setlength(yl);
	for (int xi = 0; xi < xl; xi++){
		for (int zi = 0; zi < zl; zi++){
			for (int yi = 0; yi < yl; yi++){
				x[yi] = y_mi + yi*ydist;
				y[yi] = Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dxdy from spline interpolation
			for (int yi = 0; yi < yl; yi++)
				Tab12[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[yi];
		}
	}

	x.setlength(zl);
	y.setlength(zl);
	diff.setlength(zl);
	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dxdz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab13[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab2[INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dydz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab23[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab12[INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d3F/dxdydz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab123[INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

};


void TabField3::PreInterpol(std::vector<std::vector<double> > &coeff, const std::vector<double> &Tab) const{
	std::vector<double> Tab1(xl*yl*zl), Tab2(xl*yl*zl), Tab3(xl*yl*zl);
	std::vector<double> Tab12(xl*yl*zl), Tab13(xl*yl*zl), Tab23(xl*yl*zl), Tab123(xl*yl*zl);
	CalcDerivs(Tab,Tab1,Tab2,Tab3,Tab12,Tab13,Tab23,Tab123);

	// allocating space for the preinterpolation
	// The B*c are 5D arrays with xl x yl x zl x 4 x 4 fields
	coeff.resize((xl-1)*(yl-1)*(zl-1));

	int indx, indy, indz;
	double yyy[8], yyy1[8], yyy2[8], yyy3[8], yyy12[8], yyy13[8], yyy23[8], yyy123[8]; //rectangle with values at the 4 corners for interpolation

	for (indx=0; indx < xl-1; indx++){
		for (indy=0; indy < yl-1; indy++){
			for (indz=0; indz < zl-1; indz++){
				// fill cube with values and derivatives
				// order: see Lekien, Marsden: "Tricubic interpolation in three dimensions"
				int i3[8] ={INDEX_3D(indx  , indy  , indz  , xl, yl, zl),
							INDEX_3D(indx+1, indy  , indz  , xl, yl, zl),
							INDEX_3D(indx  , indy+1, indz  , xl, yl, zl),
							INDEX_3D(indx+1, indy+1, indz  , xl, yl, zl),
							INDEX_3D(indx  , indy  , indz+1, xl, yl, zl),
							INDEX_3D(indx+1, indy  , indz+1, xl, yl, zl),
							INDEX_3D(indx  , indy+1, indz+1, xl, yl, zl),
							INDEX_3D(indx+1, indy+1, indz+1, xl, yl, zl)};

				for (int i = 0; i < 8; i++){
					yyy[i]    = Tab[i3[i]];
					yyy1[i]   = Tab1[i3[i]]*xdist;
					yyy2[i]   = Tab2[i3[i]]*ydist;
					yyy3[i]   = Tab3[i3[i]]*zdist;
					yyy12[i]  = Tab12[i3[i]]*xdist*ydist;
					yyy13[i]  = Tab13[i3[i]]*xdist*zdist;
					yyy23[i]  = Tab23[i3[i]]*ydist*zdist;
					yyy123[i] = Tab123[i3[i]]*xdist*ydist*zdist;
				}

				// determine coefficients of interpolation
				int ind3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
				coeff[ind3].resize(64);
				tricubic_get_coeff(&coeff[ind3][0], yyy, yyy1, yyy2, yyy3, yyy12, yyy13, yyy23, yyy123);
			}
		}
	}
};


TabField3::TabField3(const std::string &tabfile, const std::string &Bscale, const std::string &Escale,
		const double aBoundaryWidth, const double alengthconv, const double aBconv)
	: TField(Bscale, Escale){
	BoundaryWidth = aBoundaryWidth;
	lengthconv = alengthconv;
	Bconv = aBconv;


	std::vector<double> BxTab, ByTab, BzTab;	// Bx/By/Bz values
	std::vector<double> VTab; // potential values
	ReadTabFile(tabfile,BxTab,ByTab,BzTab,VTab); // open tabfile and read values into arrays

	CheckTab(BxTab,ByTab,BzTab,VTab); // print some info

	std::cout << "Starting Preinterpolation ... ";
	float size = 0;
	if (BxTab.size() > 0){
		std::cout << "Bx ... ";
		std::cout.flush();
		PreInterpol(Bxc,BxTab); // precalculate interpolation coefficients for B field
		size += float(Bxc.size()*64*sizeof(double)/1024/1024);
		BxTab.clear();
	}
	if (ByTab.size() > 0){
		std::cout << "By ... ";
		std::cout.flush();
		PreInterpol(Byc,ByTab);
		size += float(Byc.size()*64*sizeof(double)/1024/1024);
		ByTab.clear();
	}
	if (BzTab.size() > 0){
		std::cout << "Bz ... ";
		std::cout.flush();
		PreInterpol(Bzc,BzTab);
		size += float(Bzc.size()*64*sizeof(double)/1024/1024);
		BzTab.clear();
	}
	if (VTab.size() > 0){
		std::cout << "V ... ";
		std::cout.flush();
		PreInterpol(Vc,VTab);
		size += float(Vc.size()*64*sizeof(double)/1024/1024);
		VTab.clear();
	}
	std::cout << "Done (" << size << " MB)\n";
}


void TabField3::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	double Bscale = BScaling(t);
	// get coordinate index
	int indx = (int)floor((x - x_mi)/xdist);
	int indy = (int)floor((y - y_mi)/ydist);
	int indz = (int)floor((z - z_mi)/zdist);
	if (Bscale != 0 && indx >= 0 && indx < xl - 1 && indy >= 0 && indy < yl - 1 && indz >= 0 && indz < zl - 1){
		// scale coordinates to unit cube
		double xu = (x - x_mi - indx*xdist)/xdist;
		double yu = (y - y_mi - indy*ydist)/ydist;
		double zu = (z - z_mi - indz*zdist)/zdist;
		int i3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
		// tricubic interpolation
		if (Bxc.size() > 0){
			double *coeff = const_cast<double*>(&Bxc[i3][0]);
			B[0] = Bscale*tricubic_eval(coeff, xu, yu, zu);
			if (dBidxj != NULL){
				dBidxj[0][0] = Bscale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 0)/xdist;
				dBidxj[0][1] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 0)/ydist;
				dBidxj[0][2] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 0, 1)/zdist;
			}
		}
		if (Byc.size() > 0){
			double *coeff = const_cast<double*>(&Byc[i3][0]);
			B[1] = Bscale*tricubic_eval(coeff, xu, yu, zu);
			if (dBidxj != NULL){
				dBidxj[1][0] = Bscale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 0)/xdist;
				dBidxj[1][1] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 0)/ydist;
				dBidxj[1][2] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 0, 1)/zdist;
			}
		}
		if (Bzc.size() > 0){
			double *coeff = const_cast<double*>(&Bzc[i3][0]);
			B[2] = Bscale*tricubic_eval(coeff, xu, yu, zu);
			if (dBidxj != NULL){
				dBidxj[2][0] = Bscale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 0)/xdist;
				dBidxj[2][1] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 0)/ydist;
				dBidxj[2][2] = Bscale*tricubic_eval(coeff, xu, yu, zu, 0, 0, 1)/zdist;
			}
		}

		for (int i = 0; i < 3; i++){
			if (dBidxj != NULL)
				FieldSmthr(x, y, z, B[i], dBidxj[i]); // apply field smoothing to each component
			else
				FieldSmthr(x, y, z, B[i], NULL);
		}

	}
}

void TabField3::FieldSmthr(const double x, const double y, const double z, double &F, double dFdxi[3]) const{
	if (BoundaryWidth != 0 && F != 0){ // skip, if BoundaryWidth is set to zero
		double dxlo = (x - x_mi)/BoundaryWidth; // calculate distance to edges in units of BoundaryWidth
		double dxhi = (x_mi + xdist*(xl - 1) - x)/BoundaryWidth;
		double dylo = (y - y_mi)/BoundaryWidth;
		double dyhi = (y_mi + ydist*(yl - 1) - y)/BoundaryWidth;
		double dzlo = (z - z_mi)/BoundaryWidth;
		double dzhi = (z_mi + zdist*(zl - 1) - z)/BoundaryWidth;

		// F'(x,y,z) = F(x,y,z)*f(x)*f(y)*f(z)
		// dF'/dx = dF/dx*f(x)*f(y)*f(z) + F*df(x)/dx*f(y)*f(z) --> similar for dF'/dy and dF'/dz

		double Fscale = 1; // f(x)*f(y)*f(z)
		double dFadd[3] = {0, 0, 0}; // F*df(x_i)/dx_i / f(x_i)

		if (dxhi < 1) { // if point in upper x boundary (distance smaller than BoundaryWidth)
			Fscale *= SmthrStp(dxhi); // scale field by value of smoother function
			if (dFdxi != NULL)
				dFadd[0] = -F*SmthrStpDer(dxhi)/BoundaryWidth/SmthrStp(dxhi); // add derivative of smoother function according to product rule
		}
		if (dyhi < 1){ // if point in upper y boundary
			Fscale *= SmthrStp(dyhi);
			if (dFdxi != NULL)
				dFadd[1] = -F*SmthrStpDer(dyhi)/BoundaryWidth/SmthrStp(dyhi);
		}
		if (dzhi < 1){ // if point in upper z boundary
			Fscale *= SmthrStp(dzhi);
			if (dFdxi != NULL)
				dFadd[2] = -F*SmthrStpDer(dzhi)/BoundaryWidth/SmthrStp(dzhi);
		}

		if (dxlo < 1){ // if point in lower x boundary
			Fscale *= SmthrStp(dxlo);
			if (dFdxi != NULL)
				dFadd[0] = F*SmthrStpDer(dxlo)/BoundaryWidth/SmthrStp(dxlo);
		}
		if (dylo < 1){ // if point in lower y boundary
			Fscale *= SmthrStp(dylo);
			if (dFdxi != NULL)
				dFadd[1] = F*SmthrStpDer(dylo)/BoundaryWidth/SmthrStp(dylo);
		}
		if (dzlo < 1){ // if point in lower z boundary
			Fscale *= SmthrStp(dzlo);
			if (dFdxi != NULL)
				dFadd[2] = F*SmthrStpDer(dzlo)/BoundaryWidth/SmthrStp(dzlo);
		}

		if (Fscale != 1){
			F *= Fscale; // scale field value
			if (dFdxi != NULL){
				for (int i = 0; i < 3; i++){
					dFdxi[i] = dFdxi[i]*Fscale + dFadd[i]*Fscale; // scale derivatives according to product rule
				}
			}
		}

	}
}

double TabField3::SmthrStp(const double x) const{
	return 6*pow(x, 5) - 15*pow(x, 4) + 10*pow(x, 3);
}

double TabField3::SmthrStpDer(const double x) const{
	return 30*pow(x, 4) - 60*pow(x, 3) + 30*pow(x,2);
}

void TabField3::EField(const double x, const double y, const double z, const double t,
		double &V, double Ei[3], double dEidxj[3][3]) const{
	double Escale = EScaling(t);
	if (Escale != 0 && Vc.size() > 0 &&
		(x - x_mi)/xdist > 0 && (x - x_mi - xl*xdist)/xdist < 0 &&
		(y - y_mi)/ydist > 0 && (y - y_mi - yl*ydist)/ydist < 0 &&
		(z - z_mi)/zdist > 0 && (z - z_mi - zl*zdist)/zdist < 0){
		// get coordinate index
		int indx = (int)floor((x - x_mi)/xdist);
		int indy = (int)floor((y - y_mi)/ydist);
		int indz = (int)floor((z - z_mi)/zdist);
		// scale coordinates to unit cube
		double xu = (x - x_mi - indx*xdist)/xdist;
		double yu = (y - y_mi - indy*ydist)/ydist;
		double zu = (z - z_mi - indz*zdist)/zdist;
		int i3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
		// tricubic interpolation
		double *coeff = const_cast<double*>(&Vc[i3][0]);
		V = tricubic_eval(coeff, xu, yu, zu);
		Ei[0] = -Escale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 0)/xdist;
		Ei[1] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 0)/ydist;
		Ei[2] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 0, 1)/zdist;
		if (dEidxj != NULL){ // calculate higher derivatives
			dEidxj[0][0] = -Escale*tricubic_eval(coeff, xu, yu, zu, 2, 0, 0)/xdist/xdist;
			dEidxj[0][1] = -Escale*tricubic_eval(coeff, xu, yu, zu, 1, 1, 0)/xdist/ydist;
			dEidxj[0][2] = -Escale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 1)/xdist/zdist;
			dEidxj[1][0] = -Escale*tricubic_eval(coeff, xu, yu, zu, 1, 1, 0)/ydist/xdist;
			dEidxj[1][1] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 2, 0)/ydist/ydist;
			dEidxj[1][2] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 1)/ydist/zdist;
			dEidxj[2][0] = -Escale*tricubic_eval(coeff, xu, yu, zu, 1, 0, 1)/zdist/xdist;
			dEidxj[2][1] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 1, 1)/zdist/ydist;
			dEidxj[2][2] = -Escale*tricubic_eval(coeff, xu, yu, zu, 0, 0, 2)/zdist/zdist;
		}

		for (int i = 0; i < 3; i++)
			FieldSmthr(x, y, z, Ei[i], dEidxj[i]); // apply field smoothing to each component
	}
}
