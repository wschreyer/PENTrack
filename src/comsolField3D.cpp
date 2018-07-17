/**
 * \file
 * Tricubic interpolation of a COMSOL arrow volume plot
 * 
 */


#include "comsolField3D.h"

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <set>
#include <string>
#include <fstream>

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
int COMSOL_INDEX_3D(int xi, int yi, int zi, int xsize, int ysize, int zsize){
	return zi*xsize*ysize + yi*xsize + xi;
}


void comsolField3::ReadTabFile(const std::string &tabfile,
		std::vector<double> &BxTab, std::vector<double> &ByTab, std::vector<double> &BzTab){

	std::string line;
  std::vector<std::string> line_parts;
  std::vector<double> xind, yind, zind;
  std::vector<double> x, y, z;
  std::vector<double> bx, by, bz;
  std::set<double> setX, setY, setZ;

  std::ifstream FIN(tabfile, std::ifstream::in);
  if (!FIN.is_open()){
		std::cout << "\nCould not open " << tabfile << "!\n";
		exit(-1);
	}
	std::cout << "\nReading " << tabfile << "\n";

  // Read in file data
  int lineNum = 0;

  while (getline(FIN,line)){
    lineNum++;
    if (line.substr(0,1) == "%" || line.substr(0,1) == "#") continue;     // Skip commented lines
    boost::split(line_parts, line, boost::is_any_of("\t, "), boost::token_compress_on); //Delineate tab, space, commas

    if (line_parts.size() != 6){
      std::cout << "\nError reading line " << lineNum << " of file "<< tabfile << "\n...quitting\n";
      exit(-1);
    }

		// String to double conversion with error handling
		try {
	    x.push_back( boost::lexical_cast<double>(line_parts[0]) * lengthconv);
	    y.push_back( boost::lexical_cast<double>(line_parts[1]) * lengthconv);
	    z.push_back( boost::lexical_cast<double>(line_parts[2]) * lengthconv);
	    bx.push_back( boost::lexical_cast<double>(line_parts[3]) * Bconv);
	    by.push_back( boost::lexical_cast<double>(line_parts[4]) * Bconv);
	    bz.push_back( boost::lexical_cast<double>(line_parts[5]) * Bconv);
		}
		catch(boost::bad_lexical_cast const& e)
		{
			std::cout << "\nError reading line " << lineNum << " of file "<< tabfile << std::endl;
			std::cout << e.what() << std::endl;
			exit(-1);
		}

    // Put x, y, z values into ordered set for later use
    setX.insert( x.back() );
    setY.insert( y.back() );
    setZ.insert( z.back() );
  }

  FIN.close();
  if (x.empty() || y.empty() || z.empty() || bx.empty() || by.empty()|| bz.empty() ) {
    std::cout<< "Error: no data read in from " << tabfile << "\n...quitting\n";
    exit(-1);
  }

  // We need to find xl, yl, zl. This is the size of the table in each
  // respective direction (aka the # of x, y, z values if we ignore duplicates)
  xl = setX.size();
  yl = setY.size();
  zl = setZ.size();

  // Fill in vectors xind, yind, zind, which contain the afore mentioned values
  // of non-duplicated x, y, z
  // (thanks to set, these are now sorted in ascending order)
  std::set<double>::iterator itX, itY, itZ;
  for (itX = setX.begin(), itY = setY.begin(), itZ = setZ.begin();
        (itX != setX.end()) && (itY != setY.end()) && (itZ != setZ.end());
        ++itX, ++itY, ++itZ)
  {
    xind.push_back( *itX );
    yind.push_back( *itY );
    zind.push_back( *itZ );
  }

	// Resize and initially fill in B field vectors
  BxTab.resize(xl*yl*zl);
  ByTab.resize(xl*yl*zl);
  BzTab.resize(xl*yl*zl);

  // xi, yi, zi are the indices of the values in xind, yind, zind
  // i3 is some magical index defined by COMSOL_INDEX_3D used for interpolation
	int xi = 0, yi = 0, zi = 0;
  int i3 = 0;

	// Now fill in BxTab, ByTab, BzTab with the correct indexing for interpolation
  for (unsigned int j = 0; j < x.size(); j++){

    xi = std::distance( setX.begin(), setX.find(x[j]) ) ;
    yi = std::distance( setY.begin(), setY.find(y[j]) ) ;
    zi = std::distance ( setZ.begin(), setZ.find(z[j]) ) ;

    i3 = COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl);

    BxTab[i3] = bx[j];
    ByTab[i3] = by[j];
    BzTab[i3] = bz[j];

	}

	xdist = xind[1] - xind[0];
	ydist = yind[1] - yind[0];
	zdist = zind[1] - zind[0];
	x_mi = xind[0];
	y_mi = yind[0];
	z_mi = zind[0];

};


void comsolField3::CheckTab(const std::vector<double> &BxTab, const std::vector<double> &ByTab,
		const std::vector<double> &BzTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
	std::cout << "The arrays are " << xl << " by " << yl << " by " << zl << ".\n";
	std::cout << "The x values go from " << x_mi << " to " << x_mi + xdist*(xl-1) << "[m]\n";
	std::cout << "The y values go from " << y_mi << " to " << y_mi + ydist*(yl-1) << "[m]\n";
	std::cout << "The z values go from " << z_mi << " to " << z_mi + zdist*(zl-1) << "[m]\n";
	std::cout << "xdist = " << xdist << " ydist = " << ydist << " zdist = " << zdist << "\n";

	double Babsmax = 0, Babsmin = 9e99, Babs;
	for (int j=0; j < xl; j++)
	{
		for (int k=0; k < yl; k++)
		{
			for (int l=0; l < zl; l++)
			{
				Babs = 0;
				int i3 = COMSOL_INDEX_3D(j, k, l, xl, yl, zl);
				if (BxTab.size() > 0)
					Babs += BxTab[i3] * BxTab[i3];
				if (ByTab.size() > 0)
					Babs += ByTab[i3] * ByTab[i3];
				if (BzTab.size() > 0)
					Babs += BzTab[i3] * BzTab[i3];
				Babsmax = std::max(sqrt(Babs),Babsmax);
				Babsmin = std::min(sqrt(Babs),Babsmin);

			}
		}
	}

	std::cout << "The input table file has values of |B| from " << Babsmin << " T to " << Babsmax << " T\n";
};


void comsolField3::CalcDerivs(const std::vector<double> &Tab, std::vector<double> &Tab1, std::vector<double> &Tab2, std::vector<double> &Tab3,
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
				y[xi] = Tab[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)]; // get derivatives dF/dx from spline interpolation
			}
			alglib::spline1dgriddiffcubic(x, y, diff);
			for (int xi = 0; xi < xl; xi++)
				Tab1[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[xi];
		}
	}

	x.setlength(yl);
	y.setlength(yl);
	diff.setlength(yl);
	for (int xi = 0; xi < xl; xi++){
		for (int zi = 0; zi < zl; zi++){
			for (int yi = 0; yi < yl; yi++){
				x[yi] = y_mi + yi*ydist;
				y[yi] = Tab[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)]; // get derivatives dF/dy from spline interpolation
			}
			alglib::spline1dgriddiffcubic(x, y, diff);
			for (int yi = 0; yi < yl; yi++)
				Tab2[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[yi];
		}
	}

	x.setlength(zl);
	y.setlength(zl);
	diff.setlength(zl);
	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get derivatives dF/dz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab3[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	x.setlength(yl);
	y.setlength(yl);
	diff.setlength(yl);
	for (int xi = 0; xi < xl; xi++){
		for (int zi = 0; zi < zl; zi++){
			for (int yi = 0; yi < yl; yi++){
				x[yi] = y_mi + yi*ydist;
				y[yi] = Tab1[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dxdy from spline interpolation
			for (int yi = 0; yi < yl; yi++)
				Tab12[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[yi];
		}
	}

	x.setlength(zl);
	y.setlength(zl);
	diff.setlength(zl);
	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab1[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dxdz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab13[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab2[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d2F/dydz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab23[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

	for (int xi = 0; xi < xl; xi++){
		for (int yi = 0; yi < yl; yi++){
			for (int zi = 0; zi < zl; zi++){
				x[zi] = z_mi + zi*zdist;
				y[zi] = Tab12[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)];
			}
			alglib::spline1dgriddiffcubic(x, y, diff); // get cross derivatives d3F/dxdydz from spline interpolation
			for (int zi = 0; zi < zl; zi++)
				Tab123[COMSOL_INDEX_3D(xi, yi, zi, xl, yl, zl)] = diff[zi];
		}
	}

};


void comsolField3::PreInterpol(std::vector<std::vector<double> > &coeff, const std::vector<double> &Tab) const{
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
				int i3[8] ={COMSOL_INDEX_3D(indx  , indy  , indz  , xl, yl, zl),
							COMSOL_INDEX_3D(indx+1, indy  , indz  , xl, yl, zl),
							COMSOL_INDEX_3D(indx  , indy+1, indz  , xl, yl, zl),
							COMSOL_INDEX_3D(indx+1, indy+1, indz  , xl, yl, zl),
							COMSOL_INDEX_3D(indx  , indy  , indz+1, xl, yl, zl),
							COMSOL_INDEX_3D(indx+1, indy  , indz+1, xl, yl, zl),
							COMSOL_INDEX_3D(indx  , indy+1, indz+1, xl, yl, zl),
							COMSOL_INDEX_3D(indx+1, indy+1, indz+1, xl, yl, zl)};

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
				int ind3 = COMSOL_INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
				coeff[ind3].resize(64);
				tricubic_get_coeff(&coeff[ind3][0], yyy, yyy1, yyy2, yyy3, yyy12, yyy13, yyy23, yyy123);
			}
		}
	}
};


comsolField3::comsolField3(const std::string &tabfile, const std::string &Bscale,
		const double aBoundaryWidth, const double alengthconv, const double aBconv)
	: TField(Bscale, "0"){
	BoundaryWidth = aBoundaryWidth;
	lengthconv = alengthconv;
	Bconv = aBconv;


	std::vector<double> BxTab, ByTab, BzTab;	// Bx/By/Bz values
	ReadTabFile(tabfile,BxTab,ByTab,BzTab); // open tabfile and read values into arrays

	CheckTab(BxTab,ByTab,BzTab); // print some info

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
	std::cout << "Done (" << size << " MB)\n";
}


void comsolField3::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
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
		int i3 = COMSOL_INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
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

void comsolField3::FieldSmthr(const double x, const double y, const double z, double &F, double dFdxi[3]) const{
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

double comsolField3::SmthrStp(const double x) const{
	return 6*pow(x, 5) - 15*pow(x, 4) + 10*pow(x, 3);
}

double comsolField3::SmthrStpDer(const double x) const{
	return 30*pow(x, 4) - 60*pow(x, 3) + 30*pow(x,2);
}
