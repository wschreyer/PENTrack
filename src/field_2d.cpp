/**
 * \file
 * Bicubic interpolation of axisymmetric field tables.
 */

#include "field_2d.h"

#include <string>
#include <iostream>
#include <fstream>

#include "boost/format.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
void CylToCart(const double v_r, const double v_phi, const double phi, double &v_x, double &v_y){
	v_x = v_r*cos(phi) - v_phi*sin(phi);
	v_y = v_r*sin(phi) + v_phi*cos(phi);
}


TFieldContainer ReadOperaField2(const std::string &params, const std::map<std::string, std::string> &formulas){
    std::istringstream ss(params);
    boost::filesystem::path ft;
    std::string fieldtype, Bscale, Escale;
    double lengthconv;
    ss >> fieldtype >> ft >> Bscale >> Escale;
	Bscale = ResolveFormula(Bscale, formulas);
	Escale = ResolveFormula(Escale, formulas);
    if (fieldtype == "2Dtable"){
        std::cout << "Field type " << fieldtype << " is deprecated. Consider using the new OPERA2D format. I'm assuming that file " << ft << " is using centimeters, Gauss, Volt/centimeter, and Volts as units.\n";
        Bscale = "(" + Bscale + ")*0.0001"; // scale magnetic field to Tesla
        Escale = "(" + Escale + ")*100"; // scale electric field to Volt/meter
        lengthconv = 0.01;
    }
    else if (fieldtype == "OPERA2D") {
        ss >> lengthconv; // read fieldtype, tablefilename, and rest of parameters
    }
    else{
        throw std::runtime_error("Tried to load 2D table file for unknown field type " + fieldtype + "!\n");
    }
    if (!ss){
        throw std::runtime_error((boost::format("Could not read all required parameters for field %1%!") % fieldtype).str());
    }

    return TFieldContainer(std::unique_ptr<TabField>(new TabField(boost::filesystem::absolute(ft, configpath.parent_path()).string(), lengthconv)), Bscale, Escale);
}


void TabField::ReadTabFile(const std::string &tabfile, const double lengthconv, alglib::real_1d_array &rind, alglib::real_1d_array &zind,
		alglib::real_1d_array BTabs[3], alglib::real_1d_array ETabs[3], alglib::real_1d_array &VTab){
	
	boost::filesystem::path filepath = boost::filesystem::absolute(tabfile, configpath.parent_path());
	std::ifstream FINstream(filepath.string(), std::ifstream::in | std::ifstream::binary);
	boost::iostreams::filtering_istream FIN;
	if (filepath.extension() == ".bz2"){
		FIN.push(boost::iostreams::bzip2_decompressor());
	}
	else if (filepath.extension() == ".gz"){
		FIN.push(boost::iostreams::gzip_decompressor());
	}
	FIN.push(FINstream);
	if (!FINstream.is_open() or !FIN.is_complete()){
		throw std::runtime_error("Could not open " + filepath.string());
	}

	cout << "\nReading " << filepath << "!\n";
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
		std::cout << filepath << " not found or corrupt! Exiting...\n";
		exit(-1);
	}

	int ri = 0,zi = -1;
	double r, z, val;
	progress_display progress(n*m, std::cout);
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
		++progress;


		rind[ri] = r;
		zind[zi] = z;
		int i2 = zi * m + ri;
		if (BTabs[0].length() > 0){
			FIN >> val;
			BTabs[0][i2] = val;
		}
		if (BTabs[1].length() > 0){
			FIN >> val;
			BTabs[1][i2] = val;
		}
		if (BTabs[2].length() > 0){
			FIN >> val;
			BTabs[2][i2] = val;
		}
		if (ETabs[0].length() > 0){
			FIN >> val;
			ETabs[0][i2] = val;
		}
		if (ETabs[1].length() > 0){
			FIN >> val;
			ETabs[1][i2] = val;
		}
		if (ETabs[2].length() > 0){
			FIN >> val;
			ETabs[2][i2] = val;
		}
		if (VTab.length() > 0){
			FIN >> val;
			VTab[i2] = val;
		}
		FIN >> ws;
	}

	std::cout << "\n";
	if (ri+1 != (int)m || zi+1 != (int)n){
		std::cout << "The header says the size is " << m << " by " << n << ", actually it is " << ri+1 << " by " << zi+1 << "! Exiting...\n";
		exit(-1);
	}
}


void TabField::CheckTab(const alglib::real_1d_array &rind, const alglib::real_1d_array &zind,
		const alglib::real_1d_array BTabs[3], const alglib::real_1d_array ETabs[3], const alglib::real_1d_array &VTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
	r_mi = rind[0];
	z_mi = zind[0];
	rdist = rind[1] - rind[0];
	zdist = zind[1] - zind[0];
	std::cout << "The arrays are " << m << " by " << n << "\n";
	std::cout << "The r values go from " << rind[0] << " to " << rind[m-1] << "\n";
	std::cout << "The z values go from " << zind[0] << " to " << zind[n-1] << ".\n";
	std::cout << "rdist = " << rdist << " zdist = " << zdist << "\n";

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

	std::cout << "The input table file has values of magnetic field |B| from " << Babsmin << " to " << Babsmax << " and values of electric potential from " << Vmin << " to " << Vmax << "\n";
}

TabField::TabField(const std::string &tabfile, const double alengthconv){
	alglib::real_1d_array rind, zind, BTabs[3], ETabs[3], VTab;

	ReadTabFile(tabfile, alengthconv, rind, zind, BTabs, ETabs, VTab); // open tabfile and read values into arrays

	CheckTab(rind, zind, BTabs, ETabs, VTab); // print some info

	std::cout << "Starting Preinterpolation ... ";
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
	cout << "Done\n";
}

void TabField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	double r = sqrt(x*x+y*y);
	if (r >= r_mi && r < r_mi + rdist*(m - 1) && z >= z_mi && z < z_mi + zdist*(n - 1)){
		// bicubic interpolation
		double Br = 0, Bphi = 0;
		double Bx = 0, By = 0, Bz = 0;
		double phi = atan2(y,x);
		if (dBidxj != NULL){
			double dBrdr = 0, dBrdz = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0, dBxdz = 0, dBydz = 0, dBzdz = 0;
			// double dummy;
			if (fBrc){
				alglib::spline2ddiff(Brc, r, z, Br, dBrdr, dBrdz);
			}
			if (fBphic){
				alglib::spline2ddiff(Bphic, r, z, Bphi, dBphidr, dBphidz);
			}
			if (r > 0){
				dBidxj[0][0] = dBrdr*cos(phi)*cos(phi) - dBphidr*cos(phi)*sin(phi) + (Br*sin(phi)*sin(phi) + Bphi*cos(phi)*sin(phi))/r;
				dBidxj[0][1] = dBrdr*cos(phi)*sin(phi) - dBphidr*sin(phi)*sin(phi) - (Br*cos(phi)*sin(phi) + Bphi*cos(phi)*cos(phi))/r;
				dBidxj[1][0] = dBrdr*cos(phi)*sin(phi) + dBphidr*cos(phi)*cos(phi) - (Br*cos(phi)*sin(phi) - Bphi*sin(phi)*sin(phi))/r;
				dBidxj[1][1] = dBrdr*sin(phi)*sin(phi) + dBphidr*cos(phi)*sin(phi) + (Br*cos(phi)*cos(phi) - Bphi*cos(phi)*sin(phi))/r;
			}
			CylToCart(dBrdz,dBphidz,phi,dBxdz,dBydz);
			dBidxj[0][2] = dBxdz;
			dBidxj[1][2] = dBydz;
			if (fBzc){
				alglib::spline2ddiff(Bzc, r, z, Bz, dBzdr, dBzdz);
				dBidxj[2][0] = dBzdr*cos(phi);
				dBidxj[2][1] = dBzdr*sin(phi);
				dBidxj[2][2] = dBzdz;
			}
		}
		else{
			if (fBrc)
				Br = alglib::spline2dcalc(Brc, r, z);
			if (fBphic)
				Bphi = alglib::spline2dcalc(Bphic, r, z);
			if (fBzc)
				Bz = alglib::spline2dcalc(Bzc, r, z);
		}
		CylToCart(Br,Bphi,phi,Bx,By);
		B[0] = Bx;
		B[1] = By;
		B[2] = Bz;
	}
}


void TabField::EField(const double x, const double y, const double z, const double t,
		double &V, double Ei[3]) const{
	double r = sqrt(x*x+y*y);
	if (r >= r_mi && r < r_mi + rdist*(m - 1) && z >= z_mi && z < z_mi + zdist*(n - 1)){
		if (fErc || fEphic || fEzc){ // prefer E-field interpolation over potential interpolation
			double phi = atan2(y,x);
//			if (dEidxj == nullptr){
				alglib::real_1d_array Er("[0]"), Ephi("[0]"), Ez("[0]");
				double Ex, Ey;
				if (fErc)
					alglib::spline2dcalcvbuf(Erc, r, z, Er); // use cheaper interpolation if derivatives not required
				if (fEphic)
					alglib::spline2dcalcvbuf(Ephic, r, z, Ephi);
				if (fEzc)
					alglib::spline2dcalcvbuf(Ezc, r, z, Ez);
				CylToCart(Er[0], Ephi[0], phi, Ex, Ey); // convert r,phi components to x,y components
				Ei[0] = Ex; // set electric field
				Ei[1] = Ey;
				Ei[2] = Ez[0];
/*			}
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
			}*/
		}

        else if (fVc){
	  double dVdrj[3]; //, dummy;
            // bicubic interpolation
            alglib::spline2ddiff(Vc, r, z, V, dVdrj[0], dVdrj[2]);
	    double phi = atan2(y,x);
            Ei[0] = -dVdrj[0]*cos(phi);
            Ei[1] = -dVdrj[0]*sin(phi);
            Ei[2] = -dVdrj[2];
        }



    }
}
