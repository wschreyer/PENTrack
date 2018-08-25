/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */


#include "field_3d.h"

#include <fstream>
#include <iostream>

#include "interpolation.h"
#include "boost/format.hpp"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "tricubic.h"
#include "globals.h"


std::unique_ptr<TabField3> ReadComsolField(const std::string &params){
  std::istringstream ss(params);
  boost::filesystem::path ft;
  std::string fieldtype, Bscale;
  double BoundaryWidth, lengthconv;
  ss >> fieldtype >> ft >> Bscale >> BoundaryWidth >> lengthconv; // read fieldtype, tablefilename, and rest of parameters
  if (!ss)
      throw std::runtime_error("Could not read all required parameters for field COMSOL!");

  std::string line;
  std::vector<std::string> line_parts;
  std::vector<double> x, y, z;
  std::vector<double> bx, by, bz;

  std::ifstream FIN(boost::filesystem::absolute(ft, configpath.parent_path()).string(), std::ifstream::in);
  if (!FIN.is_open()){
    throw std::runtime_error("Could not open " + ft.string());
  }
  std::cout << "\nReading " << ft << "\n";

  // Read in file data
  int lineNum = 0;

  while (getline(FIN,line)){
    lineNum++;
    if (line.substr(0,1) == "%" || line.substr(0,1) == "#") continue;     // Skip commented lines
    boost::split(line_parts, line, boost::is_any_of("\t, "), boost::token_compress_on); //Delineate tab, space, commas

    if (line_parts.size() != 6){
      throw std::runtime_error((boost::format("Error reading line %1% of file %2%") % lineNum % ft.string()).str());
    }

    x.push_back( std::stod(line_parts[0], nullptr) * lengthconv);
    y.push_back( std::stod(line_parts[1], nullptr) * lengthconv);
    z.push_back( std::stod(line_parts[2], nullptr) * lengthconv);
    bx.push_back( std::stod(line_parts[3], nullptr));
    by.push_back( std::stod(line_parts[4], nullptr));
    bz.push_back( std::stod(line_parts[5], nullptr));
  }

  FIN.close();
  if (x.empty() || y.empty() || z.empty() || bx.empty() || by.empty()|| bz.empty() ) {
    throw std::runtime_error("No data read from " + ft.string());
  }

  return std::unique_ptr<TabField3>(new TabField3({x,y,z}, {bx,by,bz}, std::vector<double>(), Bscale, "0", BoundaryWidth));
}


std::unique_ptr<TabField3> ReadOperaField3(const std::string &params){
    std::istringstream ss(params);
    boost::filesystem::path ft;
    std::string fieldtype, Bscale, Escale;
    double BoundaryWidth, lengthconv;
    ss >> fieldtype >> ft >> Bscale >> Escale >> BoundaryWidth >> lengthconv;

    std::ifstream FIN(boost::filesystem::absolute(ft, configpath.parent_path()).string(), std::ifstream::in);
    if (!FIN.is_open())
        throw std::runtime_error((boost::format("\nCould not open %1%!\n") % ft).str());
    std::cout << "\nReading " << ft << " ";
	std::string line;
    int xl, yl, zl;
	FIN >> xl >> yl >> zl;

	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);

    std::array<std::vector<double>, 3> xyzTab, BTab;
    std::vector<double> VTab;
    for (auto &xi: xyzTab)
        xi.resize(xl*yl*zl);

	if (line.find("BX") != std::string::npos){
        BTab[0].resize(xl*yl*zl);
		getline(FIN,line);
	}
	if (line.find("BY") != std::string::npos){
        BTab[1].resize(xl*yl*zl);
		getline(FIN,line);
	}
	if (line.find("BZ") != std::string::npos){
        BTab[2].resize(xl*yl*zl);
		getline(FIN,line);
	}

	if (line.find("V") != std::string::npos){	// file contains potential?
		VTab.resize(xl*yl*zl);
		getline(FIN,line);
	}

	if (!FIN || line.substr(0,2) != " 0"){
        std::cout << ft << " not found or corrupt! Exiting...\n";
		exit(-1);
	}

    int perc = 0, i = 0;
	while (FIN.good()){
        double x, y, z, val;
        FIN >> x;
		FIN >> y;
		FIN >> z;
		if (!FIN) break;
        x *= lengthconv;
        y *= lengthconv;
        z *= lengthconv;

		// status if read is displayed
        PrintPercent((float)i/(xl*yl*zl), perc);

        xyzTab[0][i] = x;
        xyzTab[1][i] = y;
        xyzTab[2][i] = z;
        if (BTab[0].size() > 0){
			FIN >> val;
            BTab[0][i] = val;
		}
        if (BTab[1].size() > 0){
			FIN >> val;
            BTab[1][i] = val;
		}
        if (BTab[2].size() > 0){
			FIN >> val;
            BTab[2][i] = val;
		}
		if (VTab.size() > 0){
			FIN >> val;
            VTab[i] = val;
		}
		FIN >> std::ws;
        ++i;
    }

	std::cout << "\n";
    if (i != xl*yl*zl)
        throw std::runtime_error((boost::format("The header says the size is %1%, actually it is %2%! Exiting...\n") % (xl*yl*zl) % i).str());
	FIN.close();

    return std::unique_ptr<TabField3>(new TabField3(xyzTab, BTab, VTab, Bscale, Escale, BoundaryWidth));
}


void TabField3::CheckTab(const std::array<std::vector<double>, 3> &B, const std::vector<double> &V){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
    std::cout << "The arrays are " << xyz[0].size() << " by " << xyz[1].size() << " by " << xyz[2].size()
              << " (" << B[0].size() << ", " << B[1].size() << ", " << B[2].size() << ", " << V.size() << " field components).\n";
    std::cout << "The x values go from " << xyz[0].front() << " to " << xyz[0].back() << "\n";
    std::cout << "The y values go from " << xyz[1].front() << " to " << xyz[1].back() << "\n";
    std::cout << "The z values go from " << xyz[2].front() << " to " << xyz[2].back() << ".\n";

    std::vector<double> Babs;
    std::transform( boost::make_zip_iterator(boost::make_tuple(B[0].begin(), B[1].begin(), B[2].begin())),
                    boost::make_zip_iterator(boost::make_tuple(B[0].end(), B[1].end(), B[2].end())), std::back_inserter(Babs),
                    [](const boost::tuple<double,double,double> &Bs){ return std::sqrt(std::pow(Bs.get<0>(), 2) + std::pow(Bs.get<1>(), 2) + std::pow(Bs.get<2>(), 2)); });
    std::copy(V.begin(), V.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout.flush();
    std::cout << "The input table file has values of |B| from " << *std::min_element(Babs.begin(), Babs.end()) << " to " << *std::max_element(Babs.begin(), Babs.end());
    if (V.size() > 0)
        std::cout << " and values of V from " << *std::min_element(V.begin(), V.end()) << " V to " << *std::max_element(V.begin(), V.end()) << " V";
    std::cout << std::endl;

    unsigned long l = xyz[0].size() * xyz[1].size() * xyz[2].size();
    if ((B[0].size() > 0 && B[0].size() != l) || (B[1].size() > 0 && B[1].size() != l) || (B[2].size() > 0 && B[2].size() != l) || (V.size() > 0 && V.size() != l))
        std::cout << "Warning: Number of field samples does not match number of grid points! Missing points will default to zero field\n";
}


void TabField3::CalcDerivs(const array3D &Tab, const unsigned long diff_dim, array3D &DiffTab) const
{
    alglib::real_1d_array x, y, diff;
    std::array<unsigned long, 3> index;
    DiffTab.resize(boost::extents[Tab.shape()[0]][Tab.shape()[1]][Tab.shape()[2]]);
    unsigned long len = xyz[diff_dim].size();
    x.setcontent(len, &xyz[diff_dim][0]);
    y.setlength(len);
    diff.setlength(len);
//    std::cout << "x: " << x.tostring(2) << std::endl;

    unsigned long dim1 = (diff_dim + 1) % 3;
    unsigned long dim2 = (diff_dim + 2) % 3;
    for (index[dim1] = 0; index[dim1] < xyz[dim1].size(); index[dim1] = index[dim1] + 1){ // iterate over all indices in both other dimensions
        for (index[dim2] = 0; index[dim2] < xyz[dim2].size(); index[dim2] = index[dim2] + 1){
//            std::copy(index.begin(), index.end(), std::ostream_iterator<unsigned long>(std::cout, " "));
            for (index[diff_dim] = 0; index[diff_dim] < len; index[diff_dim] = index[diff_dim] + 1){
                y[index[diff_dim]] = Tab(index);
            }
//            std::cout << "y: " << y.tostring(2) << std::endl;
            alglib::spline1dgriddiffcubic(x,y,diff); // get derivatives with respect to coordinate dimension diff_dim
//            std::cout << " diff: " << diff.tostring(2) << std::endl;
            for (index[diff_dim] = 0; index[diff_dim] < len; index[diff_dim] = index[diff_dim] + 1){
                DiffTab(index) = diff[index[diff_dim]];
            }
//            std::cout << " diff: " << diff.tostring(2) << std::endl;
        }
    }
}


void TabField3::PreInterpol(const array3D &Tab, field_type &coeff) const{
    std::array<unsigned long, 3> len;
    for (unsigned long i = 0; i < 3; ++i){
        len[i] = xyz[i].size() - 1;
    }
    coeff.resize(len);

    array3D dFdx, dFdy, dFdz, dFdxdy, dFdxdz, dFdydz, dFdxdydz; // derivatives with respect to x, y, z, xy, xz, yz, xyz
    CalcDerivs(Tab, 0, dFdx); // dF/dx
    CalcDerivs(Tab, 1, dFdy); // dF/dy
    CalcDerivs(Tab, 2, dFdz); // dF/dz
    CalcDerivs(dFdx, 1, dFdxdy); // d2F/dxdy
    CalcDerivs(dFdx, 2, dFdxdz); // d2F/dxdz
    CalcDerivs(dFdy, 2, dFdydz); // d2F/dydz
    CalcDerivs(dFdxdy, 2, dFdxdydz); // d3F/dxdydz

    for (unsigned long ix = 0; ix < len[0]; ++ix){
        for (unsigned long iy = 0; iy < len[1]; ++iy){
            for (unsigned long iz = 0; iz < len[2]; ++iz){
                std::array<std::array<unsigned long, 3>, 8> indices;
                indices[0] = {ix  ,iy  ,iz  }; // collect indices of corners of each grid cell
                indices[1] = {ix+1,iy  ,iz  }; // order according to tricubic manual
                indices[2] = {ix  ,iy+1,iz  };
                indices[3] = {ix+1,iy+1,iz  };
                indices[4] = {ix  ,iy  ,iz+1};
                indices[5] = {ix+1,iy  ,iz+1};
                indices[6] = {ix  ,iy+1,iz+1};
                indices[7] = {ix+1,iy+1,iz+1};

                double cellx = xyz[0][ix+1] - xyz[0][ix];
                double celly = xyz[1][iy+1] - xyz[1][iy];
                double cellz = xyz[2][iz+1] - xyz[2][iz];
                std::array<std::array<double, 8>, 8> yyy;
                for (unsigned i = 0; i < 8; ++i){
                    yyy[0][i] = Tab(indices[i]); // get values and derivatives at each corner of grid cell
                    yyy[1][i] = dFdx(indices[i])*cellx;
                    yyy[2][i] = dFdy(indices[i])*celly;
                    yyy[3][i] = dFdz(indices[i])*cellz;
                    yyy[4][i] = dFdxdy(indices[i])*cellx*celly;
                    yyy[5][i] = dFdxdz(indices[i])*cellx*cellz;
                    yyy[6][i] = dFdydz(indices[i])*celly*cellz;
                    yyy[7][i] = dFdxdydz(indices[i])*cellx*celly*cellz;
                }
                tricubic_get_coeff(&coeff(indices[0])[0], &yyy[0][0], &yyy[1][0], &yyy[2][0], &yyy[3][0], &yyy[4][0], &yyy[5][0], &yyy[6][0], &yyy[7][0]); // calculate tricubic interpolation coefficients and store in coeff
            }
        }
    }
}


TabField3::TabField3(const std::array<std::vector<double>, 3> &xyzTab, const std::array<std::vector<double>, 3> &BTab, const std::vector<double> &VTab,
                     const std::string &Bscale, const std::string &Escale, const double aBoundaryWidth)
	: TField(Bscale, Escale){
	BoundaryWidth = aBoundaryWidth;

    for (unsigned i = 0; i < 3; ++i){
        std::unique_copy(xyzTab[i].begin(), xyzTab[i].end(), std::back_inserter(xyz[i])); // get list of unique x, y, and z coordinates
        std::sort(xyz[i].begin(), xyz[i].end());
        auto last = std::unique(xyz[i].begin(), xyz[i].end());
        xyz[i].erase(last, xyz[i].end());
    }
    CheckTab(BTab,VTab); // print some info


    std::array<array3D, 3> B;
    array3D V;
    for (unsigned i = 0; i < 3; ++i){
        if (BTab[i].size() > 0){
            B[i].resize(boost::extents[xyz[0].size()][xyz[1].size()][xyz[2].size()]);
            std::fill_n(B[i].data(), B[i].num_elements(), 0.);
        }
    }
    if (VTab.size() > 0){
        V.resize(boost::extents[xyz[0].size()][xyz[1].size()][xyz[2].size()]);
        std::fill_n(V.data(), V.num_elements(), 0.);
    }

    for (unsigned int i = 0; i < xyzTab[0].size(); ++i){
        std::array<long, 3> index;
        for (unsigned j = 0; j < 3; ++j){
            auto found = std::lower_bound(xyz[j].begin(), xyz[j].end(), xyzTab[j][i]);
            assert(found != xyz[j].end());
            index[j] = std::distance(xyz[j].begin(), found);
        }
        for (unsigned j = 0; j < 3; ++j){
            if (BTab[j].size() > 0){
                B[j](index) = BTab[j][i];
            }
        }
        if (VTab.size() > 0)
            V(index) = VTab[i];
    }

	std::cout << "Starting Preinterpolation ... ";
	float size = 0;
    if (BTab[0].size() > 0){
		std::cout << "Bx ... ";
		std::cout.flush();
        PreInterpol(B[0], Bc[0]); // precalculate interpolation coefficients for B field
        size += float(Bc[0].num_elements()*64*sizeof(double)/1024/1024);
    }
    if (BTab[1].size() > 0){
		std::cout << "By ... ";
		std::cout.flush();
        PreInterpol(B[1], Bc[1]);
        size += float(Bc[1].num_elements()*64*sizeof(double)/1024/1024);
	}
    if (BTab[2].size() > 0){
		std::cout << "Bz ... ";
		std::cout.flush();
        PreInterpol(B[2], Bc[2]);
        size += float(Bc[2].num_elements()*64*sizeof(double)/1024/1024);
	}
	if (VTab.size() > 0){
		std::cout << "V ... ";
		std::cout.flush();
        PreInterpol(V, Vc);
        size += float(Vc.num_elements()*64*sizeof(double)/1024/1024);
	}
	std::cout << "Done (" << size << " MB)\n";
}


void TabField3::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
    double Bscale = BScaling(t);
    if (Bscale == 0)
        return;

    std::array<double, 3> r = {x, y, z}, dist;
    // get coordinate index
    std::array<long, 3> index;
    for (unsigned i = 0; i < 3; ++i){
      auto it = std::upper_bound(xyz[i].begin(), xyz[i].end(), r[i]); // find first coordinate larger than x/y/z
      if (it == xyz[i].end() || it == xyz[i].begin()){ // if x,y,z are outside bounds of field
          return;
      }
      auto low = it - 1;
      index[i] = std::distance(xyz[i].begin(), low);
      dist[i] = *it - *low;
      r[i] = (r[i] - *low)/dist[i]; // scale coordinates to unit cube
    }

    // tricubic interpolation
    for (unsigned i = 0; i < 3; ++i){
        if (Bc[i].size() > 0){
            double *coeff = const_cast<double*>(&Bc[i](index)[0]);
            B[i] = Bscale*tricubic_eval(coeff, r[0], r[1], r[2]);
            if (dBidxj != nullptr){
                dBidxj[i][0] = Bscale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 0, 0)/dist[0];
                dBidxj[i][1] = Bscale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 1, 0)/dist[1];
                dBidxj[i][2] = Bscale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 0, 1)/dist[2];
            }
        }
    }
    for (int i = 0; i < 3; i++){
        if (dBidxj != nullptr)
            FieldSmthr(x, y, z, B[i], dBidxj[i]); // apply field smoothing to each component
        else
            FieldSmthr(x, y, z, B[i], nullptr);
    }
}

void TabField3::FieldSmthr(const double x, const double y, const double z, double &F, double dFdxi[3]) const{
    if (BoundaryWidth != 0 && F != 0){ // skip, if BoundaryWidth is set to zero
        double dxlo = (x - xyz[0].front())/BoundaryWidth; // calculate distance to edges in units of BoundaryWidth
        double dxhi = (xyz[0].back() - x)/BoundaryWidth;
        double dylo = (y - xyz[1].front())/BoundaryWidth;
        double dyhi = (xyz[1].back() - y)/BoundaryWidth;
        double dzlo = (z - xyz[2].front())/BoundaryWidth;
        double dzhi = (xyz[2].back() - z)/BoundaryWidth;

		// F'(x,y,z) = F(x,y,z)*f(x)*f(y)*f(z)
		// dF'/dx = dF/dx*f(x)*f(y)*f(z) + F*df(x)/dx*f(y)*f(z) --> similar for dF'/dy and dF'/dz

		double Fscale = 1; // f(x)*f(y)*f(z)
		double dFadd[3] = {0, 0, 0}; // F*df(x_i)/dx_i / f(x_i)

		if (dxhi < 1) { // if point in upper x boundary (distance smaller than BoundaryWidth)
			Fscale *= SmthrStp(dxhi); // scale field by value of smoother function
            if (dFdxi != nullptr)
				dFadd[0] = -F*SmthrStpDer(dxhi)/BoundaryWidth/SmthrStp(dxhi); // add derivative of smoother function according to product rule
		}
		if (dyhi < 1){ // if point in upper y boundary
			Fscale *= SmthrStp(dyhi);
            if (dFdxi != nullptr)
				dFadd[1] = -F*SmthrStpDer(dyhi)/BoundaryWidth/SmthrStp(dyhi);
		}
		if (dzhi < 1){ // if point in upper z boundary
			Fscale *= SmthrStp(dzhi);
            if (dFdxi != nullptr)
				dFadd[2] = -F*SmthrStpDer(dzhi)/BoundaryWidth/SmthrStp(dzhi);
		}

		if (dxlo < 1){ // if point in lower x boundary
			Fscale *= SmthrStp(dxlo);
            if (dFdxi != nullptr)
				dFadd[0] = F*SmthrStpDer(dxlo)/BoundaryWidth/SmthrStp(dxlo);
		}
		if (dylo < 1){ // if point in lower y boundary
			Fscale *= SmthrStp(dylo);
            if (dFdxi != nullptr)
				dFadd[1] = F*SmthrStpDer(dylo)/BoundaryWidth/SmthrStp(dylo);
		}
		if (dzlo < 1){ // if point in lower z boundary
			Fscale *= SmthrStp(dzlo);
            if (dFdxi != nullptr)
				dFadd[2] = F*SmthrStpDer(dzlo)/BoundaryWidth/SmthrStp(dzlo);
		}

		if (Fscale != 1){
			F *= Fscale; // scale field value
            if (dFdxi != nullptr){
				for (int i = 0; i < 3; i++){
					dFdxi[i] = dFdxi[i]*Fscale + dFadd[i]*Fscale; // scale derivatives according to product rule
				}
			}
		}

	}
}

double TabField3::SmthrStp(const double x) const{
    return 6*std::pow(x, 5) - 15*std::pow(x, 4) + 10*std::pow(x, 3);
}

double TabField3::SmthrStpDer(const double x) const{
    return 30*std::pow(x, 4) - 60*std::pow(x, 3) + 30*std::pow(x,2);
}

void TabField3::EField(const double x, const double y, const double z, const double t,
		double &V, double Ei[3], double dEidxj[3][3]) const{
    double Escale = EScaling(t);
    if (Escale != 0 && Vc.size() > 0){
        std::array<double, 3> r = {x, y, z}, dist;
        // get coordinate index
        std::array<int, 3> index;
        for (int i = 0; i < 3; ++i){
          auto it = std::upper_bound(xyz[i].begin(), xyz[i].end(), r[i]); // find first coordinate larger than x/y/z
          if (it == xyz[i].end() || it == xyz[i].begin()) // if x,y,z are outside bounds of field
              return;
          auto low = it - 1;
          index[i] = std::distance(low, xyz[i].begin());
          dist[i] = *it - *low;
          r[i] = (r[i] - *low)/dist[i]; // scale coordinates to unit cube
        }

        double *coeff = const_cast<double*>(&Vc[index[0]][index[1]][index[2]][0]);
        V = tricubic_eval(coeff, r[0], r[1], r[2]);
        Ei[0] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 0, 0)/dist[0];
        Ei[1] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 1, 0)/dist[1];
        Ei[2] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 0, 1)/dist[2];
        if (dEidxj != nullptr){ // calculate higher derivatives
            dEidxj[0][0] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 2, 0, 0)/dist[0]/dist[0];
            dEidxj[0][1] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 1, 0)/dist[0]/dist[1];
            dEidxj[0][2] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 0, 1)/dist[0]/dist[2];
            dEidxj[1][0] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 1, 0)/dist[1]/dist[0];
            dEidxj[1][1] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 2, 0)/dist[1]/dist[1];
            dEidxj[1][2] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 1, 1)/dist[1]/dist[2];
            dEidxj[2][0] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 1, 0, 1)/dist[2]/dist[0];
            dEidxj[2][1] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 1, 1)/dist[2]/dist[1];
            dEidxj[2][2] = -Escale*tricubic_eval(coeff, r[0], r[1], r[2], 0, 0, 2)/dist[2]/dist[2];
		}

        for (int i = 0; i < 3; i++){
            if (dEidxj != nullptr)
                FieldSmthr(x, y, z, Ei[i], dEidxj[i]); // apply field smoothing to each component
            else
                FieldSmthr(x, y, z, Ei[i], nullptr);
        }
	}
}
