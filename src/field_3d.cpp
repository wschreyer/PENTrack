/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */


#include "field_3d.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "interpolation.h"
#include "boost/format.hpp"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "tricubic.h"
#include "globals.h"

/**
 * A faster implementation of the tricubic_eval function from libtricubic.
 * 
 * Uses Horner's scheme.
 * 
 * @param a Interpolation parameters (64 doubles)
 * @param x X coordinate of point field should be evaluated at
 * @param y Y coordinate
 * @param z Z coordinate
 * 
 * @return Interpolated value at point
 */
inline double tricubic_eval_fast(double a[64], double x, double y, double z) {
  /*	double result = 0.0;
	for (int i = 0; i < 4; ++i) {
	for (int j = 0; j < 4; ++j) {
	for (int k = 0; k < 4; ++k) {
	result += a[i + 4 * j + 16 * k] * std::pow(x, i) * std::pow(y, j) * std::pow(z, k);
	}
	}
	}
	return result;
  */
  double rx = 0.;
  for (int i = 3; i >= 0; --i){
    double ry = 0.;
    for (int j = 3; j >= 0; --j){
      double rz = 0.;
      for (int k = 3; k >= 0; --k){
	rz += a[i + 4*j + 16*k];
	if (k > 0) rz *= z;
      }
      ry += rz;
      if (j > 0) ry *= y;
    }
    rx += ry;
    if (i > 0) rx *= x;
  }

  return rx;
}


/**
 * A faster implementation of the tricubic_eval function from libtricubic.
 * 
 * Uses Horner's scheme.
 * 
 * @param a Interpolation parameters (64 doubles)
 * @param x X coordinate of point field should be evaluated at
 * @param y Y coordinate
 * @param z Z coordinate
 * @param derx Return derx-th derivative with respect to x
 * @param dery Return dery-th derivative with respect to y
 * @param derz Return derz-th derivative with respect to z
 * 
 * @return Interpolated value at point
 */
inline double tricubic_eval_fast(double a[64], double x, double y, double z, int derx, int dery, int derz) {
  // fact[i][n] is the coefficient in front of the nth derivative of x^i.
  // For instance, d2/dx2 x^3 = 6 x, so fact[3][2] = 6.
  static double fact[4][4] {
    { 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 0.0, 0.0 },
    { 1.0, 2.0, 2.0, 0.0 },
    { 1.0, 3.0, 6.0, 6.0 },
  };
  /*
    double result = 0.0;
    for (int i = derx; i <= 3; ++i) {
    for (int j = dery; j <= 3; ++j) {
    for (int k = derz; k <= 3; ++k) {
    double deriv_factor = fact[i][derx] * fact[j][dery] * fact[k][derz];
    double pow_factor = std::pow(x, i - derx) * std::pow(y, j - dery) * std::pow(z, k - derz);
    result += a[i + 4 * j + 16 * k] * pow_factor * deriv_factor;
    }
    }
    }
    return result;
  */
  double rx = 0.;
  for (int i = 3; i >= derx; --i){
    double ry = 0.;
    for (int j = 3; j >= dery; --j){
      double rz = 0.;
      for (int k = 3; k >= derz; --k){
	rz += fact[k][derz] * a[i + 4*j + 16*k];
	if (k > derz) rz *= z;
      }
      ry += fact[j][dery] * rz;
      if (j > dery) ry *= y;
    }
    rx += fact[i][derx] * ry;
    if (i > derx) rx *= x;
  }

  return rx;
}


TFieldContainer ReadComsolField(const std::string &params, const std::map<std::string, std::string> &formulas){
  std::istringstream ss(params);
  boost::filesystem::path ft;
  std::string fieldtype, Bscale;
  double BoundaryWidth, lengthconv;
  ss >> fieldtype >> ft >> Bscale >> BoundaryWidth >> lengthconv; // read fieldtype, tablefilename, and rest of parameters
  Bscale = ResolveFormula(Bscale, formulas);
  if (!ss){
    throw std::runtime_error((boost::format("Could not read all required parameters for field %1%!") % fieldtype).str());
  }
  ft = boost::filesystem::absolute(ft, configpath.parent_path());

  std::string line;
  std::vector<std::string> line_parts;
  std::vector<double> x, y, z;
  std::vector<double> bx, by, bz;

  std::ifstream FINstream(ft.string(), std::ifstream::in);
  boost::iostreams::filtering_istream FIN;
  if (ft.extension() == ".bz2"){
    FIN.push(boost::iostreams::bzip2_decompressor());
  }
  else if (ft.extension() == ".gz"){
    FIN.push(boost::iostreams::gzip_decompressor());
  }
  FIN.push(FINstream);
  if (!FINstream.is_open() or !FIN.is_complete()){
    throw std::runtime_error("Could not open " + ft.string());
  }
  std::cout << "\nReading " << ft << "\n";

  // Read in file data
  int lineNum = 0;

  while (getline(FIN,line)){
    lineNum++;
    if (line.substr(0,1) == "%" || line.substr(0,1) == "#") continue;     // Skip commented lines
    boost::split(line_parts, line, boost::is_any_of("\t, "), boost::token_compress_on); //Delineate tab, space, commas

    if (line_parts.size() < 6){
      throw std::runtime_error((boost::format("Error reading line %1% of file %2%") % lineNum % ft.string()).str());
    }

    x.push_back( std::stod(line_parts[0], nullptr) * lengthconv);
    y.push_back( std::stod(line_parts[1], nullptr) * lengthconv);
    z.push_back( std::stod(line_parts[2], nullptr) * lengthconv);
    bx.push_back( std::stod(line_parts[3], nullptr));
    by.push_back( std::stod(line_parts[4], nullptr));
    bz.push_back( std::stod(line_parts[5], nullptr));
  }

  if (x.empty() || y.empty() || z.empty() || bx.empty() || by.empty()|| bz.empty() ) {
    throw std::runtime_error("No data read from " + ft.string());
  }

  auto xminmax = std::minmax_element(x.begin(), x.end());
  auto yminmax = std::minmax_element(y.begin(), y.end());
  auto zminmax = std::minmax_element(z.begin(), z.end());
  return TFieldContainer(std::unique_ptr<TabField3>(new TabField3({x,y,z}, {bx,by,bz}, std::vector<double>(), std::array<std::array<std::vector<double>, 7>, 3>())), Bscale, "0", *xminmax.second, *xminmax.first, *yminmax.second, *yminmax.first, *zminmax.second, *zminmax.first, BoundaryWidth);
}

TFieldContainer ReadOperaField3(const std::string &params, const std::map<std::string, std::string> &formulas){
  std::istringstream ss(params);
  boost::filesystem::path ft;
  std::string fieldtype, Bscale, Escale;
  double BoundaryWidth, lengthconv;
  ss >> fieldtype >> ft >> Bscale >> Escale >> BoundaryWidth;
  Bscale = ResolveFormula(Bscale, formulas);
  Escale = ResolveFormula(Escale, formulas);
  if (fieldtype == "3Dtable"){
    std::cout << "Field type " << fieldtype << " is deprecated. Consider using the new OPERA3D format. I'm assuming that file " << ft << " is using centimeters, Gauss, Volt/centimeter, and Volts as units.\n";
    Bscale = "(" + Bscale + ")*0.0001"; // scale magnetic field to Tesla
    Escale = "(" + Escale + ")*100"; // scale electric field to Volt/meter
    lengthconv = 0.01;
  }
  else if (fieldtype == "OPERA3D"){
    ss >> lengthconv;
  }
  else{
    throw std::runtime_error("Tried to load 3D table file for unknown field type " + fieldtype + "!\n");
  }
  if (!ss){
    throw std::runtime_error((boost::format("Could not read all required parameters for field %1%!") % fieldtype).str());
  }

  ft = boost::filesystem::absolute(ft, configpath.parent_path());
  std::ifstream FINstream(ft.string(), std::ifstream::in);
  boost::iostreams::filtering_istream FIN;
  if (ft.extension() == ".bz2"){
    FIN.push(boost::iostreams::bzip2_decompressor());
  }
  else if (ft.extension() == ".gz"){
    FIN.push(boost::iostreams::gzip_decompressor());
  }
  FIN.push(FINstream);
  if (!FINstream.is_open() or !FIN.is_complete()){
    throw std::runtime_error("Could not open " + ft.string());
  }

  
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
  for (auto &xi: xyzTab){
    xi.resize(xl * yl * zl);
  }
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

  std::vector<double> xind(xl), yind(yl), zind(zl);
  progress_display progress(xl*yl*zl, std::cout);
  int i = 0;
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
    ++progress;

    xyzTab[0][i] = x;
    xyzTab[1][i] = y;
    xyzTab[2][i] = z;
    if (not BTab[0].empty()){
      FIN >> val;
      BTab[0][i] = val;
    }
    if (not BTab[1].empty()){
      FIN >> val;
      BTab[1][i] = val;
    }
    if (not BTab[2].empty()){
      FIN >> val;
      BTab[2][i] = val;
    }
    if (not VTab.empty()){
      FIN >> val;
      VTab[i] = val;
    }
    FIN >> std::ws;
    ++i;
  }

  std::cout << "\n";
  if (i != xl*yl*zl){
    throw std::runtime_error((boost::format("The header says the size is %1%, actually it is %2%! Exiting...\n") % (xl*yl*zl) % i).str());
  }

  auto xminmax = std::minmax_element(xyzTab[0].begin(), xyzTab[0].end());
  auto yminmax = std::minmax_element(xyzTab[1].begin(), xyzTab[1].end());
  auto zminmax = std::minmax_element(xyzTab[2].begin(), xyzTab[2].end());
  return TFieldContainer(std::unique_ptr<TabField3>(new TabField3(xyzTab, BTab, VTab, std::array<std::array<std::vector<double>, 7>, 3>())), Bscale, Escale, *xminmax.second, *xminmax.first, *yminmax.second, *yminmax.first, *zminmax.second, *zminmax.first, BoundaryWidth);
}



TFieldContainer ReadTricubicField(const std::string &params, const std::map<std::string, std::string> &formulas){
  std::istringstream ss(params);
  boost::filesystem::path ft;
  std::string fieldtype, Bscale;
  double BoundaryWidth, lengthconv;
  ss >> fieldtype >> ft >> Bscale >> BoundaryWidth >> lengthconv; // read fieldtype, tablefilename, and rest of parameters
  Bscale = ResolveFormula(Bscale, formulas);
  if (!ss){
    throw std::runtime_error((boost::format("Could not read all required parameters for field %1%!") % fieldtype).str());
  }
  ft = boost::filesystem::absolute(ft, configpath.parent_path());

  std::string line;
  std::vector<std::string> line_parts;
  std::vector<double> x, y, z;
  std::vector<double> bx, by, bz;
  std::vector<double> dBx_dx, dBx_dy, dBx_dz, dBx_dxdy, dBx_dxdz, dBx_dydz, dBx_dxdydz;
  std::vector<double> dBy_dx, dBy_dy, dBy_dz, dBy_dxdy, dBy_dxdz, dBy_dydz, dBy_dxdydz;
  std::vector<double> dBz_dx, dBz_dy, dBz_dz, dBz_dxdy, dBz_dxdz, dBz_dydz, dBz_dxdydz;
  
  std::ifstream FINstream(ft.string(), std::ifstream::in);
  boost::iostreams::filtering_istream FIN;
  if (ft.extension() == ".bz2"){
    FIN.push(boost::iostreams::bzip2_decompressor());
  }
  else if (ft.extension() == ".gz"){
    FIN.push(boost::iostreams::gzip_decompressor());
  }
  FIN.push(FINstream);
  if (!FINstream.is_open() or !FIN.is_complete()){
    throw std::runtime_error("Could not open " + ft.string());
  }
  std::cout << "\nReading " << ft << "\n";

  // Read in file data
  int lineNum = 0;

  while (getline(FIN,line)){
    lineNum++;
    if (line.substr(0,1) == "%" || line.substr(0,1) == "#") continue;     // Skip commented lines
    boost::split(line_parts, line, boost::is_any_of("\t, "), boost::token_compress_on); //Delineate tab, space, commas

    if (line_parts.size() != 21){
      throw std::runtime_error((boost::format("Error reading line %1% of file %2%") % lineNum % ft.string()).str());
    }

    x.push_back( std::stod(line_parts[0], nullptr) * lengthconv);
    y.push_back( std::stod(line_parts[1], nullptr) * lengthconv);
    z.push_back( std::stod(line_parts[2], nullptr) * lengthconv);

    bx.push_back( std::stod(line_parts[3], nullptr));
    by.push_back( std::stod(line_parts[4], nullptr));
    bz.push_back( std::stod(line_parts[5], nullptr));

    dBx_dx.push_back( std::stod(line_parts[6], nullptr));
    dBx_dy.push_back( std::stod(line_parts[7], nullptr));
    dBx_dz.push_back( std::stod(line_parts[8], nullptr));
    dBx_dxdy.push_back( std::stod(line_parts[9], nullptr));
    dBx_dxdz.push_back( std::stod(line_parts[10], nullptr));
    dBx_dydz.push_back( std::stod(line_parts[11], nullptr));
    dBx_dxdydz.push_back( std::stod(line_parts[12], nullptr));
    
    dBy_dx.push_back( std::stod(line_parts[7], nullptr)); // = dBx_dy
    dBy_dy.push_back( std::stod(line_parts[13], nullptr));
    dBy_dz.push_back( std::stod(line_parts[14], nullptr));
    dBy_dxdy.push_back( std::stod(line_parts[15], nullptr));
    dBy_dxdz.push_back( std::stod(line_parts[11], nullptr)); // = dBx_dydz
    dBy_dydz.push_back( std::stod(line_parts[16], nullptr));
    dBy_dxdydz.push_back( std::stod(line_parts[17], nullptr));
    
    dBz_dx.push_back( std::stod(line_parts[8], nullptr)); // = dBx_dz
    dBz_dy.push_back( std::stod(line_parts[14], nullptr)); // = dBy_dz
    dBz_dz.push_back( -std::stod(line_parts[6], nullptr) -std::stod(line_parts[13], nullptr)); // = - dBx_dx - dBy_dy
    dBz_dxdy.push_back( std::stod(line_parts[11], nullptr)); // = dBx_dydz
    dBz_dxdz.push_back( std::stod(line_parts[18], nullptr));
    dBz_dydz.push_back( std::stod(line_parts[19], nullptr));
    dBz_dxdydz.push_back( std::stod(line_parts[20], nullptr));
    
  }
  
  if (x.empty() || y.empty() || z.empty() || bx.empty() || by.empty()|| bz.empty() ) {
    throw std::runtime_error("No data read from " + ft.string());
  }
  
  const std::array<std::vector<double>, 7> dBxTab={dBx_dx, dBx_dy, dBx_dz, dBx_dxdy, dBx_dxdz, dBx_dydz, dBx_dxdydz};
  const std::array<std::vector<double>, 7> dByTab={dBy_dx, dBy_dy, dBy_dz, dBy_dxdy, dBy_dxdz, dBy_dydz, dBy_dxdydz};
  const std::array<std::vector<double>, 7> dBzTab={dBz_dx, dBz_dy, dBz_dz, dBz_dxdy, dBz_dxdz, dBz_dydz, dBz_dxdydz};
    
  auto xminmax = std::minmax_element(x.begin(), x.end());
  auto yminmax = std::minmax_element(y.begin(), y.end());
  auto zminmax = std::minmax_element(z.begin(), z.end());
  return TFieldContainer(std::unique_ptr<TabField3>(new TabField3({x,y,z}, {bx,by,bz}, std::vector<double>(), {dBxTab, dByTab, dBzTab})), Bscale, "0", *xminmax.second, *xminmax.first, *yminmax.second, *yminmax.first, *zminmax.second, *zminmax.first, BoundaryWidth);
}




void TabField3::CheckTab(const std::array<std::vector<double>, 3> &B, const std::vector<double> &V){
  //  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
  std::cout << "The arrays are " << xyz[0].size() << " by " << xyz[1].size() << " by " << xyz[2].size()
	    << " (" << B[0].size() << ", " << B[1].size() << ", " << B[2].size() << ", " << V.size() << " field components).\n";
  std::cout << "The x values go from " << xyz[0].front() << " to " << xyz[0].back() << " [m]\n";
  std::cout << "The y values go from " << xyz[1].front() << " to " << xyz[1].back() << " [m]\n";
  std::cout << "The z values go from " << xyz[2].front() << " to " << xyz[2].back() << " [m].\n";

  std::vector<double> Babs;
  std::transform( boost::make_zip_iterator(boost::make_tuple(B[0].begin(), B[1].begin(), B[2].begin())),
		  boost::make_zip_iterator(boost::make_tuple(B[0].end(), B[1].end(), B[2].end())), std::back_inserter(Babs),
		  [](const boost::tuple<double,double,double> &Bs){ return std::sqrt(std::pow(Bs.get<0>(), 2) + std::pow(Bs.get<1>(), 2) + std::pow(Bs.get<2>(), 2)); });
  std::cout << "The input table file has values of magnetic field |B| from " << *std::min_element(Babs.begin(), Babs.end()) << " to " << *std::max_element(Babs.begin(), Babs.end());
  if (not V.empty())
    std::cout << " and values of electric potential from " << *std::min_element(V.begin(), V.end()) << " to " << *std::max_element(V.begin(), V.end());
  std::cout << std::endl;

  unsigned long l = xyz[0].size() * xyz[1].size() * xyz[2].size();
  if ((not B[0].empty() && B[0].size() != l) || (not B[1].empty() && B[1].size() != l) || (not B[2].empty() && B[2].size() != l) || (not V.empty() && V.size() != l))
    std::cout << "Warning: Number of field samples does not match number of grid points! Missing points will default to zero field\n";
}


void TabField3::CalcDerivs(const array3D &Tab, const unsigned long diff_dim, array3D &DiffTab) const
{
  alglib::real_1d_array x, y, diff;
  std::array<unsigned long, 3> index;
  DiffTab.resize(boost::extents[Tab.shape()[0]][Tab.shape()[1]][Tab.shape()[2]]);

  // std::cout << Tab.shape()[0] << ", " << Tab.shape()[1] << ", " << Tab.shape()[2] << std::endl;
  
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


void TabField3::PreInterpol(const array3D &Tab, const std::array<array3D, 7> &dTab, field_type &coeff) const{
  std::array<unsigned long, 3> len;
  for (unsigned long i = 0; i < 3; ++i){
    len[i] = xyz[i].size() - 1;
  }
  coeff.resize(len);

  array3D dFdx, dFdy, dFdz, dFdxdy, dFdxdz, dFdydz, dFdxdydz; // derivatives with respect to x, y, z, xy, xz, yz, xyz

  if (dTab[0].empty()){    	  	
    CalcDerivs(Tab, 0, dFdx); // dF/dx
    CalcDerivs(Tab, 1, dFdy); // dF/dy
    CalcDerivs(Tab, 2, dFdz); // dF/dz
    CalcDerivs(dFdx, 1, dFdxdy); // d2F/dxdy
    CalcDerivs(dFdx, 2, dFdxdz); // d2F/dxdz
    CalcDerivs(dFdy, 2, dFdydz); // d2F/dydz
    CalcDerivs(dFdxdy, 2, dFdxdydz); // d3F/dxdydz
  }
  // else {
  //   dFdx = dTab[0];
  //   dFdy = dTab[1];
  //   dFdz = dTab[2];
  //   dFdxdy = dTab[3];
  //   dFdxdz = dTab[4];
  //   dFdydz = dTab[5];
  //   dFdxdydz = dTab[6];
  // }

  // if (!dTab[0].empty()){
  //   // // print partial derivatives infos (sly)
  //   // long unsigned int ixx=100, iyy=20, izz=25;
  //   long unsigned int ixx=30, iyy=18, izz=18;
  //   std::array<unsigned long, 3> myindices = {ixx, iyy, izz}; // collect indices of corners of each grid cell

  //   double mycellx = xyz[0][ixx+1] - xyz[0][ixx];
  //   double mycelly = xyz[1][iyy+1] - xyz[1][iyy];
  //   double mycellz = xyz[2][izz+1] - xyz[2][izz];

  //   std::cout << "\nx = " <<  xyz[0][ixx] <<" y = " << xyz[1][iyy] << " z = " << xyz[2][izz] << std::endl;
  //   std::cout << "Tab = " <<  Tab(myindices) << std::endl;
  //   std::cout << "dFdx = " <<  dFdx(myindices)*mycellx << " | dTab = "<< dTab[0](myindices)*mycellx << std::endl;
  //   std::cout << "dFdy = "   <<  dFdy(myindices)*mycelly << " | dTab = "<< dTab[1](myindices)*mycelly << std::endl;
  //   std::cout << "dFdz = "   <<  dFdz(myindices)*mycellz <<  " | dTab = "<< dTab[2](myindices)*mycellz << std::endl;
  //   std::cout << "dFdxdy = " <<  dFdxdy(myindices)*mycellx*mycelly << " | dTab = "<< dTab[3](myindices)*mycellx*mycelly << std::endl;
  //   std::cout << "dFdxdz = " <<  dFdxdz(myindices)*mycellx*mycellz << " | dTab = "<< dTab[4](myindices)*mycellx*mycellz << std::endl;
  //   std::cout << "dFdydz = " <<  dFdydz(myindices)*mycelly*mycellz << " | dTab = "<< dTab[5](myindices)*mycelly*mycellz << std::endl;
  //   std::cout << "dFdxdydz = " <<  dFdxdydz(myindices)*mycellx*mycelly*mycellz << " | dTab = "<< dTab[6](myindices)*mycellx*mycelly*mycellz << std::endl;
  // }

  
  // // Pre interpol
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

	  if (dTab[0].empty()){
	    // std::cout << "dTab empty" << std::endl;
	    yyy[0][i] = Tab(indices[i]); // get values and derivatives at each corner of grid cell
	    yyy[1][i] = dFdx(indices[i])*cellx;
	    yyy[2][i] = dFdy(indices[i])*celly;
	    yyy[3][i] = dFdz(indices[i])*cellz;
	    yyy[4][i] = dFdxdy(indices[i])*cellx*celly;
	    yyy[5][i] = dFdxdz(indices[i])*cellx*cellz;
	    yyy[6][i] = dFdydz(indices[i])*celly*cellz;
	    yyy[7][i] = dFdxdydz(indices[i])*cellx*celly*cellz;
	  }
	  else{
	    yyy[0][i] = Tab(indices[i]); // get values and derivatives at each corner of grid cel
	    yyy[1][i] = dTab[0](indices[i])*cellx;
	    yyy[2][i] = dTab[1](indices[i])*celly;
	    yyy[3][i] = dTab[2](indices[i])*cellz;
	    yyy[4][i] = dTab[3](indices[i])*cellx*celly;
	    yyy[5][i] = dTab[4](indices[i])*cellx*cellz;
	    yyy[6][i] = dTab[5](indices[i])*celly*cellz;
	    yyy[7][i] = dTab[6](indices[i])*cellx*celly*cellz;
	  }
	  // std::cout << "\ndFdydz -> " << dFdxdy(indices[i])*celly*cellz << " dTab[3] -> " << dTab[3](indices[i]) << std::endl; 
	  // std::cout << "dFdz -> " << dFdz(indices[i])*cellz << " dTab[2] -> " <<  dTab[2](indices[i]) <<  std::endl; 
	}		
	tricubic_get_coeff(&coeff(indices[0])[0], &yyy[0][0], &yyy[1][0], &yyy[2][0], &yyy[3][0], &yyy[4][0], &yyy[5][0], &yyy[6][0], &yyy[7][0]); // calculate tricubic interpolation coefficients and store in coeff
      }
    }
  }

}


TabField3::TabField3(const std::array<std::vector<double>, 3> &xyzTab, const std::array<std::vector<double>, 3> &BTab, const std::vector<double> &VTab, const std::array<std::array<std::vector<double>, 7>, 3> &dBTab){


  for (unsigned i = 0; i < 3; ++i){
    std::unique_copy(xyzTab[i].begin(), xyzTab[i].end(), std::back_inserter(xyz[i])); // get list of unique x, y, and z coordinates
    std::sort(xyz[i].begin(), xyz[i].end());
    auto last = std::unique(xyz[i].begin(), xyz[i].end());
    xyz[i].erase(last, xyz[i].end());
  }
  CheckTab(BTab,VTab); // print some info


  std::array<array3D, 3> B;
  std::array<std::array<array3D, 7>, 3> dB;
  array3D V;

  std::cout << " Initializing BTab  " << std::endl;

  for (unsigned j = 0; j < 3; ++j){
    if (not BTab[j].empty()){
      B[j].resize(boost::extents[xyz[0].size()][xyz[1].size()][xyz[2].size()]);
      std::fill_n(B[j].data(), B[j].num_elements(), 0.);
      for (unsigned k = 0; k < 7; ++k){
	if (not dBTab[j][k].empty()){
	  dB[j][k].resize(boost::extents[xyz[0].size()][xyz[1].size()][xyz[2].size()]);
	  // std::cout << " dBTab not empty  " << std::endl;
	  std::fill_n(dB[j][k].data(), dB[j][k].num_elements(), 0.);
	}
      }
    }
  }
  if (not VTab.empty()){
    V.resize(boost::extents[xyz[0].size()][xyz[1].size()][xyz[2].size()]);
    std::fill_n(V.data(), V.num_elements(), 0.);
  }

  std::cout << " Indexing BTab  " << std::endl;

  for (unsigned long i = 0; i < xyzTab[0].size(); ++i){
    std::array<unsigned long, 3> index;
    for (unsigned j = 0; j < 3; ++j){
      auto found = std::lower_bound(xyz[j].begin(), xyz[j].end(), xyzTab[j][i]);
      assert(found != xyz[j].end());
      index[j] = std::distance(xyz[j].begin(), found);
    }
    for (unsigned j = 0; j < 3; ++j){
      if (not BTab[j].empty()){
	B[j](index) = BTab[j][i];
	for (unsigned k = 0; k < 7; ++k){
	  if (not dBTab[j][k].empty()){
	    dB[j][k](index) = dBTab[j][k][i]; // dB
	    // std::cout << " dB[j][k](index)  " << dB[j][k](index) << std::endl;
	    // std::cout << "x = " << xyz[0][index[0]] << ", y = " << xyz[1][index[1]] << ", z = " << xyz[2][index[2]] << " dB[j][k](index)  " << dB[j][k](index) << std::endl;
	  }
	}
	// std::cout << "xyzTab[" << j << "][" << i << "]=" << xyzTab[j][i] << " and " << "BTab[" << j<< "][" << i << "]=" << BTab[j][i] <<  std::endl;
      }
    }
    // std::cout << "xyzTab[" << 0 << "][" << i << "]=" << xyzTab[0][i] << " and " << "BTab[" << 0 << "][" << i << "]=" << BTab[0][i] <<  std::endl;

    if (not VTab.empty())
      V(index) = VTab[i];
  }

  // print some info on dB table
  // long unsigned int ixx=10, iyy=2, izz=5;
  // std::array<unsigned long, 3> myindices = {ixx, iyy, izz}; // collect indices of corners of each grid cell

  // std::cout<<"x= "<<xyz[0][ixx]<<", y = "<<xyz[1][iyy]<<", z = "<<xyz[2][izz]<<" dB[j][k](index)="<<dB[0][0](myindices)<<std::endl;

  // ixx=0, iyy=0, izz=0;
  // myindices = {ixx, iyy, izz}; // collect indices of corners of each grid cell

  // std::cout<<"x= "<<xyz[0][ixx]<<", y = "<<xyz[1][iyy]<<", z = "<<xyz[2][izz]<<" dB[j][k](index)="<<dB[0][0](myindices)<<std::endl;


  // for(unsigned j = 0; j < 3; ++j){
  //   std::cout << "\nx = " <<  xyz[0][ixx] <<" y = " << xyz[1][iyy] << " z = " << xyz[2][izz] << std::endl;
  //   std::cout << "Tab = " <<  B[j](myindices) << std::endl;
  //   std::cout << "dB[j][0] = "<< dB[j][0](myindices) << std::endl;
  //   std::cout << "dB[j][1] = "<< dB[j][1](myindices) << std::endl;
  //   std::cout << "dB[j][2] = "<< dB[j][2](myindices) << std::endl;
  //   std::cout << "dB[j][3] = "<< dB[j][3](myindices) << std::endl;
  //   std::cout << "dB[j][4] = "<< dB[j][4](myindices) << std::endl;
  //   std::cout << "dB[j][5] = "<< dB[j][5](myindices) << std::endl;
  //   std::cout << "dB[j][6] = "<< dB[j][6](myindices) << std::endl;
  // }
  
  
  std::cout << "Starting Preinterpolation ... ";
  float size = 0;
  if (not BTab[0].empty()){
    std::cout << "Bx ... ";
    std::cout.flush();
    PreInterpol(B[0], dB[0], Bc[0]); // precalculate interpolation coefficients for B field
    size += float(Bc[0].num_elements()*64*sizeof(double)/1024/1024);
  }
  if (not BTab[1].empty()){
    std::cout << "By ... ";
    std::cout.flush();
    PreInterpol(B[1], dB[1], Bc[1]);
    size += float(Bc[1].num_elements()*64*sizeof(double)/1024/1024);
  }
  if (not BTab[2].empty()){
    std::cout << "Bz ... ";
    std::cout.flush();
    PreInterpol(B[2], dB[2], Bc[2]);
    size += float(Bc[2].num_elements()*64*sizeof(double)/1024/1024);
  }
  if (not VTab.empty()){
    std::cout << "V ... ";
    std::cout.flush();
    PreInterpol(V, std::array<array3D, 7>(), Vc);
    size += float(Vc.num_elements()*64*sizeof(double)/1024/1024);
  }
  std::cout << "Done (" << size << " MB)\n";
}


void TabField3::Interpolate(const double x, const double y, const double z,
                            const field_type& coeffs, double &F, double dFdxi[3]) const {
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
  if (not coeffs.empty()){
    double *coeff = const_cast<double*>(&coeffs(index)[0]);
    F = tricubic_eval_fast(coeff, r[0], r[1], r[2]);
    if (dFdxi != nullptr){
      dFdxi[0] = tricubic_eval_fast(coeff, r[0], r[1], r[2], 1, 0, 0)/dist[0];
      dFdxi[1] = tricubic_eval_fast(coeff, r[0], r[1], r[2], 0, 1, 0)/dist[1];
      dFdxi[2] = tricubic_eval_fast(coeff, r[0], r[1], r[2], 0, 0, 1)/dist[2];
    }
  }
}


void TabField3::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  for (unsigned i = 0; i < 3; ++i){
    Interpolate(x, y, z, Bc[i], B[i], dBidxj == nullptr ? nullptr : dBidxj[i]);
  }

  // if(dBidxj != nullptr){
  //   std::cout << " \ndBidxj" << std::endl;
  //   for(int i=0; i<3;++i){
  //     for(int j=0; j<3; ++j){
  // 	std::cout << dBidxj[i][j] << ", ";
  //     }
  //     std::cout << std::endl;
  //   }  
  // }
  
}

void TabField3::EField(const double x, const double y, const double z, const double t,
		       double &V, double Ei[3]) const{
  double dVdxi[3];
  Interpolate(x, y, z, Vc, V, dVdxi);
  for (int i = 0; i < 3; i++){
    Ei[i] = -dVdxi[i]; // Ei = -dV/dxi
  }
}
