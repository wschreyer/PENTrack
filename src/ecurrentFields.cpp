#include <cmath>
#include <fstream>
#include <iostream>

#include "boost/format.hpp"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "boost/multi_array.hpp"

#include "globals.h"
#include "ecurrentFields.h"


TECurrentField::TECurrentField(const std::string sft, const std::string &_It) {

  // TCustomBField::TCustomBField(const std::string &_Bx, const std::string &_By, const std::string &_Bz){
  
  tvar = std::unique_ptr<double>(new double(0.0));
  exprtk::symbol_table<double> symbol_table;
  symbol_table.add_variable("t",*tvar);
  symbol_table.add_constants();
  exprtk::parser<double> parser;
  std::string expr{_It,};
  Iexpr.register_symbol_table(symbol_table);
  if (not parser.compile(expr, Iexpr)){
    throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing Custom I(t) Current formula '" + expr + "': " + parser.get_error(0).diagnostic);
  }

  
  boost::filesystem::path ft(sft);
  ft = boost::filesystem::absolute(ft, configpath.parent_path());

  std::string line;
  std::vector<std::string> line_parts;
  WireSegment segment;

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
    
    if (line_parts.size() < 3){
      throw std::runtime_error((boost::format("Error reading line %1% of file %2%") % lineNum % ft.string()).str());
    }
    
    segment.x1 = std::stod(line_parts[0], nullptr);
    segment.y1 = std::stod(line_parts[1], nullptr);
    segment.z1 = std::stod(line_parts[2], nullptr);
    segment.x2 = std::stod(line_parts[3], nullptr);
    segment.y2 = std::stod(line_parts[4], nullptr);
    segment.z2 = std::stod(line_parts[5], nullptr);
    
    wireSegments.push_back(segment);    
  }

  if (wireSegments.empty()) {
    throw std::runtime_error("No data read from " + ft.string());
  } 
}



void TECurrentField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {

  *tvar = t;
  double I =  Iexpr.value();
  
  for(int i=0; i<3; ++i){
    B[i] = 0.0;
    if (dBidxj != nullptr){
      for(int j=0; j<3; ++j){
	dBidxj[i][j] = 0.0;
      }
    }
  }
  
  double vorfaktor = mu0 * I / (4 * pi);
		
  for (const auto& segment : wireSegments) {

    double SW1x = segment.x1; ///< x coordinate of start point
    double SW1y = segment.y1; ///< y coordinate of start point
    double SW1z = segment.z1; ///< z coordinate of start point
    double SW2x = segment.x2; ///< x coordinate of end point
    double SW2y = segment.y2; ///< y coordinate of end point
    double SW2z = segment.z2; ///< z coordinate of end point

    double t1 = SW2z * SW2z;
    double t2 = z * z;
    double t3 = SW2y * SW2y;
    double t4 = SW2z * z;
    double t5 = 2 * t4;
    double t6 = y * y;
    double t7 = SW2y * y;
    double t8 = 2 * t7;
    double t10 = SW1x * SW1x;
    double t12 = -SW2y + y;
    double t13 = -SW2x + x;
    double t17 = -SW2z + z;
    double t23 = t17 * SW2z;
    double t29 = x * x;
    double t30 = SW2x * x;
    double t31 = 2 * t30;
    double t32 = SW2x * SW2x;
    double t34 = SW1y * SW1y;
    double t40 = 2 * y * t32;
    double t41 = y + SW2y;
    double t47 = SW2z * y;
    double t53 = SW1z * SW1z;
    double t55 = z * t32;
    double t56 = z + SW2z;
    double t59 = z * t3;
    double t62 = t29 + t6;
    double t68 = t4 + t7;
    double t74 = SW2y * z;
    double t78 = t10 * (t1 + t2 + t3 - t5 + t6 - t8) +
      SW1x * (-2 * SW1y * t13 * t12 - 2 * SW1z * t13 * t17 + 2 * SW2x * (t4 - t6 + t7 - t2) + 2 * (-t3 + t7 + t23) * x) +
      t34 * (t29 + t1 + t2 - t5 - t31 + t32) +
      SW1y * (-2 * SW1z * t12 * t17 - t40 + 2 * SW2x * t41 * x + 2 * SW2y * (-t2 + t4 - t29) + 2 * t17 * t47) +
      t53 * (t29 - t8 + t32 + t6 - t31 + t3) +
      2 * SW1z * (-t55 + SW2x * t56 * x - t59 + SW2y * t56 * y - t62 * SW2z) +
      t32 * (t2 + t6) - 2 * SW2x * t68 * x + t3 * (t2 + t29) - 2 * t74 * t47 + t62 * t1;
    double t79 = sqrt(t78);
    double t80 = 0.1e1 / t79;
    double t81 = t80 * vorfaktor;
    double t82 = SW1y - y;
    double t90 = sqrt(t32 - 2 * SW2x * SW1x + t10 + t3 - 2 * SW2y * SW1y + t34 + t1 - 2 * SW2z * SW1z + t53);
    double t91 = 0.1e1 / t90;
    double t92 = t91 * t82;
    double t93 = SW2z - SW1z;
    double t95 = SW1z - z;
    double t96 = t91 * t95;
    double t97 = SW2y - SW1y;
    double t99 = t93 * t92 - t97 * t96;
    double t100 = t99 * t99;
    double t101 = SW2x - SW1x;
    double t103 = SW1x - x;
    double t104 = t91 * t103;
    double t106 = t101 * t96 - t93 * t104;
    double t107 = t106 * t106;
    double t110 = t97 * t104 - t101 * t92;
    double t111 = t110 * t110;
    double t112 = t100 + t107 + t111;
    double t113 = sqrt(t112);
    double t114 = 0.1e1 / t113;
    double t115 = t32 - t31 + t29 + t3 - t8 + t6 + t1 - t5 + t2;
    double t116 = sqrt(t115);
    double t117 = 0.1e1 / t116;
    double t119 = t101 * t91;
    double t122 = t97 * t91;
    double t125 = t93 * t91;
    double t133 = t10 - 2 * SW1x * x + t29 + t34 - 2 * SW1y * y + t6 + t53 - 2 * SW1z * z + t2;
    double t134 = sqrt(t133);
    double t135 = 0.1e1 / t134;
    double t142 = -t119 * t13 * t117 - t122 * t12 * t117 - t125 * t17 * t117 - t119 * t103 * t135 - t122 * t82 * t135 - t125 * t95 * t135;
    double t143 = fabs(t142);
    double t144 = t143 * t114;
    double t146 = SW1z * t12;
    double t147 = -SW1y * t17 + t146 + t74 - t47;
    double t162 = x * SW2y;
    double t169 = SW2z * x;
    double t234 = SW1y * t13;
    double t296 = SW1z * t13;
    double t359 = (-SW1x * t17 + t296 + SW2x * z - t169) * vorfaktor;
    double t361 = t143 * t114 * t80;
    double t398 = -SW1x * t12 + t234 + y * SW2x - t162;
    B[0] += t147 * t144 * t81;
    B[1] += -t361 * t359;
    B[2] += t398 * t144 * t81;

    if (dBidxj != nullptr){
      double t151 = 0.1e1 / t79 / t78;
      double t153 = t114 * t151 * vorfaktor;
      double t154 = t147 * t143;
      double t156 = SW1z * t17;
      double t179 = 2 * SW1x * (-SW1y * t12 - t156 - t3 + t7 + t23) + 2 * t34 * t13 + SW1y * (2 * SW2x * t41 - 4 * t162) + 2 * t53 * t13 + SW1z * (2 * SW2x * t56 - 4 * t169) - 2 * SW2x * t68 + 2 * x * t3 + 2 * x * t1;
      double t184 = 0.1e1 / t113 / t112;
      double t185 = t184 * t81;
      double t186 = t91 * t106;
      double t188 = t91 * t110;
      double t190 = t93 * t186 - t97 * t188;
      double t194 = t114 * t81;
      double t195 = fabs(t142) / t142;
      double t197 = 0.1e1 / t116 / t115;
      double t198 = -t13 * t197;
      double t202 = t91 * t117;
      double t204 = -t12 * t197;
      double t208 = -t17 * t197;
      double t213 = 0.1e1 / t134 / t133;
      double t214 = t103 * t213;
      double t218 = t91 * t135;
      double t220 = t82 * t213;
      double t224 = t95 * t213;
      double t228 = -t13 * t119 * t198 - t101 * t202 - t13 * t122 * t204 - t13 * t125 * t208 - t103 * t119 * t214 + t101 * t218 - t103 * t122 * t220 - t103 * t125 * t224;
      double t229 = t228 * t195;
      double t257 = 2 * t10 * t12 + SW1x * (-2 * t234 + (-4 * y + 2 * SW2y) * SW2x + 2 * t162) + 2 * SW1y * (-t156 - t32 + t30 + t23) + 2 * t53 * t12 + SW1z * (2 * SW2y * t56 - 4 * t47) + t40 - 2 * SW2x * t162 - 2 * SW2y * t4 + 2 * y * t1;
      double t261 = t91 * t99;
      double t264 = -t93 * t261 + t101 * t188;
      double t288 = -t12 * t119 * t198 - t12 * t122 * t204 - t97 * t202 - t12 * t125 * t208 - t82 * t119 * t214 - t82 * t122 * t220 + t97 * t218 - t82 * t125 * t224;
      double t289 = t288 * t195;
      double t293 = -t93 * t144 * t81;
      double t300 = 2 * SW2z - 4 * z;
      double t319 = 2 * t10 * t17 + SW1x * (-2 * t296 + SW2x * t300 + 2 * t169) + 2 * t34 * t17 + SW1y * (-2 * t146 + SW2y * t300 + 2 * t47) + 2 * SW1z * (-t32 + t30 - t3 + t7) + 2 * t55 - 2 * SW2z * t30 + 2 * t59 - 2 * SW2y * t47;
      double t325 = t97 * t261 - t101 * t186;
      double t349 = -t17 * t119 * t198 - t17 * t122 * t204 - t17 * t125 * t208 - t93 * t202 - t95 * t119 * t214 - t95 * t122 * t220 - t95 * t125 * t224 + t93 * t218;
      double t350 = t349 * t195;
      double t363 = t151 * t359;
      double t367 = t80 * t359;
      double t368 = t143 * t184;
      double t372 = t195 * t114;
      double t386 = t361 * t101 * vorfaktor;
      double t401 = t398 * t143;
      dBidxj[0][0] += -t179 * t154 * t153 / 2 - t190 * t154 * t185 + t147 * t229 * t194;
      dBidxj[0][1] += -t257 * t154 * t153 / 2 - t264 * t154 * t185 + t147 * t289 * t194 + t293;
      dBidxj[0][2] += -t319 * t154 * t153 / 2 - t325 * t154 * t185 + t147 * t350 * t194 + t97 * t144 * t81;
      dBidxj[1][0] += -t293 + t179 * t144 * t363 / 2 + t190 * t368 * t367 - t228 * t372 * t367;
      dBidxj[1][1] += t257 * t144 * t363 / 2 + t264 * t368 * t367 - t288 * t372 * t367;
      dBidxj[1][2] += -t386 + t319 * t144 * t363 / 2 + t325 * t368 * t367 - t349 * t372 * t367;
      dBidxj[2][0] += -t179 * t401 * t153 / 2 - t190 * t401 * t185 + t398 * t229 * t194 - t97 * t144 * t81;
      dBidxj[2][1] += -t257 * t401 * t153 / 2 - t264 * t401 * t185 + t398 * t289 * t194 + t386;
      dBidxj[2][2] += -t319 * t401 * t153 / 2 - t325 * t401 * t185 + t398 * t350 * t194;
    }
  }
}



// void TECurrentField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {
  
//   *tvar = t;
//   double I =  Iexpr.value();
  
//   for(int i=0; i<3; ++i){
//     B[i] = 0.0;
//     if (dBidxj != nullptr){
//       for(int j=0; j<3; ++j){
// 	dBidxj[i][j] = 0.0;
//       }
//     }
//   }
 
//   for (const auto& segment : wireSegments) {
//     double x1 = segment.x1;
//     double y1 = segment.y1;
//     double z1 = segment.z1;
//     double x2 = segment.x2;
//     double y2 = segment.y2;
//     double z2 = segment.z2;

//     // Calculate the differences in coordinates
//     double dx1 = x - x1;
//     double dy1 = y - y1;
//     double dz1 = z - z1;

//     double dx2 = x - x2;
//     double dy2 = y - y2;
//     double dz2 = z - z2;

//     // Calculate the magnitudes of the position vectors
//     double r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
//     double r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

//     // Calculate the cross products and their magnitudes
//     double cross1 = (y2 - y1) * dz1 - (z2 - z1) * dy1;
//     double cross2 = (y2 - y1) * dz2 - (z2 - z1) * dy2;

//     double cross1_mag = std::sqrt(cross1 * cross1);
//     double cross2_mag = std::sqrt(cross2 * cross2);

//     // Calculate the magnetic field components
//     B[0] += (mu0 * I / (4.0 * M_PI)) * (((y2 - y1) * dz1 - (z2 - z1) * dy1) / (r1 * r1 * r1) - (mu0 * I / (4.0 * M_PI)) * ((y2 - y1) * dz2 - (z2 - z1) * dy2) / (r2 * r2 * r2));

//     B[1] += (mu0 * I / (4.0 * M_PI)) * (((z2 - z1) * dx1 - (x2 - x1) * dz1) / (r1 * r1 * r1) - ((z2 - z1) * dx2 - (x2 - x1) * dz2) / (r2 * r2 * r2));

//     B[2] += (mu0 * I / (4.0 * M_PI)) * (((x2 - x1) * dy1 - (y2 - y1) * dx1) / (r1 * r1 * r1) - ((x2 - x1) * dy2 - (y2 - y1) * dx2) / (r2 * r2 * r2));


//     // First derivatives
//     if (dBidxj != nullptr){
//       // Calculate the partial derivatives
//       for (int i = 0; i < 3; ++i) {
//         for (int j = 0; j < 3; ++j) {
// 	  double dridxj = 0.0;
// 	  if (i == 0 && j == 0)
// 	    dridxj = dx1;
// 	  else if (i == 0 && j == 1)
// 	    dridxj = dy1;
// 	  else if (i == 0 && j == 2)
// 	    dridxj = dz1;
// 	  else if (i == 1 && j == 0)
// 	    dridxj = dx2;
// 	  else if (i == 1 && j == 1)
// 	    dridxj = dy2;
// 	  else if (i == 1 && j == 2)
// 	    dridxj = dz2;

// 	  double dBidxj_val = (mu0 * I / (4.0 * M_PI)) * ((cross1_mag * dridxj / (r1 * r1 * r1 * r1 * r1 * r1)) - (cross2_mag * dridxj / (r2 * r2 * r2 * r2 * r2 * r2)));

// 	  dBidxj[i][j] += dBidxj_val;
//         }
//       }
//     }
//   }
// }




// void TECurrentField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {
//   *tvar = t;
//   double current =  Iexpr.value();

//   double bx = 0.0;
//   double by = 0.0;
//   double bz = 0.0;

//   double dbx_dxi[3] = {0.0};
//   double dby_dxi[3] = {0.0};
//   double dbz_dxi[3] = {0.0};  

//   for (const auto& segment : wireSegments) {
//     double x1 = segment.x1;
//     double y1 = segment.y1;
//     double z1 = segment.z1;
//     double x2 = segment.x2;
//     double y2 = segment.y2;
//     double z2 = segment.z2;

//     double dx = x2 - x1;
//     double dy = y2 - y1;
//     double dz = z2 - z1;

//     double r = std::sqrt(dx * dx + dy * dy + dz * dz);
//     double r3 = r * r * r;

//     double dlx = mu0 / (4 * pi) * current * (dy * (z - z1) - dz * (y - y1)) / r3;
//     double dly = mu0 / (4 * pi) * current * (dz * (x - x1) - dx * (z - z1)) / r3;
//     double dlz = mu0 / (4 * pi) * current * (dx * (y - y1) - dy * (x - x1)) / r3;

//     bx += dlx;
//     by += dly;
//     bz += dlz;

//     // First derivatives
//     if (dBidxj != nullptr){
//       double dr_dxi[3] = {0.0};
//       dr_dxi[0] = (x - x1) / r;
//       dr_dxi[1] = (y - y1) / r;
//       dr_dxi[2] = (z - z1) / r;

//       for (int i = 0; i < 3; i++) {
// 	dbx_dxi[i] += mu0 / (4 * pi) * current * (dy * dr_dxi[i] - dz * dr_dxi[(i + 1) % 3]) / r3;
// 	dby_dxi[i] += mu0 / (4 * pi) * current * (dz * dr_dxi[i] - dx * dr_dxi[(i + 1) % 3]) / r3;
// 	dbz_dxi[i] += mu0 / (4 * pi) * current * (dx * dr_dxi[i] - dy * dr_dxi[(i + 1) % 3]) / r3;
//       }
//     }
//   }

//   B[0] = bx;
//   B[1] = by;
//   B[2] = bz;

//   if (dBidxj != nullptr){

//     dBidxj[0][0] = dbx_dxi[0];
//     dBidxj[0][1] = dbx_dxi[1];
//     dBidxj[0][2] = dbx_dxi[2];

//     dBidxj[1][0] = dby_dxi[0];
//     dBidxj[1][1] = dby_dxi[1];
//     dBidxj[1][2] = dby_dxi[2];

//     dBidxj[2][0] = dbz_dxi[0];
//     dBidxj[2][1] = dbz_dxi[1];
//     dBidxj[2][2] = dbz_dxi[2];
//   }
// }



// void TECurrentField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {

//   *tvar = t;
//   double I =  Iexpr.value();
  
//   for(int i=0; i<3; ++i){
//     B[i] = 0;
//     if (dBidxj != nullptr){
//       for(int j=0; j<3; ++j){
// 	dBidxj[i][j] = 0;
//       }
//     }
//   }
  
//   // double vorfaktor = mu0 * I / (4 * pi);

//   for (const auto& segment : wireSegments) {

//     // Compute wire segment properties
//     double dx = segment.x2 - segment.x1;
//     double dy = segment.y2 - segment.y1;
//     double dz = segment.z2 - segment.z1;

//     double length = std::sqrt(dx * dx + dy * dy + dz * dz);
//     double invLength = 1.0 / length;

//     // Calculate wire segment orientation
//     double ux = dx * invLength;
//     double uy = dy * invLength;
//     double uz = dz * invLength;

//     // Calculate vector r from wire segment to observation point (x, y, z)
//     double rx = x - segment.x1;
//     double ry = y - segment.y1;
//     double rz = z - segment.z1;

//     // Calculate scalar projection of r onto wire segment
//     double projection = (rx * ux + ry * uy + rz * uz);

//     // Calculate vector projection of r onto wire segment
//     double px = projection * ux;
//     double py = projection * uy;
//     double pz = projection * uz;

//     // Calculate vector projection of r perpendicular to wire segment
//     double qx = rx - px;
//     double qy = ry - py;
//     double qz = rz - pz;

//     // Calculate distance from wire segment to observation point
//     double r = std::sqrt(qx * qx + qy * qy + qz * qz);

//     // Calculate magnetic field components
//     double factor = mu0 * I / (4.0 * pi * r * r * r);

//     B[0] += factor * (3.0 * qx * px - r * r * ux);
//     B[1] += factor * (3.0 * qy * py - r * r * uy);
//     B[2] += factor * (3.0 * qz * pz - r * r * uz);

//     // Calculate derivatives of B with respect to x, y, and z
//     if (dBidxj != nullptr){
//       double factor2 = mu0 * I / (4.0 * pi * r * r * r * r * r);

//       dBidxj[0][0] = factor2 * (3.0 * px * px - r * r + 2.0 * qx * qx);
//       dBidxj[0][1] = factor2 * (3.0 * px * py + 2.0 * qx * qy);
//       dBidxj[0][2] = factor2 * (3.0 * px * pz + 2.0 * qx * qz);
//       dBidxj[1][0] = factor2 * (3.0 * py * px + 2.0 * qy * qx);
//       dBidxj[1][1] = factor2 * (3.0 * py * py - r * r + 2.0 * qy * qy);
//       dBidxj[1][2] = factor2 * (3.0 * py * pz + 2.0 * qy * qz);
//       dBidxj[2][0] = factor2 * (3.0 * pz * px + 2.0 * qz * qx);
//       dBidxj[2][1] = factor2 * (3.0 * pz * py + 2.0 * qz * qy);
//       dBidxj[2][2] = factor2 * (3.0 * pz * pz - r * r + 2.0 * qz * qz);
//     }
//   }
// }

 
