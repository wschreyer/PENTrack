#include "mc.h"

#include <chrono>
#include <iostream>
#include "exprtk.hpp"

#include "globals.h"

const int PIECEWISE_LINEAR_DIST_INTERVALS = 1000;

//double TMCGenerator::NeutronSpectrum() const{
//		return SqrtDist(0, 300e-9);


////////////////// Todo : sly



// // Function to compute the cumulative distribution function (CDF) from a histogram
// std::vector<double> computeCDF(const Histogram& histogram) {
//   std::vector<double> cdf(histogram.counts.size());
//   std::partial_sum(histogram.counts.begin(), histogram.counts.end(), cdf.begin());
//   const double total = cdf.back();
//   std::transform(cdf.begin(), cdf.end(), cdf.begin(), [total](double count) {
//     return count / total;
//   });
//   return cdf;
// }

	
/////////////////////


	
template<typename UnaryFunction>
std::piecewise_linear_distribution<double> parse_distribution(UnaryFunction f, const double range_min, const double range_max){
  if (range_min > range_max)
    throw std::runtime_error("Invalid range used to parse distribution");
  double rmin = range_min;
  double rmax = range_max;
  int nw = PIECEWISE_LINEAR_DIST_INTERVALS;
  if (range_min == range_max){
    rmax = std::nextafter(range_max, std::numeric_limits<double>::max()); // use slightly larger double value for range_max, else distribution will default to range 0..1
    nw = 1;
  }

  return std::piecewise_linear_distribution<double>(nw, rmin, rmax, f);
}


std::piecewise_linear_distribution<double> parse_distribution(const std::string &func, const double range_min, const double range_max){
  double x;
  exprtk::symbol_table<double> symbol_table;
  symbol_table.add_variable("x", x);
  symbol_table.add_function("ProtonBetaSpectrum", ProtonBetaSpectrum);
  symbol_table.add_function("ElectronBetaSpectrum", ElectronBetaSpectrum);
  symbol_table.add_function("MaxwellBoltzSpectrum", MaxwellBoltzSpectrum);
  symbol_table.add_function("CustomSpectra", CustomSpectra); //for PSI spectrum // Utkarsh
  symbol_table.add_constants();

  exprtk::expression<double> expression;
  expression.register_symbol_table(symbol_table);

  exprtk::parser<double> parser;
  if (not parser.compile(func, expression)){
    throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing formula '" + func + "': " + parser.get_error(0).diagnostic);
  }
  return parse_distribution(
			    [&x, &expression](const double px){
			      x = px;
			      return expression.value();
			    },
			    range_min,
			    range_max
			    );
}

std::piecewise_linear_distribution<double> proton_beta_distribution = parse_distribution(ProtonBetaSpectrum, 0., 750.);
std::piecewise_linear_distribution<double> electron_beta_distribution = parse_distribution(ElectronBetaSpectrum, 0., 782000.);
