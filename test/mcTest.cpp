/**
 * This file contains unit tests for Monte Carlo methods
 */

#include <numeric>

#include <boost/test/unit_test.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include "mc.h"

BOOST_AUTO_TEST_CASE(adaptive_interpolation_test){
    auto f = [](double x){ return std::exp(-std::pow(x-50,2)); };
    auto dist = interpolate_distribution(f, 0., 100.);
    auto dx = dist.intervals();
    std::vector<double> y2;
    std::transform(dx.begin(), dx.end(), std::back_inserter(y2), f);
    std::adjacent_difference(dx.begin(), dx.end(), dx.begin());
//    std::copy(dx.begin(), dx.end(), std::ostream_iterator<double>(std::cout, ","));
//    std::cout << '\n';
//    std::copy(y2.begin(), y2.end(), std::ostream_iterator<double>(std::cout, ","));
//    std::cout << '\n';
    std::adjacent_difference(y2.begin(), y2.end(), y2.begin(), std::plus<> {});
    auto integral = 0.5*std::inner_product(dx.begin() + 1, dx.end(), y2.begin() + 1, 0.);
    BOOST_CHECK_SMALL(integral - std::sqrt(M_PI), 1e-12);
}