/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <random>

#include <boost/math/special_functions/erf.hpp>

#include "exprtk.hpp"

#include "globals.h"

typedef std::mt19937_64 TMCGenerator; ///< typedef to default random-number generator

/**
 * Generate polarization +1 or -1 weighted by the projection of the spin-vector onto an axis.
 *
 * Satisfies the concept RandomNumberDistribution of STL
 */
template<typename T>
class polarization_distribution{
public:
	typedef T result_type; ///< type returned by operator()
	typedef T param_type; ///< type of parameter (projection of spin onto axis)
private:
	param_type _p; ///< member containing parameter
public:
	polarization_distribution(){ reset(); } ///< empty constructor, calls reset()
	polarization_distribution(const param_type &p){ param(p); } ///< construct with specific projection parameter
	void reset(){ _p = static_cast<param_type>(0); } ///< reset to default state (projection == 0, +1 and -1 equally likely)
	param_type param() const { return _p; } ///< returns stored parameter
	void param(const param_type &p){ _p = p; } ///< set stored parameter

	/**
	 * Return polarization +1/-1 weighted by projection parameter p
	 */
	template<class Random> result_type operator()(Random &r, const param_type &p) const{
		std::bernoulli_distribution d(static_cast<param_type>(0.5)*(p + static_cast<param_type>(1)));
		return d(r) ? static_cast<result_type>(1) : static_cast<result_type>(-1); // return +1/-1 according to Bernoulli distribution with parameter 0.5*(p + 1)
	}
	template<class Random> result_type operator()(Random &r) const { return operator()(r, _p); } ///< Return polarization +1/-1 weighted by internally stored projection parameter
	result_type min() const { return static_cast<result_type>(-1); } ///< return min random value (-1)
	result_type max() const { return static_cast<result_type>(1); } ///< return max random value (+1)
	bool operator==(const polarization_distribution<T> &rhs) const { return _p == rhs._p; } ///< equality operator (compares internal parameters)
	bool operator!=(const polarization_distribution<T> &rhs) const { return !(operator==(rhs)); } ///< inequality operator (compares internal parameters)
};


/**
 * Sample microfacet slope in x direction from Beckmann distribution with given width,
 * taking into account visibility of microfacet normals to incoming direction [-sin(theta), 0, cos(theta)].
 * See E. Heitz, E. d'Eon: "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
 * https://dx.doi.org/10.1111/cgf.12417
*/
template<typename T> class beckmann_visible_x_slope_distribution{
public:
	typedef T result_type; ///< type returned by operator()
	typedef std::tuple<T, T> param_type; ///< type of parameter (pair of polar angle of incoming direction and Beckmann distribution width)
private:
	param_type _p; ///< member containing parameter
public:
	beckmann_visible_x_slope_distribution(){ reset(); } ///< empty constructor, calls reset()

	/**
	 * Constructor.
	 * 
	 * @param p Pair of parameters: polar angle theta of incoming direction and Beckmann distribution width
	*/
	beckmann_visible_x_slope_distribution(const param_type &p){ param(p); }
	beckmann_visible_x_slope_distribution(const T incomingAngle, const T beckmannWidth){ param(std::make_tuple(incomingAngle, beckmannWidth)); }
	void reset (){ _p = std::tuple<T, T>(0, 1); } ///< Reset to default state (incoming angle 0, distribution width 1)
	param_type param() const { return _p; } ///< Return stored parameters
	void param(const param_type &p){ _p = p; } ///< Set stored parameters

	/**
	 * Sample microfacet slope in x direction from Beckmann distribution with given width,
	 * taking into account visibility of microfacet normals to incoming direction [-sin(theta), 0, cos(theta)].
	 * See E. Heitz, E. d'Eon: "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
	 * https://dx.doi.org/10.1111/cgf.12417
	 * 
	 * @param rng Random number generator
	 * @param p Pair of parameters: polar angle theta of incoming direction and Beckmann distribution width
	 * @return Returns slope of microfacet in x direction
	*/
	template<class Random> result_type operator()(Random &rng, const param_type &p) const {
		T width = std::get<1>(p);
		if (width == 0){
			return 0;
		}

		T slope_x = 0;
		T incoming_polar_angle = std::get<0>(p);
		if (incoming_polar_angle == 0){
			// if incoming polar angle is 0 sample from normal distribution with standard deviation 1
			slope_x = std::normal_distribution<T>(0, M_SQRT1_2)(rng);
		}
		else{
			// by rescaling the incoming polar angle by the width we can sample from a Beckmann distribution with width 1
			incoming_polar_angle = std::atan(std::tan(incoming_polar_angle)*width);

			T a = static_cast<T>(1)/std::tan(incoming_polar_angle);
			T erf_a = std::erf(a);
			T exp_a2 = std::exp(-a*a);
			T G1 = static_cast<T>(2)/(1 + erf_a + exp_a2*M_2_SQRTPI/2/a);
			T C = 1 - G1 * erf_a;
			
			if (std::generate_canonical<T>(rng) < C){
				T p = exp_a2*M_2_SQRTPI/2 / (exp_a2*M_2_SQRTPI/2 + a*(1 - erf_a));

				if (std::generate_canonical<T>(rng) < p){
					slope_x = -std::sqrt(-std::log(std::generate_canonical<T>(rng)*exp_a2));
				}
				else{
					slope_x = boost::math::erf_inv(std::generate_canonical<T>(rng)*(1 - erf_a) - 1);
				}
			}
			else{
				slope_x = boost::math::erf_inv((2*std::generate_canonical<T>(rng) - 1)*erf_a);
				T p = (-slope_x/a + 1)/2;
				if (std::generate_canonical<T>(rng) > p){
					slope_x = -slope_x;
				}
			}
		}
		return slope_x*width; // rescale slope with width
	}

	template<class Random> result_type operator()(Random &rng, const T incomingAngle){
		return operator()(rng, std::make_tuple(incomingAngle, std::get<1>(_p)));
	}

	/**
	 * Sample microfacet slope in x direction from Beckmann distribution with stored width,
	 * taking into account visibility of microfacet normals to incoming direction [-sin(theta), 0, cos(theta)].
	 * See E. Heitz, E. d'Eon: "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
	 * https://dx.doi.org/10.1111/cgf.12417
	 * 
	 * @param rng Random number generator
	 * @return Returns slope of microfacet in x direction
	*/
	template<class Random> result_type operator()(Random &rng) const { return operator()(rng, _p); }

	result_type min() const { return static_cast<result_type>(0); } ///< Return minimum random value
	result_type max() const { return std::numeric_limits<result_type>::infinity; } ///< Return maximum random value
	bool operator==(const beckmann_visible_x_slope_distribution<T> &rhs) const { return _p == rhs._p; } ///< equality operator (compares internal parameters)
	bool operator!=(const beckmann_visible_x_slope_distribution<T> &rhs) const { return !(operator==(rhs)); } ///< inequality operator (compares internal parameters)
};


/**
 * Generate distribution of polar angles theta of microfacets following (Beckmann distribution)*cos(theta), where Beckmann distribution has a width << 1
 * 
 * See B. Walter, S.R. Marschner, H. Li, K.E. Torrance: "Microfacet Models for Refraction through Rough Surfaces", Eurographics Symposium on Rendering (2007)
 * https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
 * Satisfies the concept RandomNumberDistribution of STL
 */
template<typename T>
class beckmann_cos_distribution{
public:
	typedef T result_type; ///< type returned by operator()
	typedef T param_type; ///< type of parameter (width of the distribution)
private:
	param_type _p; ///< member containing parameter
public:
	beckmann_cos_distribution(){ reset(); } ///< empty constructor, calls reset()
	beckmann_cos_distribution(const param_type &p){ param(p); } ///< construct with specific width parameter
	void reset(){ _p = static_cast<param_type>(0); } ///< reset to default state (0 width)
	param_type param() const { return _p; } ///< returns stored parameter
	void param(const param_type &p){ _p = p; } ///< set stored parameter

	/**
	 * Generate distribution of polar angles theta of microfacets following (Beckmann distribution)*cos(theta), where Beckmann distribution has a width << 1
	 * 
	 * See B. Walter, S.R. Marschner, H. Li, K.E. Torrance: "Microfacet Models for Refraction through Rough Surfaces", Eurographics Symposium on Rendering (2007)
	 * https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
	 */
	template<class Random> result_type operator()(Random &r, const param_type &p) const{
		result_type u = std::generate_canonical<result_type>(r); // generate random number between 0 and 1
		return std::atan(std::sqrt(-p*p * std::log(1 - u))); // x = atan(sqrt(-p^2 log(1 - u)))
	}
	template<class Random> result_type operator()(Random &r) const {return operator()(r, _p); } ///< Return random numbers distributed according to Beckmann distribution with previously set width
	result_type min() const { return static_cast<result_type>(0); } ///< return min random value (0)
	result_type max() const { return static_cast<result_type>(std::asin(static_cast<param_type>(1))); } ///< return max random value (pi/2)
	bool operator==(const beckmann_cos_distribution<T> &rhs) const { return _p == rhs._p; } ///< equality operator (compares internal parameters)
	bool operator!=(const beckmann_cos_distribution<T> &rhs) const { return !(operator==(rhs)); } ///< inequality operator (compares internal parameters)
};


/**
 * Generate random numbers between two values using inverse transform sampling.
 *
 * Template parameter func expects function result_type(result_type x, result_type min, result_type max) returning the inverse of the cumulative distribution function.
 * Satisfies the concept RandomNumberDistribution of STL.
 */
template<typename T, typename func>
class inversetransform_sampler{
public:
	typedef T result_type; ///< type returned by operator()
	typedef std::pair<result_type, result_type > param_type; ///< type of parameter (pair of min and max values)
private:
	param_type _p; ///< internally stored parameter
public:
	inversetransform_sampler(){ reset(); } ///< default constructor, calls reset()
	inversetransform_sampler(const param_type &p){ param(p); } ///< construct using specific parameter
	inversetransform_sampler(const result_type min, const result_type max){	param(std::make_pair(min, max)); } ///< construct using min and max values
	void reset(){ param(std::make_pair(static_cast<result_type>(0), static_cast<result_type>(1))); } ///< reset to default parameter (min = 0, max = 1)
	param_type param() const { return _p; } ///< return internally stored parameter
	void param(const param_type &p){ _p = p; } ///< set internally stored parameter
	/**
	 * Return random number with between min and max values using inverse transform sampling of func
	 */
	template<class Random> result_type operator()(Random &r, const param_type &p) const {
		func f;
		result_type u = std::generate_canonical<T, std::numeric_limits<T>::digits>(r); // generate random number between 0 and 1
		return f(u, p.first, p.second); // call inverse of cumulative-distribution function with random number, min and max, creating correctly distributed numbers
	}
	template<class Random> result_type operator()(Random &r) const { return operator()(r, _p);	} ///< Return random number using internally stored min and max values
	result_type min() const { return _p.first; } ///< return internally stored min value
	result_type max() const { return _p.second; } ///< return internally stored max value
	bool operator==(const inversetransform_sampler<result_type, func> &rhs) const { return _p == rhs._p; } ///< comparison operator, compares internal min and max values
	bool operator!=(const inversetransform_sampler<result_type, func> &rhs) const { return !(operator==(rhs)); } ///< inequality operator, compares internal min and max values
};

/*
template<class CharT, class Traits, typename T, typename func>
std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits> &os, const inversetransform_sampler<T, func> &d){
	return os << "Distribution generating random numbers between " << d.min() << " and " << d.max() << " using inverse transform sampling\n";
}

template<class CharT, class Traits, typename T, typename func>
std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits> &is, inversetransform_sampler<T, func> &d){
	typename inversetransform_sampler<T, func>::param_type p;
	is >> p;
	if (is)
		d.param(p);
	return is;
}
*/

/**
 * Functor to create sine-distributed random numbers with inversetransform_sampler
 */
template<typename T>
struct inversetransform_sin{
	/**
	 * Inverse of cumulative distribution = acos(-x)
	 */
	T operator()(const T x, const T amin, const T amax) const {
		return std::acos(std::cos(amin) - x * (std::cos(amin) - std::cos(amax)));
	}
};
/**
 * Specialization of inversetransform_sampler to create sine-distributed random numbers
 */
template<typename T>
using sin_distribution = inversetransform_sampler<T, inversetransform_sin<T> >;


/**
 * Functor to create sin*cos-distributed random numbers with inversetransform_sampler
 */
template<typename T>
struct inversetransform_sincos{
	/**
	 * Inverse of cumulative distribution = acos(sqrt(x))
	 */
	T operator()(const T x, const T amin, const T amax) const {
		return std::acos(std::sqrt(x*(std::cos(amax)*std::cos(amax) - std::cos(amin)*std::cos(amin)) + std::cos(amin)*std::cos(amin)));
	}
};
/**
 * Specialization of inversetransform_sampler to create sin*cos-distributed random numbers
 */
template<typename T>
using sincos_distribution = inversetransform_sampler<T, inversetransform_sincos<T> >;


/**
 * Functor to create x^2-distributed random numbers with inversetransform_sampler
 */
template<typename T>
struct inversetransform_parabolic{
	/**
	 * Inverse of cumulative distribution = x^1/3
	 */
	T operator()(const T x, const T amin, const T amax) const {
		return std::pow(x*(std::pow(amax,3) - std::pow(amin,3)) + std::pow(amin,3),1.0/3.0);
	}
};
/**
 * Specialization of inversetransform_sampler to create x^2-distributed random numbers
 */
template<typename T>
using parabolic_distribution = inversetransform_sampler<T, inversetransform_parabolic<T> >;


/**
 * Functor to create linearly distributed random numbers with inversetransform_sampler
 */
template<typename T>
struct inversetransform_linear{
	/**
	 * Inverse of cumulative distribution = sqrt(x)
	 */
	T operator()(const T &x, const T &amin, const T &amax) const {
		return std::sqrt(x*(amax*amax - amin*amin) + amin*amin);
	}
};
/**
 * Specialization of inversetransform_sampler to create linearly distributed random numbers
 */
template<typename T>
using linear_distribution = inversetransform_sampler<T, inversetransform_linear<T> >;


/**
 * Functor to create sqrt-distributed random numbers with inversetransform_sampler
 */
template<typename T>
struct inversetransform_sqrt{
	/**
	 * Inverse of cumulative distirbution = x^2/3
	 */
	T operator()(const T x, const T amin, const T amax) const {
		return std::pow((std::pow(amax, 1.5) - std::pow(amin, 1.5))*x + std::pow(amin, 1.5), 2.0/3.0);
	}
};
/**
 * Specialization of inversetransform_sampler to create sqrt-distributed random numbers
 */
template<typename T> using sqrt_distribution = inversetransform_sampler<T, inversetransform_sqrt<T> >;


/**
 * Calculate points for linear interpolation of function that best approximates integral
 * 
 * Calculate integral of f between last value in abscissas and next_abscissa using trapezoidal rule.
 * Repeat after bisecting the integration range and compare result. If difference is above tolerance recursively bisect again.
 * 
 * @param f Function with single parameter, returning a single value of the same type
 * @param next_abscissa End of integration range
 * @param next_value Value of f at end of integration range
 * @param abscissas Container of increasing abscissas already calculated
 * @param values Container of values of f already calculated at each abscissa
 * @param tolerance Once the difference between integrals of f over full range and bisected range is lower than this tolerance the bisection is stopped
 * @param max_abscissa_difference Keep bisecting range at least until difference betwen abscissas is smaller than this value
 */
template<class RealType, class UnaryFunction, class RealTypeContainer>
void adaptive_linear_interpolation(UnaryFunction f, RealType next_abscissa, RealType next_value,
									RealTypeContainer &abscissas, RealTypeContainer &values, RealType tolerance, RealType max_abscissa_difference){
	RealType h = next_abscissa - abscissas.back();
	RealType I0 = h*(values.back() + next_value)/2.;
	RealType intermediate_abscissa = (abscissas.back() + next_abscissa)/2.;
	RealType intermediate_value = f(intermediate_abscissa);
	RealType I1 = (intermediate_abscissa - abscissas.back())*(values.back() + intermediate_value)/2.;
	RealType I2 = (next_abscissa - intermediate_abscissa)*(intermediate_value + next_value)/2.;
	if ((h > max_abscissa_difference) or 
	    ((std::abs(I1 + I2 - I0) > tolerance) and (h > max_abscissa_difference/100.))){
		adaptive_linear_interpolation(f, intermediate_abscissa, intermediate_value, abscissas, values, tolerance/2, max_abscissa_difference);
		adaptive_linear_interpolation(f, next_abscissa, next_value, abscissas, values, tolerance/2, max_abscissa_difference);
	}
	else{
		abscissas.push_back(intermediate_abscissa);
		values.push_back(intermediate_value);
		abscissas.push_back(next_abscissa);
		values.push_back(next_value);
	}
};

/**
 * Create a piecewise linear distribution from a function with single parameter
 *
 * @param f Function returning a double using a single double parameter
 * @param range_min Lower range limit of distribution
 * @param range_max Upper range limit of distribution
 * @param tolerance Bisect range until this tolerance is reached
 *
 * @return Return piecewise linear distribution
 */
template<class RealType, class UnaryFunction>
std::piecewise_linear_distribution<RealType> interpolate_distribution(UnaryFunction f, RealType range_min, RealType range_max, RealType tolerance = std::numeric_limits<RealType>::epsilon()*1000.){
	if (range_min > range_max)
		throw std::runtime_error("Invalid range used to parse distribution");
	std::vector<RealType> abscissas, values;
	abscissas.push_back(range_min);
	values.push_back(f(range_min));
	if (range_min == range_max){
		abscissas.push_back(std::nextafter(range_max, std::numeric_limits<RealType>::max())); // if range_min and range_max are equal use slightly larger double value for range_max, else distribution will default to range 0..1
		values.push_back(f(abscissas.back()));
	}
	else{
		adaptive_linear_interpolation(f, range_max, f(range_max), abscissas, values, tolerance, (range_max - range_min)/10.); // do linear interpolation of function with adaptive step size
	}

	return std::piecewise_linear_distribution<RealType>(abscissas.begin(), abscissas.end(), values.begin());

};


/**
 * Create a piecewise linear distribution from a function defined in a string.
 * The string is interpreted by the ExprTk library.
 *
 * @param func Formula string to parse
 * @param range_min Lower range limit of distribution
 * @param range_max Upper range limit of distribution
 *
 * @return Return piecewise linear distribution
 */
inline std::piecewise_linear_distribution<double> parse_distribution(std::string func, double range_min, double range_max){
	double x;
	exprtk::symbol_table<double> symbol_table;
	symbol_table.add_variable("x", x);
	symbol_table.add_function("ProtonBetaSpectrum", ProtonBetaSpectrum);
	symbol_table.add_function("ElectronBetaSpectrum", ElectronBetaSpectrum);
	symbol_table.add_function("MaxwellBoltzSpectrum", MaxwellBoltzSpectrum);
	symbol_table.add_constants();

	exprtk::expression<double> expression;
	expression.register_symbol_table(symbol_table);

	exprtk::parser<double> parser;
	if (not parser.compile(func, expression)){
	    throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing formula '" + func + "': " + parser.get_error(0).diagnostic);
	}
	return interpolate_distribution(
			[&x, &expression](double px){
				x = px;
				return expression.value();
			},
			range_min,
			range_max
	);
};


static std::piecewise_linear_distribution<double> proton_beta_distribution = interpolate_distribution(ProtonBetaSpectrum, 0., 750.); ///< Creates random numbers according to the energy spectrum of protons from free-neutron decay
static std::piecewise_linear_distribution<double> electron_beta_distribution = interpolate_distribution(ElectronBetaSpectrum, 0., 782000.); ///< Creates random numbers according to the energy spectrum of electrons from free-neutron decay

#endif /*MC_H_*/
