/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <random>

typedef std::mt19937_64 TMCGenerator; ///< typedef to default random-number generator

namespace std{

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
 * Generate random numbers between 0 and pi/2 following the Beckmann distribution with width << 1
 * 
 * See B. Walter, S.R. Marschner, H. Li, K.E. Torrance: "Microfacet Models for Refraction through Rough Surfaces", Eurographics Symposium on Rendering (2007)
 * https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
 * Satisfies the concept RandomNumberDistribution of STL
 */
template<typename T>
class beckmann_distribution{
public:
	typedef T result_type; ///< type returned by operator()
	typedef T param_type; ///< type of parameter (width of the distribution)
private:
	param_type _p; ///< member containing parameter
public:
	beckmann_distribution(){ reset(); } ///< empty constructor, calls reset()
	beckmann_distribution(const param_type &p){ param(p); } ///< construct with specific width parameter
	void reset(){ _p = static_cast<param_type>(0); } ///< reset to default state (0 width)
	param_type param() const { return _p; } ///< returns stored parameter
	void param(const param_type &p){ _p = p; } ///< set stored parameter

	/**
	 * Return random numbers distributed according to Beckmann distribution with width p << 1
	 * 
	 * See B. Walter, S.R. Marschner, H. Li, K.E. Torrance: "Microfacet Models for Refraction through Rough Surfaces", Eurographics Symposium on Rendering (2007)
	 * https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
	 */
	template<class Random> result_type operator()(Random &r, const param_type &p) const{
		std::uniform_real_distribution<param_type> d;
		return std::atan(std::sqrt(-_p*_p * std::log(1 - d()))); // x = atan(sqrt(-p^2 log(1 - u)))
	}
	template<class Random> result_type operator()(Random &r) const {return operator()(r, _p); } ///< Return random numbers distributed according to Beckmann distribution with previously set width
	result_type min() const { return static_cast<result_type>(0); } ///< return min random value (0)
	result_type max() const { return static_cast<result_type>(std::asin(1)); } ///< return max random value (pi/2)
	bool operator==(const polarization_distribution<T> &rhs) const { return _p == rhs._p; } ///< equality operator (compares internal parameters)
	bool operator!=(const polarization_distribution<T> &rhs) const { return !(operator==(rhs)); } ///< inequality operator (compares internal parameters)
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
	inversetransform_sampler(const result_type min, const result_type max){	param(make_pair(min, max)); } ///< construct using min and max values
	void reset(){ param(make_pair(static_cast<result_type>(0), static_cast<result_type>(1))); } ///< reset to default parameter (min = 0, max = 1)
	param_type param() const { return _p; } ///< return internally stored parameter
	void param(const param_type &p){ _p = p; } ///< set internally stored parameter
	/**
	 * Return random number with between min and max values using inverse transform sampling of func
	 */
	template<class Random> result_type operator()(Random &r, const param_type &p) const {
		func f;
		std::uniform_real_distribution<result_type> unidist(static_cast<result_type>(0), static_cast<result_type>(1)); // generate random numbers between 0 and 1
		return f(unidist(r), p.first, p.second); // call inverse of cumulative-distribution function with random number, min and max, creating correctly distributed numbers
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
		return acos(sqrt(x*(cos(amax)*cos(amax) - cos(amin)*cos(amin)) + cos(amin)*cos(amin)));
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
		return pow(x*(pow(amax,3) - pow(amin,3)) + pow(amin,3),1.0/3.0);
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
		return sqrt(x*(amax*amax - amin*amin) + amin*amin);
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
		return pow((pow(amax, 1.5) - pow(amin, 1.5))*x + pow(amin, 1.5), 2.0/3.0);
	}
};
/**
 * Specialization of inversetransform_sampler to create sqrt-distributed random numbers
 */
template<typename T> using sqrt_distribution = inversetransform_sampler<T, inversetransform_sqrt<T> >;

} // end namespace std

/**
 * Create a piecewise linear distribution from a function with single parameter
 *
 * @param f Function returning a double using a single double parameter
 * @param range_min Lower range limit of distribution
 * @param range_max Upper range limit of distribution
 *
 * @return Return piecewise linear distribution
 */
template<typename UnaryFunction>
std::piecewise_linear_distribution<double> parse_distribution(UnaryFunction f, const double range_min, const double range_max);


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
std::piecewise_linear_distribution<double> parse_distribution(const std::string &func, const double range_min, const double range_max);


extern std::piecewise_linear_distribution<double> proton_beta_distribution; ///< Creates random numbers according to the energy spectrum of protons from free-neutron decay
extern std::piecewise_linear_distribution<double> electron_beta_distribution; ///< Creates random numbers according to the energy spectrum of electrons from free-neutron decay

#endif /*MC_H_*/
