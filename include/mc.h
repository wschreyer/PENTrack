/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <random>

/**
 * Class to generate random numbers in several distributions.
 */
typedef std::mt19937_64 TMCGenerator;

namespace std{

template<typename T>
class polarization_distribution{
public:
	typedef T result_type;
	typedef T param_type;
private:
	param_type _p;
public:
	polarization_distribution(){ reset(); }
	polarization_distribution(const param_type &p){ param(p); }
	void reset(){ _p = static_cast<T>(0.5); }
	param_type param() const { return _p; }
	void param(const param_type &p){ _p = p; }
	template<class Random> result_type operator()(Random &r, const param_type &p) const{
		std::bernoulli_distribution d(p);
		return d(r) ? 1 : -1;
	}
	template<class Random> result_type operator()(Random &r) const { return operator()(r, _p); }
	result_type min() const { return static_cast<T>(-1); }
	result_type max() const { return static_cast<T>(1); }
	bool operator==(const polarization_distribution<T> &rhs) const { return _p == rhs._p; }
	bool operator!=(const polarization_distribution<T> &rhs) const { return !(operator==(rhs)); }
};

template<typename T, typename func>
class inversetransform_sampler{
public:
	typedef T result_type;
	typedef std::pair<result_type, result_type > param_type;
private:
	param_type _p;
public:
	inversetransform_sampler(){ reset(); }
	inversetransform_sampler(const param_type &p){ param(p); }
	inversetransform_sampler(const T min, const T max){	param(make_pair(min, max)); }
	void reset(){ param(make_pair(static_cast<T>(0), static_cast<T>(1))); };
	param_type param() const { return _p; }
	void param(const param_type &p){ _p = p; }
	template<class Random> result_type operator()(Random &r, const param_type &p) const {
		func f;
		std::uniform_real_distribution<T> unidist(static_cast<T>(0), static_cast<T>(1));
		return f(unidist(r), p.first, p.second);
	}
	template<class Random> result_type operator()(Random &r) const { return operator()(r, _p);	}
	result_type min() const { return _p.first; }
	result_type max() const { return _p.second; }
	bool operator==(const inversetransform_sampler<T, func> &rhs) const { return _p == rhs._p; }
	bool operator!=(const inversetransform_sampler<T, func> &rhs) const { return !(operator==(rhs)); }
};

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

template<typename T>
struct inversetransform_sin{
	T operator()(const T x, const T amin, const T amax) const {
		return std::acos(std::cos(amin) - x * (std::cos(amin) - std::cos(amax)));
	}
};
template<typename T>
using sin_distribution = inversetransform_sampler<T, inversetransform_sin<T> >;

template<typename T>
struct inversetransform_sincos{
	T operator()(const T x, const T amin, const T amax) const {
		return acos(sqrt(x*(cos(amax)*cos(amax) - cos(amin)*cos(amin)) + cos(amin)*cos(amin)));
	}
};
template<typename T>
using sincos_distribution = inversetransform_sampler<T, inversetransform_sincos<T> >;

template<typename T>
struct inversetransform_parabolic{
	T operator()(const T x, const T amin, const T amax) const {
		return pow(x*(pow(amax,3) - pow(amin,3)) + pow(amin,3),1.0/3.0);
	}
};
template<typename T>
using parabolic_distribution = inversetransform_sampler<T, inversetransform_parabolic<T> >;

template<typename T>
struct inversetransform_linear{
	T operator()(const T &x, const T &amin, const T &amax) const {
		return sqrt(x*(amax*amax - amin*amin) + amin*amin);
	}
};
template<typename T>
using linear_distribution = inversetransform_sampler<T, inversetransform_linear<T> >;

template<typename T>
struct inversetransform_sqrt{
	T operator()(const T x, const T amin, const T amax) const {
		return pow((pow(amax, 1.5) - pow(amin, 1.5))*x + pow(amin, 1.5), 2.0/3.0);
	}
};
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


extern std::piecewise_linear_distribution<double> proton_beta_distribution;
extern std::piecewise_linear_distribution<double> electron_beta_distribution;

#endif /*MC_H_*/
