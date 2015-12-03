/**
 * \file
 * Do "brute force" integration of the Bloch equation.
 */

#ifndef BRUTEFORCE_H_
#define BRUTEFORCE_H_

#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <boost/numeric/odeint.hpp>

/**
 * Bloch equation integrator.
 *
 * Create this class to do "brute force" tracking of your particle's spin in magnetic fields along its track.
 * It starts tracking the spin by integrating the Bloch equation when the absolute magnetic field drops below TBFIntegrator::Bmax and stops it when the field rises above this value again.
 * Then it calculates the spin flip probability after such a low-field-pass.
 */
struct TBFIntegrator{
private:
	typedef double value_type; ///< define floating point type for spin integration
	typedef std::vector<value_type> state_type; ///< define type which contains spin state vector
	typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, value_type> stepper_type; ///< define integration stepper type
	typedef	boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type;
	typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type;
	state_type I_n; ///< Spin vector
//	stepper_type stepper;

	value_type gamma; ///< Particle's gyromagnetic ration
	std::string particlename; ///< Name of particle whose spin is to be tracked, needed for logging.
	double Bmax; ///< Spin tracking is only done when absolut magnetic field drops below this value.
	std::vector<double> BFtimes; ///< Pairs of absolute time in between which spin tracking shall be done.
	double BFBminmem; ///< Stores minimum field during one spin track for information
	bool spinlog; ///< Should the tracking be logged to file?
	double spinloginterval; ///< Time interval between log file entries.
	double spinlogtimeoutinterval; ///< Prints to the spin log only after the specified spinlogtimeoutinterval has passed
	long int intsteps; ///< Count integrator steps during spin tracking for information.
	std::ofstream &fspinout; ///< file to log into
	double starttime; ///< time of last integration start
	double prevTimeOut; ///< Time when the previous write to the spinlog occurred
	
	value_type t1; ///< field interpolation start time
	value_type t2; ///< field interpolation end time
	value_type cx[4]; ///< cubic spline coefficients for magnetic field x-component
	value_type cy[4]; ///< cubic spline coefficients for magnetic field y-component
	value_type cz[4]; ///< cubic spline coefficients for magnetic field z-component

public:
	/**
	 * Constructor.
	 *
	 * Set initial values and read options from config file.
	 *
	 * @param agamma Gyromagnetic ration of particle whose spin is to be tracked.
	 * @param aparticlename Particle name.
	 * @param conf Option map containing particle specific spin tracking options.
	 * @param spinout Stream to which spin track is written
	 */
	TBFIntegrator(double agamma, std::string aparticlename, std::map<std::string, std::string> &conf, std::ofstream &spinout);
private:
	/**
	 * Do cubic spline interpolation of magnetic field components with coefficients determined in TBFderivs::TBFderivs
	 *
	 * @param t Time
	 * @param B Magnetic field components
	 */
	void Binterp(value_type t, value_type B[3]);

public:
	/**
	 *  Bloch equation integrator calls TBFIntegrator(x,y,dydx) to get derivatives
	 *
	 *  @param y Spin vector
	 *  @param dydx Returns temporal derivative of spin vector
	 *  @param x Time
	 */
	void operator()(state_type y, state_type &dydx, value_type x);

	/**
	 * Integration observer
	 *
	 * Bloch equation integrator calls TBFIntegrator(y, x) on each integration step
	 *
	 * @param y Current state vector of the ODE system
	 * @param x Current time
	 */
	void operator()(const state_type &y, value_type x);

	/**
	 * Track spin between two particle track points.
	 *
	 * Checks if spin should be tracked between given track points. Returns spin flip probability when tracking is finished.
	 *
	 * @param x1 Time at first track point.
	 * @param y1 State vector at first track point.
	 * @param dy1dx Temporal derivative of state vector at first track point.
	 * @param B1 Magnetic field at first track point.
	 * @param E1 Electric field at first track point.
	 * @param x2 Time at second track point.
	 * @param y2 State vector at second track point.
	 * @param dy2dx Temporal derivative of state vector at second track point.
	 * @param B2 Magnetic field at second track point.
	 *
	 * @return Probability, that NO spin flip occured (usually close to 1).
	 */
	long double Integrate(double x1, double y1[6], double dy1dx[6], double B1[4][4], double E1[3], 
						double x2, double y2[6], double dy2dx[6], double B2[4][4], double E2[3]);

};

#endif // BRUTEFORCE_H_
