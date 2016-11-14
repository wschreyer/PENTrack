/**
 * \file
 * Particle base class definition.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <fstream>
#include <vector>
#include <map>

#include <boost/numeric/odeint.hpp>

#include "geometry.h"
#include "mc.h"
#include "fields.h"

static const double MAX_SAMPLE_DIST = 0.01; ///< max spatial distance of reflection checks, spin flip calculation, etc; longer integration steps will be interpolated
static const int STATE_VARIABLES = 8; ///< number of variables in trajectory integration (position, velocity, proper time, polarization)
static const int SPIN_STATE_VARIABLES = 5; ///< number of variables in spin integration (spin vector, time, total phase)

/**
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor and methods TParticle::OnHit, TParticle::OnStep and TParticle::Decay.
 * Optionally, derived particles can also re-implement TParticle::Epot, TParticle::PrintStartEnd, TParticle::PrintTrack, TParticle::PrintSnapshots, TParticle::PrintHits and define its own constructors.
 */
struct TParticle{
protected:
	typedef double value_type; ///< data type used for trajectory integration
	typedef std::vector<value_type> state_type; ///< type representing current particle state (position, velocity, proper time, and polarization) or spin state (x,y,z component)
	typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, value_type> stepper_type; ///< basic integration stepper (5th-order Runge-Kutta)
	typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type; ///< integration step length controller
	typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type; ///< integration step interpolator
//	typedef boost::numeric::odeint::bulirsch_stoer_dense_out<state_type, value_type> dense_stepper_type2; ///< alternative stepper type (Bulirsch-Stoer)

	/**
	 * Enum containing all types of log files.
	 *
	 * Is passed to GetLogStream to identify, which file stream should be returned.
	 */
	enum LogStream { endLog, ///< end log
					snapshotLog, ///< snapshot log
					hitLog, ///< hit log
					trackLog, ///< track log
					spinLog ///< spin log
					};
public:
	const char *name; ///< particle name (has to be initialized in all derived classes!)
	const long double q; ///< charge [C] (has to be initialized in all derived classes!)
	const long double m; ///< mass [eV/c^2] (has to be initialized in all derived classes!)
	const long double mu; ///< magnetic moment [J/T] (has to be initialized in all derived classes!)
	const long double gamma; ///< gyromagnetic ratio [rad/(s T)] (has to be initialized in all derived classes!)
	int particlenumber; ///< particle number
	stopID ID; ///< particle fate (defined in globals.h)
	double tau; ///< particle life time
	double maxtraj; ///< max. simulated trajectory length

	/// start time
	value_type tstart;

	/// stop time
	value_type tend;

	/// state vector before integration (position, velocity, proper time, and polarization)
	state_type ystart;

	/// state vector after integration (position, velocity, proper time, and polarization)
	state_type yend;

	/// spin vector before integration
	state_type spinstart;

	/// spin vector after integration
	state_type spinend;

	/// solid in which the particle started
	solid solidstart;

	/// solid in which particle stopped
	solid solidend;

	/// total start energy
	double Hstart() const;

	/// total end energy
	double Hend() const;

	/// max total energy
	double Hmax;

	/// kinetic start energy
	double Estart() const;

	/// kinetic end energy
	double Eend() const;

	/// trajectory length
	double lend;

	/// number of material boundary hits
	int Nhit;

	/// number of spin flips
	int Nspinflip;

	/// total probability of NO spinflip calculated by spin tracking
	long double noflipprob;

	/// number of integration steps
	int Nstep;
	
	std::vector<TParticle*> secondaries; ///< list of secondary particles

	/**
	 * Constructor, initializes TParticle::type, TParticle::q, TParticle::m, TParticle::mu
	 *
	 * Has to be called by every derived class constructor with the respective values
	 *
	 * @param aname Particle name
	 * @param qq Electric charge
	 * @param mm Mass
	 * @param mumu Magnetic dipole moment
	 * @param agamma Gyromagnetic ratio
	 * @param number Particle number
	 * @param t Starting time
	 * @param x Initial x coordinate
	 * @param y Initial y coordinate
	 * @param z Initial z coordinate
	 * @param E Initial kinetic energy
	 * @param phi Azimuth of initial velocity vector
	 * @param theta Polar angle of initial velocity vector
	 * @param polarisation polarisation of particle (-1..1)
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, int number,
			double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

	/**
	 * Destructor, deletes secondaries
	 */
	virtual ~TParticle();


	/**
	 * Integrate particle trajectory.
	 *
	 * Takes inital state vector ystart and integrates the trajectory step by step.
	 * If a step is longer than MAX_SAMPLE_DIST, the step is split by interpolating intermediate points.
	 * On each step it checks for interaction with solids, prints snapshots and track into files and calls TParticle::OnStep.
	 * TParticle::StopIntegration is called if TParticle::tau or tmax are reached; or if something happens to the particle (absorption, error, ...)
	 *
	 * @param tmax Max. absolute time at which integration will be stopped
	 * @param conf Option map containing particle specific options from particle.in
	 */
	void Integrate(double tmax, std::map<std::string, std::string> &conf);

protected:
	TGeometry *geom; ///< TGeometry structure passed by "Integrate"
	TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
	TFieldManager *field; ///< TFieldManager structure passed by "Integrate"
private:
	std::map<solid, bool> currentsolids; ///< solids in which particle is currently inside

	/**
	 * Return first non-ignored solid in TParticle::currentsolids list
	 */
	const solid& GetCurrentsolid() const;


	/**
	 * Equations of motion dy/dx = f(x,y).
	 *
	 * Equations of motion (fully relativistic).
	 * Including gravitation, Lorentz-force and magnetic interaction with magnetic moment.
	 *
	 * @param y	State vector (position, velocity, proper time, and polarization)
	 * @param dydx Returns derivatives of y with respect to x
	 * @param x Time
	 */
	void derivs(const state_type &y, state_type &dydx, const value_type x) const;


	/**
	 * Equations of motion dy/dx = f(x,y).
	 *
	 * Equations of motion (fully relativistic).
	 * Including gravitation, Lorentz-force and magnetic interaction with magnetic moment.
	 *
	 * @param y	State vector (position, velocity, proper time, and polarization)
	 * @param dydx Returns derivatives of y with respect to x
	 * @param x Time
	 * @param B Magnetic field
	 * @param E Electric field
	 */
	void EquationOfMotion(const state_type &y, state_type &dydx, const value_type x, const double B[4][4], const double E[3]) const;


	/**
	 * Check if current collision is consistent with list of current solids
	 *
	 * @param hitsolid Solid that was hit.
	 * @param distnormal Distance to surface, can be positive (outgoing) or negative (incoming), depending on direction of particle velocity.
	 *
	 * @return Returns true, if an inconsistency was found.
	 */
	bool CheckHitError(const solid &hitsolid, double distnormal) const;


	/**
	 * Check, if particle hit a material boundary or was absorbed.
	 *
	 * Checks if a particle which flies from y1 to y2 in time x2-x1 hits a surface or is absorbed inside a material.
	 * If a surface is hit the routine splits the line segment y1->y2 on both sides of the collision point and calls itself recursively
	 * with the three new line segments as parameters. This is repeated until both split points are nearer than REFLECTION_TOLERANCE
	 * to the collision point (much like an bisection algorithm). For each line segment "OnStep" is called to check for scattering/absorption/etc.
	 * For each short segment crossing a collision point "OnHit" is called to check for reflection/refraction/etc.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment
	 * @param y2 End point of line segment
	 * @param stepper Trajectory integrator, used to calculate intermediate state vectors
	 * @param hitlog Should hits be logged to file?
	 * @param iteration Iteration counter (incremented by recursive calls to avoid infinite loop)
	 * @return Returns true if particle was reflected/absorbed
	 */
	bool CheckHit(const value_type x1, const state_type y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper, const bool hitlog, const int iteration = 1);


	/**
	 * Simulate spin precession
	 *
	 * Integrates general BMT equation over one time step.
	 * If the conditions given by times and Bmax are not fulfilled, the spin vector will simply be rotated along the magnetic field, keeping the spin projection onto the magnetic field constant.
	 *
	 * @param spin Spin vector, returns new spin vector after step
	 * @param stepper Trajectory integrator containing last step
	 * @param x2 Time at end of step [s]
	 * @param y2 Particle state vector at end of step (position, velocity, proper time, and polarisation)
	 * @param times Absolute time intervals in between spin integration should be carried out [s]
	 * @param interpolatefields If this is set to true, the magnetic and electric fields will be interpolated between the trajectory-step points. This will speed up spin tracking in high, static fields, but might break spin tracking in small, quickly varying fields (e.g. spin-flip pulses)
	 * @param Bmax Spin integration will only be carried out, if magnetic field is below this value [T]
	 * @param flipspin If set to true, polarisation in y2 will be randomly set when magnetic field rises above Bmax, weighted by spin projection onto the magnetic field
	 * @param spinloginterval Min. distance [s] between spin-trajectory prints
	 * @param nextspinlog Time at which the next spin-trajectory point should be written to file
	 *
	 * @return Return probability of spin flip
	 */
	double IntegrateSpin(state_type &spin, const dense_stepper_type &stepper, const double x2, state_type &y2, const std::vector<double> &times,
						const bool interpolatefields, const double Bmax, const bool flipspin, const double spinloginterval, double &nextspinlog);

private:
	/**
	 * Calculate spin precession axis.
	 *
	 * Includes relativistic distortion of magnetic and electric fields and Thomas precession.
	 *
	 * @param t Time
	 * @param stepper Trajectory integrator used to calculate position and velocity at time t
	 * @param Omega Returns precession axis as 3-vector in lab frame
	 */
	void SpinPrecessionAxis(const double t, const dense_stepper_type &stepper, double &Omegax, double &Omegay, double &Omegaz) const;

	/**
	 * Calculate spin precession axis.
	 *
	 * Includes relativistic distortion of magnetic and electric fields and Thomas precession.
	 *
	 * @param t Time
	 * @param stepper Trajectory integrator used to calculate position and velocity at time t
	 * @param B Magnetic field
	 * @param E Electric field
	 * @param dydt Right-hand side of particle's equation of motion, velocity and acceleration requried to calculate vxE effect and Thomas precession
	 * @param Omega Returns precession axis as 3-vector in lab frame
	 */
	void SpinPrecessionAxis(const double t, const double B[4][4], const double E[3], const state_type &dydt, double &Omegax, double &Omegay, double &Omegaz) const;

	/**
	 * Equations of motion of spin vector.
	 *
	 * Calculates spin-precession axis either directly or from pre-calculated splines
	 *
	 * @param y Current spin vector
	 * @param dydx Calculated time derivative of spin vector
	 * @param x Current time
	 * @param stepper Trajectory integrator used to calculate spin-precession axis
	 * @param omega Vector of three splines used to interpolate spin-precession axis (can be empty)
	 */
	void SpinDerivs(const state_type &y, state_type &dydx, const value_type x,
			const dense_stepper_type &stepper, const std::vector<alglib::spline1dinterpolant> &omega_int) const;

protected:
	/**
	 * This virtual method is executed, when a particle crosses a material boundary.
	 *
	 * Can be used to reflect/refract/etc. at material boundaries.
	 * Has to be implemented separately for each derived particle.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	virtual void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
						const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed) = 0;


	/**
	 * This virtual method is executed on each integration step.
	 *
	 * Can be used to scatter/absorb/etc. in materials.
	 * Has to be implemented separately for each derived particle.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param stepper Trajectory integrator, can be used to calculate intermediate state vectors
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle path was changed
	 */
	virtual bool OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper, const solid &currentsolid) = 0;


	/**
	 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
	 */
	virtual void Decay() = 0;


	/**
	 * Return stream for each log type.
	 *
	 * Purely virtual function, has to be implemented for each particle type.
	 *
	 * @param str Enum specifying which stream should be returned
	 *
	 * @return Returns stream corresponding to str
	 */
	virtual std::ofstream& GetLogStream(const LogStream str) const = 0;

	
	/**
	 * Print start and current values to the endLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param spin Current spin vector
	 * @param sld Solid in which the particle is currently.
	 * @param logType Select either endlog or snapshotlog
	 */
	virtual void Print(const value_type x, const state_type &y, const state_type &spin, const solid &sld, const LogStream logType) const;


	/**
	 * Print current track point  to the trackLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param spin Spin vector
	 * @param sld Solid in which the particle is currently.
	 */
	virtual void PrintTrack(const value_type x, const state_type &y, const state_type &spin, const solid &sld) const;


	/**
	 * Print material boundary hits to the hitLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Time of material hit
	 * @param y1 State vector before material hit
	 * @param y2 State vector after material hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Material which is left at this boundary
	 * @param entering Material which is entered at this boundary
	 */
	virtual void PrintHit(const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering) const;


	/**
	 * Write spin state to the spinLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x time
	 * @param spin Spin vector
	 * @param stepper Trajectory integrator used to calculate spin-precession axis at time t
	 */
	virtual void PrintSpin(const value_type x, const state_type &spin, const dense_stepper_type &stepper) const;


	/**
	 * Calculate kinetic energy.
	 *
	 * Kinetic energy is calculated by series expansion of rel. gamma factor for small velocities, else it is calculated exactly by (gamma-1)mc^2
	 *
	 * @param v Velocity vector [m/s]
	 *
	 * @return Kinetic energy [eV]
	 */
	double Ekin(const value_type v[3]) const;


	/**
	 * Calculate potential energy of particle
	 *
	 * @param t Time
	 * @param y Coordinate vector
	 * @param field Pointer to TFieldManager structure for electric and magnetic potential
	 * @param sld Solid in which the particle is currently.
	 *
	 * @return Returns potential energy [eV]
	 */
	virtual double Epot(const value_type t, const state_type &y, TFieldManager *field, const solid &sld) const;

};



#endif // PARTICLE_H_
