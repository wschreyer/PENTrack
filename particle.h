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

using namespace std;

static const double MAX_SAMPLE_DIST = 0.01; ///< max spatial distance of reflection checks, spin flip calculation, etc; longer integration steps will be interpolated


/**
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor and methods TParticle::OnHit, TParticle::OnStep and TParticle::Decay.
 * Optionally, derived particles can also re-implement TParticle::Epot, TParticle::PrintStartEnd, TParticle::PrintTrack, TParticle::PrintSnapshots, TParticle::PrintHits and define its own constructors.
 */
struct TParticle{
protected:
	typedef double value_type;
	typedef std::vector<value_type> state_type;
	typedef boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper_type;
	typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type;
	typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type;
public:
	const char *name; ///< particle name (has to be initialized in all derived classes!)
	const long double q; ///< charge [C] (has to be initialized in all derived classes!)
	const long double m; ///< mass [eV/c^2] (has to be initialized in all derived classes!)
	const long double mu; ///< magnetic moment [J/T] (has to be initialized in all derived classes!)
	const long double gamma; ///< gyromagnetic ratio [rad/(s T)] (has to be initialized in all derived classes!)
	int particlenumber; ///< particle number
	int ID; ///< particle fate (defined in globals.h)
	double tau; ///< particle life time
	double maxtraj; ///< max. simulated trajectory length
	double inttime; ///< integration computing time
	double refltime; ///< reflection computing time

	/// start time
	value_type tstart;

	/// stop time
	value_type tend;

	/// state vector before integration
	state_type ystart;

	/// state vector after integration)
	state_type yend;

	/// initial polarisation of particle (-1,0,1)
	int polstart;

	/// final polarisation of particle (-1,0,1)
	int polend;

	/// total start energy
	double Hstart();

	/// total end energy
	double Hend();

	/// max total energy
	double Hmax;

	/// kinetic start energy
	double Estart();

	/// kinetic end energy
	double Eend();

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
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, int number,
			double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

	/**
	 * Destructor, deletes secondaries
	 */
	virtual ~TParticle();


	/**
	 * Returns equations of motion.
	 *
	 * Class TParticle is given to integrator, which calls TParticle(x,y,dydx)
	 *
	 * @param x Time
	 * @param y State vector (position + velocity)
	 * @param dydx Returns derivatives of y with respect to t
	 */
	void operator()(state_type y, state_type &dydx, value_type x);


	/**
	 * Integrate particle trajectory.
	 *
	 * Takes inital state vector ystart and integrates the trajectory step by step.
	 * If a step is longer than MAX_SAMPLE_DIST, the step is split by interpolating intermediate points.
	 * On each step it checks for interaction with solids, prints snapshots and track into files and calls "OnEachStep".
	 * The Integration stops if TParticle::tau or tmax are reached; or if something happens to the particle (absorption, numerical error)
	 *
	 * @param tmax Max. absolute time at which integration will be stopped
	 * @param geometry Integrator checks for collisions with walls defined in this TGeometry structure
	 * @param mcgen Random number generator (e.g. used for reflection probability dicing)
	 * @param afield Pointer to field which acts on particle (can be NULL)
	 * @param conf Option map containing particle specific options from particle.in
	 * @param output TOutput struct containing output options and streams
	 */
	void Integrate(double tmax, TGeometry &geometry, TMCGenerator &mcgen, TFieldManager *afield, std::map<std::string, std::string> &conf);

protected:
	std::map<solid, bool> currentsolids; ///< solids in which particle is currently inside
	TGeometry *geom; ///< TGeometry structure passed by "Integrate"
	TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
	TFieldManager *field; ///< TFieldManager structure passed by "Integrate"
	dense_stepper_type stepper; ///< ODE integrator

	solid GetCurrentsolid();


	/**
	 * Equations of motion dy/dx = f(x,y).
	 *
	 * Equations of motion (fully relativistic), called by operator().
	 * Including gravitation, Lorentz-force and magnetic interaction with magnetic moment.
	 *
	 * @param x Time
	 * @param y	State vector (position and velocity)
	 * @param dydx Returns derivatives of y with respect to x
	 */
	void derivs(value_type x, state_type y, state_type &dydx);


	bool CheckHitError(solid *hitsolid, double distnormal);


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
	 * @param pol Particle polarisation
	 * @param hitlog Should hits be logged to file?
	 * @param hitout File stream to which hits are logged.
	 * @param iteration Iteration counter (incremented by recursive calls to avoid infinite loop)
	 * @return Returns true if particle was reflected/absorbed
	 */
	bool CheckHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &pol, bool hitlog, int iteration = 1);


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
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param travered Returns true if the material boundary was traversed by the particle
	 */
	virtual void OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
						const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed) = 0;


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
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle path was changed
	 */
	virtual bool OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid) = 0;


	/**
	 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
	 */
	virtual void Decay() = 0;


	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * This is a purely virtual function that has to be implemented by all derived classes.
	 * It can e.g. simply call the simple prototype TParticle::Print below.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	virtual void Print(value_type x, state_type y, int polarisation) = 0;


	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * This is a purely virtual function that has to be implemented by all derived classes.
	 * It can e.g. simply call the simple prototype TParticle::Print below.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	virtual void PrintSnapshot(value_type x, state_type y, int polarisation) = 0;


	/**
	 * Write the particle's trajectory into a file.
	 *
	 * This is a purely virtual function that has to be implemented by all derived classes.
	 * It can e.g. simply call the simple prototype TParticle::PrintTrack below.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	virtual void PrintTrack(value_type x, state_type y, int polarisation) = 0;


	/**
	 * Write the particle properties into a file, before and after it hit a material boundary.
	 *
	 * This is a purely virtual function that has to be implemented by all derived classes.
	 * It can e.g. simply call the simple prototype TParticle::PrintHit below.
	 *
	 * @param x Time of material hit
	 * @param y1 State vector before material hit
	 * @param y2 State vector after material hit
	 * @param pol1 Polarisation before material hit
	 * @param pol2 Polarisation after material hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Material which is left at this boundary
	 * @param entering Material which is entered at this boundary
	 */
	virtual void PrintHit(value_type x, state_type y1, state_type y2, int pol1, int pol2, const double *normal, solid *leaving, solid *entering) = 0;


	/**
	 * Get spin log stream.
	 *
	 * Has to be derived by all derived classes.
	 *
	 * @return Reference to spin log stream
	 */
	virtual ofstream& GetSpinOut() = 0;


	/**
	 * Print start and current values to a stream.
	 *
	 * This is a simple prototype that can be called by derived particle classes.
	 *
	 * @param file stream to print into
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 * @param filesuffix Optional suffix added to the file name (default: "end.out")
	 */
	void Print(std::ofstream &file, value_type x, state_type y, int polarisation, std::string filesuffix = "end.out");


	/**
	 * Print current track point into stream to allow visualization of the particle's trajectory.
	 *
	 * This is a simple prototype that can be called by derived particle classes.
	 *
	 * @param trackfile Stream to print into
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	void PrintTrack(std::ofstream &trackfile, value_type x, state_type y, int polarisation);


	/**
	 * Print material boundary hits into stream.
	 *
	 * This is a simple prototype that can be called by derived particle classes.
	 *
	 * @param hitfile stream to print into
	 * @param x Time of material hit
	 * @param y1 State vector before material hit
	 * @param y2 State vector after material hit
	 * @param pol1 Polarisation before material hit
	 * @param pol2 Polarisation after material hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Material which is left at this boundary
	 * @param entering Material which is entered at this boundary
	 */
	void PrintHit(std::ofstream &hitfile, value_type x, state_type y1, state_type y2, int pol1, int pol2, const double *normal, solid *leaving, solid *entering);


	/**
	 * Calculate kinetic energy.
	 *
	 * Kinetic energy is calculated by series expansion of rel. gamma factor, if v/c < ::RELATIVSTIC_THRESHOLD, else it is calculated exactly by (gamma-1)mc^2
	 *
	 * @return Kinetic energy [eV]
	 */
	double Ekin(value_type v[3]);


	/**
	 * Calculate potential energy of particle
	 *
	 * @param t Time
	 * @param y Coordinate vector
	 * @param polarisation Particle polarisation
	 * @param field Pointer to TFieldManager structure for electric and magnetic potential
	 * @param solids List of solids for optional material potential
	 *
	 * @return Returns potential energy [eV]
	 */
	virtual double Epot(value_type t, state_type y, int polarisation, TFieldManager *field, std::map<solid, bool> &solids);

};



#endif // PARTICLE_H_
