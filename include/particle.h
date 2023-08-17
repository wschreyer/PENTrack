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
#include "interpolation.h"

#include "geometry.h"
#include "mc.h"
#include "fields.h"

static const double MAX_TRACK_DEVIATION = 0.001; ///< max deviation of actual trajectory from straight line between start and end points of a step used for geometry-intersection test. If deviation is larger, the step will be split
static const int STATE_VARIABLES = 9; ///< number of variables in trajectory integration (position, velocity, proper time, polarization, path length)
static const int SPIN_STATE_VARIABLES = 5; ///< number of variables in spin integration (spin vector, time, total phase)

typedef double value_type; ///< data type used for trajectory integration
typedef std::vector<value_type> state_type; ///< type representing current particle state (position, velocity, proper time, and polarization) or spin state (x,y,z component)
typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, value_type> stepper_type; ///< basic integration stepper (5th-order Runge-Kutta)
typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type; ///< integration step length controller
typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type; ///< integration step interpolator
//	typedef boost::numeric::odeint::bulirsch_stoer_dense_out<state_type, value_type> dense_stepper_type2; ///< alternative stepper type (Bulirsch-Stoer)

/**
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor and methods TParticle::OnHit, TParticle::OnStep and TParticle::Decay.
 * Optionally, derived particles can also re-implement TParticle::Epot, TParticle::PrintStartEnd, TParticle::PrintTrack, TParticle::PrintSnapshots, TParticle::PrintHits and define its own constructors.
 */
struct TParticle{
private:
	std::string name; ///< particle name (has to be initialized in all derived classes!)
	const long double q; ///< charge [C] (has to be initialized in all derived classes!)
	const long double m; ///< mass [eV/c^2] (has to be initialized in all derived classes!)
	const long double mu; ///< magnetic moment [J/T] (has to be initialized in all derived classes!)
	const long double gamma; ///< gyromagnetic ratio [rad/(s T)] (has to be initialized in all derived classes!)
	int particlenumber; ///< particle number
	stopID ID; ///< particle fate (defined in globals.h)
	
	value_type tstart; ///< start time
	value_type tend; ///< stop time
	state_type ystart; ///< state vector before integration (position, velocity, proper time, polarization, and path length)
	state_type yend; ///< state vector after integration (position, velocity, proper time, polarization, and path length)
	state_type spinstart; ///< spin vector before integration
	state_type spinend; ///< spin vector after integration
	solid solidstart; ///< solid in which the particle started
	solid solidend; ///< solid in which particle stopped

	double Hmax; ///< max total energy
	int Nhit; ///< number of material boundary hits
	int Nspinflip; ///< number of spin flips
	long double noflipprob; ///< total probability of NO spinflip calculated by spin tracking
	int Nstep; ///< number of integration steps

	std::vector<std::unique_ptr<TParticle> > secondaries; ///< list of secondary particles
public:
	/**
	 * Return name of particle
	 *
	 * @return Name of particle
	 */
	std::string GetName() const { return name; };

	/**
	 * Return mass of particle
	 *
	 * @return Mass of particle
	 */
	double GetMass() const { return m; };

	/**
	 * Return charge of particle
	 *
	 * @return Charge of particle
	 */
	double GetCharge() const { return q; };

	/**
	 * Return magnetic moment of particle
	 *
	 * @return Magnetic moment of particle
	 */
	double GetMagneticMoment() const { return mu; };

	/**
	 * Return gyromagnetic ratio of particle
	 *
	 * @return Gyromagnetic ratio
	 */
	double GetGyromagneticRatio() const { return gamma; };

	/**
	 * Return number of particle
	 *
	 * @return Number of particle
	 */
	int GetParticleNumber() const { return particlenumber; };

	/**
	 * Return ID indicating why particle was stopped
	 *
	 * @return ID indicating why particle was stopped
	 */
	stopID GetStopID() const { return ID; };

	/**
	 * Return time when particle was created
	 *
	 * @return Time when particle was created
	 */
	value_type GetInitialTime() const { return tstart; };

	/**
	 * Return time when particle was stopped
	 *
	 * @return Time when particle was stopped
	 */
	value_type GetFinalTime() const { return tend; };

	/**
	 * Return initial state of particle (position, velocity, proper time, and polarization)
	 *
	 * @return Initial state of particle
	 */
	state_type GetInitialState() const { return ystart; };

	/**
	 * Return final state of particle (position, velocity, proper time, and polarization)
	 *
	 * @return Final state of particle
	 */
	state_type GetFinalState() const { return yend; };

	/**
	 * Return initial spin vector of particle
	 *
	 * @return Initial spin vector of particle
	 */
	state_type GetInitialSpin() const { return spinstart; };

	/**
	 * Return final spin vector of particle
	 *
	 * @return Final spin vector of particle
	 */
	state_type GetFinalSpin() const { return spinend; };

	/**
	 * Return solid in which particle was created
	 *
	 * @return Solid in which particle was created
	 */
	solid GetInitialSolid() const { return solidstart; };

	/**
	 * Return solid in which particle was stopped
	 *
	 * @return Solid in which particle stopped
	 */
	solid GetFinalSolid() const { return solidend; };

	/**
	 * Return maximal total energy on trajectory of particle
	 *
	 * @return Maximal total energy
	 */
	double GetMaxTotalEnergy() const { return Hmax; };

	/**
	 * Return total length of particle trajectory
	 *
	 * @return Length of trajectory
	 */
	double GetTrajectoryLength() const { return yend[8]; };

	/**
	 * Return number of times particle hit any part of geometry
	 *
	 * @return Number of hits
	 */
	int GetNumberOfHits() const { return Nhit; };

	/**
	 * Return number of times the particle's spin was flipped
	 * 
	 * @return Number of spin flips
	 */
	int GetNumberOfSpinflips() const { return Nspinflip; };

	/**
	 * Return probability that spin of particle was not flipped, determined with integration of BMT equation
	 *
	 * @return Probability of no spin flip
	 */
	double GetNoSpinFlipProbability() const { return noflipprob; };

	/**
	 * Return number of steps taken by integrator
	 *
	 * @return number of integration steps
	 */
	int GetNumberOfSteps() const { return Nstep; };

	/**
	 * Return initial total energy
	 *
	 * @return Initial total energy
	 */
	double GetInitialTotalEnergy(const TGeometry &geom, const TFieldManager &field) const;

	/**
	 * Return final total energy
	 *
	 * @return Final total energy
	 */
	double GetFinalTotalEnergy(const TGeometry &geom, const TFieldManager &field) const;

	/**
	 * Return initial kinetic energy
	 *
	 * @return Initial kinetic energy
	 */
	double GetInitialKineticEnergy() const;

	/**
	 * Return final kinetic energy
	 *
	 * @return Final kinetic energy
	 */
	double GetFinalKineticEnergy() const;

	/**
	 * Return list of secondary particles
	 * 
	 * @return List of particles
	 */
	std::vector<std::unique_ptr<TParticle> >& GetSecondaryParticles(){ return secondaries; }

	/**
	 * Set ID of particle
	 * 
	 * @param aID ID to set
	 */
	void SetStopID(const stopID aID){ ID = aID; }

	/**
	 * Set final state of particles
	 * 
	 * @param x Stop time
	 * @param y Stop state (position, velocity, proper time, polarization, trajectory length)
	 * @param spin 3-vector of particle spin
	 * @param sld Solid the particle was stopped in
	 */
	void SetFinalState(const value_type& x, const state_type& y, const state_type& spin, const solid& sld);

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
	 * @param polarisation polarisation of particle (-1 or 1) used for trajectory tracking in magnetic field
	 * @param spinprojection component of semi-classical (Bloch) spin vector parallel to magnetic field
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield TFieldManager containing all electromagnetic fields
	 */
	TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, const int number,
			const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation,
			const double spinprojection, TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);

	TParticle(const TParticle &p) = delete; ///< TParticle is not copyable
	TParticle& operator=(const TParticle &p) = delete; ///< TParticle is not copyable

	/**
	 * Destructor, deletes secondaries
	 */
	virtual ~TParticle(){ };


	/**
	 * Equations of motion dy/dx = f(x,y).
	 *
	 * Equations of motion (fully relativistic).
	 * Including gravitation, Lorentz-force and magnetic interaction with magnetic moment.
	 *
	 * @param y	State vector (position, velocity, proper time, and polarization)
	 * @param dydx Returns derivatives of y with respect to x
	 * @param x Time
	 * @param field TFieldManager used to calculate magnetic and electric fields
	 */
	void derivs(const state_type &y, state_type &dydx, const value_type x, const TFieldManager *field) const;


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
	 * @param dBidxj Spatial derivatives of magnetic field
	 * @param E Electric field
	 */
	void EquationOfMotion(const state_type &y, state_type &dydx, const value_type x, const double B[3], const double dBidxj[3][3], const double E[3]) const;


    /**
     * Call OnStep for particle-dependent physics processes on a step.
     *
     * Check if trajectory has been altered by physics processes, return true if it was.
     *
     * @param x1 Start time of line segment
     * @param y1 Start point of line segment
     * @param x2 End time of line segment
     * @param y2 End point of line segment
     * @param stepper Trajectory integrator, used to calculate intermediate state vectors
     * @param currentsolid Material the particle is in during this step
     * @param mc Random-number generator
     * @param field TFieldManager used to calculate electric and magnetic fields
	 * 
     * @return Returns true if trajectory was altered
     */
    void DoStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
                const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field);

    /**
     * Call OnHit to check if particle should cross material boundary.
     *
     * Update list of solids the particle is in, check for geometry-tracking errors.
     *
     * @param x1 Start time of line segment
     * @param y1 Start point of line segment
     * @param x2 End time of line segment
     * @param y2 End point of line segment
	 * @param normal Normal vector of the hit surface
	 * @param leaving Material that particle is leaving
	 * @param entering Material that particle is entering
     * @param mc Random-number generator
	 * 
     * @return Returns true if trajectory was altered
     */
    void DoHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
               const double normal[3], const solid &leaving, const solid &entering,
               TMCGenerator &mc);


    /**
     * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
     *
     * @param t Time when the particle decayed
     * @param y State (position, velocity, proper time, polarization, path length) of particle when it decayed
     * @param mc Random-number generator
     * @param geom Geometry of the simulation
     * @param field TFieldManager containing all electromagnetic fields
     */
    void DoDecay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field);
    
	/**
	 * Set polarization of particle
	 * 
	 * @param t Time
	 * @param y Returns state (position, velocity, proper time, polarization, path length) after polarization was set
	 * @param polarization Projection of spin onto polarization axis [-1..1]
	 * @param flipspin Indicates if polarization should actually be able to flip
	 * @param mc TMCGenerator random number generator
	 */
    void DoPolarization(const double t, state_type &y, const double polarization, const bool flipspin, TMCGenerator &mc);

	/**
	 * Randomly set polarization of particle
	 * 
	 * @param t Time
	 * @param y Returns state (position, velocity, proper time, polarization, path length) after polarization was set
	 * @param polarization Projection of spin onto polarization axis [-1..1]
	 * @param flipspin Indicates if polarization should actually be able to flip
	 * @param mc TMCGenerator random number generator
	 */
    void DoPolarize(const double t, state_type &y, const double polarization, const bool flipspin, TMCGenerator &mc);

    /**
	 * Calculate spin precession axis.
	 *
	 * Includes relativistic distortion of magnetic and electric fields and Thomas precession.
	 *
	 * @param t Time
	 * @param stepper Trajectory integrator used to calculate position and velocity at time t
	 * @param field TFieldManager used to calculate electric and magnetic fields
	 * @param Omegax Returns x component of precession axis in lab frame
	 * @param Omegay Returns y component of precession axis in lab frame
	 * @param Omegaz Returns z component of precession axis in lab frame
	 */
	void SpinPrecessionAxis(const double t, const dense_stepper_type &stepper, const TFieldManager &field, double &Omegax, double &Omegay, double &Omegaz) const;

	/**
	 * Calculate spin precession axis.
	 *
	 * Includes relativistic distortion of magnetic and electric fields and Thomas precession.
	 *
	 * @param t Time
	 * @param B Magnetic field
	 * @param E Electric field
	 * @param dydt Right-hand side of particle's equation of motion, velocity and acceleration required to calculate vxE effect and Thomas precession
	 * @param Omegax Returns x component of precession axis
	 * @param Omegay Returns y component of precession axis
	 * @param Omegaz Returns z component of precession axis
	 */
	void SpinPrecessionAxis(const double t, const double B[3], const double E[3], const state_type &dydt, double &Omegax, double &Omegay, double &Omegaz) const;

    /**
	 * Equations of motion of spin vector.
	 *
	 * Calculates spin-precession axis either directly or from pre-calculated splines
	 *
	 * @param y Current spin vector
	 * @param dydx Calculated time derivative of spin vector
	 * @param x Current time
	 * @param stepper Trajectory integrator used to calculate spin-precession axis
	 * @param field TFieldManager used to calculate magnetic and electric fields
	 * @param omega_int Vector of three splines used to interpolate spin-precession axis (can be empty)
	 */
	void SpinDerivs(const state_type &y, state_type &dydx, const value_type x,
			const dense_stepper_type &stepper, const TFieldManager *field, const std::vector<alglib::spline1dinterpolant> &omega_int) const;

	/**
	 * Calculate kinetic energy.
	 *
	 * Kinetic energy is calculated by series expansion of rel. gamma factor for small velocities, else it is calculated exactly by (gamma-1)mc^2
	 *
	 * @param v Velocity vector [m/s]
	 *
	 * @return Kinetic energy [eV]
	 */
	double GetKineticEnergy(const value_type v[3]) const;

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
	 * @param mc Random-number generator
	 * @param ID If particle is stopped, set this to the appropriate stopID
	 * @param secondaries Add any secondary particles produced in this interaction
	 */
	virtual void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
						const double normal[3], const solid &leaving, const solid &entering,
						TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const = 0;


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
	 * @param mc Random-number generator
	 * @param ID If particle is stopped, set this to the appropriate stopID
	 * @param secondaries Add any secondary particles produced in this interaction
	 */
	virtual void OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
			const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const = 0;


	/**
	 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
	 *
	 * @param t Time when the particle decayed
	 * @param y State (position, velocity, proper time, polarization, path length) of particle when it decayed
	 * @param mc Random-number generator
	 * @param geom Geometry of the simulation
	 * @param field TFieldManager containing all electromagnetic fields
	 * @param secondaries Add any secondary particles produced in this decay
	 */
	virtual void Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const = 0;


public:
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
	virtual double GetPotentialEnergy(const value_type t, const state_type &y, const TFieldManager &field, const solid &sld) const;
};

#endif // PARTICLE_H__
