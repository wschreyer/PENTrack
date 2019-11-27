#ifndef PENTRACK_TRACKING_H
#define PENTRACK_TRACKING_H

#include <string>
#include <memory>

#include "mc.h"
#include "geometry.h"
#include "fields.h"
#include "particle.h"
#include "logger.h"


class TTracker {
private:
    std::vector<std::pair<solid, bool> > currentsolids; ///< solids in which particle is currently inside
    std::unique_ptr<TLogger> logger;
public:
    TTracker(TConfig& config);
    /**
     * Integrate particle trajectory.
     *
     * Takes state vector yend and integrates the trajectory step by step.
     * If a step is longer than MAX_SAMPLE_DIST, the step is split by interpolating intermediate points.
     * On each step it checks for interaction with solids, prints snapshots and track into files and calls TParticle::OnStep.
     * TParticle::StopIntegration is called if TParticle::tau or tmax are reached; or if something happens to the particle (absorption, error, ...)
     *
     * @param tmax Max. absolute time at which integration will be stopped
     * @param particleconf Option map containing particle specific options from particle.in
     * @param mc Random-number generator
     * @param geom Geometry of the simulation
     * @param field TFieldManager containing all electromagnetic fields
     */
    void IntegrateParticle(std::unique_ptr<TParticle>& p, const double tmax, std::map<std::string, std::string> &particleconf,
                           TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field);
private:
    /**
     * Check if particle hit a material boundary
     *
     * Check if particle hit a material boundary. Iterate exact collision point and call DoStep and DoHit appropriately.
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
    bool CheckHit(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
             const dense_stepper_type &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field);

    /**
     * Iterate collision point
     *
     * Split trajectory step right before and after collision point and call function recursively for each segment.
     *
     * @param x1 Start time of line segment
     * @param y1 Start point of line segment
     * @param x2 End time of line segment
     * @param y2 End point of line segment
     * @param coll Collision found in this segment
     * @param stepper Trajectory integrator, used to calculate intermediate state vectors
     * @param geom Geometry
     * @param interation Increase iteration count for each recursive call to limit number of iterations
     * @return Returns true if collision point was successfully iterated
     */
    bool iterate_collision(value_type &x1, state_type &y1, value_type &x2, state_type &y2,
                           const TCollision &coll, const dense_stepper_type &stepper, const TGeometry &geom,
                           const unsigned int iteration = 0);

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
     * @param geom Geometry
     * @return Returns true if trajectory was altered
     */
    bool DoStep(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
            const dense_stepper_type &stepper, const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field);

    /**
     * Call OnHit to check if particle should cross material boundary.
     *
     * Update list of solids the particle is in, check for geometry-tracking errors.
     *
     * @param x1 Start time of line segment
     * @param y1 Start point of line segment
     * @param x2 End time of line segment
     * @param y2 End point of line segment
     * @param stepper Trajectory integrator, used to calculate intermediate state vectors
     * @param mc Random-number generator
     * @param geom Geometry
     * @param hitlog Set true if hits should be logged to hitlog.out
     * @return Returns true if trajectory was altered
     */
    bool DoHit(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
            const dense_stepper_type &stepper, TMCGenerator &mc, const TGeometry &geom);

    /**
     * Return first non-ignored solid in TParticle::currentsolids list
     */
    const solid& GetCurrentsolid() const;

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
    void IntegrateSpin(const std::unique_ptr<TParticle>& p, state_type &spin, const dense_stepper_type &stepper,
            const double x2, state_type &y2, const std::vector<double> &times, const TFieldManager &field,
            const bool interpolatefields, const double Bmax, TMCGenerator &mc, const bool flipspin) const;



};


#endif //PENTRACK_TRACKING_H
