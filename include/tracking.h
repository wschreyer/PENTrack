#ifndef PENTRACK_TRACKING_H
#define PENTRACK_TRACKING_H

#include <string>
#include <memory>

#include "mc.h"
#include "geometry.h"
#include "fields.h"
#include "particle.h"
#include "logger.h"
#include "stepper.h"

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
    void DoStep(const std::unique_ptr<TParticle>& p, TStep &stepper, const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field);

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
    void DoHit(const std::unique_ptr<TParticle>& p, TStep &stepper, TMCGenerator &mc, const TGeometry &geom);

    /**
     * Return first non-ignored solid in TParticle::currentsolids list
     */
    const solid& GetCurrentsolid() const;
};


#endif //PENTRACK_TRACKING_H
