#ifndef PENTRACK_LOGGER_H
#define PENTRACK_LOGGER_H

#include <memory>
#include <map>
#include <fstream>

#include "particle.h"
#include "geometry.h"
#include "fields.h"

class TLogger {
private:
    std::map<std::string, std::ofstream> logstreams;
	TConfig config;
	std::map<std::string, std::istringstream> snapshots;
	std::map<std::string, double> nextsnapshot;
public:
	TLogger(const TConfig& aconfig): config(aconfig){ };

	/**
	 * Print start and current values to the endLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param spin Current spin vector
	 * @param geom Geometry of the simulation
	 * @param field TFieldManager containing all electromagnetic fields
	 * @param logType Select either endlog or snapshotlog
	 */
	void Print(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y, const state_type &spin, const TGeometry &geom, const TFieldManager &field, const std::string& suffix = "end");


	/**
	 * Print start and current values to the endLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param spin Current spin vector
	 * @param geom Geometry of the simulation
	 * @param field TFieldManager containing all electromagnetic fields
	 * @param logType Select either endlog or snapshotlog
	 */
	void PrintSnapshot(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x2, const state_type &y2,
		const state_type &spin, const dense_stepper_type& stepper, const TGeometry &geom, const TFieldManager &field);


	/**
	 * Print current track point  to the trackLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param spin Spin vector
	 * @param sld Solid in which the particle is currently.
	 * @param field TFieldManager containing all electromagnetic fields
	 */
	void PrintTrack(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x, const state_type& y,
		const state_type &spin, const solid &sld, const TFieldManager &field);


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
	void PrintHit(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering);


	/**
	 * Write spin state to the spinLog returned by GetLogStream.
	 *
	 * This is a simple prototype that can be overridden by derived particle classes.
	 *
	 * @param x time
	 * @param spin Spin vector
	 * @param stepper Trajectory integrator used to calculate spin-precession axis at time t
	 * @param field TFieldManager containing all electromagnetic fields
	 */
	// virtual void PrintSpin(const value_type x, const state_type &spin, const dense_stepper_type &stepper, const TFieldManager &field) const;
	void PrintSpin(const std::unique_ptr<TParticle>& p, const value_type x, const dense_stepper_type& spinstepper,
		const dense_stepper_type &trajectory_stepper, const TFieldManager &field);
};

extern std::unique_ptr<TLogger> logger;

void CreateLogger(TConfig& config);


#endif //PENTRACK_LOGGER_H
