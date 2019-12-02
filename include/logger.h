#ifndef PENTRACK_LOGGER_H
#define PENTRACK_LOGGER_H

#include <memory>
#include <map>
#include <fstream>

#include "particle.h"
#include "geometry.h"
#include "fields.h"
#include "stepper.h"

#ifdef USEROOT
#include "TFile.h"
#include "TNtupleD.h"
#endif

class TLogger {
protected:
    TConfig config;
    virtual void Log(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) = 0;
public:
    virtual ~TLogger(){ };
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
    void Print(const std::unique_ptr<TParticle>& p, const value_type &x, const TStep &stepper, const state_type &spin,
            const TGeometry &geom, const TFieldManager &field, const std::string suffix = "end");

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
    void PrintSnapshot(const std::unique_ptr<TParticle>& p, const TStep& stepper, const state_type &spin, const TGeometry &geom, const TFieldManager &field);


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
    void PrintTrack(const std::unique_ptr<TParticle>& p, const TStep &stepper,
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
    void PrintHit(const std::unique_ptr<TParticle>& p, const TStep &stepper, const double *normal, const solid &leaving, const solid &entering);


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
    void PrintSpin(const std::unique_ptr<TParticle>& p, const dense_stepper_type& spinstepper,
                   const TStep &trajectory_stepper, const TFieldManager &field);

    virtual void Close(){ };
};

class TTextLogger: public TLogger {
private:
    std::map<std::string, std::ofstream> logstreams;
    void Log(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) final;
public:
    TTextLogger(TConfig& aconfig){ config = aconfig; };
    ~TTextLogger() final { for (auto &s: logstreams){ s.second.close(); } };
};

#ifdef USEROOT
class TROOTLogger: public TLogger {
private:
    TFile* ROOTfile;
    void Log(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) final;
public:
    TROOTLogger(TConfig &aconfig);
    ~TROOTLogger() final;
};
#endif

std::unique_ptr<TLogger> CreateLogger(TConfig& config);


#endif //PENTRACK_LOGGER_H
