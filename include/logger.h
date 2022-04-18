#ifndef PENTRACK_LOGGER_H
#define PENTRACK_LOGGER_H

#include <memory>
#include <map>
#include <fstream>

#include "particle.h"
#include "geometry.h"
#include "fields.h"

#ifdef USEROOT
#include "TFile.h"
#include "TNtupleD.h"
#endif

#ifdef USEHDF5
#include "hdf5.h"
#include "hdf5_hl.h"
#endif


/**
 * Virtual base class printing particle states, track, spin
 */
class TLogger {
protected:
    TConfig config; ///< configuration parameters read from config files

    /**
     * Evaluates the logvars and corresponding filters and formulas set in the config file and calls DoLog
     * 
     * @param particlename Name of particle being logged
     * @param suffix Indicates logging type (e.g. "end", "snapshot", "track", "spin")
     * @param variables Maps of variable names and their values used to evaluate logvars
     * @param default_titles Optional parameter containing default variables to be logged in case none are given in the config
     */
    void Log(const std::string &particlename, const std::string &suffix, const std::map<std::string, double> &variables, const std::vector<std::string> &default_titles = {});

    /**
     * Virtual function actually doing the logging. Must be implemented in all derived classes
     * 
     * @param particlename Name of particle being logged
     * @param suffix Indicates logging type (e.g. "end", "snapshot", "track", "spin"), can be used for e.g. filename
     * @param titles List of names for each variable
     * @param vars List of variables to be logged
     */
    virtual void DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) = 0;
public:
    virtual ~TLogger(){ }; ///< Virtual desctructor (empty)
    /**
     * Print start and current states of a particle
     *
     * Collects variables and passes them to the virtual Log function
     *
     * @param p Particle to be printed
     * @param x Current time
     * @param y Current state vector
     * @param spin Current spin vector
     * @param geom Geometry of the simulation
     * @param field TFieldManager containing all electromagnetic fields
     * @param suffix Indicates logging type (e.g. "end", "snapshot"), passed to the virtual Log function
     */
    void Print(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y, const state_type &spin,
            const TGeometry &geom, const TFieldManager &field, const std::string suffix = "end");

    /**
     * Print start and current values at a specific time within an integration step
     *
     * Collects variables and passes them to the virtual Log function
     *
     * @param p Particle to be printed
     * @param x1 Start time of integration step
     * @param y1 Start state of integration step
     * @param x2 End time of integration step
     * @param y2 End state of integration step
     * @param spin Current spin vector
     * @param stepper Stepper can be used to interpolate state between x1 and x2
     * @param geom Geometry of the simulation
     * @param field TFieldManager containing all electromagnetic fields
     */
    void PrintSnapshot(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x2, const state_type &y2,
                       const state_type &spin, const dense_stepper_type& stepper, const TGeometry &geom, const TFieldManager &field);


    /**
     * Print point on particle trajectory
     *
     * Collects variables and passes them to the virtual Log function
     *
     * @param p Particle to be printed
     * @param x1 Start time of integration step
     * @param y1 Start state of integration step
     * @param x Current time
     * @param y Current state
     * @param spin Spin vector
     * @param sld Solid in which the particle is currently.
     * @param field TFieldManager containing all electromagnetic fields
     */
    void PrintTrack(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x, const state_type& y,
                    const state_type &spin, const solid &sld, const TFieldManager &field);


    /**
     * Print material boundary hit
     *
     * Collects variables and passes them to the virtual Log function
     *
     * @param p Particle to be printed
     * @param x Time of material hit
     * @param y1 State vector before material hit
     * @param y2 State vector after material hit
     * @param normal Normal vector of hit surface
     * @param leaving Material which is left at this boundary
     * @param entering Material which is entered at this boundary
     */
    void PrintHit(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering);


    /**
     * Write spin state of particle
     *
     * Collects variables and passes them to the virtual Log function
     *
     * @param p Particle to be printed
     * @param x Time to print the spin state at
     * @param spinstepper Integration stepper used to interpolate spin trajectory
     * @param trajectory_stepper Trajectory integrator used to calculate spin-precession axis at time t
     * @param field TFieldManager containing all electromagnetic fields
     */
    void PrintSpin(const std::unique_ptr<TParticle>& p, const value_type x, const dense_stepper_type& spinstepper,
                   const dense_stepper_type &trajectory_stepper, const TFieldManager &field);

};

/**
 * Class to print particles states to text files
 */
class TTextLogger: public TLogger {
private:
    std::map<std::string, std::ofstream> logstreams; ///< List of file streams used for logging

    /**
     * Logs given variables to selected text file
     * 
     * @param particlename Name of particle to be printed
     * @param suffix Select file to log to (e.g. "end", "snapshot", "track", "spin")
     * @param titles List of variable names
     * @param vars List of variables to be logged
     */
    void DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) override;
public:
    /**
     * Constructor, reads relevant configuration parameters
     * 
     * @param aconfig List of configuration parameters read from config file
     */
    TTextLogger(TConfig& aconfig){ config = aconfig; };

    /**
     * Destructor, closes all opened file streams
     */
    ~TTextLogger() final { for (auto &s: logstreams){ s.second.close(); } };
};

#ifdef USEROOT
/**
 * Class to print particle states to ROOT trees. Only available if cmake found ROOT libraries.
 */
class TROOTLogger: public TLogger {
private:
    TFile* ROOTfile; ///< ROOT file to print to

    /**
     * Logs given variables to selected ROOT tree
     * 
     * @param particlename Name of particle to be printed
     * @param suffix Select file to log to (e.g. "end", "snapshot", "track", "spin")
     * @param titles List of variable names
     * @param vars List of variables to be logged
     */
    void DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) override;
public:
    /**
     * Constructor, reads relevant configuration parameters from config and opens ROOT file
     * 
     * @param config List of configuration parameters read from config file
     */
    TROOTLogger(TConfig &aconfig);

    /**
     * Destructor, writes ROOT trees to file and closes it
     */
    ~TROOTLogger() final;
};
#endif

#ifdef USEHDF5
class THDF5Logger: public TLogger {
private:
    hid_t HDF5file;

    /**
     * Logs given variables to selected HDF5 table
     * 
     * @param particlename Name of particle to be printed
     * @param suffix Select file to log to (e.g. "end", "snapshot", "track", "spin")
     * @param titles List of variable names
     * @param vars List of variables to be logged
     */
    void DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) override;
public:
    /**
     * Constructor, reads relevant configuration parameters from config and opens HDF5 file
     * 
     * @param config List of configuration parameters read from config file
     */
    THDF5Logger(TConfig &aconfig);

    /**
     * Destructor, closes HDF5 file
     */
    ~THDF5Logger() final;
};
#endif

/**
 * Instantiates one of the classes derived from TLogger, depending on configuration variables
 * 
 * @param config List of configuration variables
 * 
 * @returns Class derived from TLogger
 */
std::unique_ptr<TLogger> CreateLogger(TConfig& config);


#endif //PENTRACK_LOGGER_H
