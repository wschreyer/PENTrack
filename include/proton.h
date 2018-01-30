/**
 * \file
 * Proton class definition.
 */

#ifndef PROTON_H_
#define PROTON_H_

extern const char* NAME_PROTON; ///< name of TProton class

#include "globals.h"
#include "particle.h"

/**
 * Proton particle class.
 *
 * Simulates a proton including gravitation and Lorentz-force
 */
struct TProton: TParticle{
public:
	/**
	 * Create proton
	 *
	 * Wraps basic TParticle constructor
	 *
	 * @param number Particle number
	 * @param t Starting time
	 * @param x Initial x coordinate
	 * @param y Initial y coordinate
	 * @param z Initial z coordinate
	 * @param E Initial kinetic energy
	 * @param phi Initial azimuth of velocity vector
	 * @param theta Initial polar angle of velocity vector
	 * @param polarisation polarisation of particle (+/-1)
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TProton(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
			TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);

protected:
	static std::ofstream endout; ///< endlog file stream
	static std::ofstream snapshotout; ///< snapshot file stream
	static std::ofstream trackout; ///< tracklog file stream
	static std::ofstream hitout; ///< hitlog file stream
	static std::ofstream spinout; ///< spinlog file stream
	
	/**
	 * This method is executed, when a particle crosses a material boundary.
	 *
	 * Nothing happens to protons.
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
			const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * This method is executed on each step.
	 *
	 * Protons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
			const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Proton decay (not used)
	 *
	 * For parameter doc see TParticle::Decay
	 */
	void Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;


	/**
	 * Return this particle type's log files
	 *
	 * @param str Enum specifying, which log file should be returned.
	 *
	 * @return log file
	 */
	std::ofstream& GetLogStream(const LogStream str) const{
		switch (str){
			case endLog: return endout;
			case snapshotLog: return snapshotout;
			case hitLog: return hitout;
			case trackLog: return trackout;
			case spinLog: return spinout;
		}
		throw std::out_of_range("Unknown LogStream requested!\n");
	};
};

#endif // PROTON_H_
