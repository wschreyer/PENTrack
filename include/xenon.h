/**
 * \file
 * Xenon class definition. Xenon-129 is used as a comagnetometer used in the TRIUMF EDM experiment. 
 */

#ifndef XENON_H_
#define XENON_H_

#include "particle.h"

extern const char* NAME_XENON; ///< name of TXenon class

/**
 * Xenon particle class.
 *
 * Simulates a xenon atom including gravitation and Lorentz-force
 */
struct TXenon:TParticle{
public:
	/**
	 * Create xenon-129 atom. 
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
	 * @param polarisation polarisation (+/-1)
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TXenon (const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
			const TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);

protected:
	static std::ofstream endout; ///< endlog file stream
	static std::ofstream snapshotout; ///< snapshot file stream
	static std::ofstream trackout; ///< tracklog file stream
	static std::ofstream hitout; ///< hitlog file stream
	static std::ofstream spinout; ///< spinlog file stream
	
	/**
	 * This method is executed, when a particle encounters a material boundary.
	 * The atom is reflected specularly or diffusely depending on the material type. 
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
			const solid &leaving, const solid &entering, const TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;

	
	/**
	 * This method is executed on each step.
	 *
	 * Do nothing. 
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
			const solid &currentsolid, const TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Xenon decay (not used)
	 *
	 * For parameter doc see TParticle::Decay
	 */
	void Decay(const TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;


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

#endif // XENON_H_
