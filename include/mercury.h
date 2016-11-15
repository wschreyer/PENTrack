/**
 * \file
 * Mercury class definition. Mercury-199 is used as a comagnetometer in the EDM experiment at TRIUMF. 
 */

#ifndef MERCURY_H_
#define MERCURY_H_

#include "particle.h"

extern const char* NAME_MERCURY; ///< name of TMERCURY class

/**
 * Mercury-199 particle class.
 *
 * Simulates a mercury atom including gravitation and Lorentz-force
 */
struct TMercury: TParticle{
public:
	/**
	 * Create mercury-199 atom.
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
	TMercury(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation,
			TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

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
	bool OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
			const solid &leaving, const solid &entering, bool &traversed, stopID &ID, std::vector<TParticle*> &secondaries) const;

	/**
	 * Reflect mercury-199 from surface.
	 *
	 * Reflects or scatters the neutron according to specular or Lambert.
	 *
	 * For parameter doc see TMercury::OnHit
	 */
	bool Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
				const double normal[3], const solid &leaving, const solid &entering, bool &traversed) const;
	
	/**
	 * This method is executed on each step.
	 *
	 * Do nothing. 
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	bool OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const dense_stepper_type &stepper, const solid &currentsolid, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Mercury decay (not used)
	 *
	 * For parameter doc see TParticle::Decay
	 */
	void Decay(std::vector<TParticle*> &secondaries) const;


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

#endif // MERCURY_H_
