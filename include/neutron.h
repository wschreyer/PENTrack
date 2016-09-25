/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

#include "particle.h"
#include "microroughness.h"

extern const char* NAME_NEUTRON; ///< name of TNeutron class

/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Create neutron
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
	 * @param polarisation polarisation
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TNeutron(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);
	
protected:
	static std::ofstream endout; ///< endlog file stream
	static std::ofstream snapshotout; ///< snapshot file stream
	static std::ofstream trackout; ///< tracklog file stream
	static std::ofstream hitout; ///< hitlog file stream
	static std::ofstream spinout; ///< spinlog file stream
	
	/**
	 * Check for reflection/transmission/absorption on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Diffuse reflection can be done according to Lambert model or Micro Roughness model.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering (can be modified by method)
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);


	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts neutron velocity.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Transmit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects neutron specularly.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * Caclulates adiabacity condition (TNeutron::frac).
	 * Estimates spin flip probability according to Vladimirskii (TNeutron::vlad)
	 * and by doing a bruteforce integration of the Bloch equation describing spin precession.
	 * Additionally it can print out a spatial neutron distribution matrix.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param polarisation Polarisation of particle, may be altered
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper, const solid &currentsolid);


	/**
	 * Neutron decay.
	 *
	 * Create decay proton and electron and add them to secondaries list
	 */
	void Decay();


	/**
	 * Potential energy, other than electromagnetic or gravitational, which is added to total energy of neutron
	 *
	 * @param t Time
	 * @param y Position vector
	 * @param polarisation Particle polarisation
	 * @param field TFieldManager for electromagnetic potential
	 * @param sld Solid in which the particle is currently.
	 *
	 * @return Returns potential energy plus Fermi-Potential of solid
	 */
	value_type Epot(const value_type t, const state_type &y, TFieldManager *field, const solid &sld) const;


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


#endif // NEUTRON_H_
