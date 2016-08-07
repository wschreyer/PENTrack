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
	TNeutron(int number, double t, double x, double y, double z, double E, double phi, double theta, int polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);
	
protected:
	static ofstream endout; ///< endlog file stream
	static ofstream snapshotout; ///< snapshot file stream
	static ofstream trackout; ///< tracklog file stream
	static ofstream hitout; ///< hitlog file stream
	static ofstream spinout; ///< spinlog file stream
	static ofstream spinout2; ///< spinlog file stream for doing simultaneous anti-parallel Efield spin integration
	static TMicroRoughness MR;
	
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
	void OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);


	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts neutron velocity.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Transmit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitMR(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitLambert(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects neutron specularly.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Reflect(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectMR(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectLambert(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);

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
	bool OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation, solid currentsolid);


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
	value_type Epot(value_type t, state_type y, int polarisation, TFieldManager *field, solid sld);

	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * Calls the simple prototype TParticle::Print.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 * @param sld Solid in which the particle is currently.
	 */
	void Print(value_type x, state_type y, int polarisation, solid sld){
		TParticle::Print(endout, x, y, polarisation, sld);
	};


	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * Calls the simple prototype TParticle::Print.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 * @param sld Solid in which the particle is currently.
	 */
	virtual void PrintSnapshot(value_type x, state_type y, int polarisation, solid sld){
		TParticle::Print(snapshotout, x, y, polarisation, sld, "snapshot.out");
	};


	/**
	 * Write the particle's trajectory into a file.
	 *
	 * Calls the simple prototype TParticle::PrintTrack.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 * @param sld Solid in which the particle is currently.
	 */
	virtual void PrintTrack(value_type x, state_type y, int polarisation, solid sld){
		TParticle::PrintTrack(trackout, x, y, polarisation, sld);
	};


	/**
	 * Write the particle properties into a file, before and after it hit a material boundary.
	 *
	 * Calls the simple prototype TParticle::PrintHit.
	 *
	 * @param x Time of material hit
	 * @param y1 State vector before material hit
	 * @param y2 State vector after material hit
	 * @param pol1 Polarisation before material hit
	 * @param pol2 Polarisation after material hit
	 * @param normal Normal vector of hit surface
	 * @param leaving Material which is left at this boundary
	 * @param entering Material which is entered at this boundary
	 */
	virtual void PrintHit(value_type x, state_type y1, state_type y2, int pol1, int pol2, const double *normal, solid *leaving, solid *entering){
		TParticle::PrintHit(hitout, x, y1, y2, pol1, pol2, normal, leaving, entering);
	};


	/**
	 * Get spin log stream.
	 *
	 * @return Returns static spinout stream to use same stream for all TNeutrons
	 */
	ofstream& GetSpinOut(){
		return spinout;
	};
	
	
	/**
	 * Get secondary spin log stream used for spin integration of anti-parallel E-field.
	 *
	 * @return Returns static spinout stream to use same stream for all TNeutrons
	 */
	ofstream& GetSpinOut2(){
		return spinout2;
	};
};


#endif // NEUTRON_H_
