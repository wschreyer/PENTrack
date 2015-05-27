/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

#include "particle.h"

#include "optimization.h"

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
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TNeutron(int number, double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

protected:
	static ofstream endout; ///< endlog file stream
	static ofstream snapshotout; ///< snapshot file stream
	static ofstream trackout; ///< tracklog file stream
	static ofstream hitout; ///< hitlog file stream
	static ofstream spinout; ///< spinlog file stream

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
	 * Refracts or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Transmit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);


	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Lambert or Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Reflect(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
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
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid);


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
	 * @param solids List of solids for optional material potential
	 *
	 * @return Returns potential energy plus Fermi-Potential of solid
	 */
	value_type Epot(value_type t, state_type y, int polarisation, TFieldManager *field, map<solid, bool> &solids);

	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * Calls the simple prototype TParticle::Print.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	void Print(value_type x, state_type y, int polarisation){
		TParticle::Print(endout, x, y, polarisation);
	};


	/**
	 * Write the particle's start properties and current values into a file.
	 *
	 * Calls the simple prototype TParticle::Print.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	virtual void PrintSnapshot(value_type x, state_type y, int polarisation){
		TParticle::Print(snapshotout, x, y, polarisation, "snapshot.out");
	};


	/**
	 * Write the particle's trajectory into a file.
	 *
	 * Calls the simple prototype TParticle::PrintTrack.
	 *
	 * @param x Current time
	 * @param y Current state vector
	 * @param polarisation Current polarisation
	 */
	virtual void PrintTrack(value_type x, state_type y, int polarisation){
		TParticle::PrintTrack(trackout, x, y, polarisation);
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

private:
	/**
	 * Check if the MicroRoughness is model is applicable to the current interaction
	 *
	 * @param y State vector right before hitting the material
	 * @param normal normal vector of the material boundary
	 * @param leaving Material that the particle is leaving through the boundary
	 * @param entering Material, that the particle entering through the boundary
	 *
	 * @return Returns true, if the MicroRoughness model can be used.
	 */
	bool MRValid(state_type y, const double normal[3], solid *leaving, solid *entering);

	/**
	 * Return MicroRoughness model distribution for scattering angles
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param integral Compute total probability of diffuse scattering according to MicroRoughness model
	 * @param y State vector of neutron right before hit surface
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid, which the neutron would leave if it were transmitted
	 * @param entering Solid, which the neutron would enter if it were transmitted
	 * @param theta_r Polar angle of scattered velocity vector (0 < theta_r < pi/2), ignored if integral == true
	 * @param phi_r Azimuthal angle of scattered velocity vector (0 < phi_r < 2*pi), ignored if integral == true
	 *
	 * @return Returns probability of reflection/transmission, in direction (theta_r, phi_r) or total if integral == true
	 */
	double MRDist(bool transmit, bool integral, state_type y, const double normal[3], solid *leaving, solid *entering, double theta, double phi);

	/**
	 * Struct containing parameters for TNeutron::MRDist and TNeutron::NegMRDist wrapper functions.
	 */
	struct TMRParams{
		TNeutron *n;
		bool transmit, integral;
		state_type y;
		const double *normal;
		solid *leaving, *entering;
	};

	/**
	 * Wrapper function to allow ALGLIB integration of MicroRoughness model distribution
	 *
	 * @param theta Polar angle of reflected/transmitted velocity
	 * @param xminusa Required by ALGLIB???
	 * @param bminusx Required by ALGLIB???
	 * @param params Contains TMRParams struct
	 *
	 * @return Returns value of TNeutron::MRDist at (theta, phi = 0) with given params
	 */
	static void MRDist(double theta, double xminusa, double bminusx, double &y, void *params);

	/**
	 * Wrapper function to allow ALGLIB maximization of MicroRoughness model distribution
	 *
	 * @param x Array of x1,x2,... values. Contains only single variable theta in this case.
	 * @param f NEGATIVE value of distribution at (theta = x[0], phi = 0)
	 * @param params Contains TMRParams struct
	 */
	static void NegMRDist(const alglib::real_1d_array &x, double &f, void *params);

	/**
	 * Calculate total diffuse scattering probability according to MicroRoughness model
	 *
	 * @param transmit True if the particle is transmitted through the surface, false if it is reflected
	 * @param y State vector of neutron right before hit surface
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid, which the neutron would leave if it were transmitted
	 * @param entering Solid, which the neutron would enter if it were transmitted
	 *
	 * @return Returns probability of reflection/transmission, in direction (theta_r, phi_r) or total if integral == true
	 */
	double MRProb(bool transmit, state_type y, const double normal[3], solid *leaving, solid *entering);

	/*
	 * Calculate maximum of MicroRoughness model distribution
	 *
	 * @param transmit True, if the particle is transmitted through the material boundary
	 * @param normal Normal vector of material boundary
	 * @param leaving Material, that the particle is leaving at this material boundary
	 * @param entering Material, that the particle is entering at this material boundary
	 *
	 * @return Returns maximal value of MicroRoughness model distribution in range (theta = 0..pi/2, phi = 0..2pi)
	 */
	double MRDistMax(bool transmit, state_type y, const double normal[3], solid *leaving, solid *entering);
};


#endif // NEUTRON_H_
