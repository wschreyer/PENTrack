/**
 * \file
 * Electron class definition.
 */

#ifndef ELECTRON_H_
#define ELECTRON_H_

#include "particle.h"

extern const char* NAME_ELECTRON; ///< name of TElectron class

/**
 * Electron particle class.
 *
 * Simulates an electron including gravitation and Lorentz-force
 */
struct TElectron: TParticle{
public:
	/**
	 * Create electron
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
	TElectron(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield);

protected:
	static ofstream endout; ///< endlog file stream
	static ofstream snapshotout; ///< snapshot file stream
	static ofstream trackout; ///< tracklog file stream
	static ofstream hitout; ///< hitlog file stream
	static ofstream spinout; ///< spinlog file stream
	static ofstream spinout2; ///< spinlog file stream for doing simultaneous anti-parallel Efield spin integration
	
	/**
	 * This method is executed, when a particle crosses a material boundary.
	 *
	 * Nothing happens to electrons.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the electron is leaving
	 * @param entering Solid that the electron is entering
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2,
				const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed);


	/**
	 * This method is executed on each step.
	 *
	 * Electrons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param polarisation Polarisation of particle, may be altered
	 * @param currentsolid Solid in which the electron is at the moment
	 * @return Returns true if particle trajectory was altered
	 */
	bool OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid);


	/**
	 * Electron decay (not used)
	 */
	void Decay();

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
	void Print(value_type x, state_type y, solid sld){
		TParticle::Print(endout, x, y, sld);
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
	virtual void PrintSnapshot(value_type x, state_type y, solid sld){
		TParticle::Print(snapshotout, x, y, sld, "snapshot.out");
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
	virtual void PrintTrack(value_type x, state_type y, solid sld){
		TParticle::PrintTrack(trackout, x, y, sld);
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
	virtual void PrintHit(value_type x, state_type y1, state_type y2, const double *normal, solid *leaving, solid *entering){
		TParticle::PrintHit(hitout, x, y1, y2, normal, leaving, entering);
	};


	/**
	 * Get spin log stream.
	 *
	 * @return Returns static spinout stream to use same stream for all TElectrons
	 */
	ofstream& GetSpinOut(){
		return spinout;
	};

	
	/**
	 * Get secondary spin log stream used for spin integration of anti-parallel E-field.
	 *
	 * @return Returns static spinout stream to use same stream for all TElectrons
	 */
	ofstream& GetSpinOut2(){
		return spinout2;
	};

};

#endif // ELECTRON_H_
