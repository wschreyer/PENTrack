/**
 * \file
 * Proton class definition.
 */

#ifndef PROTON_H_
#define PROTON_H_

static const char* NAME_PROTON = "proton";

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
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 */
	TProton(int number, double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_PROTON, ele_e, m_p, 0, 0, number, t, x, y, z, E, phi, theta, amc, geometry, afield){

	};

protected:
	static ofstream endout; ///< endlog file stream
	static ofstream snapshotout; ///< snapshot file stream
	static ofstream trackout; ///< tracklog file stream
	static ofstream hitout; ///< hitlog file stream
	static ofstream spinout; ///< spinlog file stream

	/**
	 * This method is executed, when a particle crosses a material boundary.
	 *
	 * Nothing happens to protons.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the proton is leaving
	 * @param entering Solid that the proton is entering
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(long double x1, long double y1[6], long double &x2, long double y2[6], int &polarisation,
				const long double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
		traversed = true;
		trajectoryaltered = false;
	};


	/**
	 * This method is executed on each step.
	 *
	 * Protons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid in which the proton is at the moment
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(long double x1, long double y1[6], long double &x2, long double y2[6], solid currentsolid){
		if (currentsolid.ID != geom->defaultsolid.ID){
			x2 = x1;
			for (int i = 0; i < 6; i++)
				y2[i] = y1[i];
			ID = currentsolid.ID;
			printf("Absorption!\n");
			return true;
		}
		return false;
	};


	/**
	 * Proton decay (not used)
	 */
	void Decay(){

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
	void Print(long double x, long double y[6], int polarisation){
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
	virtual void PrintSnapshot(long double x, long double y[6], int polarisation){
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
	virtual void PrintTrack(long double x, long double y[6], int polarisation){
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
	virtual void PrintHit(long double x, long double *y1, long double *y2, int pol1, int pol2, const long double *normal, solid *leaving, solid *entering){
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

};

ofstream TProton::endout; ///< endlog file stream
ofstream TProton::snapshotout; ///< snapshot file stream
ofstream TProton::trackout; ///< tracklog file stream
ofstream TProton::hitout; ///< hitlog file stream
ofstream TProton::spinout; ///< spinlog file stream

#endif // PROTON_H_
