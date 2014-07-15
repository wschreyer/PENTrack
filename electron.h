/**
 * \file
 * Electron class definition.
 */

#ifndef ELECTRON_H_
#define ELECTRON_H_

static const char* NAME_ELECTRON = "electron";


#include "globals.h"
#include "particle.h"

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
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 */
	TElectron(int number, double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_ELECTRON, -ele_e, m_e, 0, 0, number, t, x, y, z, E, phi, theta, amc, geometry, afield){

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
	 * @param travered Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(long double x1, long double y1[6], long double &x2, long double y2[6], int &polarisation,
				const long double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
		traversed = true;
		trajectoryaltered = false;
	};


	/**
	 * This method is executed on each step.
	 *
	 * Electrons are immediately absorbed in solids other than TParticle::geom::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid in which the electron is at the moment
	 * @return Returns true if particle trajectory was altered
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
/*		else{
			long double v = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			long double s = sqrt(pow(y2[0]- y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));
			long double Ekin = m*c_0*c_0*(1/sqrt(1 - v*v/c_0/c_0) - 1);
			double Eloss, theta;
			int index;
			long double n = (1e-3*100)/(1.3806488e-23)/(77); // number density of rest gas (P [Pa])/(k [J/K])/(T [K])
			eH2(Ekin, s, n, &Eloss, &theta, &index); // did scattering happen?
			if (index != 0){
				printf("\nScattering %d, Eloss = %f\n",index,Eloss);
				double v_scat[3] = {y2[3], y2[4], y2[5]};
				eH2velocity(Eloss, theta, v_scat); // change velocity accordingly
				if (index == 3){ // on ionization create secondary electron
					TElectron *e = new TElectron(particlenumber, x2, y2, &y2[3], polarisation, *mc, *field);
					secondaries.push_back(e);
				}
			}
		}
*/		return false;
	};


	/**
	 * Electron decay (not used)
	 */
	void Decay(){

	}

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

ofstream TElectron::endout; ///< endlog file stream
ofstream TElectron::snapshotout; ///< snapshot file stream
ofstream TElectron::trackout; ///< tracklog file stream
ofstream TElectron::hitout; ///< hitlog file stream
ofstream TElectron::spinout; ///< spinlog file stream


#endif // ELECTRON_H_
