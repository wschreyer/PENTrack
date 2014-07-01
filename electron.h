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

};



#endif // ELECTRON_H_
