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
	 * Electron constructor, set start values randomly according to *.in files.
	 *
	 * Sets start time TParticle::tstart according to [SOURCE] in geometry.in.
	 * Sets start energy according to TMCGenerator::Spectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 *
	 * @param number Particle number
	 * @param ageometry TGeometry in which the particle will be simulated
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TField used to calculate energies (can be NULL)
	 */
	TElectron(int number, TGeometry &ageometry, TSource &src,
				TMCGenerator &mcgen, TFieldManager *afield): TParticle(NAME_ELECTRON, -ele_e, m_e, 0){
		Init(number, ageometry, src, mcgen, afield);
	};

	/**
	 * Generic electron constructor
	 *
	 * @param number Particle number
	 * @param t Start time
	 * @param atau Life time
	 * @param x Start x coordinate
	 * @param y Start y coordinate
	 * @param z Start z coordinate
	 * @param vx Velocity x component
	 * @param vy Velocity y component
	 * @param vz Velocity z component
	 * @param pol Start polarisation
	 * @param trajl Max. simulated trajectory length
	 * @param ageometry TGeometry in which the particle will be simulated
	 * @param afield Electric and magnetic fields (needed to determine total energy)
	 */
	TElectron(int number, long double t, long double atau, long double x, long double y, long double z,
			long double vx, long double vy, long double vz, int pol, long double trajl,
			TGeometry &ageometry, TFieldManager *afield): TParticle(NAME_ELECTRON, -ele_e, m_e, 0){
		InitV(number, t, atau, x, y, z, vx, vy, vz, pol, trajl, ageometry, afield);
	}

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
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the electron is leaving
	 * @param entering Solid that the electron is entering (can be modified by method)
	 * @return Returns true if particle path was changed
	 */
	bool OnHit(long double x1, long double y1[6], long double &x2, long double y2[6], long double normal[3], const solid *leaving, const solid *&entering){
		return false;
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
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(long double x1, long double y1[6], long double &x2, long double y2[6], const solid *currentsolid){
		if (currentsolid->ID != geom->defaultsolid.ID){
			x2 = x1;
			for (int i = 0; i < 6; i++)
				y2[i] = y1[i];
			ID = currentsolid->ID;
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
