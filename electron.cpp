/**
 * \file
 * Electron class definition.
 */

#include "globals.h"
#include "electron.h"

const char* NAME_ELECTRON = "electron";

ofstream TElectron::endout; ///< endlog file stream
ofstream TElectron::snapshotout; ///< snapshot file stream
ofstream TElectron::trackout; ///< tracklog file stream
ofstream TElectron::hitout; ///< hitlog file stream
ofstream TElectron::spinout; ///< spinlog file stream


TElectron::TElectron(int number, double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_ELECTRON, -ele_e, m_e, 0, 0, number, t, x, y, z, E, phi, theta, amc, geometry, afield){

}


void TElectron::OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
			const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
	traversed = true;
	trajectoryaltered = false;
}


bool TElectron::OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation, solid currentsolid){
	if (currentsolid.ID != geom->defaultsolid.ID){
		x2 = x1;
		for (int i = 0; i < 6; i++)
			y2[i] = y1[i];
		StopIntegration(ID_ABSORBED_IN_MATERIAL, x2, y2, polarisation, currentsolid);
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
	}*/
	return false;
}


void TElectron::Decay(){

}
