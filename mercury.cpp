/**
 * \file
 * Mercury class definition. Mercury-199 is used as a comagnetometer for the EDM experiment at TRIUMF. 
 */

#include "globals.h"
#include "mercury.h"

const char* NAME_MERCURY = "mercury";

ofstream TMercury::endout; ///< endlog file stream
ofstream TMercury::snapshotout; ///< snapshot file stream
ofstream TMercury::trackout; ///< tracklog file stream
ofstream TMercury::hitout; ///< hitlog file stream
ofstream TMercury::spinout; ///< spinlog file stream


TMercury::TMercury(int number, double t, double x, double y, double z, double E, double phi, double theta, int polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_MERCURY, 0,  m_hg, mu_hgSI, gamma_hg, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TMercury::OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
			const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
	traversed = true;
	trajectoryaltered = false;
}


bool TMercury::OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation, solid currentsolid){
	if (currentsolid.ID != geom->defaultsolid.ID){
		x2 = x1;
		for (int i = 0; i < 6; i++)
			y2[i] = y1[i];
		StopIntegration(ID_ABSORBED_IN_MATERIAL, x2, y2, polarisation, currentsolid);
		printf("Absorption!\n");
		return true;
	}

	return false;
}

//Mercury does not decay
void TMercury::Decay(){}
