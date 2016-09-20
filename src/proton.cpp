/**
 * \file
 * Proton class definition.
 */

#include "globals.h"
#include "proton.h"

const char* NAME_PROTON = "proton";

ofstream TProton::endout; ///< endlog file stream
ofstream TProton::snapshotout; ///< snapshot file stream
ofstream TProton::trackout; ///< tracklog file stream
ofstream TProton::hitout; ///< hitlog file stream
ofstream TProton::spinout; ///< spinlog file stream
ofstream TProton::spinout2; ///< spinlog file stream for doing simultaneous anti-parallel Efield spin integration

TProton::TProton(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: TParticle(NAME_PROTON, ele_e, m_p, 0, 0, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TProton::OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2,
					const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
	traversed = true;
	trajectoryaltered = false;
}


bool TProton::OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid){
	if (currentsolid.ID != geom->defaultsolid.ID){
		x2 = x1;
		y2 = y1;
		StopIntegration(ID_ABSORBED_IN_MATERIAL, x2, y2, currentsolid);
		printf("Absorption!\n");
		return true;
	}
	return false;
}


void TProton::Decay(){

};
