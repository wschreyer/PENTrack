/**
 * \file
 * Proton class definition.
 */

#include "globals.h"
#include "proton.h"

const char* NAME_PROTON = "proton";

std::ofstream TProton::endout; ///< endlog file stream
std::ofstream TProton::snapshotout; ///< snapshot file stream
std::ofstream TProton::trackout; ///< tracklog file stream
std::ofstream TProton::hitout; ///< hitlog file stream
std::ofstream TProton::spinout; ///< spinlog file stream

TProton::TProton(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: TParticle(NAME_PROTON, ele_e, m_p, 0, 0, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TProton::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed){
	traversed = true;
	trajectoryaltered = false;
}


bool TProton::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper, const solid &currentsolid){
	if (currentsolid.ID != geom->defaultsolid.ID){
		x2 = x1;
		y2 = y1;
		ID = ID_ABSORBED_IN_MATERIAL;
		printf("Absorption!\n");
		return true;
	}
	return false;
}


void TProton::Decay(){

};
