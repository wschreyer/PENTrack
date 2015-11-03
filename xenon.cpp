/**
 * \file
 * Xenon class definition. Xenon-129 is used as a comagnetometer used in the TRIUMF EDM experiment. 
 */

#include "globals.h"
#include "xenon.h"

const char* NAME_XENON = "xenon";

ofstream TXenon::endout; ///< endlog file stream
ofstream TXenon::snapshotout; ///< snapshot file stream
ofstream TXenon::trackout; ///< tracklog file stream
ofstream TXenon::hitout; ///< hitlog file stream
ofstream TXenon::spinout; ///< spinlog file stream


TXenon::TXenon(int number, double t, double x, double y, double z, double E, double phi, double theta, int polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_XENON, 0,  m_xe, mu_xeSI, gamma_xe, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TXenon::OnHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation,
			const double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
	traversed = true;
	trajectoryaltered = false;
}


bool TXenon::OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, int &polarisation, solid currentsolid){
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

//Xenon does not decay
void TXenon::Decay() {}
