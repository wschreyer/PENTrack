/**
 * \file
 * Proton class definition.
 */

#include "proton.h"

#include "globals.h"

const char* NAME_PROTON = "proton";

TProton::TProton(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
		: TParticle(NAME_PROTON, ele_e, m_p, 0, 0, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TProton::OnHit(TStep &stepper, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{

}


void TProton::OnStep(TStep &stepper,
		const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	if (currentsolid.ID > 1){
		stepper.SetStepEnd(stepper.GetStartTime());
		ID = ID_ABSORBED_IN_MATERIAL;
//		printf("Absorption!\n");
	}
}


void TProton::Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{

};
