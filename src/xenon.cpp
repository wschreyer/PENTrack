/**
 * \file
 * Xenon class definition. Xenon-129 is used as a comagnetometer used in the TRIUMF EDM experiment. 
 */

#include "xenon.h"

#include "globals.h"
#include "scattering.h"

const char* NAME_XENON = "xenon";

TXenon::TXenon(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
			: TParticle(NAME_XENON, 0,  m_xe, mu_xeSI, gamma_xe, number, t, x, y, z, E, phi, theta, polarisation, spinprojection, amc, geometry, afield){

}

void TXenon::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	double v1normal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	material mat = v1normal < 0 ? entering.mat : leaving.mat;

	if (std::generate_canonical<double, std::numeric_limits<double>::digits>(mc) < mat.SpinflipProb){
		y2[7] *= -1;
	}

	std::array<double, 3> v1{y1[3], y1[4], y1[5]};
	std::array<double, 3> n{normal[0], normal[1], normal[2]};
	lambert_scattering_distribution<double> scatteringDist(mat.DiffProb, 0, mat.LossPerBounce);
	std::array<double, 3> v2 = scattered_vector(v1, n, scatteringDist, mc);

	double v2normal = n[0]*v2[0] + n[1]*v2[1] + n[2]*v2[2];
	if (v1normal * v2normal <= 0){
		x2 = x1;
		y2[0] = y1[0];
		y2[1] = y1[1];
		y2[2] = y2[2];
		y2[3] = v2[0];
		y2[4] = v2[1];
		y2[5] = v2[2];
		y2[6] = y1[6];
		y2[8] = y1[8];
	}
	else{
		ID = ID_ABSORBED_IN_MATERIAL;
	}
} //end OnHit method

//do nothing for each for step
void TXenon::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
		const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{

}

//Xenon does not decay
void TXenon::Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{

}
