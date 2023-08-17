/**
 * \file
 * Xenon class definition. Xenon-129 is used as a comagnetometer used in the TRIUMF EDM experiment. 
 */

#include "xenon.h"

#include "globals.h"

const char* NAME_XENON = "xenon";

TXenon::TXenon(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
			: TParticle(NAME_XENON, 0,  m_xe, mu_xeSI, gamma_xe, number, t, x, y, z, E, phi, theta, polarisation, spinprojection, amc, geometry, afield){

}

void TXenon::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected
	std::uniform_real_distribution<double> unidist(0, 1);
	double prob = unidist(mc);
	material mat = vnormal < 0 ? entering.mat : leaving.mat;
	double diffprob = mat.DiffProb;
//	cout << "prob: " << diffprob << '\n';
	
	if (prob >= diffprob){
		//************** specular reflection **************
//				printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
		x2 = x1;
		y2 = y1;
		y2[3] -= 2*vnormal*normal[0]; // reflect velocity
		y2[4] -= 2*vnormal*normal[1];
		y2[5] -= 2*vnormal*normal[2];
	}
	else{
		//************** diffuse reflection no MR model ************
		double phi_r, theta_r;	
//		std::cout << "Diffuse reflection!" << std::endl;
		std::uniform_real_distribution<double> phidist(0., 2.*pi);
		phi_r = phidist(mc);
		std::sincos_distribution<double> sincosdist(0., 0.5*pi);
		theta_r = sincosdist(mc);

		if (vnormal > 0) theta_r = pi - theta_r; // if velocity points out of volume invert polar angle
		x2 = x1;
		y2 = y1;
		double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		y2[3] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
		y2[4] = vabs*sin(phi_r)*sin(theta_r);
		y2[5] = vabs*cos(theta_r);
		RotateVector(&y2[3], normal, &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	}

	if (unidist(mc) < entering.mat.SpinflipProb){
		y2[7] *= -1;
	}
} //end OnHit method

//do nothing for each for step
void TXenon::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
		const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{

}

//Xenon does not decay
void TXenon::Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{

}
