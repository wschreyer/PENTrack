/**
 * \file
 * Mercury class definition. Mercury-199 is used as a comagnetometer for the EDM experiment at TRIUMF. 
 */

#include "mercury.h"

#include "globals.h"

const char* NAME_MERCURY = "mercury";

TMercury::TMercury(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
			: TParticle(NAME_MERCURY, 0, m_hg, mu_hgSI, gamma_hg, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}


void TMercury::OnHit(TStep &stepper, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	std::array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected
	std::uniform_real_distribution<double> unidist(0, 1);
	double prob = unidist(mc);
	material mat = vnormal < 0 ? entering.mat : leaving.mat;
	double diffprob = mat.DiffProb;
//	cout << "prob: " << diffprob << '\n';
	
	if (prob >= diffprob){
		//************** specular reflection **************
//				printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
		stepper.SetStepEnd(stepper.GetStartTime());
		v1[0] -= 2*vnormal*normal[0]; // reflect velocity
		v1[1] -= 2*vnormal*normal[1];
		v1[2] -= 2*vnormal*normal[2];
		stepper.SetVelocity(v1);
	}
	else{
		//************** diffuse reflection no MR model ************
		double phi_r, theta_r;	
		
		std::uniform_real_distribution<double> phidist(0., 2.*pi);
		phi_r = phidist(mc);
		std::sincos_distribution<double> sincosdist(0., 0.5*pi);
		theta_r = sincosdist(mc);
		
		if (vnormal > 0) theta_r = pi - theta_r; // if velocity points out of volume invert polar angle
		stepper.SetStepEnd(stepper.GetStartTime());
		double vabs = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
		std::array<double, 3> v2 = {vabs*cos(phi_r)*sin(theta_r), vabs*sin(phi_r)*sin(theta_r), vabs*cos(theta_r)};	// new velocity with respect to z-axis
		RotateVector(&v2[0], normal, &v1[0]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	}

	if (unidist(mc) < entering.mat.SpinflipProb){
		stepper.SetPolarization(-stepper.GetPolarization());
	}
}


//do nothing for each for step
void TMercury::OnStep(TStep &stepper,
		const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{

}

//Mercury does not decay
void TMercury::Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{

}

