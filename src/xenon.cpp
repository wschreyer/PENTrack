/**
 * \file
 * Xenon class definition. Xenon-129 is used as a comagnetometer used in the TRIUMF EDM experiment. 
 */

#include "globals.h"
#include "xenon.h"

const char* NAME_XENON = "xenon";

std::ofstream TXenon::endout; ///< endlog file stream
std::ofstream TXenon::snapshotout; ///< snapshot file stream
std::ofstream TXenon::trackout; ///< tracklog file stream
std::ofstream TXenon::hitout; ///< hitlog file stream
std::ofstream TXenon::spinout; ///< spinlog file stream

TXenon::TXenon(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
			: TParticle(NAME_XENON, 0,  m_xe, mu_xeSI, gamma_xe, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}

void TXenon::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed){
	
	trajectoryaltered = false;
	traversed = true;
	
	Reflect(x1, y1, x2, y2, normal, leaving, entering, trajectoryaltered, traversed);
	
} //end OnHit method

void TXenon::Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const solid &leaving, const solid &entering, bool &trajectoryaltered, bool &traversed){
	
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected
	double prob = mc->UniformDist(0,1);
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
		phi_r = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
		theta_r = mc->SinCosDist(0, 0.5*pi);
		
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

	if (mc->UniformDist(0,1) < entering.mat.SpinflipProb){
		y2[7] *= -1;
		Nspinflip++;
	}

	traversed = false;
	trajectoryaltered = true;
} // end reflect method

//do nothing for each for step
bool TXenon::OnStep(value_type x1, state_type y1, value_type &x2, state_type &y2, solid currentsolid){ return false; }

//Xenon does not decay
void TXenon::Decay() {}
