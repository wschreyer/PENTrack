/**
 * \file
 * Neutron class definition.
 */

#include "neutron.h"
#include "globals.h"
#include "proton.h"
#include "electron.h"

using namespace std;

const char* NAME_NEUTRON = "neutron";

TNeutron::TNeutron(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
		: TParticle(NAME_NEUTRON, 0, m_n, mu_nSI, gamma_n, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}

void TNeutron::Transmit(TStep &stepper,const double normal[3], const solid &leaving, const solid &entering) const{
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane

	// specular transmission (refraction)
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
    	double k2 = sqrt(Enormal - CalcPotentialStep(leaving.mat, entering.mat, stepper)); // wavenumber in second solid (use only real part for transmission!)
	array<double, 3> v2 = stepper.GetVelocity();
	for (int i = 0; i < 3; i++){
		v2[i] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)
	}
	stepper.SetVelocity(v2);
}

void TNeutron::TransmitMR(TStep &stepper, const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const{
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	if (!MR::MRValid(&v1[0], normal, leaving, entering)){ // check if MicroRoughness model should be applied
		throw runtime_error("Tried to use micro-roughness model in invalid energy regime. That should not happen!");
	}
	
	double theta_t, phi_t;
	std::uniform_real_distribution<double> MRprobdist(0, 1.5 * MR::MRDistMax(true, &v1[0], normal, leaving, entering)); // scale up maximum to make sure it lies above all values of scattering distribution
	std::uniform_real_distribution<double> phidist(0, 2.*pi);
	std::sin_distribution<double> sindist(0, pi/2.);
	do{
		phi_t = phidist(mc);
		theta_t = sindist(mc);
	}while (MRprobdist(mc) > MR::MRDist(true, false, &v1[0], normal, leaving, entering, theta_t, phi_t));
	double Estep = CalcPotentialStep(leaving.mat, entering.mat, stepper);
	double E = GetKineticEnergy(v1) - Estep;
	double vabs = sqrt(2.*E/GetMass());

	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane
	if (vnormal < 0) theta_t = pi - theta_t; // if velocity points into volume invert polar angle
	array<double, 3> v2;
	v2[0] = vabs*cos(phi_t)*sin(theta_t);	// new velocity with respect to z-axis
	v2[1] = vabs*sin(phi_t)*sin(theta_t);
	v2[2] = vabs*cos(theta_t);
	RotateVector(&v2[0], normal, &v1[0]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
	stepper.SetVelocity(v2);
}

void TNeutron::TransmitLambert(TStep &stepper, const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const{
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - CalcPotentialStep(leaving.mat, entering.mat, stepper)); // wavenumber in second solid (use only real part for transmission!)
	array<double, 3> v2 = stepper.GetVelocity();
	for (int i = 0; i < 3; i++)
		v2[i] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)

	double theta_t, phi_t;
	std::uniform_real_distribution<double> unidist(0, 2.*pi);
	phi_t = unidist(mc);
	std::sincos_distribution<double> sincosdist(0, pi/2.);
	theta_t = sincosdist(mc);

	if (vnormal < 0) theta_t = pi - theta_t; // if velocity points into volume invert polar angle
	double vabs = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
	v2[0] = vabs*cos(phi_t)*sin(theta_t);	// new velocity with respect to z-axis
	v2[1] = vabs*sin(phi_t)*sin(theta_t);
	v2[2] = vabs*cos(theta_t);
	RotateVector(&v2[0], normal, &v1[0]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
	stepper.SetVelocity(v2);
}

void TNeutron::Reflect(TStep &stepper, const double normal[3], const solid &leaving, const solid &entering) const{
	//particle was neither transmitted nor absorbed, so it has to be reflected
	//************** specular reflection **************
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane

	stepper.SetStepEnd(stepper.GetStartTime());
	array<double, 3> v2 = stepper.GetVelocity();
	v2[0] -= 2*vnormal*normal[0]; // reflect velocity
	v2[1] -= 2*vnormal*normal[1];
	v2[2] -= 2*vnormal*normal[2];
	stepper.SetVelocity(v2);
}

void TNeutron::ReflectMR(TStep &stepper,
		const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const{
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane
	material mat = vnormal < 0 ? entering.mat : leaving.mat;
	if (!MR::MRValid(&v1[0], normal, leaving, entering)){
		throw runtime_error("Tried to use micro-roughness model in invalid energy regime. That should not happen!");
	}

	double phi_r, theta_r;
	std::uniform_real_distribution<double> MRprobdist(0, 1.5 * MR::MRDistMax(true, &v1[0], normal, leaving, entering)); // scale up maximum to make sure it lies above all values of scattering distribution
//			cout << "max: " << MRmax << '\n';
	std::uniform_real_distribution<double> unidist(0, 2.*pi);
	std::sin_distribution<double> sindist(0, pi/2.);
	do{
		phi_r = unidist(mc);
		theta_r = sindist(mc);
	}while (MRprobdist(mc) > MR::MRDist(false, false, &v1[0], normal, leaving, entering, theta_r, phi_r));

	if (vnormal > 0) theta_r = pi - theta_r; // if velocity points out of volume invert polar angle
	stepper.SetStepEnd(stepper.GetStartTime());
	array<double, 3> v2 = stepper.GetVelocity();
	double vabs = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
	v2[0] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
	v2[1] = vabs*sin(phi_r)*sin(theta_r);
	v2[2] = vabs*cos(theta_r);
	RotateVector(&v2[0], normal, &v1[0]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	stepper.SetVelocity(v2);
}

void TNeutron::ReflectLambert(TStep &stepper,
		const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const{
	//particle was neither transmitted nor absorbed, so it has to be reflected
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane

	double phi_r, theta_r;
	std::uniform_real_distribution<double> unidist(0, 2.*pi);
	phi_r = unidist(mc);
	std::sincos_distribution<double> sincosdist(0, pi/2.);
	theta_r = sincosdist(mc);

	if (vnormal > 0) theta_r = pi - theta_r; // if velocity points out of volume invert polar angle
	stepper.SetStepEnd(stepper.GetStartTime());
	array<double, 3> v2 = stepper.GetVelocity();
	double vabs = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
	v2[0] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
	v2[1] = vabs*sin(phi_r)*sin(theta_r);
	v2[2] = vabs*cos(theta_r);
	RotateVector(&v2[0], normal, &v1[0]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	stepper.SetVelocity(v2);
}

void TNeutron::OnHit(TStep &stepper, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	array<double, 3> v1 = stepper.GetVelocity(stepper.GetStartTime());
	double vnormal = v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]; // velocity normal to reflection plane
    double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
    material mat = vnormal < 0 ? entering.mat : leaving.mat; // use material properties of the solid whose surface was hit

    std::uniform_real_distribution<double> unidist(0, 1);
    if (unidist(mc) < mat.SpinflipProb){ // should spin be flipped?
		stepper.SetPolarization(-stepper.GetPolarization());
    }

    double Estep = CalcPotentialStep(leaving.mat, entering.mat, stepper);

//		cout << "Leaving " << leaving->ID << " Entering " << entering->ID << " Enormal = " << Enormal << " Estep = " << Estep;

    bool UseMRModel = MR::MRValid(&v1[0], normal, leaving, entering);
	double MRreflprob = 0, MRtransprob = 0;
	if (UseMRModel){ 	// handle MicroRoughness reflection/transmission separately
		MRreflprob = MR::MRProb(false, &v1[0], normal, leaving, entering);
		if (GetKineticEnergy(v1) > Estep) // MicroRoughness transmission can happen if neutron energy > potential step
			MRtransprob = MR::MRProb(true, &v1[0], normal, leaving, entering);
	}
	double prob = unidist(mc);
	if (UseMRModel && prob < MRreflprob){
		ReflectMR(stepper, normal, leaving, entering, mc);
	}
	else if (UseMRModel && prob < MRreflprob + MRtransprob){
		TransmitMR(stepper, normal, leaving, entering, mc);
	}

	else{
		complex<double> k1 = sqrt(complex<double>(Enormal, -leaving.mat.FermiImag*1e-9)); // wavenumber in first solid
		complex<double> k2 = sqrt(complex<double>(Enormal - Estep, -entering.mat.FermiImag*1e-9)); // wavenumber in second solid
		double reflprob = norm((k1 - k2)/(k1 + k2)); // specular reflection probability
		if (Enormal > Estep){ // transmission only possible if Enormal > Estep
			if (prob < MRreflprob + MRtransprob + reflprob*(1 - MRreflprob - MRtransprob)){ // reflection, scale down reflprob so MRreflprob + MRtransprob + reflprob + transprob = 1
				if (!UseMRModel && unidist(mc) < mat.DiffProb){
					ReflectLambert(stepper, normal, leaving, entering, mc); // Lambert reflection
				}
				else{
					Reflect(stepper, normal, leaving, entering); // specular reflection
				}
			}
			else{
				if (!UseMRModel && unidist(mc) < mat.DiffProb){
					TransmitLambert(stepper, normal, leaving, entering, mc); // Lambert transmission
				}
				else{
					Transmit(stepper, normal, leaving, entering); // specular transmission
				}
			}
		}
		else{ // total reflection (Enormal < Estep)
			double absprob = 1 - reflprob; // absorption probability

			if (UseMRModel){
				double kc = sqrt(2*m_n*Estep)*ele_e/hbar;
				double addtrans = 2*pow(entering.mat.RMSRoughness, 2)*kc*kc/(1 + 0.85*kc*entering.mat.CorrelLength + 2*kc*kc*pow(entering.mat.CorrelLength, 2));
				absprob *= sqrt(1 + addtrans); // second order correction for reflection on MicroRoughness surfaces
			}
	//			cout << " ReflProb = " << reflprob << '\n';

			if (prob < MRreflprob + MRtransprob + absprob*(1 - MRreflprob - MRtransprob)){ // -> absorption on reflection, scale down absprob so MRreflprob + MRtransprob + absprob + reflprob = 1
				ID = ID_ABSORBED_ON_SURFACE;
			}
			else{ // no absorption -> reflection
				if (!UseMRModel && unidist(mc) < mat.DiffProb){
					ReflectLambert(stepper, normal, leaving, entering, mc); // Lambert reflection
				}
				else{
					Reflect(stepper, normal, leaving, entering); // specular reflection
				}
			}
		}
	}

}


void TNeutron::OnStep(TStep &stepper,
					const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	if (currentsolid.mat.FermiImag > 0){
		complex<double> E(GetKineticEnergy(stepper.GetVelocity()), currentsolid.mat.FermiImag*1e-9); // E + i*W
		complex<double> k = sqrt(2*(double)m_n*E)*(double)ele_e/(double)hbar; // wave vector
		double l2 = stepper.GetPathLength();
		double l1 = stepper.GetPathLength(stepper.GetStartTime());
		std::exponential_distribution<double> expdist(2*imag(k));
		double abspath = expdist(mc);
		if (abspath < l2 - l1){
			stepper.SetStepEndToMatch([&](const double& t){ return stepper.GetPathLength(t); }, l1 + abspath);
			ID = ID_ABSORBED_IN_MATERIAL;
//			printf("Absorption!\n");
		}
	}
}


void TNeutron::Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{
	double E_p, E_e, phi_p, phi_e, theta_p, theta_e;
	double pol_p, pol_e;
	double m_nue = 1 / pow(c_0, 2); // [eV/c^2]

	uniform_real_distribution<double> phi_dist(0, 2*pi);
	sin_distribution<double> sin_dist(0., pi/2.);
 	double pabs, beta1, beta2, delta1, delta2; // 3-momentums (abs+directions)
 	vector<double> pp(4), pe(4); // four-momenta of proton and electron

//-------- Step 1 ----------------------------------------------------------------------------------------------------------
/*
 * p   = (E/c, p_x, p_y, p_z)
 * E   = E_kin + m*c^2
 * p^2 = (E/c)^2 - (m*c)^2
 *
 */
	// energy of proton
 	double Ekin = proton_beta_distribution(mc);
	pp[0] = Ekin/c_0 + m_p*c_0;
	// momentum norm of proton
	pabs = sqrt(pp[0]*pp[0] - m_p*m_p*c_0*c_0);

//-------- Step 2 ----------------------------------------------------------------------------------------------------------
	// isotropic emission characteristics (in the rest frame of the neutron)
	beta1 = phi_dist(mc);
	delta1 = sin_dist(mc);
	// 3-momentum of the proton
	pp[1] = pabs * sin(delta1) * cos(beta1);
	pp[2] = pabs * sin(delta1) * sin(beta1);
	pp[3] = pabs * cos(delta1);

//-------- Step 3 ----------------------------------------------------------------------------------------------------------
	// calculate intermediate virtual state via 4-momentum conservation
	double virt[4] = {(double)m_n*(double)c_0 - pp[0], -pp[1], -pp[2], -pp[3]}; // 4-momentum of intermediate virtual state (n - p)
	double m2_virt = virt[0]*virt[0] - virt[1]*virt[1] - virt[2]*virt[2] - virt[3]*virt[3]; // squared mass of virtual state


//-------- Step 4 ----------------------------------------------------------------------------------------------------------

	// energy of electron from two-body-decay of virtual state
	pe[0] = (m2_virt + m_e*m_e*c_0*c_0 - m_nue*m_nue*c_0*c_0)/2/sqrt(m2_virt);
	pabs = sqrt(pe[0]*pe[0] - m_e*m_e*c_0*c_0);

//-------- Step 5 ----------------------------------------------------------------------------------------------------------
	// isotropic emission characteristics (in the rest frame of the virtual state)
	beta2 = phi_dist(mc);
	delta2 = sin_dist(mc);
	// 3-momentum of the electron
	pe[1] = pabs * sin(delta2) * cos(beta2);
	pe[2] = pabs * sin(delta2) * sin(beta2);
	pe[3] = pabs * cos(delta2);

//-------- Step 6 ----------------------------------------------------------------------------------------------------------
	// boost e into moving frame of virtual state
	BOOST(vector<double>({virt[1]/virt[0], virt[2]/virt[0], virt[3]/virt[0]}), pe);

//-------- Step 7 ----------------------------------------------------------------------------------------------------------
	// get 4-momentum of neutrino via 4-momentum conservation
	std::vector<double> pnue({virt[0] - pe[0], virt[1] - pe[1], virt[2] - pe[2], virt[3] - pe[3]});

//-------- Step 8 ----------------------------------------------------------------------------------------------------------
	// boost p,e,nu into moving frame of neutron
	array<double, 3> v = stepper.GetVelocity();
	vector<double> beta({v[0]/static_cast<double>(c_0),
						 v[1]/static_cast<double>(c_0),
						 v[2]/static_cast<double>(c_0)});
	BOOST(beta,pe);
	BOOST(beta,pp);
	BOOST(beta,pnue);

//-------- Step 9 ----------------------------------------------------------------------------------------------------------
	// use neutrino mass for check
//		long double decayerror = m_nue*c_0*c_0 - sqrt(nue[0]*nue[0] - nue[1]*nue[1] - nue[2]*nue[2] - nue[3]*nue[3])*c_0;
//		printf("\n   +++ decay error : %LG eV +++\n", decayerror);


//-------- Finished --------------------------------------------------------------------------------------------------------
	double beta_p2 = (pp[1]*pp[1] + pp[2]*pp[2] + pp[3]*pp[3])/pp[0]/pp[0];
	E_p = m_p*c_0*c_0*beta_p2*(0.5 + beta_p2*(3./8. + beta_p2*(5./16.))); // use series expansion for low energy proton calculations
	E_e = pe[0]*c_0 - m_e*c_0*c_0; // use fully relativistic formula for high energy electron
	phi_p = atan2(pp[2], pp[1]);
	phi_e = atan2(pe[2], pe[1]);
	theta_p = acos(pp[3]/sqrt(pp[1]*pp[1] + pp[2]*pp[2] + pp[3]*pp[3]));
	theta_e = acos(pe[3]/sqrt(pe[1]*pe[1] + pe[2]*pe[2] + pe[3]*pe[3]));
	std::uniform_real_distribution<double> uni_dist(0, 1);
	pol_p = uni_dist(mc) < 0.5 ? 1 : -1;
	pol_e = uni_dist(mc) < 0.5 ? 1 : -1;

	double t = stepper.GetTime();
	array<double, 3> pos = stepper.GetPosition();
	secondaries.push_back(new TProton(GetParticleNumber(), t, pos[0], pos[1], pos[2], E_p, phi_p, theta_p, pol_p, mc, geom, field));
	secondaries.push_back(new TElectron(GetParticleNumber(), t, pos[0], pos[1], pos[2], E_e, phi_e, theta_e, pol_e, mc, geom, field));
}


double TNeutron::GetPotentialEnergy(const double t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const int pol, const TFieldManager &field, const solid &sld) const{
    return TParticle::GetPotentialEnergy(t, pos, v, pol, field, sld) + MaterialPotential(t, pos, v, pol, sld.mat);
}
