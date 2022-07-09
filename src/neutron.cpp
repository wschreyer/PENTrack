/**
 * \file
 * Neutron class definition.
 */

#include "neutron.h"
#include "globals.h"
#include "proton.h"
#include "electron.h"

#include <valarray>

using namespace std;

const char* NAME_NEUTRON = "neutron";


TNeutron::TNeutron(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
		: TParticle(NAME_NEUTRON, 0, m_n, mu_nSI, gamma_n, number, t, x, y, z, E, phi, theta, polarisation, spinprojection, amc, geometry, afield){

}


void TNeutron::Transmit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const double Estep) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane

	// specular transmission (refraction)
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
	for (int i = 0; i < 3; i++)
		y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)
}


void TNeutron::TransmitMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const{

	if (!MR::MRValid(&y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength)){ // check if MicroRoughness model should be applied
		throw runtime_error("Tried to use micro-roughness model in invalid energy regime. That should not happen!");
	}
	
	double theta_t, phi_t;
	std::uniform_real_distribution<double> MRprobdist(0, 1.5 * MR::MRDistMax(true, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength)); // scale up maximum to make sure it lies above all values of scattering distribution
	std::uniform_real_distribution<double> phidist(0, 2.*pi);
	std::sin_distribution<double> sindist(0, pi/2.);
	do{
		phi_t = phidist(mc);
		theta_t = sindist(mc);
	}while (MRprobdist(mc) > MR::MRDist(true, false, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength, theta_t, phi_t));

	double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5] - 2*Estep/m_n);
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	if (vnormal < 0) theta_t = pi - theta_t; // if velocity points into volume invert polar angle

	y2[3] = vabs*cos(phi_t)*sin(theta_t);	// new velocity with respect to z-axis
	y2[4] = vabs*sin(phi_t)*sin(theta_t);
	y2[5] = vabs*cos(theta_t);
	RotateVector(&y2[3], normal, &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
}


void TNeutron::TransmitLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
    std::valarray<double> n(normal, 3);
    if (vnormal < 0){
        vnormal *= -1;
        n *= -1;
    }
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
	std::valarray<double> specular_trans(&y2[3], 3);
	specular_trans += (k2/k1 - 1)*n*vnormal;

	double theta_t, phi_t;

	double vabs = std::sqrt((specular_trans*specular_trans).sum());
	do {
        std::uniform_real_distribution<double> unidist(0, 2.*pi);
        phi_t = unidist(mc);
        std::sincos_distribution<double> sincosdist(0, pi/2.);
        theta_t = sincosdist(mc);
        y2[3] = vabs * cos(phi_t) * sin(theta_t);    // new velocity with respect to z-axis
        y2[4] = vabs * sin(phi_t) * sin(theta_t);
        y2[5] = vabs * cos(theta_t);
        if (mat.DiffProb > 0) {
            RotateVector(&y2[3], &n[0], &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
        } else {
            RotateVector(&y2[3], &specular_trans[0], &y1[3]);
        }
    }while(y2[3]*n[0] + y2[4]*n[1] + y2[5]*n[2] <= 0); // resample until velocity points into correct hemisphere
}


void TNeutron::Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3]) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected

	//************** specular reflection **************
//				printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
	x2 = x1;
	y2[0] = y1[0];
	y2[1] = y1[1];
	y2[2] = y1[2];
	y2[3] -= 2*vnormal*normal[0]; // reflect velocity
	y2[4] -= 2*vnormal*normal[1];
	y2[5] -= 2*vnormal*normal[2];
	y2[6] = y1[6];
	y2[8] = y1[8];
}


void TNeutron::ReflectMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const{
	if (!MR::MRValid(&y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength)){
		throw runtime_error("Tried to use micro-roughness model in invalid energy regime. That should not happen!");
	}

	double phi_r, theta_r;
	std::uniform_real_distribution<double> MRprobdist(0, 1.5 * MR::MRDistMax(true, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength)); // scale up maximum to make sure it lies above all values of scattering distribution
//			cout << "max: " << MRmax << '\n';
	std::uniform_real_distribution<double> unidist(0, 2.*pi);
	std::sin_distribution<double> sindist(0, pi/2.);
	do{
		phi_r = unidist(mc);
		theta_r = sindist(mc);
	}while (MRprobdist(mc) > MR::MRDist(false, false, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength, theta_r, phi_r));

	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	if (vnormal > 0) theta_r = pi - theta_r; // if velocity points out of volume invert polar angle
	x2 = x1;
	y2[0] = y1[0];
	y2[1] = y1[1];
	y2[2] = y1[2];
	double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
	y2[3] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
	y2[4] = vabs*sin(phi_r)*sin(theta_r);
	y2[5] = vabs*cos(theta_r);
	RotateVector(&y2[3], normal, &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	y2[6] = y1[6];
	y2[8] = y1[8];
}


void TNeutron::ReflectLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const material &mat, TMCGenerator &mc) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	std::valarray<double> n(normal, 3);
	if (vnormal > 0){
	    vnormal *= -1;
	    n *= -1;
	}

	x2 = x1;
	y2[0] = y1[0];
	y2[1] = y1[1];
	y2[2] = y1[2];
	double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
	do {
        std::uniform_real_distribution<double> unidist(0, 2.*pi);
        double phi_r = unidist(mc);
        std::sincos_distribution<double> sincosdist(0, pi/2.);
        double theta_r = sincosdist(mc);

        y2[3] = vabs * cos(phi_r) * sin(theta_r);    // new velocity with respect to z-axis
        y2[4] = vabs * sin(phi_r) * sin(theta_r);
        y2[5] = vabs * cos(theta_r);
        if (mat.DiffProb > 0) {
            RotateVector(&y2[3], &n[0], &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
        } else { // if DiffProb == 0 use modified Lambert model
            std::valarray<double> specular_refl(&y1[3],3);
            specular_refl -= 2*vnormal*n;
            RotateVector(&y2[3], &specular_refl[0], &y1[3]); // rotate velocity into coordinate system defined by specularly reflected velocity vector
        }
    }while(y2[3]*n[0] + y2[4]*n[1] + y2[5]*n[2] <= 0); // resample until velocity points into correct reflection hemisphere
	y2[6] = y1[6];
	y2[8] = y1[8];
}


void TNeutron::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
		const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{

    double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
    double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
    material mat = vnormal < 0 ? entering.mat : leaving.mat; // use material properties of the solid whose surface was hit

    std::uniform_real_distribution<double> unidist(0, 1);
    if (unidist(mc) < mat.SpinflipProb){ // should spin be flipped?
        y2[7] *= -1;
    }

    double Estep = CalcPotentialStep(leaving.mat, entering.mat, y2);

//		cout << "Leaving " << leaving->ID << " Entering " << entering->ID << " Enormal = " << Enormal << " Estep = " << Estep;

    bool UseMRModel = MR::MRValid(&y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength);
	double MRreflprob = 0, MRtransprob = 0;
	if (UseMRModel){ 	// handle MicroRoughness reflection/transmission separately
		MRreflprob = MR::MRProb(false, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength);
		if (GetKineticEnergy(&y1[3]) > Estep) // MicroRoughness transmission can happen if neutron energy > potential step
			MRtransprob = MR::MRProb(true, &y1[3], normal, Estep, mat.RMSRoughness, mat.CorrelLength);
	}
	double prob = unidist(mc);
	if (UseMRModel && prob < MRreflprob){
		ReflectMR(x1, y1, x2, y2, normal, Estep, mat, mc);
	}
	else if (UseMRModel && prob < MRreflprob + MRtransprob){
		TransmitMR(x1, y1, x2, y2, normal, Estep, mat, mc);
	}

	else{
		complex<double> k1 = sqrt(complex<double>(Enormal, -leaving.mat.FermiImag*1e-9)); // wavenumber in first solid
		complex<double> k2 = sqrt(complex<double>(Enormal - Estep, -entering.mat.FermiImag*1e-9)); // wavenumber in second solid
		double reflprob = norm((k1 - k2)/(k1 + k2)); // specular reflection probability
		if (Enormal > Estep){ // transmission only possible if Enormal > Estep
			if (prob < MRreflprob + MRtransprob + reflprob*(1 - MRreflprob - MRtransprob)){ // reflection, scale down reflprob so MRreflprob + MRtransprob + reflprob + transprob = 1
				if (!UseMRModel && unidist(mc) < mat.DiffProb + mat.ModifiedLambertProb){
					ReflectLambert(x1, y1, x2, y2, normal, mat, mc); // Lambert reflection
				}
				else{
					Reflect(x1, y1, x2, y2, normal); // specular reflection
				}
			}
			else{
				if (!UseMRModel && unidist(mc) < mat.DiffProb + mat.ModifiedLambertProb){
					TransmitLambert(x1, y1, x2, y2, normal, Estep, mat, mc); // Lambert transmission
				}
				else{
					Transmit(x1, y1, x2, y2, normal, Estep); // specular transmission
				}
			}
		}
		else{ // total reflection (Enormal < Estep)
			double absprob = 1 - reflprob + mat.LossPerBounce; // absorption probability during total reflection, add loss per bounce

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
				if (!UseMRModel && unidist(mc) < mat.DiffProb + mat.ModifiedLambertProb){
					ReflectLambert(x1, y1, x2, y2, normal, mat, mc); // Lambert reflection
				}
				else{
					Reflect(x1, y1, x2, y2, normal); // specular reflection
				}
			}
		}
	}

}


void TNeutron::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
					const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	if (currentsolid.mat.FermiImag > 0){
		complex<double> E(0.5*(double)m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]), currentsolid.mat.FermiImag*1e-9); // E + i*W
		complex<double> k = sqrt(2*(double)m_n*E)*(double)ele_e/(double)hbar; // wave vector
		double l = sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2)); // travelled length
		std::exponential_distribution<double> expdist(2*imag(k));
		double abspath = expdist(mc);
		if (abspath < l){
			x2 = x1 + abspath/l*(x2 - x1); // if absorbed, interpolate stopping time and position
			stepper.calc_state(x2, y2);
			ID = ID_ABSORBED_IN_MATERIAL;
//			printf("Absorption!\n");
		}
	}
}


void TNeutron::Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const{
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
	vector<double> beta({y[3]/static_cast<double>(c_0),
						y[4]/static_cast<double>(c_0),
						y[5]/static_cast<double>(c_0)});
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

	secondaries.push_back(new TProton(GetParticleNumber(), t, y[0], y[1], y[2], E_p, phi_p, theta_p, pol_p, pol_p, mc, geom, field));
	secondaries.push_back(new TElectron(GetParticleNumber(), t, y[0], y[1], y[2], E_e, phi_e, theta_e, pol_e, pol_e, mc, geom, field));
}


double TNeutron::GetPotentialEnergy(const value_type t, const state_type &y, const TFieldManager &field, const solid &sld) const{
    return TParticle::GetPotentialEnergy(t, y, field, sld) + MaterialPotential(sld.mat, y);
}
