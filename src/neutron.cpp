/**
 * \file
 * Neutron class definition.
 */

#include "neutron.h"
#include "globals.h"
#include "proton.h"
#include "electron.h"
#include "scattering.h"
#include "vectormath.h"

using namespace std;

const char* NAME_NEUTRON = "neutron";


TNeutron::TNeutron(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
		: TParticle(NAME_NEUTRON, 0, m_n, mu_nSI, gamma_n, number, t, x, y, z, E, phi, theta, polarisation, spinprojection, amc, geometry, afield){

}


void TNeutron::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const std::array<double, 3> &normal,
		const solid &leaving, const solid &entering, const std::array<double, 3> &surfaceVelocity, TMCGenerator &mc,
		stopID &ID, std::vector<TParticle*> &secondaries) const{
	std::array<double, 3> v1{y1[3], y1[4], y1[5]};
    double v1normal = dot(v1 - surfaceVelocity, normal); // velocity normal to reflection plane
    material mat = v1normal < 0 ? entering.mat : leaving.mat; // use material properties of the solid whose surface was hit

    if (std::generate_canonical<double, std::numeric_limits<double>::digits>(mc) < mat.SpinflipProb){ // should spin be flipped?
        y2[7] *= -1;
    }

	double E = 0.5*m_n*mag_sqr(v1 - surfaceVelocity);
	std::complex<double> FermiPotential1(MaterialPotential(leaving.mat, y1), -leaving.mat.FermiImag*1e-9);
	std::complex<double> FermiPotential2(MaterialPotential(entering.mat, y2), -entering.mat.FermiImag*1e-9);
	std::array<double, 3> v2;
	if (mat.RMSRoughness != 0 and mat.CorrelLength != 0){
		microroughness_scattering_distribution MRdist(E, FermiPotential1, FermiPotential2, mat.RMSRoughness, mat.CorrelLength, mat.LossPerBounce);
		if (mat.microfacetDistributionWidth > 0){
			double lambda_c = 2.*M_PI*hbar/GetMass()/ele_e/std::abs(v1normal);
			double beckmannWidth = mat.microfacetDistributionWidth*std::pow(lambda_c*1e9, mat.microfacetDistributionWidthExponent);
			microfacet_scattering_distribution mfDist(MRdist, beckmannWidth);
			v2 = scattered_vector(v1, normal, surfaceVelocity, mfDist, mc);
		}
		else{
			v2 = scattered_vector(v1, normal, surfaceVelocity, MRdist, mc);
		}
	}
	else{
		lambert_scattering_distribution lambertDist(E, FermiPotential1, FermiPotential2, mat.DiffProb, mat.LossPerBounce);
		if (mat.microfacetDistributionWidth > 0){
			double lambda_c = 2.*M_PI*hbar/GetMass()/ele_e/std::abs(v1normal);
			double beckmannWidth = mat.microfacetDistributionWidth*std::pow(lambda_c*1e9, mat.microfacetDistributionWidthExponent);
//			std::cout << lambda_c << " " << beckmannWidth << "\n";
			microfacet_scattering_distribution mfDist(lambertDist, beckmannWidth);
			v2 = scattered_vector(v1, normal, surfaceVelocity, mfDist, mc);
		}
		else{
			v2 = scattered_vector(v1, normal, surfaceVelocity, lambertDist, mc);
		}
	}

	double v2normal = dot(v2 - surfaceVelocity, normal);
	if (v1normal*v2normal <= 0){
		x2 = x1;
		y2[0] = y1[0];
		y2[1] = y1[1];
		y2[2] = y1[2];
		y2[3] = v2[0];
		y2[4] = v2[1];
		y2[5] = v2[2];
		y2[6] = y1[6];
		y2[8] = y1[8];
	}
	else{
		if (0.5*m_n*v1normal*v1normal <= std::real(FermiPotential2 - FermiPotential1)){
			ID = ID_ABSORBED_ON_SURFACE;
		}
		else{
			y2[3] = v2[0];
			y2[4] = v2[1];
			y2[5] = v2[2];
		}
	}
}


void TNeutron::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
					const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const{
	if (currentsolid.mat.FermiImag > 0){
		complex<double> E(0.5*(double)m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]), currentsolid.mat.FermiImag*1e-9); // E + i*W
		complex<double> k = sqrt(2*(double)m_n*E)*(double)ele_e/(double)hbar; // wave vector
		double l = sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2)); // travelled length
		exponential_distribution<double> expdist(2*imag(k));
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
 	array<double, 4> pp, pe; // four-momenta of proton and electron

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
	BOOST(array<double, 3>({virt[1]/virt[0], virt[2]/virt[0], virt[3]/virt[0]}), pe);

//-------- Step 7 ----------------------------------------------------------------------------------------------------------
	// get 4-momentum of neutrino via 4-momentum conservation
	std::array<double, 4> pnue({virt[0] - pe[0], virt[1] - pe[1], virt[2] - pe[2], virt[3] - pe[3]});

//-------- Step 8 ----------------------------------------------------------------------------------------------------------
	// boost p,e,nu into moving frame of neutron
	array<double, 3> beta({y[3]/static_cast<double>(c_0),
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
