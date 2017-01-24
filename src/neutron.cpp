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

ofstream TNeutron::endout; ///< endlog file stream
ofstream TNeutron::snapshotout; ///< snapshot file stream
ofstream TNeutron::trackout; ///< tracklog file stream
ofstream TNeutron::hitout; ///< hitlog file stream
ofstream TNeutron::spinout; ///< spinlog file stream

TNeutron::TNeutron(int number, double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: TParticle(NAME_NEUTRON, 0, m_n, mu_nSI, gamma_n, number, t, x, y, z, E, phi, theta, polarisation, amc, geometry, afield){

}

void TNeutron::Transmit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane

	// specular transmission (refraction)
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double Estep = entering.mat.FermiReal*1e-9 - leaving.mat.FermiReal*1e-9;
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
	for (int i = 0; i < 3; i++)
		y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)
}

void TNeutron::TransmitMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	material mat = vnormal < 0 ? entering.mat : leaving.mat;

	if (!mat.UseMRModel || !MR::MRValid(&y1[3], normal, leaving, entering)){ // check if MicroRoughness model should be applied
		std::cout << "Tried to use micro-roughness model in invalid energy regime. That should not happen!\n";
		exit(-1);
	}
	
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double Estep = entering.mat.FermiReal*1e-9 - leaving.mat.FermiReal*1e-9;
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
	for (int i = 0; i < 3; i++)
		y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)

	double theta_t, phi_t;
	double MRmax = 1.5 * MR::MRDistMax(true, &y1[3], normal, leaving, entering); // scale up maximum to make sure it lies above all values of scattering distribution
	do{
		phi_t = GetRandomGenerator()->UniformDist(0, 2*pi);
		theta_t = GetRandomGenerator()->SinDist(0, pi/2);
	}while (GetRandomGenerator()->UniformDist(0, MRmax) > MR::MRDist(true, false, &y1[3], normal, leaving, entering, theta_t, phi_t));

	if (vnormal < 0) theta_t = pi - theta_t; // if velocity points into volume invert polar angle
	double vabs = sqrt(y2[3]*y2[3] + y2[4]*y2[4] + y2[5]*y2[5]);
	y2[3] = vabs*cos(phi_t)*sin(theta_t);	// new velocity with respect to z-axis
	y2[4] = vabs*sin(phi_t)*sin(theta_t);
	y2[5] = vabs*cos(theta_t);
	RotateVector(&y2[3], normal, &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
}

void TNeutron::TransmitLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double Estep = entering.mat.FermiReal*1e-9 - leaving.mat.FermiReal*1e-9;
	double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
	double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
	for (int i = 0; i < 3; i++)
		y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)

	double theta_t, phi_t;
	phi_t = GetRandomGenerator()->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
	theta_t = GetRandomGenerator()->SinCosDist(0, 0.5*pi);

	if (vnormal < 0) theta_t = pi - theta_t; // if velocity points into volume invert polar angle
	double vabs = sqrt(y2[3]*y2[3] + y2[4]*y2[4] + y2[5]*y2[5]);
	y2[3] = vabs*cos(phi_t)*sin(theta_t);	// new velocity with respect to z-axis
	y2[4] = vabs*sin(phi_t)*sin(theta_t);
	y2[5] = vabs*cos(theta_t);
	RotateVector(&y2[3], normal, &y1[3]); // rotate velocity into coordinate system defined by incoming velocity and plane normal
}

void TNeutron::Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected

	//************** specular reflection **************
//				printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
	x2 = x1;
	y2 = y1;
	y2[3] -= 2*vnormal*normal[0]; // reflect velocity
	y2[4] -= 2*vnormal*normal[1];
	y2[5] -= 2*vnormal*normal[2];
}

void TNeutron::ReflectMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	material mat = vnormal < 0 ? entering.mat : leaving.mat;
	if (!mat.UseMRModel || !MR::MRValid(&y1[3], normal, leaving, entering)){
		std::cout << "Tried to use micro-roughness model in invalid energy regime. That should not happen!\n";
		exit(-1);
	}

	double phi_r, theta_r;
	double MRmax = 1.5 * MR::MRDistMax(false, &y1[3], normal, leaving, entering);  // scale up maximum to make sure it lies above all values of scattering distribution
//			cout << "max: " << MRmax << '\n';
	do{
		phi_r = GetRandomGenerator()->UniformDist(0, 2*pi);
		theta_r = GetRandomGenerator()->SinDist(0, pi/2);
	}while (GetRandomGenerator()->UniformDist(0, MRmax) > MR::MRDist(false, false, &y1[3], normal, leaving, entering, theta_r, phi_r));

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

void TNeutron::ReflectLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		const double normal[3], const solid &leaving, const solid &entering) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	//particle was neither transmitted nor absorbed, so it has to be reflected

	double phi_r, theta_r;
	phi_r = GetRandomGenerator()->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
	theta_r = GetRandomGenerator()->SinCosDist(0, 0.5*pi);

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

void TNeutron::OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
		const solid &leaving, const solid &entering, stopID &ID, std::vector<TParticle*> &secondaries) const{
	double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
	double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	double Estep = entering.mat.FermiReal*1e-9 - leaving.mat.FermiReal*1e-9;

	material mat = vnormal < 0 ? entering.mat : leaving.mat; // use material properties of the solid whose surface was hit
	//material *mat = &entering->mat;

//		cout << "Leaving " << leaving->ID << " Entering " << entering->ID << " Enormal = " << Enormal << " Estep = " << Estep;

	bool UseMRModel = mat.UseMRModel && MR::MRValid(&y1[3], normal, leaving, entering);
	double MRreflprob = 0, MRtransprob = 0;
	if (UseMRModel){ 	// handle MicroRoughness reflection/transmission separately
		MRreflprob = MR::MRProb(false, &y1[3], normal, leaving, entering);
		if (0.5*m_n*(y1[0]*y1[0] + y1[1]*y1[1] + y1[2]*y1[2]) > Estep) // MicroRoughness transmission can happen if neutron energy > potential step
			MRtransprob = MR::MRProb(true, &y1[3], normal, leaving, entering);
	}
	double prob = GetRandomGenerator()->UniformDist(0, 1);
	if (UseMRModel && prob < MRreflprob){
		ReflectMR(x1, y1, x2, y2, normal, leaving, entering);
	}
	else if (UseMRModel && prob < MRreflprob + MRtransprob){
		TransmitMR(x1, y1, x2, y2, normal, leaving, entering);
	}

	else{
		complex<double> k1 = sqrt(complex<double>(Enormal, -leaving.mat.FermiImag*1e-9)); // wavenumber in first solid
		complex<double> k2 = sqrt(complex<double>(Enormal - Estep, -entering.mat.FermiImag*1e-9)); // wavenumber in second solid
		double reflprob = norm((k1 - k2)/(k1 + k2)); // specular reflection probability
		if (Enormal > Estep){ // transmission only possible if Enormal > Estep
			if (prob < MRreflprob + MRtransprob + reflprob*(1 - MRreflprob - MRtransprob)){ // reflection, scale down reflprob so MRreflprob + MRtransprob + reflprob + transprob = 1
				if (!UseMRModel && GetRandomGenerator()->UniformDist(0,1) < mat.DiffProb){
					ReflectLambert(x1, y1, x2, y2, normal, leaving, entering); // Lambert reflection
				}
				else{
					Reflect(x1, y1, x2, y2, normal, leaving, entering); // specular reflection
				}
			}
			else{
				if (!UseMRModel && GetRandomGenerator()->UniformDist(0,1) < mat.DiffProb){
					TransmitLambert(x1, y1, x2, y2, normal, leaving, entering); // Lambert transmission
				}
				else{
					Transmit(x1, y1, x2, y2, normal, leaving, entering); // specular transmission
				}
			}
		}
		else{ // total reflection (Enormal < Estep)
			double absprob = 1 - reflprob; // absorption probability

			if (entering.mat.UseMRModel){
				double kc = sqrt(2*m_n*Estep)*ele_e/hbar;
				double addtrans = 2*pow(entering.mat.RMSRoughness, 2)*kc*kc/(1 + 0.85*kc*entering.mat.CorrelLength + 2*kc*kc*pow(entering.mat.CorrelLength, 2));
				absprob *= sqrt(1 + addtrans); // second order correction for reflection on MicroRoughness surfaces
			}
	//			cout << " ReflProb = " << reflprob << '\n';

			if (prob < MRreflprob + MRtransprob + absprob*(1 - MRreflprob - MRtransprob)){ // -> absorption on reflection, scale down absprob so MRreflprob + MRtransprob + absprob + reflprob = 1
				ID = ID_ABSORBED_ON_SURFACE;
			}
			else{ // no absorption -> reflection
				if (!UseMRModel && GetRandomGenerator()->UniformDist(0,1) < mat.DiffProb){
					ReflectLambert(x1, y1, x2, y2, normal, leaving, entering); // Lambert reflection
				}
				else{
					Reflect(x1, y1, x2, y2, normal, leaving, entering); // specular reflection
				}
			}
		}
	}

	if (GetRandomGenerator()->UniformDist(0,1) < entering.mat.SpinflipProb){ // should spin be flipped?
		y2[7] *= -1;
	}
}


void TNeutron::OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
					const solid &currentsolid, stopID &ID, std::vector<TParticle*> &secondaries) const{
	if (currentsolid.mat.FermiImag > 0){
		complex<double> E(0.5*(double)m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]), currentsolid.mat.FermiImag*1e-9); // E + i*W
		complex<double> k = sqrt(2*(double)m_n*E)*(double)ele_e/(double)hbar; // wave vector
		double l = sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2)); // travelled length
		double abspath = GetRandomGenerator()->ExpDist(2*imag(k)); // choose random absorption length with exponential distribution
		if (abspath < l){
			x2 = x1 + abspath/l*(x2 - x1); // if absorbed, interpolate stopping time and position
			stepper.calc_state(x2, y2);
			ID = ID_ABSORBED_IN_MATERIAL;
			printf("Absorption!\n");
		}
	}
}


void TNeutron::Decay(std::vector<TParticle*> &secondaries) const{
	double E_p, E_e, phi_p, phi_e, theta_p, theta_e;
	double pol_p, pol_e;
	TParticle *p;
	GetRandomGenerator()->NeutronDecay(&GetFinalState()[3], E_p, E_e, phi_p, phi_e, theta_p, theta_e, pol_p, pol_e);
	p = new TProton(GetParticleNumber(), GetFinalTime(), GetFinalState()[0], GetFinalState()[1], GetFinalState()[2],
					E_p, phi_p, theta_p, pol_p, *GetRandomGenerator(), *GetGeometry(), GetFieldManager());
	secondaries.push_back(p);
	p = new TElectron(GetParticleNumber(), GetFinalTime(), GetFinalState()[0], GetFinalState()[1], GetFinalState()[2],
						E_e, phi_e, theta_e, pol_e, *GetRandomGenerator(), *GetGeometry(), GetFieldManager());
	secondaries.push_back(p);
}


double TNeutron::GetPotentialEnergy(const value_type t, const state_type &y, TFieldManager *field, const solid &sld) const{
	return TParticle::GetPotentialEnergy(t, y, field, sld) + sld.mat.FermiReal*1e-9;
}
