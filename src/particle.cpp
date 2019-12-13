/**
 * \file
 * Particle base class definition.
 */

#include <vector>
#include <algorithm>
#include <iomanip>

#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

#include "particle.h"

using namespace std;

double TParticle::GetInitialTotalEnergy(const TGeometry &geom, const TFieldManager &field) const{
	return GetInitialKineticEnergy() + GetPotentialEnergy(tstart, posstart, vstart, polstart, field, geom.GetSolid(tstart, &posstart[0]));
}


double TParticle::GetFinalTotalEnergy(const TGeometry &geom, const TFieldManager &field) const{
	return GetFinalKineticEnergy() + GetPotentialEnergy(tend, posend, vend, polend, field, geom.GetSolid(tend, &posend[0]));
}


double TParticle::GetInitialKineticEnergy() const{
	return GetKineticEnergy(vstart);
}


double TParticle::GetFinalKineticEnergy() const{
	return GetKineticEnergy(vend);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, const int number,
		const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield):
		name(aname),
		q(qq),
		m(mm),
		mu(mumu),
		gamma(agamma),
		particlenumber(number),
		ID(ID_UNKNOWN),
		tstart(t),
		tend(t),
		posstart({x, y, z}),
		posend({x, y, z}),
		pathlength(0),
		propertime(0),
		spinphase(0),
		spinintegrationtime(0),
		Nhit(0),
		Nspinflip(0),
		noflipprob(1),
		Nstep(0)
		{

	// for small velocities Ekin/m is very small and the relativstic claculation beta^2 = 1 - 1/gamma^2 gives large round-off errors
	// the round-off error can be estimated as 2*epsilon
	// if a series expansion to order O(Eoverm^5) has a smaller error then the round-off error, we will use the series expansion
	double Eoverm = E/m/c_0/c_0; // gamma - 1
	double beta2;
	//cout << pow(Eoverm, 6) << " " << ((((6.0*Eoverm - 5.0)*Eoverm + 4.0)*Eoverm - 3.0)*Eoverm + 2.0)*Eoverm << " " << 1.0 - 1.0/(Eoverm + 1)/(Eoverm + 1) << "\n";
	if (pow(Eoverm, 6) < 2*numeric_limits<double>::epsilon()) // if error in series expansion smaller than round-off error
		beta2 = ((((6.0*Eoverm - 5.0)*Eoverm + 4.0)*Eoverm - 3.0)*Eoverm + 2.0)*Eoverm; // use series expansion
	else
		beta2 = 1.0 - 1.0/(Eoverm + 1)/(Eoverm + 1); // relativstic beta^2
	double vabs = c_0*sqrt(beta2);

	vstart[0] = vabs*cos(phi)*sin(theta); // velocity
	vstart[1] = vabs*sin(phi)*sin(theta);
	vstart[2] = vabs*cos(theta);
	vend = vstart;

	if (polarisation < -1 || polarisation > 1)
		throw std::runtime_error("Polarisation has to be between -1 and 1");
	std::polarization_distribution<double> pdist(polarisation);
	polstart = polend = pdist(amc); // choose initial polarisation randomly, weighted by spin projection onto magnetic field

	double B[3];
	afield.BField(x, y, z, t, B);
	double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	if (Babs > 0){
		std::uniform_real_distribution<double> phidist(0., 2.*pi);
		double spinaz = phidist(amc); // random azimuth angle of spin
		spinstart[0] = sqrt(1 - polarisation*polarisation)*sin(spinaz); // set initial spin vector
		spinstart[1] = sqrt(1 - polarisation*polarisation)*cos(spinaz);
		spinstart[2] = polarisation;
		RotateVector(&spinstart[0], B); // rotate initial spin vector such that its z-component is parallel to magnetic field
	}

	spinend = spinstart;

	solidend = solidstart = geometry.GetSolid(t, &posstart[0]); // set to solid with highest priority
	Hmax = GetInitialTotalEnergy(geometry, afield);
}





std::array<double, 3> TParticle::EquationOfMotion(const double &t, const std::array<double, 3>& pos, const std::array<double, 3>& v, const int& polarization, const TFieldManager &field) const{
	double B[3], dBidxj[3][3], E[3], V; // magnetic/electric field and electric potential in lab frame
	if (q != 0 || (mu != 0 && polarization != 0)) // if particle has charge or magnetic moment, calculate magnetic field
		field.BField(pos[0], pos[1], pos[2], t, B, dBidxj);
 	if (q != 0) // if particle has charge caculate electric field
		field.EField(pos[0], pos[1], pos[2], t, V, E);
	return EquationOfMotion(t, pos, v, polarization, B, dBidxj, E);
}

std::array<double, 3> TParticle::EquationOfMotion(const double &t, const std::array<double, 3>& pos, const std::array<double, 3>& v, const int& polarization,
										const double B[3], const double dBidxj[3][3], const double E[3]) const{
	array<double, 3> F = {0,0,0}; // Force in lab frame
	F[2] += -gravconst*m*ele_e; // add gravitation to force
	if (q != 0){
		F[0] += q*(E[0] + v[1]*B[2] - v[2]*B[1]); // add Lorentz-force
		F[1] += q*(E[1] + v[2]*B[0] - v[0]*B[2]);
		F[2] += q*(E[2] + v[0]*B[1] - v[1]*B[0]);
	}
	if (mu != 0 && polarization != 0 && (B[0] != 0 || B[1] != 0 || B[2] != 0)){
		double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
		double dBdxi[3] = {	(B[0]*dBidxj[0][0] + B[1]*dBidxj[1][0] + B[2]*dBidxj[2][0])/Babs,
							(B[0]*dBidxj[0][1] + B[1]*dBidxj[1][1] + B[2]*dBidxj[2][1])/Babs,
							(B[0]*dBidxj[0][2] + B[1]*dBidxj[1][2] + B[2]*dBidxj[2][2])/Babs}; // derivatives of |B|
		F[0] += polarization*mu*dBdxi[0]; // add force on magnetic dipole moment
		F[1] += polarization*mu*dBdxi[1];
		F[2] += polarization*mu*dBdxi[2];
	}
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	double inversegamma = sqrt(1. - v2/(c_0*c_0)); // relativstic factor 1/gamma
	array<double, 3> a;
	a[0] = inversegamma/m/ele_e*(F[0] - (v[0]*v[0]*F[0] + v[0]*v[1]*F[1] + v[0]*v[2]*F[2])/c_0/c_0); // general relativstic equation of motion
	a[1] = inversegamma/m/ele_e*(F[1] - (v[1]*v[0]*F[0] + v[1]*v[1]*F[1] + v[1]*v[2]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
	a[2] = inversegamma/m/ele_e*(F[2] - (v[2]*v[0]*F[0] + v[2]*v[1]*F[1] + v[2]*v[2]*F[2])/c_0/c_0);
	return a;
}





std::array<double, 3> TParticle::SpinPrecessionAxis(const double t, const TStep &stepper, const TFieldManager &field) const{
	double B[3], dBidxj[3][3], V, E[3];
	auto pos = stepper.GetPosition(t);
	field.BField(pos[0], pos[1], pos[2], t, B, dBidxj);
	field.EField(pos[0], pos[1], pos[2], t, V, E);
	auto v = stepper.GetVelocity(t);
	auto a = EquationOfMotion(t, pos, v, stepper.GetPolarization(t), B, dBidxj, E); // calculate velocity and acceleration required for vxE effect and Thomas precession
	return SpinPrecessionAxis(t, B, E, v, a); // calculate precession axis
}

std::array<double, 3> TParticle::SpinPrecessionAxis(const double t, const double B[3], const double E[3], const std::array<double, 3> &v, const std::array<double, 3> &a) const{
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	double gamma_rel = 1./sqrt(1. - v2/c_0/c_0);
	double Bdotv = B[0]*v[0] + B[1]*v[1] + B[2]*v[2];
	double Bparallel[3];
	for (int j = 0; j < 3; j++)
		Bparallel[j] = Bdotv*v[j]/v2; // magnetic-field component parallel to velocity

	// spin precession axis due to relativistically distorted magnetic field, omega_B = -gyro/gamma * ( (1 - gamma)*(v.B)*v/v^2 + gamma*B - gamma*(v x E)/c^2 )
	double OmegaB[3];
	OmegaB[0] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[0] + gamma_rel*B[0] - gamma_rel*(v[1]*E[2] - v[2]*E[1])/c_0/c_0);
	OmegaB[1] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[1] + gamma_rel*B[1] - gamma_rel*(v[2]*E[0] - v[0]*E[2])/c_0/c_0);
	OmegaB[2] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[2] + gamma_rel*B[2] - gamma_rel*(v[0]*E[1] - v[1]*E[0])/c_0/c_0);

	// Thomas precession in lab frame, omega_T = gamma^2/(gamma + 1)/c^2*(dv/dt x v)
	double OmegaT[3] = {0,0,0};
	OmegaT[0] = gamma_rel*gamma_rel/(gamma_rel + 1)*(a[1]*v[2] - a[2]*v[1])/c_0/c_0;
	OmegaT[1] = gamma_rel*gamma_rel/(gamma_rel + 1)*(a[2]*v[0] - a[0]*v[2])/c_0/c_0;
	OmegaT[2] = gamma_rel*gamma_rel/(gamma_rel + 1)*(a[0]*v[1] - a[1]*v[0])/c_0/c_0;

	// Total spin precession is sum of magnetic-field precession and Thomas precession
	return {OmegaB[0] + OmegaT[0], OmegaB[1] + OmegaT[1], OmegaB[2] + OmegaT[2]};
}


std::array<double, 3> TParticle::SpinEquationOfMotion(const double &t, const std::array<double, 3> &spin, const std::array<double, 3> &omega) const{
/*	array<double, 3> W(omega);
	if (omega.size() == 3){ // if interpolator exists, use it
		W[0] = alglib::spline1dcalc(omega_int[0], t);
		W[1] = alglib::spline1dcalc(omega_int[1], t);
		W[2] = alglib::spline1dcalc(omega_int[2], t);
	}
*/
	return {omega[1]*spin[2] - omega[2]*spin[1], omega[2]*spin[0] - omega[0]*spin[2], omega[0]*spin[1] - omega[1]*spin[0]}; // dS/dt = W x S
}


void TParticle::DoStep(TStep &stepper, const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field){
    int prevpol = stepper.GetPolarization();
    vector<TParticle*> secs;
    OnStep(stepper, currentsolid, mc, ID, secs);
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
    Hmax = max(GetKineticEnergy(stepper.GetVelocity()) + GetPotentialEnergy(stepper.GetTime(), stepper.GetPosition(), stepper.GetVelocity(), stepper.GetPolarization(), field, currentsolid), Hmax);
    if (prevpol != stepper.GetPolarization())
        Nspinflip++;
    Nstep++;
}

void TParticle::DoHit(TStep &stepper, const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc){
    int prevpol = stepper.GetPolarization();
    vector<TParticle*> secs;
    OnHit(stepper, normal, leaving, entering, mc, ID, secs); // do particle specific things
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
    if (prevpol != stepper.GetPolarization())
        Nspinflip++;
    Nhit++;
}

void TParticle::DoDecay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field){
    vector<TParticle*> secs;
    Decay(stepper, mc, geom, field, secs);
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
}

void TParticle::DoPolarize(TStep &stepper, const TFieldManager &field, const bool flipspin, TMCGenerator &mc){
	int prevpol = stepper.GetPolarization();
	double B[3];
	auto pos = stepper.GetPosition();
	field.BField(pos[0], pos[1], pos[2], stepper.GetTime(), B);
	double polarization = stepper.GetSpinPolarization({B[0], B[1], B[2]});
	double flipprob = 0.5*(1 - prevpol*polarization);
	cout << " polarizing: " << polarization << ", " << flipprob << '\n';
	noflipprob *= 1. - flipprob;
	if (flipspin){
		stepper.SetPolarization(polarization_distribution<double>(polarization)(mc));
		if (stepper.GetPolarization() != prevpol)
			++Nspinflip;
	}
}


double TParticle::GetKineticEnergy(const std::array<double, 3> &v) const{
	double beta2 = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/c_0/c_0;
	double gammaminusone;

	// for small velocities the relativistic gamma factor = 1/sqrt(1 - beta^2) is very close to one, givin large round-off errors when calculatin E = m*c^2*(gamma - 1)
	// calculating gamma - 1 has a round-off error of epsilon
	// if a series expansion to order O(beta^8) has a smaller error than epsilon, we will use the series expansion

	//cout << pow(beta2, 5) << " " << (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2 << " " << 1.0/sqrt(1.0 - beta2) - 1.0 << "\n";
	if (pow(beta2, 5) < numeric_limits<double>::epsilon()) // if error in series expansion O(beta^10) is smaller than rounding error in gamma factor
		gammaminusone = (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2; // use series expansion for energy calculation with small beta
	else
		gammaminusone = 1.0/sqrt(1.0 - beta2) - 1.0; // use relativistic formula for larger beta
	return c_0*c_0*m*gammaminusone;
}


double TParticle::GetPotentialEnergy(const double t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const int pol, const TFieldManager &field, const solid &sld) const{
	double result = 0;
	if (q != 0 || mu != 0){
		double B[3], E[3], V;
		if (mu != 0){
			field.BField(pos[0], pos[1], pos[2], t, B);
			result += -pol*mu/ele_e*sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
		}
		if (q != 0){
			field.EField(pos[0], pos[1], pos[2],t,V,E);
			result += q/ele_e*V;
		}
	}
	result += m*gravconst*pos[2];
	return result;
}

void TParticle::SetFinalState(const double& t, const std::array<double, 3>& pos, const std::array<double, 3>& v, const int& polarization, const double& propert, const double& path,
					const std::array<double, 3>& spin, const double& phase, const double& spintime, const solid& sld) {
    tend = t;
    posend = pos;
	vend = v;
	polend = polarization;
	propertime = propert;
	pathlength = path;
    spinend = spin;
	spinphase = phase;
	spinintegrationtime = spintime;
    solidend = sld;
}
