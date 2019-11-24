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
	return GetKineticEnergy(&ystart[3]) + GetPotentialEnergy(tstart, ystart, field, geom.GetSolid(tstart, &ystart[0]));
}


double TParticle::GetFinalTotalEnergy(const TGeometry &geom, const TFieldManager &field) const{
	return GetKineticEnergy(&yend[3]) + GetPotentialEnergy(tend, yend, field, geom.GetSolid(tend, &yend[0]));
}


double TParticle::GetInitialKineticEnergy() const{
	return GetKineticEnergy(&ystart[3]);
}


double TParticle::GetFinalKineticEnergy() const{
	return GetKineticEnergy(&yend[3]);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, const int number,
		const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const double polarisation,
		TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
		: name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN),
		  tstart(t), tend(t), Hmax(0), Nhit(0), Nspinflip(0), noflipprob(1), Nstep(0){

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
	double vstart = c_0*sqrt(beta2);

	if (polarisation < -1 || polarisation > 1)
		throw std::runtime_error("Polarisation has to be between -1 and 1");

	ystart.resize(STATE_VARIABLES, 0);
	ystart[0] = x; // position
	ystart[1] = y;
	ystart[2] = z;
	ystart[3] = vstart*cos(phi)*sin(theta); // velocity
	ystart[4] = vstart*sin(phi)*sin(theta);
	ystart[5] = vstart*cos(theta);
	ystart[6] = 0; // proper time
	std::polarization_distribution<double> pdist(polarisation);
	ystart[7] = pdist(amc); // choose initial polarisation randomly, weighted by spin projection onto magnetic field
	ystart[8] = 0;
	yend = ystart;

	spinstart.resize(SPIN_STATE_VARIABLES, 0);

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

		spinstart[3] = 0; // initial integration time is zero
		spinstart[4] = 0; // initial phase is zero
	}

	spinend = spinstart;

	solidend = solidstart = geometry.GetSolid(t, &ystart[0]); // set to solid with highest priority
	Hmax = GetInitialTotalEnergy(geometry, afield);
}





void TParticle::derivs(const state_type &y, state_type &dydx, const value_type x, const TFieldManager *field) const{
	double B[3], dBidxj[3][3], E[3], V; // magnetic/electric field and electric potential in lab frame
	if (q != 0 || (mu != 0 && y[7] != 0)) // if particle has charge or magnetic moment, calculate magnetic field
		field->BField(y[0],y[1],y[2], x, B, dBidxj);
 	if (q != 0) // if particle has charge caculate electric field
		field->EField(y[0],y[1],y[2], x, V, E);
	EquationOfMotion(y, dydx, x, B, dBidxj, E);
}

void TParticle::EquationOfMotion(const state_type &y, state_type &dydx, const value_type x, const double B[3], const double dBidxj[3][3], const double E[3]) const{
	dydx[0] = y[3]; // time derivatives of position = velocity
	dydx[1] = y[4];
	dydx[2] = y[5];

	value_type F[3] = {0,0,0}; // Force in lab frame
	F[2] += -gravconst*m*ele_e; // add gravitation to force
	if (q != 0){
		F[0] += q*(E[0] + y[4]*B[2] - y[5]*B[1]); // add Lorentz-force
		F[1] += q*(E[1] + y[5]*B[0] - y[3]*B[2]);
		F[2] += q*(E[2] + y[3]*B[1] - y[4]*B[0]);
	}
	if (mu != 0 && y[7] != 0 && (B[0] != 0 || B[1] != 0 || B[2] != 0)){
		double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
		double dBdxi[3] = {	(B[0]*dBidxj[0][0] + B[1]*dBidxj[1][0] + B[2]*dBidxj[2][0])/Babs,
							(B[0]*dBidxj[0][1] + B[1]*dBidxj[1][1] + B[2]*dBidxj[2][1])/Babs,
							(B[0]*dBidxj[0][2] + B[1]*dBidxj[1][2] + B[2]*dBidxj[2][2])/Babs}; // derivatives of |B|
		F[0] += y[7]*mu*dBdxi[0]; // add force on magnetic dipole moment
		F[1] += y[7]*mu*dBdxi[1];
		F[2] += y[7]*mu*dBdxi[2];
	}
	double v2 = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
	value_type inversegamma = sqrt(1 - v2/(c_0*c_0)); // relativstic factor 1/gamma
	dydx[3] = inversegamma/m/ele_e*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
	dydx[4] = inversegamma/m/ele_e*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
	dydx[5] = inversegamma/m/ele_e*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);

	dydx[6] = inversegamma; // derivative of proper time is 1/gamma
	dydx[7] = 0; // polarisaton does not change
	dydx[8] = sqrt(v2); // derivative of path length is abs(velocity)
}





void TParticle::SpinPrecessionAxis(const double t, const dense_stepper_type &stepper, const TFieldManager &field, double &Omegax, double &Omegay, double &Omegaz) const{
	double B[3], dBidxj[3][3], V, E[3];
	state_type y(STATE_VARIABLES), dydt(STATE_VARIABLES);
	stepper.calc_state(t, y); // calculate particle state at time t
	field.BField(y[0], y[1], y[2], t, B, dBidxj);
	field.EField(y[0], y[1], y[2], t, V, E);
	EquationOfMotion(y, dydt, t, B, dBidxj, E); // calculate velocity and acceleration required for vxE effect and Thomas precession
	SpinPrecessionAxis(t, B, E, dydt, Omegax, Omegay, Omegaz); // calculate precession axis
}

void TParticle::SpinPrecessionAxis(const double t, const double B[3], const double E[3], const state_type &dydt, double &Omegax, double &Omegay, double &Omegaz) const{
	double v2 = dydt[0]*dydt[0] + dydt[1]*dydt[1] + dydt[2]*dydt[2];
	double gamma_rel = 1./sqrt(1. - v2/c_0/c_0);
	double Bdotv = B[0]*dydt[0] + B[1]*dydt[1] + B[2]*dydt[2];
	double Bparallel[3];
	for (int j = 0; j < 3; j++)
		Bparallel[j] = Bdotv*dydt[j]/v2; // magnetic-field component parallel to velocity

	// spin precession axis due to relativistically distorted magnetic field, omega_B = -gyro/gamma * ( (1 - gamma)*(v.B)*v/v^2 + gamma*B - gamma*(v x E)/c^2 )
	double OmegaB[3];
	OmegaB[0] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[0] + gamma_rel*B[0] - gamma_rel*(dydt[1]*E[2] - dydt[2]*E[1])/c_0/c_0);
	OmegaB[1] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[1] + gamma_rel*B[1] - gamma_rel*(dydt[2]*E[0] - dydt[0]*E[2])/c_0/c_0);
	OmegaB[2] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[2] + gamma_rel*B[2] - gamma_rel*(dydt[0]*E[1] - dydt[1]*E[0])/c_0/c_0);

	// Thomas precession in lab frame, omega_T = gamma^2/(gamma + 1)/c^2*(dv/dt x v)
	double OmegaT[3] = {0,0,0};
	OmegaT[0] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[4]*dydt[2] - dydt[5]*dydt[1])/c_0/c_0;
	OmegaT[1] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[5]*dydt[0] - dydt[3]*dydt[2])/c_0/c_0;
	OmegaT[2] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[3]*dydt[1] - dydt[4]*dydt[0])/c_0/c_0;

	// Total spin precession is sum of magnetic-field precession and Thomas precession
	Omegax = OmegaB[0] + OmegaT[0];
	Omegay = OmegaB[1] + OmegaT[1];
	Omegaz = OmegaB[2] + OmegaT[2];
}


void TParticle::SpinDerivs(const state_type &y, state_type &dydx, const value_type x, const dense_stepper_type &stepper, const TFieldManager *field, const std::vector<alglib::spline1dinterpolant> &omega) const{
	double omegax, omegay, omegaz;
	if (omega.size() == 3){ // if interpolator exists, use it
		omegax = alglib::spline1dcalc(omega[0], x);
		omegay = alglib::spline1dcalc(omega[1], x);
		omegaz = alglib::spline1dcalc(omega[2], x);
	}
	else
		SpinPrecessionAxis(x, stepper, *field, omegax, omegay, omegaz); // else calculate precession axis directly

	dydx[0] = omegay*y[2] - omegaz*y[1]; // dS/dt = W x S
	dydx[1] = omegaz*y[0] - omegax*y[2];
	dydx[2] = omegax*y[1] - omegay*y[0];
	dydx[3] = 1.; // integrate time
	dydx[4] = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz); // integrate precession phase
}


void TParticle::DoStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
                       const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field){
    state_type y2temp = y2;
    vector<TParticle*> secs;
    OnStep(x1, y1, x2, y2, stepper, currentsolid, mc, ID, secs);
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
    Hmax = max(GetKineticEnergy(&y2[3]) + GetPotentialEnergy(x2, y2, field, currentsolid), Hmax);
    if (y2temp[7] != y2[7])
        Nspinflip++;
    Nstep++;
}

void TParticle::DoHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
                      const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc){
    state_type y2temp = y2;
    vector<TParticle*> secs;
    OnHit(x1, y1, x2, y2, normal, leaving, entering, mc, ID, secs); // do particle specific things
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
    if (y2temp[7] != y2[7])
        Nspinflip++;
    Nhit++;
}

void TParticle::DoDecay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field){
    vector<TParticle*> secs;
    Decay(t, y, mc, geom, field, secs);
    for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
}

void TParticle::DoPolarize(const double t, state_type &y, const double polarization, TMCGenerator &mc){
    double prevpol = y[7];
    polarization_distribution<double> pdist(polarization);
    y[7] = pdist(mc);
    if (y[7] != prevpol)
        ++Nspinflip;
}


void TParticle::Print(const value_type x, const state_type &y, const state_type &spin, const TGeometry &geom, const TFieldManager &field, const LogStream logType) const{
	ofstream &file = GetLogStream(logType);
	if (!file.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << name;
		if (logType == endLog)
			filename << "end.out";
		else if (logType == snapshotLog)
			filename << "snapshot.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		file.open(outfile.c_str());
		if (!file.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		file <<	"jobnumber particle "
					"tstart xstart ystart zstart "
					"vxstart vystart vzstart polstart "
					"Sxstart Systart Szstart "
					"Hstart Estart Bstart Ustart solidstart "
					"tend xend yend zend "
					"vxend vyend vzend polend "
					"Sxend Syend Szend "
					"Hend Eend Bend Uend solidend "
					"stopID Nspinflip spinflipprob "
					"Nhit Nstep trajlength Hmax wL\n";
		file << std::setprecision(std::numeric_limits<double>::digits10); // need maximum precision for wL and delwL 
	}
//	cout << "Printing status\n";

	value_type E = GetKineticEnergy(&y[3]);
	double B[3], Ei[3], V;

	field.BField(ystart[0], ystart[1], ystart[2], tstart, B);
	field.EField(ystart[0], ystart[1], ystart[2], tstart, V, Ei);

	file	<< jobnumber << " " << particlenumber << " "
			<< tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
			<< ystart[3] << " " << ystart[4] << " " << ystart[5] << " " << ystart[7] << " "
			<< spinstart[0] << " " << spinstart[1] << " " << spinstart[2] << " " << GetInitialTotalEnergy(geom, field) << " " << GetInitialKineticEnergy() << " "
			<< sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) << " " << V << " " << solidstart.ID << " ";

	double H;
	solid sld = geom.GetSolid(x, &y[0]);
	H = E + GetPotentialEnergy(x, y, field, sld);

	field.BField(y[0], y[1], y[2], x, B);
	field.EField(y[0], y[1], y[2], x, V, Ei);

	double wL = 0;
	if (spin[3] > 0)
		wL = spin[4]/spin[3];

	file	<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< y[3] << " " << y[4] << " " << y[5] << " " << y[7] << " "
			<< spin[0] << " " << spin[1] << " " << spin[2] << " " << H << " " << E << " "
			<< sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) << " " << V << " " << sld.ID << " "
			<< ID << " " << Nspinflip << " " << 1 - noflipprob << " "
			<< Nhit << " " << Nstep << " " << y[8] << " " << Hmax << " " << wL << '\n';
}


void TParticle::PrintTrack(const value_type x, const state_type &y, const state_type &spin, const solid &sld, const TFieldManager &field) const{
	ofstream &trackfile = GetLogStream(trackLog);
	if (!trackfile.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << name << "track.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		trackfile.open(outfile.c_str());
		if (!trackfile.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		trackfile << 	"jobnumber particle polarisation "
						"t x y z vx vy vz "
						"H E Bx dBxdx dBxdy dBxdz By dBydx "
						"dBydy dBydz Bz dBzdx dBzdy dBzdz Ex Ey Ez V\n";
		trackfile.precision(10);
	}

//	cout << "-";
	double B[3] = {0,0,0};
	double dBidxj[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double E[3] = {0,0,0};
	double V = 0;
	field.BField(y[0],y[1],y[2],x,B, dBidxj);
	field.EField(y[0],y[1],y[2],x,V,E);
	value_type Ek = GetKineticEnergy(&y[3]);
	value_type H = Ek + GetPotentialEnergy(x, y, field, sld);

	trackfile << jobnumber << " " << particlenumber << " " << y[7] << " "
				<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
				<< H << " " << Ek << " ";
	for (int i = 0; i < 3; i++){
		trackfile << B[i] << " ";
		for (int j = 0; j < 3; j++)
			trackfile << dBidxj[i][j] << " ";
	}
	trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
}


void TParticle::PrintHit(const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering) const{
	ofstream &hitfile = GetLogStream(hitLog);
	if (!hitfile.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << name << "hit.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		hitfile.open(outfile.c_str());
		if (!hitfile.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		hitfile << "jobnumber particle "
					"t x y z v1x v1y v1z pol1 "
					"v2x v2y v2z pol2 "
					"nx ny nz solid1 solid2\n";
		hitfile.precision(10);
	}

//	cout << ":";
	hitfile << jobnumber << " " << particlenumber << " "
			<< x << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << y1[7] << " "
			<< y2[3] << " " << y2[4] << " " << y2[5] << " " << y2[7] << " "
			<< normal[0] << " " << normal[1] << " " << normal[2] << " " << leaving.ID << " " << entering.ID << '\n';
}


// void TParticle::PrintSpin(const value_type x, const state_type &spin, const dense_stepper_type &stepper, const TFieldManager &field) const{
void TParticle::PrintSpin(const value_type x, const state_type &y, const state_type &spin, const dense_stepper_type &stepper, const TFieldManager &field) const{
	ofstream &spinfile = GetLogStream(spinLog);
	double B[3] = {0,0,0};
	field.BField(y[0],y[1],y[2],x,B);
	if (!spinfile.is_open()){
		std::ostringstream filename;
		filename << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << name << "spin.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		std::cout << "Creating " << outfile << '\n';
		spinfile.open(outfile.c_str());
		if(!spinfile.is_open())
		{
			throw std::runtime_error("Could not open " + outfile.native());
		}

		//need the maximum accuracy in spinoutlog for the larmor frequency to see any difference
		spinfile << std::setprecision(std::numeric_limits<double>::digits10);
		// spinfile << "jobnumber particle t Sx Sy Sz Wx Wy Wz\n";
		spinfile << "jobnumber particle t x y z Sx Sy Sz Wx Wy Wz Bx By Bz\n";
	}
//	std::cout << "/";
	double Omega[3];
	SpinPrecessionAxis(x, stepper, field, Omega[0], Omega[1], Omega[2]);

	// spinfile << jobnumber << " " << particlenumber << " "
	// 		<< x << " " << spin[0] << " " << spin[1] << " " << spin[2] << " "
	// 		<< Omega[0] << " " << Omega[1] << " " << Omega[2] << "\n";
	spinfile << jobnumber << " " << particlenumber << " "
			<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< spin[0] << " " << spin[1] << " " << spin[2] << " "
			<< Omega[0] << " " << Omega[1] << " " << Omega[2] << " "
			<< B[0] << " " << B[1] << " " << B[2] << "\n";
}


double TParticle::GetKineticEnergy(const value_type v[3]) const{
	double beta2 = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/c_0/c_0;
	double gammaminusone;

	// for small velocities the relativistic gamma factor = 1/sqrt(1 - beta^2) is very close to one, givin large round-off errors when calculatin E = m*c^2*(gamma - 1)
	// calculating gamma - 1 has a round-off error of epsilon
	// if a series expansion to order O(beta^8) has a smaller error than epsilon, we will use the series expansion

	//cout << pow(beta2, 5) << " " << (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2 << " " << 1.0/sqrt(1.0 - beta2) - 1.0 << "\n";
	if (pow(beta2, 5) < numeric_limits<value_type>::epsilon()) // if error in series expansion O(beta^10) is smaller than rounding error in gamma factor
		gammaminusone = (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2; // use series expansion for energy calculation with small beta
	else
		gammaminusone = 1.0/sqrt(1.0 - beta2) - 1.0; // use relativistic formula for larger beta
	return c_0*c_0*m*gammaminusone;
}


double TParticle::GetPotentialEnergy(const value_type t, const state_type &y, const TFieldManager &field, const solid &sld) const{
	value_type result = 0;
	if (q != 0 || mu != 0){
		double B[3], E[3], V;
		if (mu != 0){
			field.BField(y[0],y[1],y[2],t,B);
			result += -y[7]*mu/ele_e*sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
		}
		if (q != 0){
			field.EField(y[0],y[1],y[2],t,V,E);
			result += q/ele_e*V;
		}
	}
	result += m*gravconst*y[2];
	return result;
}

void TParticle::SetFinalState(const value_type& x, const state_type& y, const state_type& spin,
        const double anoflipprob, const solid& sld) {
    tend = x;
    yend = y;
    spinend = spin;
    noflipprob *= anoflipprob;
    solidend = sld;
}
