/**
 * \file
 * Particle base class definition.
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <boost/numeric/odeint.hpp>

#include "particle.h"

using namespace std;

double TParticle::GetInitialTotalEnergy() const{
	return GetKineticEnergy(&ystart[3]) + GetPotentialEnergy(tstart, ystart, field, geom->GetSolid(tstart, &ystart[0]));
}


double TParticle::GetFinalTotalEnergy() const{
	return GetKineticEnergy(&yend[3]) + GetPotentialEnergy(tend, yend, field, geom->GetSolid(tend, &yend[0]));
}


double TParticle::GetInitialKineticEnergy() const{
	return GetKineticEnergy(&ystart[3]);
}


double TParticle::GetFinalKineticEnergy() const{
	return GetKineticEnergy(&yend[3]);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, int number,
		double t, double x, double y, double z, double E, double phi, double theta, double polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN),
		  tstart(t), tend(t), Hmax(0), lend(0), Nhit(0), Nspinflip(0), noflipprob(1), Nstep(0),
		  geom(&geometry), mc(&amc), field(afield){
	value_type Eoverm = E/m/c_0/c_0; // gamma - 1
	value_type beta = sqrt(1-(1/((Eoverm + 1)*(Eoverm + 1)))); // beta = sqrt(1 - 1/gamma^2)
	value_type vstart;
//			cout << "(E/m)^3.5 = " << pow(Eoverm, 3.5L) << "; err_beta = " << numeric_limits<value_tyoe>::epsilon()/(1 - beta) << '\n';
	if (pow(Eoverm, 3.5L) < numeric_limits<value_type>::epsilon()/(1 - beta)) // if error in series O((E/m)^3.5) is smaller than rounding error in beta
		vstart = c_0*(sqrt(2*Eoverm) - 3*pow(Eoverm, 1.5L)/2/sqrt(2) + 23*pow(Eoverm, 2.5L)/16/sqrt(2)); // use series expansion
	else
		vstart = c_0 * beta;

	if (polarisation < -1 || polarisation > 1){
		std::cout << "Polarisation is " << polarisation << "! Has to be between -1 and 1!\n";
		exit(-1);
	}

	ystart.resize(STATE_VARIABLES, 0);
	ystart[0] = x; // position
	ystart[1] = y;
	ystart[2] = z;
	ystart[3] = vstart*cos(phi)*sin(theta); // velocity
	ystart[4] = vstart*sin(phi)*sin(theta);
	ystart[5] = vstart*cos(theta);
	ystart[6] = 0; // proper time
	ystart[7] = amc.DicePolarisation(polarisation); // choose initial polarisation randomly, weighted by spin projection onto magnetic field
	yend = ystart;

	spinstart.resize(SPIN_STATE_VARIABLES, 0);
	if (field){
		double B[4][4];
		field->BField(x, y, z, t, B);
		if (B[3][0] > 0){
			double spinaz = mc->UniformDist(0, 2*pi); // random azimuth angle of spin
			spinstart[0] = sqrt(1 - polarisation*polarisation)*sin(spinaz); // set initial spin vector
			spinstart[1] = sqrt(1 - polarisation*polarisation)*cos(spinaz);
			spinstart[2] = polarisation;
			double B0[3] = {B[0][0], B[1][0], B[2][0]};
			RotateVector(&spinstart[0], B0); // rotate initial spin vector such that its z-component is parallel to magnetic field

			spinstart[3] = 0; // initial integration time is zero
			spinstart[4] = 0; // initial phase is zero
		}
	}
	spinend = spinstart;

	Hmax = GetInitialTotalEnergy();
	maxtraj = mc->MaxTrajLength(name); // max. trajectory length
	tau = mc->LifeTime(name); // lifetime either exponentially distributed or fixed
	geom->GetSolids(t, &ystart[0], currentsolids); // detect solids that surround the particle
	solidend = solidstart = GetCurrentsolid(); // set to solid with highest priority
}


TParticle::~TParticle(){
	for (vector<TParticle*>::reverse_iterator i = secondaries.rbegin(); i != secondaries.rend(); i++)
		delete *i;
}


void TParticle::Integrate(double tmax, std::map<std::string, std::string> &conf){
	if (currentsolids.empty())
		geom->GetSolids(tend, &yend[0], currentsolids);

	int perc = 0;
	cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
	cout << "x: " << yend[0] << "m y: " << yend[1] << "m z: " << yend[2]
		 << "m E: " << GetFinalKineticEnergy() << "eV t: " << tend << "s tau: " << tau << "s lmax: " << maxtraj << "m\n";

	// set initial values for integrator
	value_type x = tend;
	state_type y = yend; 

	bool resetintegration = true;

	float nextsnapshot = -1;
	bool snapshotlog = false;
	istringstream(conf["snapshotlog"]) >> snapshotlog;
	istringstream snapshots(conf["snapshots"]);
	if (snapshotlog){
		do{
			snapshots >> nextsnapshot;
		}while (snapshots.good() && nextsnapshot < x); // find first snapshot time
	}

	bool tracklog = false;
	istringstream(conf["tracklog"]) >> tracklog;
	if (tracklog)
		PrintTrack(tend, yend, spinend, solidend);
	double trackloginterval = 1e-3;
	istringstream(conf["trackloginterval"]) >> trackloginterval;
	value_type lastsave = x;

	bool hitlog = false;
	istringstream(conf["hitlog"]) >> hitlog;

	bool flipspin = false;
	istringstream(conf["flipspin"]) >> flipspin;
	
	bool spininterpolatefields = false;
	istringstream(conf["interpolatefields"]) >> spininterpolatefields;

	int spinlog = false;
	double SpinBmax = 0, spinloginterval = 0, nextspinlog = std::numeric_limits<double>::infinity();
	vector<double> SpinTimes;
	std::istringstream(conf["Bmax"]) >> SpinBmax;
	std::istringstream SpinTimess(conf["spintimes"]);
	do{
		double t;
		SpinTimess >> t;
		if (SpinTimess.good())
			SpinTimes.push_back(t);
	}while(SpinTimess.good());
	std::istringstream(conf["spinlog"]) >> spinlog;
	std::istringstream(conf["spinloginterval"]) >> spinloginterval;
	if (spinlog)
		nextspinlog = 0;
	state_type spin = spinend;

	dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type());

	while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
		if (resetintegration){
			stepper.initialize(y, x, MAX_SAMPLE_DIST/sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])); // (re-)start integration with first guess of step size
		}
		value_type x1 = x; // save point before next step
		state_type y1 = y;

		try{
			stepper.do_step(std::bind(&TParticle::derivs, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
			x = stepper.current_time();
			y = stepper.current_state();
			Nstep++;
		}
		catch(...){ // catch Exceptions thrown by numerical recipes routines
			ID = ID_ODEINT_ERROR;
		}

		if (y[6] > tau){ // if proper time is larger than lifetime
			x = x1 + (x - x1)*(tau - y1[6])/(y[6] - y1[6]); // interpolate decay time in lab frame
			stepper.calc_state(x, y);
		}
		if (x > tmax){	//If stepsize overshot max simulation time
			x = tmax;
			stepper.calc_state(x, y);
		}

		while (x1 < x){ // split integration step in pieces (x1,y1->x2,y2) with spatial length SAMPLE_DIST, go through all pieces
			value_type v1 = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			value_type x2 = x1 + MAX_SAMPLE_DIST/v1; // time length = spatial length/velocity
			state_type y2(STATE_VARIABLES);
			if (x2 >= x){
				x2 = x;
				y2 = y;
			}
			else{
				stepper.calc_state(x2, y2);
			}

			resetintegration = CheckHit(x1, y1, x2, y2, stepper, hitlog); // check if particle hit a material boundary or was absorbed between y1 and y2
			if (resetintegration){
				x = x2; // if particle path was changed: reset integration end point
				y = y2;
			}

			lend += sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));

			x1 = x2;
			y1 = y2;
		}

		// take snapshots at certain times
		if (snapshotlog && snapshots.good()){
			if (stepper.previous_time() <= nextsnapshot && x > nextsnapshot){
				state_type ysnap(STATE_VARIABLES);
				stepper.calc_state(nextsnapshot, ysnap);
				cout << "\n Snapshot at " << nextsnapshot << " s \n";

				Print(nextsnapshot, ysnap, spin, geom->GetSolid(nextsnapshot, &ysnap[0], currentsolids), snapshotLog);
				snapshots >> nextsnapshot;
			}
		}

		double prevpol = y[7];
		noflipprob *= 1 - IntegrateSpin(spin, stepper, x, y, SpinTimes, spininterpolatefields, SpinBmax, flipspin, spinloginterval, nextspinlog); // calculate spin precession and spin-flip probability
		if (y[7] != prevpol)
			Nspinflip++;

		Hmax = max(GetKineticEnergy(&y[3]) + GetPotentialEnergy(x, y, field, GetCurrentsolid()), Hmax);
		if (tracklog && x - lastsave > trackloginterval/sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])){
			PrintTrack(x, y, spin, GetCurrentsolid());
			lastsave = x;
		}

		PrintPercent(max(y[6]/tau, max((x - tstart)/(tmax - tstart), lend/maxtraj)), perc);
		
		if (ID == ID_UNKNOWN && y[6] >= tau) // proper time >= tau?
			ID = ID_DECAYED;
		else if (ID == ID_UNKNOWN && (x >= tmax || lend >= maxtraj)) // time > tmax or trajectory length > max length?
			ID = ID_NOT_FINISH;
	}

	tend = x;
	yend = y;
	spinend = spin;
	solidend = GetCurrentsolid();
	Print(tend, yend, spinend, solidend, endLog);

	if (ID == ID_DECAYED){ // if particle reached its lifetime call TParticle::Decay
		cout << "Decayed!\n";
		Decay(secondaries);
	}

	cout << "x: " << yend[0];
	cout << " y: " << yend[1];
	cout << " z: " << yend[2];
	cout << " E: " << GetFinalKineticEnergy();
	cout << " Code: " << ID;
	cout << " t: " << tend;
	cout << " l: " << lend;
	cout << " hits: " << Nhit;
	cout << " spinflips: " << Nspinflip << '\n';
	cout << "Computation took " << Nstep << " steps\n";
	cout << "Done!!\n\n";
//			cout.flush();
}


void TParticle::derivs(const state_type &y, state_type &dydx, const value_type x) const{
	double B[4][4], E[3], V; // magnetic/electric field and electric potential in lab frame
	if (field){
		if (q != 0 || (mu != 0 && y[7] != 0)) // if particle has charge or magnetic moment, calculate magnetic field
			field->BField(y[0],y[1],y[2], x, B);
		if (q != 0) // if particle has charge caculate electric field
			field->EField(y[0],y[1],y[2], x, V, E);

	}
	EquationOfMotion(y, dydx, x, B, E);
}

void TParticle::EquationOfMotion(const state_type &y, state_type &dydx, const value_type x, const double B[4][4], const double E[3]) const{
	dydx[0] = y[3]; // time derivatives of position = velocity
	dydx[1] = y[4];
	dydx[2] = y[5];

	value_type F[3] = {0,0,0}; // Force in lab frame
	F[2] += -gravconst*m*ele_e; // add gravitation to force
	if (field){
		if (q != 0){
			F[0] += q*(E[0] + y[4]*B[2][0] - y[5]*B[1][0]); // add Lorentz-force
			F[1] += q*(E[1] + y[5]*B[0][0] - y[3]*B[2][0]);
			F[2] += q*(E[2] + y[3]*B[1][0] - y[4]*B[0][0]);
		}
		if (mu != 0 && y[7] != 0){
			F[0] += y[7]*mu*B[3][1]; // add force on magnetic dipole moment
			F[1] += y[7]*mu*B[3][2];
			F[2] += y[7]*mu*B[3][3];
		}
	}
	value_type inversegamma = sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0)); // relativstic factor 1/gamma
	dydx[3] = inversegamma/m/ele_e*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
	dydx[4] = inversegamma/m/ele_e*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
	dydx[5] = inversegamma/m/ele_e*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);

	dydx[6] = inversegamma; // derivative of proper time is 1/gamma
	dydx[7] = 0; // polarisaton does not change
}


double TParticle::IntegrateSpin(state_type &spin, const dense_stepper_type &stepper, const double x2, state_type &y2, const std::vector<double> &times,
								const bool interpolatefields, const double Bmax, const bool flipspin, const double spinloginterval, double &nextspinlog) const{
	value_type x1 = stepper.previous_time();
	if (gamma == 0 || x1 == x2 || !field)
		return 0;

	state_type y1 = stepper.previous_state();
	double B1[4][4], B2[4][4], polarisation;
	field->BField(y1[0], y1[1], y1[2], x1, B1);
	field->BField(y2[0], y2[1], y2[2], x2, B2);

	if (B2[3][0] == 0){ // if there's no magnetic field or particle is leaving magnetic field, do nothing
		std::fill(spin.begin(), spin.end(), 0); // set all spin components to zero
		return 0;
	}
	else if (B1[3][0] == 0 && B2[3][0] > 0){ // if particle enters magnetic field
		std::cout << "Entering magnetic field\n";
		if (!flipspin) // if polarization should stay constant
			polarisation = y1[7]; // take polarization from state vector
		else{
			if (mc->UniformDist(0,1) < 0.5) // if spin flips are allowed, choose random polarization
				polarisation = 1;
			else
				polarisation = -1;
			y2[7] = polarisation;
		}
		spin[0] = polarisation*B2[0][0]/B2[3][0]; // set spin (anti-)parallel to field
		spin[1] = polarisation*B2[1][0]/B2[3][0];
		spin[2] = polarisation*B2[2][0]/B2[3][0];
		return 0.5;
	}
	else{ // if particle already is in magnetic field, calculate spin-projection on magnetic field
		polarisation = (spin[0]*B1[0][0] + spin[1]*B1[1][0] + spin[2]*B1[2][0])/B1[3][0]/sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
	}

	bool integrate1 = false, integrate2 = false;;
	for (unsigned int i = 0; i < times.size(); i += 2){
		integrate1 |= (x1 >= times[i] && x1 < times[i+1]);
		integrate2 |= (x2 >= times[i] && x2 < times[i+1]);
	}

	if ((integrate1 || integrate2) && (B1[3][0] < Bmax || B2[3][0] < Bmax)){ // do spin integration only, if time is in specified range and field is smaller than Bmax
		if ((!integrate1 && integrate2) || (B1[3][0] > Bmax && B2[3][0] < Bmax))
			std::cout << x1 << "s " << y1[7] - polarisation << " ";

		vector<alglib::spline1dinterpolant> omega_int;
		if (interpolatefields){
			alglib::real_1d_array ts; // set up values for interpolation of precession axis
			vector<alglib::real_1d_array> omega(3);
			omega_int.resize(3);
			const int int_points = 10;
			ts.setlength(int_points + 1);
			for (int i = 0; i < 3; i++)
				omega[i].setlength(int_points + 1);

			for (int i = 0; i <= int_points; i++){ // calculate precession axis at several points along trajectory step
				double t = x1 + i*(x2 - x1)/int_points;
				ts[i] = t;
				SpinPrecessionAxis(t, stepper, omega[0][i], omega[1][i], omega[2][i]);
			}

			for (int i = 0; i < 3; i++)
				alglib::spline1dbuildcubic(ts, omega[i], omega_int[i]); // interpolate all three components of precession axis
		}

		if (x1 >= nextspinlog){
			PrintSpin(x1, spin, stepper);
			nextspinlog += spinloginterval;
		}

		dense_stepper_type spinstepper = boost::numeric::odeint::make_dense_output(1e-12, 1e-12, stepper_type());
		spinstepper.initialize(spin, x1, std::abs(pi/gamma/B1[3][0])); // initialize integrator with step size = half rotation
		unsigned int steps = 0;
		while (true){
			// take an integration step, SpinDerivs contains right-hand side of equation of motion
			spinstepper.do_step(std::bind(&TParticle::SpinDerivs, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, stepper, omega_int));
			steps++;
			double t = spinstepper.current_time();
			if (t > x2){ // if stepper overshot, calculate end point and stop
				t = x2;
				spinstepper.calc_state(t, spin);
			}
			else
				spin = spinstepper.current_state();

			if (t >= nextspinlog){
				PrintSpin(t, spin, stepper);
				nextspinlog += spinloginterval;
			}
			if (t >= x2)
				break;
		}

		// calculate new spin projection
		polarisation = (spin[0]*B2[0][0] + spin[1]*B2[1][0] + spin[2]*B2[2][0])/B2[3][0]/sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
	}
	else if ((B1[3][0] < Bmax || B2[3][0] < Bmax)){ // if time outside selected ranges, parallel-transport spin along magnetic field
		if (polarisation*polarisation >= 1){ // catch rounding errors
			spin[0] = 0;
			spin[1] = 0;
		}
		else{
			double spinaz = mc->UniformDist(0, 2*pi); // random azimuth
			spin[0] = sqrt(1 - polarisation*polarisation)*sin(spinaz);
			spin[1] = sqrt(1 - polarisation*polarisation)*cos(spinaz);
		}
		spin[2] = polarisation;
		double B0[3] = {B2[0][0], B2[1][0], B2[2][0]};
		RotateVector(&spin[0], B0); // rotate spin vector onto magnetic field vector
	}


	if (B2[3][0] > Bmax){ // if magnetic field grows above Bmax, collapse spin state to one of the two polarisation states
		double flipprob = 0.5*(1 - y2[7]*polarisation);
		if ((integrate1 || integrate2) && B1[3][0] < Bmax)
			std::cout << x2 << "s flipprob " << flipprob << "\n";
		if (flipspin)
			y2[7] = mc->DicePolarisation(polarisation); // set new polarisation weighted by spin-projection on magnetic field
		spin[0] = B2[0][0]*y2[7]/B2[3][0];
		spin[1] = B2[1][0]*y2[7]/B2[3][0];
		spin[2] = B2[2][0]*y2[7]/B2[3][0];
		return flipprob;
	}
	return 0;
}


void TParticle::SpinPrecessionAxis(const double t, const dense_stepper_type &stepper, double &Omegax, double &Omegay, double &Omegaz) const{
	double B[4][4], V, E[3];
	state_type y(STATE_VARIABLES), dydt(STATE_VARIABLES);
	stepper.calc_state(t, y); // calculate particle state at time t
	field->BField(y[0], y[1], y[2], t, B);
	field->EField(y[0], y[1], y[2], t, V, E);
	EquationOfMotion(y, dydt, t, B, E); // calculate velocity and acceleration required for vxE effect and Thomas precession
	SpinPrecessionAxis(t, B, E, dydt, Omegax, Omegay, Omegaz); // calculate precession axis
}

void TParticle::SpinPrecessionAxis(const double t, const double B[4][4], const double E[3], const state_type &dydt, double &Omegax, double &Omegay, double &Omegaz) const{
	double v2 = dydt[0]*dydt[0] + dydt[1]*dydt[1] + dydt[2]*dydt[2];
	double gamma_rel = 1./sqrt(1. - v2/c_0/c_0);
	double Bparallel[3];
	for (int j = 0; j < 3; j++)
		Bparallel[j] = (B[0][0]*dydt[0] + B[1][0]*dydt[1] + B[2][0]*dydt[2])*dydt[j]/v2; // magnetic-field component parallel to velocity

	// spin precession axis due to relativistically distorted magnetic field, omega_B = -gyro/gamma * ( (1 - gamma)*(v.B)*v/v^2 + gamma*B - gamma*(v x E)/c^2 )
	double OmegaB[3];
	OmegaB[0] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[0] + gamma_rel*B[0][0] - gamma_rel*(dydt[1]*E[2] - dydt[2]*E[1])/c_0/c_0);
	OmegaB[1] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[1] + gamma_rel*B[1][0] - gamma_rel*(dydt[2]*E[0] - dydt[0]*E[2])/c_0/c_0);
	OmegaB[2] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[2] + gamma_rel*B[2][0] - gamma_rel*(dydt[0]*E[1] - dydt[1]*E[0])/c_0/c_0);

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


void TParticle::SpinDerivs(const state_type &y, state_type &dydx, const value_type x, const dense_stepper_type &stepper, const std::vector<alglib::spline1dinterpolant> &omega) const{
	double omegax, omegay, omegaz;
	if (omega.size() == 3){ // if interpolator exists, use it
		omegax = alglib::spline1dcalc(omega[0], x);
		omegay = alglib::spline1dcalc(omega[1], x);
		omegaz = alglib::spline1dcalc(omega[2], x);
	}
	else
		SpinPrecessionAxis(x, stepper, omegax, omegay, omegaz); // else calculate precession axis directly

	dydx[0] = omegay*y[2] - omegaz*y[1]; // dS/dt = W x S
	dydx[1] = omegaz*y[0] - omegax*y[2];
	dydx[2] = omegax*y[1] - omegay*y[0];
	dydx[3] = 1.; // integrate time
	dydx[4] = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz); // integrate precession phase
}


const solid& TParticle::GetCurrentsolid() const{
	map<solid, bool>::const_iterator it = currentsolids.begin();
	while (it->second) // skip over ignored solids
		it++;
	return it->first; // return first non-ignored solid in list
}


bool TParticle::CheckHitError(const solid &hitsolid, double distnormal) const{
	if (distnormal < 0){ // particle is entering solid
		if (currentsolids.find(hitsolid) != currentsolids.end()){ // if solid already is in currentsolids list something went wrong
			cout << "Particle entering " << hitsolid.name << " which it did enter before! Stopping it! Did you define overlapping solids with equal priorities?\n";
			return true;
		}
	}
	else if (distnormal > 0){ // particle is leaving solid
		if (currentsolids.find(hitsolid) == currentsolids.end()){ // if solid is not in currentsolids list something went wrong
			cout << "Particle inside '" << hitsolid.name << "' which it did not enter before! Stopping it!\n";
			return true;
		}
	}
	else{
		cout << "Particle is crossing material boundary with parallel track! Stopping it!\n";
		return true;
	}
	return false;
}


bool TParticle::CheckHit(const value_type x1, const state_type y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper, const bool hitlog, const int iteration){
	solid currentsolid = GetCurrentsolid();
	if (!geom->CheckSegment(&y1[0], &y2[0])){ // check if start point is inside bounding box of the simulation geometry
		printf("\nParticle has hit outer boundaries: Stopping it! t=%g x=%g y=%g z=%g\n",x2,y2[0],y2[1],y2[2]);
		ID = ID_HIT_BOUNDARIES;
		return true;
	}

	map<TCollision, bool> colls;
	bool collfound = false;
	try{
		collfound = geom->GetCollisions(x1, &y1[0], x2, &y2[0], colls);
	}
	catch(...){
		ID = ID_CGAL_ERROR;
		return true;
	}

	if (collfound){	// if there is a collision with a wall
		TCollision coll = colls.begin()->first;
		if ((abs(coll.s*coll.distnormal) < REFLECT_TOLERANCE
			&& (1 - coll.s)*abs(coll.distnormal) < REFLECT_TOLERANCE) // if first collision is closer to y1 and y2 than REFLECT_TOLERANCE
			|| iteration > 99) // or if iteration counter reached a certain maximum value
		{
			bool trajectoryaltered = false, traversed = true;

			map<solid, bool> newsolids = currentsolids;
			for (map<TCollision, bool>::iterator it = colls.begin(); it != colls.end(); it++){ // go through list of collisions
				solid sld = geom->solids[it->first.sldindex];
				if (CheckHitError(sld, it->first.distnormal)){ // check all hits for errors
					ID = ID_GEOMETRY_ERROR; // stop trajectory integration if error found
					return true;
				}

				if (it->first.distnormal < 0){ // if particle is entering hit surface
//					cout << "-->" << sld.name << ' ';
					newsolids[sld] = it->second; // add new solid to list
				}
				else{ // if particle is leaving hit surface
//					cout << sld.name << "--> ";
					newsolids.erase(sld); // remove solid from list of new solids
				}

				if (sld.ID > geom->solids[coll.sldindex].ID)
					coll = it->first; // use geometry from collision with highest priority
}
//			cout << '\n';

			solid leaving = currentsolid; // particle can only leave highest-priority solid
			solid entering = geom->defaultsolid;
			for (map<solid, bool>::iterator it = newsolids.begin(); it != newsolids.end(); it++){ // go through list of new solids
				if (!it->second && it->first.ID > entering.ID){ // highest-priority solid in newsolids list will be entered in collisions
					entering = it->first;
				}
			}
//			cout << "Leaving " << leaving.name << ", entering " << entering.name << ", currently " << currentsolid.name << '\n';
			if (leaving.ID != entering.ID){ // if the particle actually traversed a material interface
				trajectoryaltered = OnHit(x1, y1, x2, y2, coll.normal, leaving, entering, traversed, ID, secondaries); // do particle specific things
				if (hitlog)
					PrintHit(x1, y1, y2, coll.normal, leaving, entering); // print collision to file if requested
				Nhit++;
			}

			if (traversed){
				currentsolids = newsolids; // if surface was traversed (even if it was  physically ignored) replace current solids with list of new solids
			}

			if (OnStep(x1, y1, x2, y2, stepper, GetCurrentsolid(), ID, secondaries)){ // check for absorption/scattering
				if (y2[7] != y1[7])
					Nspinflip++;
				return true;
			}
			if (trajectoryaltered)
				return true;

		}
		else{
			// else cut integration step right before and after first collision point
			// and call CheckHit again for each smaller step (quite similar to bisection algorithm)
			value_type xnew, xbisect1 = x1, xbisect2 = x1;
			state_type ybisect1(STATE_VARIABLES);
			state_type ybisect2(STATE_VARIABLES);
			ybisect1 = ybisect2 = y1;

			xnew = x1 + (x2 - x1)*(coll.s - 0.01*iteration); // cut integration right before collision point
			if (xnew > x1 && xnew < x2){ // check that new line segment is in correct time interval
				xbisect1 = xbisect2 = xnew;
				stepper.calc_state(xbisect1, ybisect1);
				ybisect2 = ybisect1;
				if (CheckHit(x1, y1, xbisect1, ybisect1, stepper, hitlog, iteration+1)){ // recursive call for step before coll. point
					x2 = xbisect1;
					y2 = ybisect1;
					return true;
				}
			}

			xnew = x1 + (x2 - x1)*(coll.s + 0.01*iteration); // cut integration right after collision point
			if (xnew > xbisect1 && xnew < x2){ // check that new line segment does not overlap with previous one
				xbisect2 = xnew;
				stepper.calc_state(xbisect2, ybisect2);
				if (CheckHit(xbisect1, ybisect1, xbisect2, ybisect2, stepper, hitlog, iteration+1)){ // recursive call for step over coll. point
					x2 = xbisect2;
					y2 = ybisect2;
					return true;
				}
			}

			if (CheckHit(xbisect2, ybisect2, x2, y2, stepper, hitlog, iteration+1)) // recursive call for step after coll. point
				return true;
		}
	}
	else if (OnStep(x1, y1, x2, y2, stepper, GetCurrentsolid(), ID, secondaries)){ // if there was no collision: just check for absorption in solid with highest priority
		if (y2[7] != y1[7])
			Nspinflip++;
		return true;
	}
	return false;
}


void TParticle::Print(const value_type x, const state_type &y, const state_type &spin, const solid &sld, const LogStream logType) const{
	ofstream &file = GetLogStream(logType);
	if (!file.is_open()){
		ostringstream filename;
		filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name;
		if (logType == endLog)
			filename << "end.out";
		else if (logType == snapshotLog)
			filename << "snapshot.out";
		cout << "Creating " << filename.str() << '\n';
		file.open(filename.str().c_str());
		if (!file.is_open()){
			cout << "Could not create" << filename.str() << '\n';
			exit(-1);
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
		file << std::setprecision(std::numeric_limits<double>::digits); // need maximum precision for wL and delwL 
	}
	cout << "Printing status\n";

	value_type E = GetKineticEnergy(&y[3]);
	double B[4][4], Ei[3], V;

	field->BField(ystart[0], ystart[1], ystart[2], tstart, B);
	field->EField(ystart[0], ystart[1], ystart[2], tstart, V, Ei);

	file	<< jobnumber << " " << particlenumber << " "
			<< tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
			<< ystart[3] << " " << ystart[4] << " " << ystart[5] << " " << ystart[7] << " "
			<< spinstart[0] << " " << spinstart[1] << " " << spinstart[2] << " " << GetInitialTotalEnergy() << " " << GetInitialKineticEnergy() << " "
			<< B[3][0] << " " << V << " " << solidstart.ID << " ";

	double H;
	if (ID == ID_ABSORBED_ON_SURFACE){
		std::map<solid, bool>::const_iterator nextsolid = currentsolids.begin();
		do{
			nextsolid++;
		}while (nextsolid->second);
		H = E + GetPotentialEnergy(x, y, field, nextsolid->first); // if particle was absorbed on surface, use previous solid to calculate potential energy
	}
	else
		H = E + GetPotentialEnergy(x, y, field, GetCurrentsolid());

	field->BField(y[0], y[1], y[2], x, B);
	field->EField(y[0], y[1], y[2], x, V, Ei);

	double wL = 0;
	if (spin[3] > 0)
		wL = spin[4]/spin[3];

	file	<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< y[3] << " " << y[4] << " " << y[5] << " " << y[7] << " "
			<< spin[0] << " " << spin[1] << " " << spin[2] << " " << H << " " << E << " "
			<< B[3][0] << " " << V << " " << sld.ID << " "
			<< ID << " " << Nspinflip << " " << 1 - noflipprob << " "
			<< Nhit << " " << Nstep << " " << lend << " " << Hmax << " " << wL << '\n';
}


void TParticle::PrintTrack(const value_type x, const state_type &y, const state_type &spin, const solid &sld) const{
	ofstream &trackfile = GetLogStream(trackLog);
	if (!trackfile.is_open()){
		ostringstream filename;
		filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << "track.out";
		cout << "Creating " << filename.str() << '\n';
		trackfile.open(filename.str().c_str());
		if (!trackfile.is_open()){
			cout << "Could not create" << filename.str() << '\n';
			exit(-1);
		}
		trackfile << 	"jobnumber particle polarisation "
						"t x y z vx vy vz "
						"H E Bx dBxdx dBxdy dBxdz By dBydx "
						"dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V\n";
		trackfile.precision(10);
	}

	cout << "-";
	double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double E[3] = {0,0,0};
	double V = 0;
	if (field){
		field->BField(y[0],y[1],y[2],x,B);
		field->EField(y[0],y[1],y[2],x,V,E);
	}
	value_type Ek = GetKineticEnergy(&y[3]);
	value_type H = Ek + GetPotentialEnergy(x, y, field, sld);

	trackfile << jobnumber << " " << particlenumber << " " << y[7] << " "
				<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
				<< H << " " << Ek << " ";
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			trackfile << B[i][j] << " ";
	trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
}


void TParticle::PrintHit(const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering) const{
	ofstream &hitfile = GetLogStream(hitLog);
	if (!hitfile.is_open()){
		ostringstream filename;
		filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << "hit.out";
		cout << "Creating " << filename.str() << '\n';
		hitfile.open(filename.str().c_str());
		if (!hitfile.is_open()){
			cout << "Could not create" << filename.str() << '\n';
			exit(-1);
		}
		hitfile << "jobnumber particle "
					"t x y z v1x v1y v1z pol1 "
					"v2x v2y v2z pol2 "
					"nx ny nz solid1 solid2\n";
		hitfile.precision(10);
	}

	cout << ":";
	hitfile << jobnumber << " " << particlenumber << " "
			<< x << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << y1[7] << " "
			<< y2[3] << " " << y2[4] << " " << y2[5] << " " << y2[7] << " "
			<< normal[0] << " " << normal[1] << " " << normal[2] << " " << leaving.ID << " " << entering.ID << '\n';
}


void TParticle::PrintSpin(const value_type x, const state_type &spin, const dense_stepper_type &stepper) const{
	ofstream &spinfile = GetLogStream(spinLog);
	if (!spinfile.is_open()){
		std::ostringstream filename;
		filename << outpath << "/" << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << name << "spin.out";
		std::cout << "Creating " << filename.str() << '\n';
		spinfile.open(filename.str().c_str());
		if(!spinfile.is_open())
		{
			std::cout << "Could not open " << filename.str() << '\n';
			exit(-1);
		}

		//need the maximum accuracy in spinoutlog for the larmor frequency to see any difference
		spinfile << std::setprecision(std::numeric_limits<double>::digits);
		spinfile << "jobnumber particle t Sx Sy Sz Wx Wy Wz\n";
	}
	std::cout << "/";
	double Omega[3];
	SpinPrecessionAxis(x, stepper, Omega[0], Omega[1], Omega[2]);

	spinfile << jobnumber << " " << particlenumber << " "
			<< x << " " << spin[0] << " " << spin[1] << " " << spin[2] << " "
			<< Omega[0] << " " << Omega[1] << " " << Omega[2] << "\n";
}


double TParticle::GetKineticEnergy(const value_type v[3]) const{
	value_type v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	value_type beta2 = v2/c_0/c_0;
	value_type gammarel = 1/sqrt(1 - beta2); // use relativistic formula for larger beta
//			cout << "beta^8 = " << beta2*beta2*beta2*beta2 << "; err_gamma = " << numeric_limits<value_type>::epsilon()/(gammarel - 1) << '\n';
	if (beta2*beta2*beta2*beta2 < numeric_limits<value_type>::epsilon()/(gammarel - 1)) // if error in series expansion O(beta^8) is smaller than rounding error in gamma factor
		return 0.5*m*v2 + (3.0/8.0*m + (5.0/16.0*m + 35.0/128.0*beta2)*beta2)*beta2*v2; // use series expansion for energy calculation with small beta
	else{
		return c_0*c_0*m*(gammarel - 1); // else use fully relativstic formula
	}
}


double TParticle::GetPotentialEnergy(const value_type t, const state_type &y, TFieldManager *field, const solid &sld) const{
	value_type result = 0;
	if ((q != 0 || mu != 0) && field){
		double B[4][4], E[3], V;
		if (mu != 0){
			field->BField(y[0],y[1],y[2],t,B);
			result += -y[7]*mu/ele_e*B[3][0];
		}
		if (q != 0){
			field->EField(y[0],y[1],y[2],t,V,E);
			result += q/ele_e*V;
		}
	}
	result += m*gravconst*y[2];
	return result;
}
