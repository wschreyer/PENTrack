/**
 * \file
 * Particle base class definition.
 */

#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sys/time.h>

#include <boost/numeric/odeint.hpp>

#include "particle.h"
#include "bruteforce.h"


double TParticle::Hstart(){
	map<solid, bool> solids;
	geom->GetSolids(tstart, &ystart[0], solids);
	return Ekin(&ystart[3]) + Epot(tstart, ystart, polstart, field, solids);
}


double TParticle::Hend(){
	return Ekin(&yend[3]) + Epot(tend, yend, polend, field, currentsolids);
}


double TParticle::Estart(){
	return Ekin(&ystart[3]);
}


double TParticle::Eend(){
	return Ekin(&yend[3]);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, int number,
		double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN), inttime(0), refltime(0),
		  tstart(t), tend(t), Hmax(0), lend(0), Nhit(0), Nspinflip(0), noflipprob(1), Nstep(0),
		  geom(&geometry), mc(&amc), field(afield){
	value_type Eoverm = E/m/c_0/c_0; // gamma - 1
	value_type beta = sqrt(1-(1/((Eoverm + 1)*(Eoverm + 1)))); // beta = sqrt(1 - 1/gamma^2)
	value_type vstart;
//			cout << "(E/m)^3.5 = " << pow(Eoverm, 3.5L) << "; err_beta = " << numeric_limits<typeof(beta)>::epsilon()/(1 - beta) << '\n';
	if (pow(Eoverm, 3.5L) < numeric_limits<typeof(beta)>::epsilon()/(1 - beta)) // if error in series O((E/m)^3.5) is smaller than rounding error in beta
		vstart = c_0*(sqrt(2*Eoverm) - 3*pow(Eoverm, 1.5L)/2/sqrt(2) + 23*pow(Eoverm, 2.5L)/16/sqrt(2)); // use series expansion
	else
		vstart = c_0 * beta;

	ystart.resize(6);
	yend.resize(6);
	yend[0] = ystart[0] = x;
	yend[1] = ystart[1] = y;
	yend[2] = ystart[2] = z;
	yend[3] = ystart[3] = vstart*cos(phi)*sin(theta);
	yend[4] = ystart[4] = vstart*sin(phi)*sin(theta);
	yend[5] = ystart[5] = vstart*cos(theta);
	polstart = polend = mc->DicePolarisation(name);
	Hmax = Hstart();
	maxtraj = mc->MaxTrajLength(name);
	tau = mc->LifeTime(name);
	geom->GetSolids(t, &ystart[0], currentsolids);
}


TParticle::~TParticle(){
	for (vector<TParticle*>::reverse_iterator i = secondaries.rbegin(); i != secondaries.rend(); i++)
		delete *i;
}


void TParticle::operator()(state_type y, state_type &dydx, value_type x){
	derivs(x,y,dydx);
}


void TParticle::Integrate(double tmax, TGeometry &geometry, TMCGenerator &mcgen, TFieldManager *afield, map<string, string> &conf){
	geom = &geometry;
	if (currentsolids.empty())
		geom->GetSolids(tend, &yend[0], currentsolids);
	mc = &mcgen;
	field = afield;
	timespec clock_start, clock_end, refl_start, refl_end;

	int perc = 0;
	cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
	cout << "x: " << yend[0] << "m y: " << yend[1] << "m z: " << yend[2]
		 << "m E: " << Eend() << "eV t: " << tend << "s tau: " << tau << "s lmax: " << maxtraj << "m\n";

	// set initial values for integrator
	value_type x = tend, x1, x2;
	state_type y = yend, dydx;
	state_type y1(6), y2(6);
	int polarisation = polend;
	value_type h = 0.001/sqrt(yend[3]*yend[3] + yend[4]*yend[4] + yend[5]*yend[5]); // first guess for stepsize

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
		PrintTrack(tend, yend, polend);
	double trackloginterval = 1e-3;
	istringstream(conf["trackloginterval"]) >> trackloginterval;
	value_type lastsave = x;

	bool hitlog = false;
	istringstream(conf["hitlog"]) >> hitlog;

	bool flipspin;
	istringstream(conf["flipspin"]) >> flipspin;
	TBFIntegrator BFint(gamma, name, conf, GetSpinOut());

	stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type());

	while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
		if (resetintegration){
			stepper.initialize(y, x, h);
		}
		x1 = x; // save point before next step
		y1 = y;

		clock_gettime(CLOCK_REALTIME, &clock_start); // start computing time measure
		try{
			stepper.do_step(boost::ref(*this));
			x = stepper.current_time();
			y = stepper.current_state();
			h = stepper.current_time_step();
			Nstep++;
		}
		catch(...){ // catch Exceptions thrown by numerical recipes routines
			ID = ID_NRERROR;
		}
		clock_gettime(CLOCK_REALTIME, &clock_end);
		inttime += clock_end.tv_sec - clock_start.tv_sec + (double)(clock_end.tv_nsec - clock_start.tv_nsec)/1e9;

		if (stepper.current_time() > tstart + tau)
			stepper.calc_state(tstart + tau, y);	//If stepsize overshot, decrease.
		if (stepper.current_time() > tmax)
			stepper.calc_state(tmax, y);

		while (x1 < x){ // split integration step in pieces (x1,y1->x2,y2) with spatial length SAMPLE_DIST, go through all pieces
			value_type v1 = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			x2 = x1 + MAX_SAMPLE_DIST/v1; // time length = spatial length/velocity
			if (x2 >= x){
				x2 = x;
				y2 = y;
			}
			else{
				stepper.calc_state(x2, y2);
			}

			clock_gettime(CLOCK_REALTIME, &refl_start);
			resetintegration = CheckHit(x1, y1, x2, y2, polarisation, hitlog); // check if particle hit a material boundary or was absorbed between y1 and y2
			if (resetintegration){
				x = x2; // if particle path was changed: reset integration end point
				y = y2;
			}
			clock_gettime(CLOCK_REALTIME, &refl_end);
			refltime += refl_end.tv_sec - refl_start.tv_sec + (double)(refl_end.tv_nsec - refl_start.tv_nsec)/1e9;

			lend += sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));
			Hmax = max(Ekin(&y2[3]) + Epot(x, y2, polarisation, field, currentsolids), Hmax);

			// take snapshots at certain times
			if (snapshotlog && snapshots.good()){
				if (x1 <= nextsnapshot && x2 > nextsnapshot){
					state_type ysnap(6);
					stepper.calc_state(nextsnapshot, ysnap);
					cout << "\n Snapshot at " << nextsnapshot << " s \n";

					PrintSnapshot(nextsnapshot, ysnap, polarisation);
					snapshots >> nextsnapshot;
				}
			}

			if (tracklog && x2 - lastsave > trackloginterval/v1){
				PrintTrack(x2, y2, polarisation);
				lastsave = x2;
			}

			if (field){
				double B1[4][4], B2[4][4];
				field->BField(y1[0], y1[1], y1[2], x1, B1);
				field->BField(y2[0], y2[1], y2[2], x2, B2);
				long double noflip = BFint.Integrate(x1, &y1[0], B1, x2, &y2[0], B2);
				if (mc->UniformDist(0,1) > noflip)
					polarisation *= -1;
				noflipprob *= noflip; // accumulate no-spin-flip probability
			}

			x1 = x2;
			for (int i = 0; i < 6; i++)
				y1[i] = y2[i];
			polend = polarisation;
		}

		PrintPercent(max((x - tstart)/tau, max((x - tstart)/(tmax - tstart), lend/maxtraj)), perc);

		// x >= tstart + tau?
		if (ID == ID_UNKNOWN && x >= tstart + tau)
			ID = ID_DECAYED;
		else if (ID == ID_UNKNOWN && (x >= tmax || lend >= maxtraj))
			ID = ID_NOT_FINISH;
	}

	tend = x;
	for (int i = 0; i < 6; i++)
		yend[i] = y[i];
	polend = polarisation;
	Print(tend, yend, polend);

	if (ID == ID_DECAYED){ // if particle reached its lifetime call TParticle::Decay
		cout << "Decayed!\n";
		Decay();
	}

	cout << "x: " << yend[0];
	cout << " y: " << yend[1];
	cout << " z: " << yend[2];
	cout << " E: " << Eend();
	cout << " Code: " << ID;
	cout << " t: " << tend;
	cout << " l: " << lend;
	cout << " hits: " << Nhit;
	cout << " spinflips: " << Nspinflip << '\n';
	cout << "Computation took " << Nstep << " steps, " << inttime << " s for integration and " << refltime << " s for geometry checks\n";
	cout << "Done!!\n\n";
//			cout.flush();
}


solid TParticle::GetCurrentsolid(){
	map<solid, bool>::iterator it = currentsolids.begin();
	while (it->second)
		it++;
	return it->first;
}


void TParticle::derivs(value_type x, state_type y, state_type &dydx){
	dydx[0] = y[3]; // time derivatives of position = velocity
	dydx[1] = y[4];
	dydx[2] = y[5];
	value_type F[3] = {0,0,0}; // Force in lab frame
	if (m != 0)
		F[2] += -gravconst*m*ele_e; // add gravitation to force
	if (field){
		double B[4][4], E[3], V; // magnetic/electric field and electric potential in lab frame
		if (q != 0 || (mu != 0 && polend != 0))
			field->BField(y[0],y[1],y[2], x, B); // if particle has charge or magnetic moment, calculate magnetic field
		if (q != 0)
			field->EField(y[0],y[1],y[2], x, V, E); // if particle has charge caculate electric field+potential
		if (q != 0){
			F[0] += q*(E[0] + y[4]*B[2][0] - y[5]*B[1][0]); // add Lorentz-force
			F[1] += q*(E[1] + y[5]*B[0][0] - y[3]*B[2][0]);
			F[2] += q*(E[2] + y[3]*B[1][0] - y[4]*B[0][0]);
		}
		if (mu != 0 && polend != 0){
			F[0] += polend*mu*B[3][1]; // add force on magnetic dipole moment
			F[1] += polend*mu*B[3][2];
			F[2] += polend*mu*B[3][3];
		}
	}
	value_type rel = sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0))/m/ele_e; // relativstic factor 1/gamma/m
	dydx[3] = rel*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
	dydx[4] = rel*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
	dydx[5] = rel*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);
}


bool TParticle::CheckHitError(solid *hitsolid, double distnormal){
	if (distnormal < 0){ // particle is entering solid
		if (currentsolids.find(*hitsolid) != currentsolids.end()){ // if solid already is in currentsolids list something went wrong
			cout << "Particle entering " << hitsolid->name.c_str() << " which it did enter before! Stopping it! Did you define overlapping solids with equal priorities?\n";
			return true;
		}
	}
	else if (distnormal > 0){ // particle is leaving solid
		if (currentsolids.find(*hitsolid) == currentsolids.end()){ // if solid is not in currentsolids list something went wrong
			cout << "Particle inside '" << hitsolid->name << "' which it did not enter before! Stopping it!\n";
			return true;
		}
	}
	else{
		cout << "Particle is crossing material boundary with parallel track! Stopping it!\n";
		return true;
	}
	return false;
}


bool TParticle::CheckHit(value_type x1, state_type y1, value_type &x2, state_type &y2, int &pol, bool hitlog, int iteration){
	if (!geom->CheckSegment(&y1[0], &y2[0])){ // check if start point is inside bounding box of the simulation geometry
		printf("\nParticle has hit outer boundaries: Stopping it! t=%g x=%g y=%g z=%g\n",x2,y2[0],y2[1],y2[2]);
		ID = ID_HIT_BOUNDARIES;
		return true;
	}

	map<TCollision, bool> colls;
	if (geom->GetCollisions(x1, &y1[0], x2, &y2[0], colls)){	// if there is a collision with a wall
		map<TCollision, bool>::iterator coll = colls.begin();
		if ((abs(coll->first.s*coll->first.distnormal) < REFLECT_TOLERANCE
			&& (1 - coll->first.s)*abs(coll->first.distnormal) < REFLECT_TOLERANCE) // if first collision is closer to y1 and y2 than REFLECT_TOLERANCE
			|| iteration > 99) // or if iteration counter reached a certain maximum value
		{
			int prevpol = pol;
			solid currentsolid = GetCurrentsolid();
			solid *hitsolid = &geom->defaultsolid;
			coll = colls.end();
			for (map<TCollision, bool>::iterator it = colls.begin(); it != colls.end(); it++){
				if (CheckHitError(&geom->solids[it->first.sldindex], it->first.distnormal)){ // check all hits for errors
					x2 = x1;
					for (int i = 0; i < 6; i++)
						y2[i] = y1[i];
					ID = ID_NRERROR;
					return true;
				}
				if (!it->second && (geom->solids[it->first.sldindex].ID > hitsolid->ID)){ // find non-ignored collision with highest priority
					coll = it;
					hitsolid = &geom->solids[it->first.sldindex];
				}
			}

			bool trajectoryaltered = false, traversed = true;

			map<solid, bool> newsolids = currentsolids;
			for (map<TCollision, bool>::iterator it = colls.begin(); it != colls.end(); it++){ // compile list of solids after material boundary
				if (it->first.distnormal < 0){
					newsolids[geom->solids[it->first.sldindex]] = it->second;
				}else{
					newsolids.erase(geom->solids[it->first.sldindex]);
				}
			}

			if (coll != colls.end()){ // if a non-ignored collision was found
				solid leaving, entering;
				bool hit = false;
				if (coll->first.distnormal < 0){ // particle is entering solid
					leaving = currentsolid;
					entering = *hitsolid;
					if (entering.ID > leaving.ID){ // check for reflection only if priority of entered solid is larger than that of current solid
						hit = true;
						OnHit(x1, y1, x2, y2, pol, coll->first.normal, &leaving, &entering, trajectoryaltered, traversed); // do particle specific things
					}
				}
				else if (coll->first.distnormal > 0){ // particle is leaving solid
					leaving = *hitsolid;
					map<solid, bool>::iterator it = newsolids.begin();
					while (it->second) // find first non-ignored solid in list of new solids
						it++;
					entering = it->first;
					if (leaving.ID == currentsolid.ID){ // check for reflection only, if left solid is the solid with highest priority
						hit = true;
						OnHit(x1, y1, x2, y2, pol, coll->first.normal, &leaving, &entering, trajectoryaltered, traversed); // do particle specific things
					}
				}

				if (hit){
					if (hitlog)
						PrintHit(x1, y1, y2, prevpol, pol, coll->first.normal, &leaving, &entering);
					Nhit++;
				}
			}

			if (traversed){
				currentsolids = newsolids; // if material boundary was traversed, replace current solids with list of new solids
			}

			if (trajectoryaltered)
				return true;
			else if (OnStep(x1, y1, x2, y2, GetCurrentsolid())){ // check for absorption
				printf("Absorption!\n");
				return true;
			}

		}
		else{
			// else cut integration step right before and after first collision point
			// and call ReflectOrAbsorb again for each smaller step (quite similar to bisection algorithm)
			value_type xnew, xbisect1 = x1, xbisect2 = x1;
			state_type ybisect1(6);
			state_type ybisect2(6);
			ybisect1 = ybisect2 = y1;

			xnew = x1 + (x2 - x1)*(coll->first.s - 0.01*iteration); // cut integration right before collision point
			if (xnew > x1 && xnew < x2){ // check that new line segment is in correct time interval
				xbisect1 = xbisect2 = xnew;
				stepper.calc_state(xbisect1, ybisect1);
				ybisect2 = ybisect1;
				if (CheckHit(x1, y1, xbisect1, ybisect1, pol, hitlog, iteration+1)){ // recursive call for step before coll. point
					x2 = xbisect1;
					y2 = ybisect1;
					return true;
				}
			}

			xnew = x1 + (x2 - x1)*(coll->first.s + 0.01*iteration); // cut integration right after collision point
			if (xnew > xbisect1 && xnew < x2){ // check that new line segment does not overlap with previous one
				xbisect2 = xnew;
				stepper.calc_state(xbisect2, ybisect2);
				if (CheckHit(xbisect1, ybisect1, xbisect2, ybisect2, pol, hitlog, iteration+1)){ // recursive call for step over coll. point
					x2 = xbisect2;
					y2 = ybisect2;
					return true;
				}
			}

			if (CheckHit(xbisect2, ybisect2, x2, y2, pol, hitlog, iteration+1)) // recursive call for step after coll. point
				return true;
		}
	}
	else if (OnStep(x1, y1, x2, y2, GetCurrentsolid())){ // if there was no collision: just check for absorption in solid with highest priority
		return true;
	}
	return false;
}


void TParticle::Print(ofstream &file, value_type x, state_type y, int polarisation, string filesuffix){
	if (!file.is_open()){
		ostringstream filename;
		filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << filesuffix;
		cout << "Creating " << filename.str() << '\n';
		file.open(filename.str().c_str());
		if (!file.is_open()){
			cout << "Could not create" << filename.str() << '\n';
			exit(-1);
		}
		file <<	"jobnumber particle "
					"tstart xstart ystart zstart "
					"vxstart vystart vzstart "
					"polstart Hstart Estart "
					"tend xend yend zend "
					"vxend vyend vzend "
					"polend Hend Eend stopID Nspinflip spinflipprob "
					"ComputingTime Nhit Nstep trajlength Hmax\n";
		file.precision(10);
	}
	cout << "Printing status\n";
	value_type E = Ekin(&y[3]);
	file	<<	jobnumber << " " << particlenumber << " "
			<<	tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
			<<	ystart[3] << " " << ystart[4] << " " << ystart[5] << " "
			<< polstart << " " << Hstart() << " " << Estart() << " "
			<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< y[3] << " " << y[4] << " " << y[5] << " "
			<< polarisation << " " << E + Epot(x, y, polarisation, field, currentsolids) << " " << E << " " << ID << " " << Nspinflip << " " << 1 - noflipprob << " "
			<< inttime << " " << Nhit << " " << Nstep << " " << lend << " " << Hmax << '\n';
}


void TParticle::PrintTrack(ofstream &trackfile, value_type x, state_type y, int polarisation){
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
	value_type Ek = Ekin(&y[3]);
	value_type H = Ek + Epot(x, y, polarisation, field, currentsolids);

	trackfile << jobnumber << " " << particlenumber << " " << polarisation << " "
				<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
				<< H << " " << Ek << " ";
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			trackfile << B[i][j] << " ";
	trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
}


void TParticle::PrintHit(ofstream &hitfile, value_type x, state_type y1, state_type y2, int pol1, int pol2, const double *normal, solid *leaving, solid *entering){
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
			<< x << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << pol1 << " "
			<< y2[3] << " " << y2[4] << " " << y2[5] << " " << pol2 << " "
			<< normal[0] << " " << normal[1] << " " << normal[2] << " " << leaving->ID << " " << entering->ID << '\n';
}


double TParticle::Ekin(value_type v[3]){
	value_type v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	value_type beta2 = v2/c_0/c_0;
	value_type gammarel = 1/sqrt(1 - beta2); // use relativistic formula for larger beta
//			cout << "beta^8 = " << beta2*beta2*beta2*beta2 << "; err_gamma = " << numeric_limits<typeof(gammarel)>::epsilon()/(gammarel - 1) << '\n';
	if (beta2*beta2*beta2*beta2 < numeric_limits<typeof(gammarel)>::epsilon()/(gammarel - 1)) // if error in series expansion O(beta^8) is smaller than rounding error in gamma factor
		return 0.5*m*v2 + (3.0/8.0*m + (5.0/16.0*m + 35.0/128.0*beta2)*beta2)*beta2*v2; // use series expansion for energy calculation with small beta
	else{
		return c_0*c_0*m*(gammarel - 1);
	}
}


double TParticle::Epot(value_type t, state_type y, int polarisation, TFieldManager *field, map<solid, bool> &solids){
	value_type result = 0;
	if ((q != 0 || mu != 0) && field){
		double B[4][4], E[3], V;
		if (mu != 0){
			field->BField(y[0],y[1],y[2],t,B);
			result += -polarisation*mu/ele_e*B[3][0];
		}
		if (q != 0){
			field->EField(y[0],y[1],y[2],t,V,E);
			result += q/ele_e*V;
		}
	}
	result += m*gravconst*y[2];
	return result;
}
