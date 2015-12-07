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
	return Ekin(&ystart[3]) + Epot(tstart, ystart, polstart, field, geom->GetSolid(tstart, &ystart[0]));
}


double TParticle::Hend(){
	return Ekin(&yend[3]) + Epot(tend, yend, polend, field, geom->GetSolid(tend, &yend[0]));
}


double TParticle::Estart(){
	return Ekin(&ystart[3]);
}


double TParticle::Eend(){
	return Ekin(&yend[3]);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, int number,
		double t, double x, double y, double z, double E, double phi, double theta, int polarisation, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
		: name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN),
		  tstart(t), tend(t), polstart(polarisation), polend(polarisation), Hmax(0), lend(0), Nhit(0), Nspinflip(0), noflipprob(1), Nstep(0),
		  geom(&geometry), mc(&amc), field(afield){
	value_type Eoverm = E/m/c_0/c_0; // gamma - 1
	value_type beta = sqrt(1-(1/((Eoverm + 1)*(Eoverm + 1)))); // beta = sqrt(1 - 1/gamma^2)
	value_type vstart;
//			cout << "(E/m)^3.5 = " << pow(Eoverm, 3.5L) << "; err_beta = " << numeric_limits<value_tyoe>::epsilon()/(1 - beta) << '\n';
	if (pow(Eoverm, 3.5L) < numeric_limits<value_type>::epsilon()/(1 - beta)) // if error in series O((E/m)^3.5) is smaller than rounding error in beta
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
	Hmax = Hstart();
	maxtraj = mc->MaxTrajLength(name);
	tau = mc->LifeTime(name);
	geom->GetSolids(t, &ystart[0], currentsolids);
	solidend = solidstart = GetCurrentsolid();
}


TParticle::~TParticle(){
	for (vector<TParticle*>::reverse_iterator i = secondaries.rbegin(); i != secondaries.rend(); i++)
		delete *i;
}


void TParticle::operator()(state_type y, state_type &dydx, value_type x){
	derivs(x,y,dydx);
}


void TParticle::Integrate(double tmax, map<string, string> &conf){
	if (currentsolids.empty())
		geom->GetSolids(tend, &yend[0], currentsolids);

	int perc = 0;
	cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
	cout << "x: " << yend[0] << "m y: " << yend[1] << "m z: " << yend[2]
		 << "m E: " << Eend() << "eV t: " << tend << "s tau: " << tau << "s lmax: " << maxtraj << "m\n";

	// set initial values for integrator
	value_type x = tend, x1, x2;
	state_type y = yend; 
	state_type y1(6), y2(6), dy1dx(6), dy2dx(6); //acceleration is required for finding temporal derivatives of vxE
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
		PrintTrack(tend, yend, polend, solidend);
	double trackloginterval = 1e-3;
	istringstream(conf["trackloginterval"]) >> trackloginterval;
	value_type lastsave = x;

	bool hitlog = false;
	istringstream(conf["hitlog"]) >> hitlog;

	bool simulEFieldSpinInteg = false; 
	istringstream(conf["simulEFieldSpinInteg"]) >> simulEFieldSpinInteg;

	bool flipspin;
	istringstream(conf["flipspin"]) >> flipspin;
	
	TBFIntegrator *BFint1, *BFint2;
	
	if (simulEFieldSpinInteg) {
		//Different names so that two different spin logs would be produced if so required
		string eupName = string(name) + "Eup"; 
		string edownName = string(name) + "Edown";
		//One of the BFintup object will integrate the spin with the E field given 
		//and the BFintdown object will integrate the spin with the reverse of the E field given
		BFint1 = new TBFIntegrator(gamma, eupName, conf, GetSpinOut()); 
		BFint2 = new TBFIntegrator(gamma, edownName, conf, GetSpinOut2()); 
	} else
		BFint1 = new TBFIntegrator(gamma, name, conf, GetSpinOut());	
	
	stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type());

	while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
		if (resetintegration){
			stepper.initialize(y, x, h);
		}
		x1 = x; // save point before next step
		y1 = y;

		try{
			stepper.do_step(boost::ref(*this));
			x = stepper.current_time();
			y = stepper.current_state();
			h = stepper.current_time_step();
			Nstep++;
		}
		catch(...){ // catch Exceptions thrown by numerical recipes routines
			StopIntegration(ID_ODEINT_ERROR, x, y, polarisation, GetCurrentsolid());
		}

		if (stepper.current_time() > tstart + tau){
			x = tstart + tau;
			stepper.calc_state(tstart + tau, y);	//If stepsize overshot, decrease.
		}
		if (stepper.current_time() > tmax){
			x = tmax;
			stepper.calc_state(tmax, y);
		}

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

			resetintegration = CheckHit(x1, y1, x2, y2, polarisation, hitlog); // check if particle hit a material boundary or was absorbed between y1 and y2
			if (resetintegration){
				x = x2; // if particle path was changed: reset integration end point
				y = y2;
			}

			lend += sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));
			Hmax = max(Ekin(&y2[3]) + Epot(x, y2, polarisation, field, GetCurrentsolid()), Hmax);

			// take snapshots at certain times
			if (snapshotlog && snapshots.good()){
				if (x1 <= nextsnapshot && x2 > nextsnapshot){
					state_type ysnap(6);
					stepper.calc_state(nextsnapshot, ysnap);
					cout << "\n Snapshot at " << nextsnapshot << " s \n";

					PrintSnapshot(nextsnapshot, ysnap, polarisation, GetCurrentsolid());
					snapshots >> nextsnapshot;
				}
			}

			if (tracklog && x2 - lastsave > trackloginterval/v1){
				PrintTrack(x2, y2, polarisation, GetCurrentsolid());
				lastsave = x2;
			}

			if (field){
				double B1[4][4], B2[4][4], E1[3], E2[3], E1down[3], E2down[3], v_pot1, v_pot2;
				field->BField(y1[0], y1[1], y1[2], x1, B1);
				field->BField(y2[0], y2[1], y2[2], x2, B2);
				
				//Need the E field to include the vxE effect in spin tracking
				field->EField( y1[0], y1[1], y1[2], x1, v_pot1, E1 ); 
				field->EField( y2[0], y2[1], y2[2], x2, v_pot2, E2 );		
				
				//Values of dy1dx and dy2dx are updated because they are passed by reference
				derivs(x1, y1, dy1dx);
				derivs(x2, y2, dy2dx);	
				
				long double noflip;
				
				if ( simulEFieldSpinInteg ) { 
					//generate the E-field that is anti-parallel to specified E-field 
					for ( unsigned i = 0; i < sizeof(E1down)/sizeof(double); ++i) {
						E1down[i] = -E1[i];
						E2down[i] = -E2[i];
					}					
					
					noflip = BFint1->Integrate(x1, &y1[0], &dy1dx[0], B1, E1down, x2, &y2[0], &dy2dx[0], B2, E2down ); 
					noflip = BFint2->Integrate(x1, &y1[0], &dy1dx[0], B1, E1down, x2, &y2[0], &dy2dx[0], B2, E2down ); 
				} else
					noflip = BFint1->Integrate(x1, &y1[0], &dy1dx[0], B1, E1, x2, &y2[0], &dy2dx[0], B2, E2);
				
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
			StopIntegration(ID_DECAYED, x, y, polarisation, GetCurrentsolid());
		else if (ID == ID_UNKNOWN && (x >= tmax || lend >= maxtraj))
			StopIntegration(ID_NOT_FINISH, x, y, polarisation, GetCurrentsolid());
	}

	Print(tend, yend, polend, solidend);

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
	cout << "Computation took " << Nstep << " steps\n";
	cout << "Done!!\n\n";
//			cout.flush();
	
	//delete the BFIntegrator objects after the integration is done
	simulEFieldSpinInteg ? delete BFint1, delete BFint2 : delete BFint1; 
}


solid TParticle::GetCurrentsolid(){
	map<solid, bool>::iterator it = currentsolids.begin();
	while (it->second) // skip over ignored solids
		it++;
	return it->first; // return first non-ignored solid in list
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
	solid currentsolid = GetCurrentsolid();
	if (!geom->CheckSegment(&y1[0], &y2[0])){ // check if start point is inside bounding box of the simulation geometry
		printf("\nParticle has hit outer boundaries: Stopping it! t=%g x=%g y=%g z=%g\n",x2,y2[0],y2[1],y2[2]);
		StopIntegration(ID_HIT_BOUNDARIES, x2, y2, pol, currentsolid);
		return true;
	}

	map<TCollision, bool> colls;
	bool collfound = false;
	try{
		collfound = geom->GetCollisions(x1, &y1[0], x2, &y2[0], colls);
	}
	catch(...){
		StopIntegration(ID_CGAL_ERROR, x2, y2, pol, currentsolid);
		return true;
	}

	if (collfound){	// if there is a collision with a wall
		TCollision coll = colls.begin()->first;
		if ((abs(coll.s*coll.distnormal) < REFLECT_TOLERANCE
			&& (1 - coll.s)*abs(coll.distnormal) < REFLECT_TOLERANCE) // if first collision is closer to y1 and y2 than REFLECT_TOLERANCE
			|| iteration > 99) // or if iteration counter reached a certain maximum value
		{
			int prevpol = pol;
			bool trajectoryaltered = false, traversed = true;

			map<solid, bool> newsolids = currentsolids;
			for (map<TCollision, bool>::iterator it = colls.begin(); it != colls.end(); it++){ // go through list of collisions
				solid sld = geom->solids[it->first.sldindex];
				if (CheckHitError(&sld, it->first.distnormal)){ // check all hits for errors
					x2 = x1;
					for (int i = 0; i < 6; i++)
						y2[i] = y1[i];
					StopIntegration(ID_GEOMETRY_ERROR, x2, y2, pol, currentsolid); // stop trajectory integration if error found
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
				OnHit(x1, y1, x2, y2, pol, coll.normal, &leaving, &entering, trajectoryaltered, traversed); // do particle specific things
				if (hitlog)
					PrintHit(x1, y1, y2, prevpol, pol, coll.normal, &leaving, &entering); // print collision to file if requested
				Nhit++;
			}

			if (traversed){
				currentsolids = newsolids; // if surface was traversed (even if it was  physically ignored) replace current solids with list of new solids
			}

			if (trajectoryaltered)
				return true;
			else if (OnStep(x1, y1, x2, y2, pol, GetCurrentsolid())){ // check for absorption
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

			xnew = x1 + (x2 - x1)*(coll.s - 0.01*iteration); // cut integration right before collision point
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

			xnew = x1 + (x2 - x1)*(coll.s + 0.01*iteration); // cut integration right after collision point
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
	else if (OnStep(x1, y1, x2, y2, pol, GetCurrentsolid())){ // if there was no collision: just check for absorption in solid with highest priority
		return true;
	}
	return false;
}


void TParticle::StopIntegration(int aID, value_type x, state_type y, int polarisation, solid sld){
	ID = aID;
	tend = x;
	for (int i = 0; i < 6; i++)
		yend[i] = y[i];
	polend = polarisation;
	solidend = sld;
}


void TParticle::Print(std::ofstream &file, value_type x, state_type y, int polarisation, solid sld, std::string filesuffix){
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
					"vxstart vystart vzstart polstart "
					"Hstart Estart Bstart Ustart solidstart "
					"tend xend yend zend "
					"vxend vyend vzend polend "
					"Hend Eend Bend Uend solidend "
					"stopID Nspinflip spinflipprob "
					"Nhit Nstep trajlength Hmax\n";
		file.precision(10);
	}
	cout << "Printing status\n";

	value_type E = Ekin(&y[3]);
	double B[4][4], Ei[3], V;

	field->BField(ystart[0], ystart[1], ystart[2], tstart, B);
	field->EField(ystart[0], ystart[1], ystart[2], tstart, V, Ei);

	file	<< jobnumber << " " << particlenumber << " "
			<< tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
			<< ystart[3] << " " << ystart[4] << " " << ystart[5] << " "
			<< polstart << " " << Hstart() << " " << Estart() << " "
			<< B[3][3] << " " << V << " " << solidstart.ID << " ";

	field->BField(y[0], y[1], y[2], x, B);
	field->EField(y[0], y[1], y[2], x, V, Ei);

	file	<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< y[3] << " " << y[4] << " " << y[5] << " "
			<< polarisation << " " << E + Epot(x, y, polarisation, field, GetCurrentsolid()) << " " << E << " " // use GetCurrentsolid() for Epot, since particle may not actually have entered sld
			<< B[3][3] << " " << V << " " << sld.ID << " "
			<< ID << " " << Nspinflip << " " << 1 - noflipprob << " "
			<< Nhit << " " << Nstep << " " << lend << " " << Hmax << '\n';
}


void TParticle::PrintTrack(std::ofstream &trackfile, value_type x, state_type y, int polarisation, solid sld){
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
	value_type H = Ek + Epot(x, y, polarisation, field, sld);

	trackfile << jobnumber << " " << particlenumber << " " << polarisation << " "
				<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
				<< H << " " << Ek << " ";
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			trackfile << B[i][j] << " ";
	trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
}


void TParticle::PrintHit(std::ofstream &hitfile, value_type x, state_type y1, state_type y2, int pol1, int pol2, const double *normal, solid *leaving, solid *entering){
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
//			cout << "beta^8 = " << beta2*beta2*beta2*beta2 << "; err_gamma = " << numeric_limits<value_type>::epsilon()/(gammarel - 1) << '\n';
	if (beta2*beta2*beta2*beta2 < numeric_limits<value_type>::epsilon()/(gammarel - 1)) // if error in series expansion O(beta^8) is smaller than rounding error in gamma factor
		return 0.5*m*v2 + (3.0/8.0*m + (5.0/16.0*m + 35.0/128.0*beta2)*beta2)*beta2*v2; // use series expansion for energy calculation with small beta
	else{
		return c_0*c_0*m*(gammarel - 1); // else use fully relativstic formula
	}
}


double TParticle::Epot(value_type t, state_type y, int polarisation, TFieldManager *field, solid sld){
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
