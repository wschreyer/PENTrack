/**
 * \file
 * Particle base class definition.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_


#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sys/time.h>

static const double MAX_SAMPLE_DIST = 0.01; ///< max spatial distance of reflection checks, spin flip calculation, etc; longer integration steps will be interpolated
static const double RELATIVISTIC_THRESHOLD = 0.01; ///< threshold (v/c) above which kinetic energy is calculated relativisticly

#include "globals.h"
#include "fields.h"
#include "mc.h"
#include "geometry.h"
#include "bruteforce.h"
#include "ndist.h"
#include "adiabacity.h"


/**
 * Struct containing all output streams
 */
struct TOutput{
	ofstream *endout; ///< endlog file stream
	ofstream *snapshotout; ///< snapshot file stream
	ofstream *trackout; ///< tracklog file stream
	ofstream *hitout; ///< hitlog file stream
	ofstream *spinout; ///< spinlog file stream
	
	/**
	 * Constructor: initializes empty streams
	 */
	TOutput(): endout(NULL), snapshotout(NULL), trackout(NULL), hitout(NULL), spinout(NULL) { };

	/**
	 * Destructor, close open streams.
	 */
	~TOutput(){
		if (endout){
			endout->close();
			delete endout;
		}
		if (snapshotout){
			snapshotout->close();
			delete snapshotout;
		}
		if (trackout){
			trackout->close();
			delete trackout;
		}
		if (hitout){
			hitout->close();
			delete hitout;
		}
		if (spinout){
			spinout->close();
			delete spinout;
		}
	}
};


/**
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor and methods TParticle::OnHit, TParticle::OnStep and TParticle::Decay.
 * Optionally, derived particles can also re-implement TParticle::Epot, TParticle::PrintStartEnd, TParticle::PrintTrack, TParticle::PrintSnapshots, TParticle::PrintHits and define its own constructors.
 */
struct TParticle{
	public:
		const char *name; ///< particle name (has to be initialized in all derived classes!)
		const long double q; ///< charge [C] (has to be initialized in all derived classes!)
		const long double m; ///< mass [eV/c^2] (has to be initialized in all derived classes!)
		const long double mu; ///< magnetic moment [J/T] (has to be initialized in all derived classes!)
		const long double gamma; ///< gyromagnetic ratio [rad/(s T)] (has to be initialized in all derived classes!)
		int ID; ///< particle fate (0: not classified; >0: Absorption in solid with according ID; -1: survived until end of simulation; -2: hit outer boundaries; -3: numerical error; -4: decayed; -5: did not find initial position)
		int particlenumber; ///< particle number
		long double tau; ///< particle life time
		long double maxtraj; ///< max. simulated trajectory length
		long double inttime; ///< integration computing time
		long double refltime; ///< reflection computing time

		/// start time
		long double tstart;

		/// stop time
		long double tend;

		/// state vector before integration
		long double ystart[6];

		/// state vector after integration)
		long double yend[6];

		/// initial polarisation of particle (-1,0,1)
		int polstart;

		/// final polarisation of particle (-1,0,1)
		int polend;

		/// total start energy
		long double Hstart(){
			map<solid, bool> solids;
			geom->GetSolids(tstart, ystart, solids);
			return Ekin(&ystart[3]) + Epot(tstart, ystart, polstart, field, solids);
		};

		/// total end energy
		long double Hend(){ return Ekin(&yend[3]) + Epot(tend, yend, polend, field, currentsolids); };

		/// max total energy
		long double Hmax;

		/// kinetic start energy
		long double Estart(){ return Ekin(&ystart[3]); };

		/// kinetic end energy
		long double Eend(){ return Ekin(&yend[3]); };

		/// trajectory length
		long double lend;

		/// number of material boundary hits
		int Nhit;

		/// number of spin flips
		int Nspinflip;

		/// total probability of NO spinflip calculated by spin tracking
		long double noflipprob;

		/// number of integration steps
		int Nstep;

		vector<TParticle*> secondaries; ///< list of secondary particles


		/**
		 * Constructor, initializes TParticle::type, TParticle::q, TParticle::m, TParticle::mu
		 *
		 * Has to be called by every derived class constructor with the respective values
		 * 
		 * @param aname Particle name
		 * @param qq Electric charge
		 * @param mm Mass
		 * @param mumu Magnetic dipole moment
		 * @param agamma Gyromagnetic ratio
		 * @param number Particle number
		 * @param t Starting time
		 * @param x Initial x coordinate
		 * @param y Initial y coordinate
		 * @param z Initial z coordinate
		 * @param E Initial kinetic energy
		 * @param phi Azimuth of initial velocity vector
		 * @param theta Polar angle of initial velocity vector
		 * @param amc Random number generator
		 * @param geometry Experiment geometry
		 * @param afield Optional fields (can be NULL)
		 */
		TParticle(const char *aname, const long double qq, const long double mm, const long double mumu, const long double agamma, int number,
				double t, double x, double y, double z, double E, double phi, double theta, TMCGenerator &amc, TGeometry &geometry, TFieldManager *afield)
				: name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN), tstart(t), tend(t),
				  Nhit(0), Hmax(0), Nstep(0), refltime(0), stepper(NULL), inttime(0), lend(0), polstart(0), Nspinflip(0), noflipprob(1),
				  field(afield), mc(&amc), geom(&geometry){
			long double Eoverm = E/m/c_0/c_0;
			long double vstart;
			if (Eoverm < 0.0001)
				vstart = c_0*(sqrt(2*Eoverm) - 3*pow(Eoverm, 1.5L)/2/sqrt(2) + 23*pow(Eoverm, 2.5L)/16/sqrt(2));
			else
				vstart = c_0 * sqrt(1-(1/((Eoverm + 1)*(Eoverm + 1))));

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
			geom->GetSolids(t, ystart, currentsolids);
		};

		/**
		 * Destructor, deletes secondaries
		 */
		virtual ~TParticle(){
			for (vector<TParticle*>::reverse_iterator i = secondaries.rbegin(); i != secondaries.rend(); i++)
				delete *i;
		};


		/**
		 * Returns equations of motion.
		 *
		 * Class TParticle is given to integrator, which calls TParticle(x,y,dydx)
		 *
		 * @param x Time
		 * @param y State vector (position + velocity)
		 * @param dydx Returns derivatives of y with respect to t
		 */
		void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx){
			derivs(x,y,dydx);
		};
		

		/**
		 * Integrate particle trajectory.
		 *
		 * Takes inital state vector ystart and integrates the trajectory step by step.
		 * If a step is longer than MAX_SAMPLE_DIST, the step is split by interpolating intermediate points.
		 * On each step it checks for interaction with solids, prints snapshots and track into files and calls "OnEachStep".
		 * The Integration stops if TParticle::tau or tmax are reached; or if something happens to the particle (absorption, numerical error)
		 *
		 * @param tmax Max. absolute time at which integration will be stopped
		 * @param geometry Integrator checks for collisions with walls defined in this TGeometry structure
		 * @param mcgen Random number generator (e.g. used for reflection probability dicing)
		 * @param afield Pointer to field which acts on particle (can be NULL)
		 * @param conf Option map containing particle specific options from particle.in
		 * @param output TOutput struct containing output options and streams
		 */
		void Integrate(long double tmax, TGeometry &geometry, TMCGenerator &mcgen, TFieldManager *afield, map<string, string> &conf, TOutput &output){
			geom = &geometry;
			if (currentsolids.empty())
				geom->GetSolids(tend, yend, currentsolids);
			mc = &mcgen;
			field = afield;
			timespec clock_start, clock_end, refl_start, refl_end;

			int perc = 0;
			cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
			cout << "x: " << yend[0] << " y: " << yend[1] << " z: " << yend[2]
			     << " E: " << Eend() << " t: " << tend << " tau: " << tau << " lmax: " << maxtraj << '\n';
		
			// set initial values for integrator
			long double x = tend, x1, x2;
			VecDoub y(6, yend), dydx(6);
			long double y1[6], y2[6];
			int polarisation = polend;
			long double h = 0.001/sqrt(yend[3]*yend[3] + yend[4]*yend[4] + yend[5]*yend[5]); // first guess for stepsize

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
				PrintTrack(output.trackout, tend, yend, polend);
			double trackloginterval = 1e-3;
			istringstream(conf["trackloginterval"]) >> trackloginterval;
			long double lastsave = x;

			bool hitlog = false;
			istringstream(conf["hitlog"]) >> hitlog;

			bool flipspin;
			istringstream(conf["flipspin"]) >> flipspin;
			TBFIntegrator BFint(gamma, name, conf);

			while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
				if (resetintegration){
					if (stepper != NULL)
						delete stepper;
					derivs(x,y,dydx);
					// create integrator class (stepperdopr853 = 8th order Runge Kutta)
					stepper = new StepperDopr853<TParticle>(y, dydx, x, 1e-13, 0, true);// y, dydx, x, atol, rtol, dense output
				}
				x1 = x; // save point before next step
				for (int i = 0; i < 6; i++)
					y1[i] = y[i];
				
				if (x + h > tstart + tau)
					h = tstart + tau - x;	//If stepsize can overshoot, decrease.
				if (x + h > tmax)
					h = tmax - x;
				
				clock_gettime(CLOCK_REALTIME, &clock_start); // start computing time measure
				try{
					stepper->step(h, *this); // integrate one step
					h = stepper->hnext; // set suggested stepsize for next step
					Nstep++;
				}
				catch(...){ // catch Exceptions thrown by numerical recipes routines
					ID = ID_NRERROR;
				}
				clock_gettime(CLOCK_REALTIME, &clock_end);
				inttime += clock_end.tv_sec - clock_start.tv_sec + (long double)(clock_end.tv_nsec - clock_start.tv_nsec)/1e9;
				
				while (x1 < x){ // split integration step in pieces (x1,y1->x2,y2) with spatial length SAMPLE_DIST, go through all pieces
					long double v1 = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
					x2 = x1 + MAX_SAMPLE_DIST/v1; // time length = spatial length/velocity
					if (x2 >= x){
						x2 = x;
						for (int i = 0; i < 6; i++)
							y2[i] = y[i];
					}
					else{
						for (int i = 0; i < 6; i++)
							y2[i] = stepper->dense_out(i, x2, stepper->hdid); // interpolate y2
					}
					
					clock_gettime(CLOCK_REALTIME, &refl_start);
					resetintegration = CheckHit(x1, y1, x2, y2, polarisation, hitlog, output.hitout); // check if particle hit a material boundary or was absorbed between y1 and y2
					if (resetintegration){
						x = x2; // if particle path was changed: reset integration end point
						for (int i = 0; i < 6; i++)
							y[i] = y2[i];
					}
					clock_gettime(CLOCK_REALTIME, &refl_end);
					refltime += refl_end.tv_sec - refl_start.tv_sec + (long double)(refl_end.tv_nsec - refl_start.tv_nsec)/1e9;
				
					lend += sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));
					Hmax = max(Ekin(&y2[3]) + Epot(x, y2, polarisation, field, currentsolids), Hmax);

					// take snapshots at certain times
					if (snapshotlog && snapshots.good()){
						if (x1 <= nextsnapshot && x2 > nextsnapshot){
							long double ysnap[6];
							for (int i = 0; i < 6; i++)
								ysnap[i] = stepper->dense_out(i, nextsnapshot, stepper->hdid);
							cout << "\n Snapshot at " << nextsnapshot << " s \n";

							Print(output.snapshotout, nextsnapshot, ysnap, polarisation, "snapshot.out");
							snapshots >> nextsnapshot;
						}
					}

					if (tracklog && x2 - lastsave > trackloginterval/v1){
						PrintTrack(output.trackout, x2, y2, polarisation);
						lastsave = x2;
					}

					long double noflip = BFint.Integrate(x1, y1, x2, y2, field, output.spinout);
					if (mc->UniformDist(0,1) > noflip)
						polarisation *= -1;
					noflipprob *= noflip; // accumulate no-spin-flip probability

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
			Print(output.endout, tend, yend, polend);

			delete stepper;
			stepper = NULL;
			
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
		};
	

	protected:
		map<solid, bool> currentsolids; ///< solids in which particle is currently inside
		TGeometry *geom; ///< TGeometry structure passed by "Integrate"
		TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
		TFieldManager *field; ///< TFieldManager structure passed by "Integrate"
		StepperDopr853<TParticle> *stepper; ///< ODE integrator

		solid GetCurrentsolid(){
			map<solid, bool>::iterator it = currentsolids.begin();
			while (it->second)
				it++;
			return it->first;
		}


		/**
		 * Equations of motion dy/dx = f(x,y).
		 *
		 * Equations of motion (fully relativistic), called by operator().
		 * Including gravitation, Lorentz-force and magnetic interaction with magnetic moment.
		 *
		 * @param x Time
		 * @param y	State vector (position and velocity)
		 * @param dydx Returns derivatives of y with respect to x
		 */
		void derivs(const Doub x, VecDoub_I &y, VecDoub_O &dydx){
			dydx[0] = y[3]; // time derivatives of position = velocity
			dydx[1] = y[4];
			dydx[2] = y[5];
			long double F[3] = {0,0,0}; // Force in lab frame
			if (m != 0)
				F[2] += -gravconst*m*ele_e; // add gravitation to force
			if (field){
				long double B[4][4], E[3], V; // magnetic/electric field and electric potential in lab frame
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
			long double rel = sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0))/m/ele_e; // relativstic factor 1/gamma/m
			dydx[3] = rel*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
			dydx[4] = rel*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
			dydx[5] = rel*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);
		};
		

		bool CheckHitError(solid *hitsolid, long double distnormal){
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


		/**
		 * Check, if particle hit a material boundary or was absorbed.
		 *
		 * Checks if a particle which flies from y1 to y2 in time x2-x1 hits a surface or is absorbed inside a material.
		 * If a surface is hit the routine splits the line segment y1->y2 on both sides of the collision point and calls itself recursively
		 * with the three new line segments as parameters. This is repeated until both split points are nearer than REFLECTION_TOLERANCE
		 * to the collision point (much like an bisection algorithm). For each line segment "OnStep" is called to check for scattering/absorption/etc.
		 * For each short segment crossing a collision point "OnHit" is called to check for reflection/refraction/etc.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment
		 * @param y2 End point of line segment
		 * @param pol Particle polarisation
		 * @param hitlog Should hits be logged to file?
		 * @param hitout File stream to which hits are logged.
		 * @param iteration Iteration counter (incremented by recursive calls to avoid infinite loop)
		 * @return Returns true if particle was reflected/absorbed
		 */
		bool CheckHit(long double &x1, long double y1[6], long double &x2, long double y2[6], int &pol, bool hitlog, ofstream *&hitout, int iteration = 1){
			if (!geom->CheckSegment(y1,y2)){ // check if start point is inside bounding box of the simulation geometry
				printf("\nParticle has hit outer boundaries: Stopping it! t=%LG x=%LG y=%LG z=%LG\n",x2,y2[0],y2[1],y2[2]);
				ID = ID_HIT_BOUNDARIES;
				return true;
			}

			map<TCollision, bool> colls;
			if (geom->GetCollisions(x1, y1, stepper->hdid, y2, colls)){	// if there is a collision with a wall
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
								PrintHit(hitout, x1, y1, y2, prevpol, pol, coll->first.normal, &leaving, &entering);
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
					long double xnew, xbisect1 = x1, xbisect2 = x1;
					long double ybisect1[6];
					long double ybisect2[6];
					for (int i = 0; i < 6; i++)
						ybisect1[i] = ybisect2[i] = y1[i];

					xnew = x1 + (x2 - x1)*(coll->first.s - 0.01*iteration); // cut integration right before collision point
					if (xnew > x1 && xnew < x2){ // check that new line segment is in correct time interval
						xbisect1 = xbisect2 = xnew;
						for (int i = 0; i < 6; i++)
							ybisect1[i] = ybisect2[i] = stepper->dense_out(i, xbisect1, stepper->hdid); // get point at cut time
						if (CheckHit(x1, y1, xbisect1, ybisect1, pol, hitlog, hitout, iteration+1)){ // recursive call for step before coll. point
							x2 = xbisect1;
							for (int i = 0; i < 6; i++)
								y2[i] = ybisect1[i];
							return true;
						}
					}

					xnew = x1 + (x2 - x1)*(coll->first.s + 0.01*iteration); // cut integration right after collision point
					if (xnew > xbisect1 && xnew < x2){ // check that new line segment does not overlap with previous one
						xbisect2 = xnew;
						for (int i = 0; i < 6; i++)
							ybisect2[i] = stepper->dense_out(i, xbisect2, stepper->hdid); // get point at cut time
						if (CheckHit(xbisect1, ybisect1, xbisect2, ybisect2, pol, hitlog, hitout, iteration+1)){ // recursive call for step over coll. point
							x2 = xbisect2;
							for (int i = 0; i < 6; i++)
								y2[i] = ybisect2[i];
							return true;
						}
					}

					if (CheckHit(xbisect2, ybisect2, x2, y2, pol, hitlog, hitout, iteration+1)) // recursive call for step after coll. point
						return true;
				}
			}
			else if (OnStep(x1, y1, x2, y2, GetCurrentsolid())){ // if there was no collision: just check for absorption in solid with highest priority
				return true;
			}
			return false;
		};


		/**
		 * This virtual method is executed, when a particle crosses a material boundary.
		 *
		 * Can be used to reflect/refract/etc. at material boundaries.
		 * Has to be implemented separately for each derived particle.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment, may be altered
		 * @param y2 End point of line segment, may be altered
		 * @param polarisation Polarisation of particle, may be altered
		 * @param normal Normal vector of hit surface
		 * @param leaving Solid that the particle is leaving
		 * @param entering Solid that the particle is entering
		 * @param trajectoryaltered Returns true if the particle trajectory was altered
		 * @param travered Returns true if the material boundary was traversed by the particle
		 */
		virtual void OnHit(long double x1, long double y1[6], long double &x2, long double y2[6], int &polarisation,
							const long double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed) = 0;


		/**
		 * This virtual method is executed on each integration step.
		 *
		 * Can be used to scatter/absorb/etc. in materials.
		 * Has to be implemented separately for each derived particle.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment, may be altered
		 * @param y2 End point of line segment, may be altered
		 * @param currentsolid Solid through which the particle is moving
		 * @return Returns true if particle path was changed
		 */
		virtual bool OnStep(long double x1, long double y1[6], long double &x2, long double y2[6], solid currentsolid) = 0;


		/**
		 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
		 */
		virtual void Decay() = 0;



		/**
		 * Print start and current values to a stream.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param file stream to print into
		 * @param x Current time
		 * @param y Current state vector
		 * @param polarisation Current polarisation
		 * @param filesuffix Optional suffix added to the file name (default: "end.out")
		 */
		virtual void Print(ofstream *&file, long double x, long double y[6], int polarisation, string filesuffix = "end.out"){
			if (!file){
				ostringstream filename;
				filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << filesuffix;
				cout << "Creating " << filename.str() << '\n';
				file = new ofstream(filename.str().c_str());
				if (!file || !file->is_open()){
					cout << "Could not create" << filename.str() << '\n';
					exit(-1);
				}
				*file <<	"jobnumber particle "
	                		"tstart xstart ystart zstart "
	                		"vxstart vystart vzstart "
	                		"polstart Hstart Estart "
	                		"tend xend yend zend "
	                		"vxend vyend vzend "
	                		"polend Hend Eend stopID Nspinflip spinflipprob "
	                		"ComputingTime Nhit Nstep trajlength Hmax\n";
				file->precision(10);
			}
			cout << "Printing status\n";
			long double E = Ekin(&y[3]);
			*file	<<	jobnumber << " " << particlenumber << " "
					<<	tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
					<<	ystart[3] << " " << ystart[4] << " " << ystart[5] << " "
					<< polstart << " " << Hstart() << " " << Estart() << " "
					<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
					<< y[3] << " " << y[4] << " " << y[5] << " "
					<< polarisation << " " << E + Epot(x, y, polarisation, field, currentsolids) << " " << E << " " << ID << " " << Nspinflip << " " << 1 - noflipprob << " "
					<< inttime << " " << Nhit << " " << Nstep << " " << lend << " " << Hmax << '\n';
		};


		/**
		 * Print current track point into stream to allow visualization of the particle's trajectory.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param trackfile Stream to print into
		 * @param x Current time
		 * @param y Current state vector
		 * @param polarisation Current polarisation
		 */
		virtual void PrintTrack(ofstream *&trackfile, long double x, long double y[6], int polarisation){
			if (!trackfile){
				ostringstream filename;
				filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << "track.out";
				cout << "Creating " << filename.str() << '\n';
				trackfile = new ofstream(filename.str().c_str());
				if (!trackfile || !trackfile->is_open()){
					cout << "Could not create" << filename.str() << '\n';
					exit(-1);
				}
				*trackfile << 	"jobnumber particle polarisation "
								"t x y z vx vy vz "
								"H E Bx dBxdx dBxdy dBxdz By dBydx "
								"dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V\n";
				trackfile->precision(10);
			}

			cout << "-";
			long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
			long double E[3] = {0,0,0};
			long double V = 0;
			if (field){
				field->BField(y[0],y[1],y[2],x,B);
				field->EField(y[0],y[1],y[2],x,V,E);
			}
			long double Ek = Ekin(&y[3]);
			long double H = Ek + Epot(x, y, polarisation, field, currentsolids);

			*trackfile << jobnumber << " " << particlenumber << " " << polarisation << " "
						<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
						<< H << " " << Ek << " ";
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					*trackfile << B[i][j] << " ";
			*trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
		};


		/**
		 * Print material boundary hits into stream.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param hitfile stream to print into
		 * @param x Time of material hit
		 * @param y1 State vector before material hit
		 * @param y2 State vector after material hit
		 * @param pol1 Polarisation before material hit
		 * @param pol2 Polarisation after material hit
		 * @param normal Normal vector of hit surface
		 * @param leaving Material which is left at this boundary
		 * @param entering Material which is entered at this boundary
		 */
		virtual void PrintHit(ofstream *&hitfile, long double x, long double *y1, long double *y2, int pol1, int pol2, const long double *normal, solid *leaving, solid *entering){
			if (!hitfile){
				ostringstream filename;
				filename << outpath << '/' << setw(12) << setfill('0') << jobnumber << name << "hit.out";
				cout << "Creating " << filename.str() << '\n';
				hitfile = new ofstream(filename.str().c_str());
				if (!hitfile || !hitfile->is_open()){
					cout << "Could not create" << filename.str() << '\n';
					exit(-1);
				}
				*hitfile << "jobnumber particle "
							"t x y z v1x v1y v1z pol1 "
							"v2x v2y v2z pol2 "
							"nx ny nz solid1 solid2\n";
				hitfile->precision(10);
			}

			cout << ":";
			*hitfile << jobnumber << " " << particlenumber << " "
					<< x << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << pol1 << " "
					<< y2[3] << " " << y2[4] << " " << y2[5] << " " << pol2 << " "
					<< normal[0] << " " << normal[1] << " " << normal[2] << " " << leaving->ID << " " << entering->ID << '\n';
		}


		/**
		 * Calculate kinetic energy.
		 *
		 * Kinetic energy is calculated by series expansion of rel. gamma factor, if v/c < ::RELATIVSTIC_THRESHOLD, else it is calculated exactly by (gamma-1)mc^2
		 *
		 * @return Kinetic energy [eV]
		 */
		long double Ekin(const long double v[3]){
			long double vend = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			long double beta = vend/c_0;
			if (beta < RELATIVISTIC_THRESHOLD) // use series expansion for energy calculation with small beta
				return 0.5*m*vend*vend + (3/8*m + 5/16*m*beta*beta + 35/128*beta*beta*beta*beta)*vend*vend*beta*beta;
			else{
				long double gammarel = 1/sqrt(1 - beta*beta); // use relativistic formula for larger beta
				return c_0*c_0*m*(gammarel - 1);
			}
		}

		/**
		 * Calculate potential energy of particle
		 *
		 * @param t Time
		 * @param y Coordinate vector
		 * @param polarisation Particle polarisation
		 * @param field Pointer to TFieldManager structure for electric and magnetic potential
		 * @param solids List of solids for optional material potential
		 *
		 * @return Returns potential energy [eV]
		 */
		virtual long double Epot(const long double t, const long double y[3], int polarisation, TFieldManager *field, map<solid, bool> &solids){
			long double result = 0;
			if ((q != 0 || mu != 0) && field){
				long double B[4][4], E[3], V;
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
		};

};



#endif // PARTICLE_H_
