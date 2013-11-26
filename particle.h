/**
 * \file
 * Particle base class definition.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_


#include <cmath>
#include <vector>
#include <algorithm>

#define MAX_SAMPLE_DIST 0.01 ///< max spatial distance of reflection checks, spin flip calculation, etc; longer integration steps will be interpolated
#define MIN_SAMPLE_DIST 0.005 ///< min spatial distance of track-entries

#include "globals.h"
#include "fields.h"
#include "mc.h"
#include "geometry.h"
#include "bruteforce.h"
#include "ndist.h"
#include "adiabacity.h"


/**
 * Track points are saved in an array of this structure
 */
struct TTrackEntry{
	long double t; ///< time
	VecDoub y; ///< state vector (x,y,z,vx,vy,vz)
	int polarisation; ///< polarisation
};


/**
 * Material boundary hits are saved in an array of this structure
 */
struct THitEntry{
	vector<TTrackEntry>::iterator point1; ///< iterator pointing to track point in front of material boundary
	vector<TTrackEntry>::iterator point2; ///< iterator pointing to track point beyond material boundary
	solid leaving; ///< solid, which the particle left at this material boundary
	solid entering; ///< solid, which the particle entered at this material boundary
	long double normal[3]; ///< normal of the hit surface
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
		int ID; ///< particle fate (0: not classified; >0: Absorption in solid with according ID; -1: survived until end of simulation; -2: hit outer boundaries; -3: numerical error; -4: decayed; -5: did not find initial position)
		int particlenumber; ///< particle number
		long double tau; ///< particle life time
		long double maxtraj; ///< max. simulated trajectory length
		long double inttime; ///< integration computing time
		long double refltime; ///< reflection computing time

		/// start time
		long double tstart(){ return track.begin()->t; };

		/// stop time
		long double tend(){ return (track.end() - 1)->t; };

		/// state vector before integration
		VecDoub ystart(){ return track.begin()->y; };

		/// state vector after integration)
		VecDoub yend(){ return (track.end() - 1)->y; };

		/// total start energy
		long double Hstart(){ return Estart() + Epot(tstart(), &ystart()[0], polstart(), field); };

		/// total end energy
		long double Hend(){ return Eend() + Epot(tend(), &yend()[0], polend(), field); };

		/// max total energy
		long double Hmax(){
			long double result = -numeric_limits<long double>::infinity();
			for (vector<TTrackEntry>::iterator i = track.begin(); i != track.end(); i++){
				long double H = Ekin(&i->y[3]) + Epot(i->t, &i->y[0], i->polarisation, field);
				result = max(result, H);
			}
			return result;
		};

		/// kinetic start energy
		long double Estart(){ return Ekin(&track.begin()->y[3]); };

		/// kinetic end energy
		long double Eend(){ return Ekin(&(track.end() - 1)->y[3]); };

		/// trajectory length
		long double lend(){
			long double result = 0;
			for (vector<TTrackEntry>::iterator i = ++track.begin(); i != track.end(); i++){
				result += sqrt(pow(i->y[0] - (i-1)->y[0], 2) + pow(i->y[1] - (i-1)->y[1], 2) + pow(i->y[2] - (i-1)->y[2], 2));
			}
			return result;
		};

		/// initial polarisation of particle (-1,0,1)
		int polstart(){ return track.begin()->polarisation; };

		/// final polarisation of particle (-1,0,1)
		int polend(){ return (track.end() - 1)->polarisation; };

		/// number of material boundary hits
		int Nhit(){ return hits.size(); };

		/// number of spin flips
		int Nspinflip(){
			int result = 0;
			for (vector<TTrackEntry>::iterator i = ++track.begin(); i != track.end(); i++){
				if (i->polarisation != (i-1)->polarisation)
					result++;
			}
			return result;
		};

		/// number of integration steps
		int Nsteps(){ return track.size(); };

		vector<TParticle*> secondaries; ///< list of secondary particles
		vector<TTrackEntry> track; ///< list of track entries
		vector<THitEntry> hits; ///< list of material boundary hits
		vector<vector<TTrackEntry>::iterator> snapshots;


		/**
		 * Constructor, initializes TParticle::type, TParticle::q, TParticle::m, TParticle::mu
		 *
		 * Has to be called by every derived class constructor with the respective values
		 */
		TParticle(const char *aname, const long double qq, const long double mm, const long double mumu)
				: name(aname), q(qq), m(mm), mu(mumu){ };

		/**
		 * Destructor, deletes secondaries and snapshots
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
		 * On each step it checks for interaction with solids, saves snapshots and track into vectors.
		 * The Integration stops if TParticle::tau or tmax are reached; or if something happens to the particle (absorption, numerical error)
		 *
		 * @param tmax Max. absolute time at which integration will be stopped
		 * @param geometry Integrator checks for collisions with walls defined in this TGeometry structure
		 * @param mcgen Random number generator (e.g. used for reflection probability dicing)
		 * @param afield Pointer to field which acts on particle (can be NULL)
		 * @param snapshottime List of times at which the track should be saved
		 */
		void Integrate(long double tmax, TGeometry &geometry, TMCGenerator &mcgen, TFieldManager *afield, set<float> &snapshottimes){
			geom = &geometry;
			if (currentsolids.empty())
				geom->GetSolids(tend(), &yend()[0], currentsolids);
			mc = &mcgen;
			field = afield;
			timespec clock_start, clock_end, refl_start, refl_end;

			int perc = 0;
			cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
			cout << "x: " << yend()[0] << " y: " << yend()[1] << " z: " << yend()[2]
			     << " E: " << Eend() << " t: " << tend() << " tau: " << tau << " lmax: " << maxtraj << '\n';
		
			// set initial values for integrator
			long double x = tend(), x1, x2;
			VecDoub y(yend()), dydx(6), y1(6), y2(6);
			long double l = lend();
			int polarisation = polend();
			long double h = 0.001/sqrt(yend()[3]*yend()[3] + yend()[4]*yend()[4] + yend()[5]*yend()[5]); // first guess for stepsize
			derivs(x,y,dydx);

			// create integrator class (stepperdopr853 = 8th order Runge Kutta)
			stepper = new StepperDopr853<TParticle>(y, dydx, x, 1e-13, 0, true);// y, dydx, x, atol, rtol, dense output

			set<float>::iterator nextsnapshot = snapshottimes.begin();
			while (nextsnapshot != snapshottimes.end() && *nextsnapshot < x) // find first snapshot time
				nextsnapshot++;
			long double lastsave = x;

			while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
				x1 = x; // save point before next step
				y1 = y;
				
				if (x + h > tstart() + tau)
					h = tstart() + tau - x;	//If stepsize can overshoot, decrease.
				if (x + h > tmax)
					h = tmax - x;
				
				clock_gettime(CLOCK_REALTIME, &clock_start); // start computing time measure
				try{
					stepper->step(h, *this); // integrate one step
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
						y2 = y;
					}
					else{
						for (int i = 0; i < 6; i++)
							y2[i] = stepper->dense_out(i, x2, stepper->hdid); // interpolate y2
					}
					
					clock_gettime(CLOCK_REALTIME, &refl_start);
					if (CheckHit(x1, y1, x2, y2, polarisation)){ // check if particle hit a material boundary or was absorbed between y1 and y2
						x = x2; // if particle path was changed: reset integration end point
						y = y2;
					}
					clock_gettime(CLOCK_REALTIME, &refl_end);
					refltime += refl_end.tv_sec - refl_start.tv_sec + (long double)(refl_end.tv_nsec - refl_start.tv_nsec)/1e9;
				
					l += sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));

					// take snapshots at certain times
					if (nextsnapshot != snapshottimes.end()){
						long double xsnap = *nextsnapshot;
						if (x1 <= xsnap && x2 > xsnap){
							VecDoub ysnap(6);
							for (int i = 0; i < 6; i++)
								ysnap[i] = stepper->dense_out(i, xsnap, stepper->hdid);
							cout << "\n Snapshot at " << xsnap << " s \n";

							TTrackEntry point = {xsnap, ysnap, polarisation};
							track.push_back(point);
							snapshots.push_back(track.end() - 1);
							nextsnapshot++;
						}
					}

					TTrackEntry trackentry = {x2, y2, polarisation};
					track.push_back(trackentry);

					x1 = x2;
					y1 = y2;
				}		
				
				PrintPercent(max((x - tstart())/tau, max((x - tstart())/(tmax - tstart()), l/maxtraj)), perc);

				// x >= tstart + tau?
				if (ID == ID_UNKNOWN && x >= tstart() + tau)
					ID = ID_DECAYED;
				else if (ID == ID_UNKNOWN && (x >= tmax || l >= maxtraj))
					ID = ID_NOT_FINISH;
			
				h = stepper->hnext; // set suggested stepsize for next step
			}
			
			delete stepper;
			
			if (ID == ID_DECAYED){ // if particle reached its lifetime call TParticle::Decay
				cout << "Decayed!\n";
				Decay();
			}

			cout << "Done!!\n";
			cout << "x: " << yend()[0];
			cout << " y: " << yend()[1];
			cout << " z: " << yend()[2];
			cout << " E: " << Eend();
			cout << " Code: " << ID;
			cout << " t: " << tend();
			cout << " l: " << lend();
			cout << " hits: " << Nhit();
			cout << " steps: " << Nsteps() << "\n";
			cout.flush();
		};
	

		/**
		 * Print *start and *end values to a stream.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param endfile stream to print into
		 */
		virtual void PrintStartEnd(ostream &endfile){
			if (endfile.tellp() == 0){
				endfile <<	"jobnumber particle "
	                		"tstart xstart ystart zstart "
	                		"vxstart vystart vzstart "
	                		"polstart Hstart Estart "
	                		"tend xend yend zend "
	                		"vxend vyend vzend "
	                		"polend Hend Eend stopID Nspinflip ComputingTime "
	                		"Nhit Nstep trajlength Hmax\n";
			}
			cout << "Printing status\n";
			endfile	<<	jobnumber << " " << particlenumber << " "
					<<	tstart() << " " << ystart()[0] << " " << ystart()[1] << " " << ystart()[2] << " "
					<<	ystart()[3] << " " << ystart()[4] << " " << ystart()[5] << " "
					<< polstart() << " " <<	Hstart() << " " << Estart() << " "
					<< tend() << " " << yend()[0] << " " << yend()[1] << " " << yend()[2] << " "
					<< yend()[3] << " " << yend()[4] << " " << yend()[5] << " "
					<< polend() << " " << Hend() << " " << Eend() << " " << ID << " " << Nspinflip() << " " << inttime << " "
					<< Nhit() << " " << Nsteps() << " " << lend() << " " << Hmax() << '\n';
		};


		/**
		 * Print snapshots to a stream.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param snapfile stream to print into
		 * @param snapshots list of snapshot times
		 */
		virtual void PrintSnapshots(ostream &snapfile){
			if (snapfile.tellp() == 0){
				snapfile <<	"jobnumber particle polarisation "
	                		"tstart xstart ystart zstart "
	                		"vxstart vystart vzstart "
	                		"Hstart Estart "
	                		"tend xend yend zend "
	                		"vxend vyend vzend "
	                		"Hend Eend stopID Nspinflip "
	                		"Nhit Nstep trajlength Hmax\n";
			}

			vector<TTrackEntry>::iterator te = track.begin();
			vector<THitEntry>::iterator he = hits.begin();
			vector<vector<TTrackEntry>::iterator>::iterator ss = snapshots.begin();
			VecDoub y;
			int Nh = 0, Ns = 0, Nf = 0;
			long double H, E, l = 0, maxH = -numeric_limits<long double>::infinity();
			while (ss != snapshots.end() && te != track.end()){
				Ns++;
				y = te->y;
				E = Ekin(&y[3]);
				H = E + Epot(te->t, &y[0], te->polarisation, field);
				if (H > maxH) maxH = H;
				if (te != track.begin()){
					if (te->polarisation != (te + 1)->polarisation)
						Nf++;
					l += sqrt(pow(te->y[0] - (te-1)->y[0], 2) + pow(te->y[1] - (te-1)->y[1], 2) + pow(te->y[2] - (te-1)->y[2], 2));
				}

				if (he->point2->t < (*ss)->t && he != hits.end()){
					Nh++;
					he++;
				}
				if ((*ss)->t == te->t){
					cout << "Printing snapshot at t = " << (*ss)->t << " s\n";
					snapfile<< jobnumber << " " << particlenumber << " " << te->polarisation << " "
							<< tstart() << " " << ystart()[0] << " " << ystart()[1] << " " << ystart()[2] << " "
							<< ystart()[3] << " " << ystart()[4] << " " << ystart()[5] << " "
							<< Hstart() << " " << Estart() << " "
							<< te->t << " " << y[0] << " " << y[1] << " " << y[2] << " "
							<< y[3] << " " << y[4] << " " << y[5] << " "
							<< H << " " << E << " " << ID << " " << Nf << " "
							<< Nh << " " << Ns << " " << l << " " << maxH << '\n';
					ss++;
				}
				te++;
			}

		};


		/**
		 * Print particle track into stream to allow visualization of the particle's trajectory.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param trackfile stream to print into
		 * @param minPointDist if points are closer than this distance, only one is printed
		 */
		virtual void PrintTrack(ostream &trackfile, long double minPointDist = 0){
			if (trackfile.tellp() == 0){
				trackfile << 	"particle polarisation "
								"t x y z vx vy vz "
								"H E Bx dBxdx dBxdy dBxdz By dBydx "
								"dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V\n";
			}

			int percent = 0;
			cout << "Printing track ";
			for (vector<TTrackEntry>::iterator i = track.begin(); i != track.end(); i++){
				long double x = i->t;
				VecDoub y = i->y;
				if (i != track.begin() && minPointDist > 0){
					VecDoub yprev = (i - 1)->y;
					if (pow(y[0] - yprev[0], 2) + pow(y[1] - yprev[1], 2) + pow(y[2] - yprev[2], 2) < minPointDist*minPointDist)
						continue;
				}
				long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
				long double E[3] = {0,0,0};
				long double V = 0;
				if (field){
					field->BField(y[0],y[1],y[2],x,B);
					field->EField(y[0],y[1],y[2],x,V,E);
				}
				long double Ek = Ekin(&y[3]);
				long double H = Ek + Epot(x, &y[0], i->polarisation, field);

				trackfile 	<< particlenumber << " " << i->polarisation << " "
							<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
							<< H << " " << Ek << " ";
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++)
						trackfile << B[i][j] << " ";
				trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';

				PrintPercent((x - tstart())/(tend() - tstart()), percent);
			}
			cout << '\n';
		};

		/**
		 * Print material boundary hits into stream.
		 *
		 * This is a virtual function and can be overwritten by derived particle classes.
		 *
		 * @param hitfile stream to print into
		 */
		virtual void PrintHits(ostream &hitfile){
			if (hitfile.tellp() == 0){
				hitfile << "jobnumber particle "
							"t x y z v1x v1y v1z pol1"
							"v2x v2y v2z pol2 "
							"nx ny nz solid1 solid2\n";
			}

			cout << "Printing hits ...\n";
			for (vector<THitEntry>::iterator i = hits.begin(); i != hits.end(); i++){
				cout << i->entering.ID;
				cout.flush();
				long double x1 = i->point1->t, x2 = i->point2->t;
				VecDoub y1 = i->point1->y, y2 = i->point2->y;
				hitfile << jobnumber << " " << particlenumber << " "
						<< x1 << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << i->point1->polarisation << " "
						<< y2[3] << " " << y2[4] << " " << y2[5] << " " << i->point2->polarisation << " "
						<< i->normal[0] << " " << i->normal[1] << " " << i->normal[2] << " " << i->leaving.ID << " " << i->entering.ID << '\n';
			}
		}

	protected:
		std::set<solid> currentsolids; ///< solids in which particle is currently inside
		TGeometry *geom; ///< TGeometry structure passed by "Integrate"
		TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
		TFieldManager *field; ///< TFieldManager structure passed by "Integrate"
		StepperDopr853<TParticle> *stepper; ///< ODE integrator


		/**
		 * Initialize particle, set start values randomly according to *.in files.
		 *
		 * Sets start time TParticle::tstart according to [SOURCE] in geometry.in.
		 * Sets start energy according to TMCGenerator::Spectrum.
		 * Then tries to find a position according to [SOURCE] in geometry.in.
		 *
		 * @param number Particle number
		 * @param ageometry TGeometry in which the particle will be simulated
		 * @param src TSource in which particle should be generated
		 * @param mcgen TMCGenerator used to dice inital values
		 * @param afield TFieldManager used to calculate energies (can be NULL)
		 */
		void Init(int number, TGeometry &ageometry, TSource &src,
					TMCGenerator &mcgen, TFieldManager *afield){
			long double t = mcgen.UniformDist(0,src.ActiveTime);
			long double E = mcgen.Spectrum(name);
			long double phi, theta;
			long double p[3];
			src.RandomPointInSourceVolume(mcgen, t, E, p[0], p[1], p[2], phi, theta);

			InitE(number, t, mcgen.LifeTime(name), p[0], p[1], p[2],
				E, phi, theta, mcgen.DicePolarisation(name), mcgen.MaxTrajLength(name), ageometry, afield);
		};


		/**
		 * Initialze all values
		 *
		 * Initialize particle using start energy, calls TParticle::InitV
		 *
		 * @param number Particle number
		 * @param t Start time
		 * @param atau Life time
		 * @param x Start x coordinate
		 * @param y Start y coordinate
		 * @param z Start z coordinate
		 * @param Ekin Kinetic start energy
		 * @param phi Azimuth of velocity vector
		 * @param theta Polar angle of velocity vector
		 * @param pol Start polarisation
		 * @param trajl Max. simulated trajectory length
		 * @param ageometry TGeometry in which the particle will be simulated
		 * @param afield Electric and magnetic fields (needed to determine total energy)
		 */
		void InitE(int number, long double t, long double atau, long double x, long double y, long double z,
				long double Ekin, long double phi, long double theta, int pol, long double trajl,
				TGeometry &ageometry, TFieldManager *afield){
			long double gammarel = Ekin/m/c_0/c_0 + 1;
			long double vstart;
			if (gammarel < 1.0001)
				vstart = sqrt(2*Ekin/m);
			else
				vstart = c_0 * sqrt(1-(1/(gammarel*gammarel)));
			InitV(number, t, atau, x, y, z,
					vstart*cos(phi)*sin(theta), vstart*sin(phi)*sin(theta), vstart*cos(theta),
					pol, trajl, ageometry, afield);
		}

		/**
		 * Initialize all values
		 *
		 * This should be called by every constructor
		 *
		 * @param number Particle number
		 * @param t Start time
		 * @param atau Life time
		 * @param x Start x coordinate
		 * @param y Start y coordinate
		 * @param z Start z coordinate
		 * @param vx Velocity x component
		 * @param vy Velocity y component
		 * @param vz Velocity z component
		 * @param pol Start polarisation
		 * @param trajl Max. simulated trajectory length
		 * @param ageometry TGeometry in which the particle will be simulated
		 * @param afield Electric and magnetic fields (needed to determine total energy)
		 */
		void InitV(int number, long double t, long double atau, long double x, long double y, long double z,
				long double vx, long double vy, long double vz, int pol, long double trajl,
				TGeometry &ageometry, TFieldManager *afield){
			ID = ID_UNKNOWN;
			inttime = refltime = 0;
			geom = &ageometry;
			mc = NULL;
			field = afield;
			stepper = NULL;

			particlenumber = number;
			tau = atau;
			maxtraj = trajl;
			long double ys[6] = {x, y, z, vx, vy, vz};
			TTrackEntry trackentry = {t, VecDoub(6, ys), pol}; // save first track point
			track.push_back(trackentry);

			if (currentsolids.empty())
				geom->GetSolids(t, ys, currentsolids);
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
				if (q != 0 || (mu != 0 && polend() != 0))
					field->BField(y[0],y[1],y[2], x, B); // if particle has charge or magnetic moment, calculate magnetic field
				if (q != 0)
					field->EField(y[0],y[1],y[2], x, V, E); // if particle has charge caculate electric field+potential
				if (q != 0){
					F[0] += q*(E[0] + y[4]*B[2][0] - y[5]*B[1][0]); // add Lorentz-force
					F[1] += q*(E[1] + y[5]*B[0][0] - y[3]*B[2][0]);
					F[2] += q*(E[2] + y[3]*B[1][0] - y[4]*B[0][0]);
				}
				if (mu != 0 && polend() != 0){
					F[0] += polend()*mu*B[3][1]; // add force on magnetic dipole moment
					F[1] += polend()*mu*B[3][2];
					F[2] += polend()*mu*B[3][3];
				}
			}
			long double rel = sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0))/m/ele_e; // relativstic factor 1/gamma/m
			dydx[3] = rel*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
			dydx[4] = rel*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
			dydx[5] = rel*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);
		};
		

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
		 * @param iteration Iteration counter (incremented by recursive calls to avoid infinite loop)
		 * @return Returns true if particle was reflected/absorbed
		 */
		bool CheckHit(long double &x1, VecDoub_IO &y1, long double &x2, VecDoub_IO &y2, int &pol, int iteration = 1){
			if (!geom->CheckSegment(&y1[0],&y2[0])){ // check if start point is inside bounding box of the simulation geometry
				printf("\nParticle has hit outer boundaries: Stopping it! t=%LG x=%LG y=%LG z=%LG\n",x2,y2[0],y2[1],y2[2]);
				ID = ID_HIT_BOUNDARIES;
				return true;
			}
			set<TCollision> colls;
			if (geom->GetCollisions(x1, &y1[0], stepper->hdid, &y2[0], colls)){	// if there is a collision with a wall
				TCollision coll = *colls.begin();
				long double u[3] = {y2[0] - y1[0], y2[1] - y1[1], y2[2] - y1[2]};
				long double distnormal = u[0]*coll.normal[0] + u[1]*coll.normal[1] + u[2]*coll.normal[2];
				if ((colls.size() == 1
					&& coll.s*abs(distnormal) < REFLECT_TOLERANCE
					&& (1 - coll.s)*abs(distnormal) < REFLECT_TOLERANCE) // if there is only one collision which is closer to y1 and y2 than REFLECT_TOLERANCE
					|| iteration > 99) // or if iteration counter reached a certain maximum value
				{
					TTrackEntry track1 = {x1,y1,pol};
					track.push_back(track1);
					solid *hitsolid = &geom->solids[coll.ID];
					if (colls.size() > 1 && coll.ID != (++colls.begin())->ID){
						// if particle hits two or more different surfaces at once, something went wrong
						printf("Reflection from two surfaces (%s & %s) at once!\n", hitsolid->name.c_str(), geom->solids[(++colls.begin())->ID].name.c_str());
						printf("Check geometry for touching surfaces!");
						x2 = x1;
						y2 = y1;
						ID = ID_NRERROR;
						return true;
					}
//					cout << "Hit " << sld->ID << " - current solid " << currentsolids.rbegin()->ID << '\n';
					const solid *leaving, *entering;
					bool resetintegration = false, hit = false;
					if (distnormal < 0){ // particle is entering solid
						if (currentsolids.find(*hitsolid) != currentsolids.end()){ // if solid already is in currentsolids list something went wrong
							printf("Particle entering '%s' which it did enter before! Stopping it! Did you define overlapping solids with equal priorities?\n", geom->solids[coll.ID].name.c_str());
							x2 = x1;
							y2 = y1;
							ID = ID_NRERROR;
							return true;
						}
						leaving = &*currentsolids.rbegin();
						entering = hitsolid;
						if (entering->ID > leaving->ID){ // check for reflection only, if priority of entered solid is larger than that of current solid
							hit = true;
							resetintegration = OnHit(x1, y1, x2, y2, coll.normal, leaving, entering); // do particle specific things
							if (entering != leaving)
								currentsolids.insert(*entering);
							entering = hitsolid;
						}
					}
					else if (distnormal > 0){ // particle is leaving solid
						if (currentsolids.find(*hitsolid) == currentsolids.end()){ // if solid is not in currentsolids list something went wrong
							cout << "Particle inside '" << hitsolid->name << "' which it did not enter before! Stopping it!\n";
							x2 = x1;
							y2 = y1;
							ID = ID_NRERROR;
							return true;
						}
						leaving = hitsolid;
						if (leaving->ID == currentsolids.rbegin()->ID){ // check for reflection only, if left solids is the solid with highest priority
							hit = true;
							entering = &*++currentsolids.rbegin();
							resetintegration = OnHit(x1, y1, x2, y2, coll.normal, leaving, entering); // do particle specific things
							if (entering != leaving)
								currentsolids.erase(*leaving);
							else
								entering = &*++currentsolids.rbegin();
						}
					}
					else{
						cout << "Particle is crossing material boundary with parallel track! Stopping it!\n";
						x2 = x1;
						y2 = y1;
						ID = ID_NRERROR;
						return true;
					}

					TTrackEntry track2 = {x2,y2,pol};
					track.push_back(track2);
					if (hit){
						THitEntry hitentry = {--track.end(), ----track.end(), *leaving, *entering, {coll.normal[0], coll.normal[1], coll.normal[2]}};
						hits.push_back(hitentry);
					}
					if (resetintegration)
						return true;
					else if (OnStep(x1, y1, x2, y2, &*currentsolids.rbegin())){ // check for absorption
						printf("Absorption!\n");
						return true;
					}
				}
				else{
					// else cut integration step right before and after first collision point
					// and call ReflectOrAbsorb again for each smaller step (quite similar to bisection algorithm)
					const TCollision *c = &*colls.begin();
					long double xnew, xbisect1 = x1, xbisect2 = x1;
					VecDoub ybisect1 = y1, ybisect2 = y1;

					xnew = x1 + (x2 - x1)*c->s*(1 - 0.01*iteration); // cut integration right before collision point
					if (xnew > x1 && xnew < x2){ // check that new line segment is in correct time interval
						xbisect1 = xbisect2 = xnew;
						for (int i = 0; i < 6; i++)
							ybisect1[i] = ybisect2[i] = stepper->dense_out(i, xbisect1, stepper->hdid); // get point at cut time
						if (CheckHit(x1, y1, xbisect1, ybisect1, pol, iteration+1)){ // recursive call for step before coll. point
							x2 = xbisect1;
							y2 = ybisect1;
							return true;
						}
					}

					xnew = x1 + (x2 - x1)*c->s*(1 + 0.01*iteration); // cut integration right after collision point
					if (xnew > xbisect1 && xnew < x2){ // check that new line segment does not overlap with previous one
						xbisect2 = xnew;
						for (int i = 0; i < 6; i++)
							ybisect2[i] = stepper->dense_out(i, xbisect2, stepper->hdid); // get point at cut time
						if (CheckHit(xbisect1, ybisect1, xbisect2, ybisect2, pol, iteration+1)){ // recursive call for step over coll. point
							x2 = xbisect2;
							y2 = ybisect2;
							return true;
						}
					}

					if (CheckHit(xbisect2, ybisect2, x2, y2, pol, iteration+1)) // recursive call for step after coll. point
						return true;
				}
			}
			else if (OnStep(x1, y1, x2, y2, &*currentsolids.rbegin())){ // if there was no collision: just check for absorption in solid with highest priority
				printf("Absorption!\n");
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
		 * @param normal Normal vector of hit surface
		 * @param leaving Solid that the particle is leaving
		 * @param entering Solid that the particle is entering (can be modified by method)
		 * @return Returns true if particle path was changed
		 */
		virtual bool OnHit(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], const solid *leaving, const solid *&entering) = 0;


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
		virtual bool OnStep(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, const solid *currentsolid) = 0;


		/**
		 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
		 */
		virtual void Decay() = 0;


		/**
		 * Calculate kinetic energy.
		 *
		 * Kinetic energy is calculated by 0.5mv^2, if v/c < ::RELATIVSTIC_THRESHOLD, else it is calculated by (gamma-1)mc^2
		 *
		 * @return Kinetic energy [eV]
		 */
		long double Ekin(const long double v[3]){
			long double vend = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			if (vend/c_0 < RELATIVISTIC_THRESHOLD)
				return 0.5*m*vend*vend;
			else{
				long double gammarel = 1/sqrt(1 - vend*vend/(c_0*c_0));
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
		 *
		 * @return Returns potential energy [eV]
		 */
		virtual long double Epot(const long double t, const long double y[3], int polarisation, TFieldManager *field = NULL){
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
