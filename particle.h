/**
 * \file
 * Particle base class definition.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_


#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/time.h>

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
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor, "Reflect", "Absorb" and "ForEachStep"
 */
struct TParticle{
	public:
		const int type; ///< particle type
		const long double q; ///< charge [C]
		const long double m; ///< mass [eV/c^2]
		const long double mu; ///< magnetic moment [J/T]
		int ID; ///< particle fate (0: not classified; >0: Absorption in solid with according ID; -1: survived until end of simulation; -2: hit outer boundaries; -3: numerical error; -4: decayed; -5: did not find initial position)
		int particlenumber; ///< particle number
		long double tstart; ///< start time
		long double tend; ///< stop time
		long double tau; ///< particle life time
		long double ystart[6]; ///< state vector before integration
		long double yend[6]; ///< state vector after integration
		long double Hstart; ///< total start energy
		long double Hend; ///< total end energy
		long double Hmax; ///< max total energy
		long double Estart; ///< kinetic start energy
		long double Eend; ///< kinetic end energy
		long double trajlength; ///< trajectory length
		long double maxtraj; ///< max. simulated trajectory length
		long double comptime; ///< computing time
		long double refltime; ///< reflection computing time
		int polarisation; ///< polarisation of particle (-1,0,1)
		int Nhit; ///< number of material boundary hits
		int Nspinflip; ///< number of spin flips
		int Nsteps; ///< number of integration steps
		std::vector<TParticle*> secondaries; ///< list of secondary particles

		/**
		 * Constructor, initializes TParticle::type, TParticle::q, TParticle::m, TParticle::mu
		 *
		 * Has to be called by every derived class constructor with the respective values
		 */
		TParticle(int pp, long double qq, long double mm, long double mumu): type(pp), q(qq), m(mm), mu(mumu){ };

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
		 * On each step it checks for reflection or absorption by solids, prints snapshots and track into files and calls "OnEachStep".
		 * The Integration stops if xend or SimulationTime are reached; or if something happens to the particle (absorption, numerical error)
		 *
		 * @param tmax Max. absolute time at which integration will be stopped
		 * @param geometry Integrator checks for collisions with walls defined in this TGeometry structure
		 * @param mcgen Random number generator (e.g. used for reflection probability dicing)
		 * @param afield Pointer to field which acts on particle (can be NULL)
		 * @param ENDLOG FILE-pointer into which particle is logged after integration
		 * @param trackfile FILE-pointer into which particle trajectory is logged
		 * @param SNAP FILE-pointer into which snapshots are logged
		 * @param snapshots Vector containing snapshot times
		 * @param areflectlog Should reflections (1) or transmissions (2) or both (3) be logged?
		 * @param aREFLECTLOG FILE-pointer into which reflection and transmission are logged
		 */
		void Integrate(long double tmax, TGeometry &geometry, TMCGenerator &mcgen, TField *afield, FILE *ENDLOG, FILE *trackfile,
						FILE *SNAP = NULL, vector<float> *snapshots = NULL, int areflectlog = 0, FILE *aREFLECTLOG = NULL){
			geom = &geometry;
			if (currentsolids.empty())
				geom->GetSolids(tstart, ystart, currentsolids);
			mc = &mcgen;
			field = afield;
			reflectlog = areflectlog;
			REFLECTLOG = aREFLECTLOG;
			timeval clock_start, clock_end, refl_start, refl_end;
			gettimeofday(&clock_start, NULL); // start computing time measure
			unsigned int nextsnapshot = 0;
			while (snapshots && nextsnapshot < snapshots->size() && (*snapshots)[nextsnapshot] < tstart)
				nextsnapshot++;
			int perc = 0;
			printf("Particle no.: %i, particle type: %i\n",particlenumber,type);
			printf("x: %LG y: %LG z: %LG E: %LG t: %LG tau: %LG lmax: %LG\n",
				ystart[0], ystart[1], ystart[2], Estart, tstart, tau, maxtraj);
		
			// set initial values for integrator
			long double x = tstart, x1, x2;
			VecDoub y(6,ystart), dydx(6), y1(6), y2(6);
			long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
			long double E[3] = {0,0,0};
			long double V = 0;
			long double h = 0.001/sqrt(ystart[3]*ystart[3] + ystart[4]*ystart[4] + ystart[5]*ystart[5]); // first guess for stepsize
			derivs(x,y,dydx);
			if (trackfile)
				PrintIntegrationStep(trackfile, x, y, h); // print start into track file
			long double lastsave = x;
			// create integrator class (stepperdopr853 = 8th order Runge Kutta)
			stepper = new StepperDopr853<TParticle>(y, dydx, x, 1e-13, 0, true);// y, dydx, x, atol, rtol, dense output
			
			while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
				x1 = x; // save point before next step
				y1 = y;
				
				if (x + h > tstart + tau)
					h = tstart + tend - x;	//If stepsize can overshoot, decrease.
				if (x + h > tmax)
					h = tmax - x;
				
				try{
					stepper->step(h, *this); // integrate one step
					Nsteps++;
				}
				catch(...){ // catch Exceptions thrown by numerical recipes routines
					ID = ID_NRERROR;
				}
				
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
					
					gettimeofday(&refl_start, NULL);
					if (ReflectOrAbsorb(x1, y1, x2, y2)){ // check for reflection or absorption between y1 and y2
						x = x2; // if there was reflection or absorption: reset integration end point
						y = y2;
					}
					gettimeofday(&refl_end, NULL);
					refltime += refl_end.tv_sec - refl_start.tv_sec + (long double)(refl_end.tv_usec - refl_start.tv_usec)/1e6;
				
					trajlength += sqrt(pow(y2[0] - y1[0],2) + pow(y2[1] - y1[1],2) + pow(y2[2] - y1[2],2)); // Trajectory length calculation		

					// take snapshots at certain times
					if(SNAP && snapshots && nextsnapshot < snapshots->size()){
						long double xsnap = (*snapshots)[nextsnapshot];
						if (x1 <= xsnap && x2 > xsnap){
							VecDoub ysnap(6);
							for (int i = 0; i < 6; i++)
								ysnap[i] = stepper->dense_out(i, xsnap, stepper->hdid);
//							printf("\n Snapshot at %LG s \n", xsnap);
							SetEndValues(xsnap,ysnap,field);
							Print(SNAP);
							nextsnapshot++;
						}
					}

					if (trackfile && x2 - lastsave > MIN_SAMPLE_DIST/v1){
						PrintIntegrationStep(trackfile, x2, y2, stepper->hdid); // print integration step in track file
						lastsave = x1;
					}

					x1 = x2;
					y1 = y2;
				}		
				
				if (trajlength/maxtraj < (x - tstart)/tau)
					percent(x, tstart, tstart + tau, perc); // print percent of calculation
				else
					percent(trajlength, 0, maxtraj, perc); // print percent of calculation

				// x >= xstart + xend?
				if (ID == ID_UNKNOWN && x >= tstart + tau)
					ID = ID_DECAYED;
				else if (ID == ID_UNKNOWN && (x >= tmax || trajlength >= maxtraj))
					ID = ID_NOT_FINISH;
			
				h = stepper->hnext; // set suggested stepsize for next step
			}
			
			gettimeofday(&clock_end, NULL);
			comptime = clock_end.tv_sec - clock_start.tv_sec + (long double)(clock_end.tv_usec - clock_start.tv_usec)/1e6;
			delete stepper;
			
			// Endwerte schreiben
			SetEndValues(x,y,field);
			Print(ENDLOG);

			if (ID == ID_DECAYED)
				Decay();

			printf("Done!!\nx: %LG y: %LG z: %LG E: %LG Code: %i t: %LG l: %LG hits: %i steps: %i\n\n",
				yend[0], yend[1], yend[2], Eend, ID, tend, trajlength, Nhit, Nsteps);
		};
	
	protected:
		std::set<solid> currentsolids; ///< solids in which particle is currently inside
		TGeometry *geom; ///< TGeometry structure passed by "Integrate"
		TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
		TField *field; ///< TField structure passed by "Integrate"
		StepperDopr853<TParticle> *stepper; ///< ODE integrator
		int reflectlog; ///< Should reflections (1) or transmissions (2) or both (3) be logged? Passed by "Integrate"
		FILE *REFLECTLOG; ///< reflection log file passed by "Integrate"


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
		 * @param afield TField used to calculate energies (can be NULL)
		 */
		void Init(int number, TGeometry &ageometry, TSource &src,
					TMCGenerator &mcgen, TField *afield){
			tstart = mcgen.UniformDist(0,src.ActiveTime);
			polarisation = mcgen.DicePolarisation(type);
			Estart = mcgen.Spectrum(type);
			long double phi, theta;
			src.RandomPointInSourceVolume(mcgen, tstart, Estart, ystart[0], ystart[1], ystart[2], phi, theta);

			InitE(number, tstart, mcgen.LifeTime(type), ystart[0], ystart[1], ystart[2],
				Estart, phi, theta, polarisation, mcgen.MaxTrajLength(type), ageometry, afield);
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
				TGeometry &ageometry, TField *afield){
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
		 * Initialze all values
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
				TGeometry &ageometry, TField *afield){
			ID = ID_UNKNOWN;
			trajlength = comptime = refltime = 0;
			Nhit = Nspinflip = Nsteps = 0;
			geom = &ageometry;
			mc = NULL;
			field = afield;
			stepper = NULL;
			REFLECTLOG = NULL;

			particlenumber = number;
			polarisation = pol;
			tstart = t;
			tau = atau;
			maxtraj = trajl;
			ystart[0] = yend[0] = x;
			ystart[1] = yend[1] = y;
			ystart[2] = yend[2] = z;
			ystart[3] = yend[3] = vx;
			ystart[4] = yend[4] = vy;
			ystart[5] = yend[5] = vz;

			if (currentsolids.empty())
				geom->GetSolids(tstart, ystart, currentsolids);
			long double vstart = sqrt(vx*vx + vy*vy + vz*vz);
			long double gammarel = 1/sqrt(1 - vstart*vstart/c_0/c_0);
			if (gammarel < 1.0001)
				Estart = 0.5*m*vstart*vstart;
			else
				Estart = (gammarel - 1)*m*c_0*c_0;
			Estart = Eend;
			Hstart = Hend = Estart + Epot(t, ystart, afield);
			Hmax = Hstart;
		}

		/**
		 * Equations of motion dy/dx = f(x,y).
		 *
		 * Equations of motion (relativistic and non-relativstic), called by operator().
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
				if (q != 0 || (mu != 0 && polarisation != 0))
					field->BFeld(y[0],y[1],y[2], x, B); // if particle has charge or magnetic moment, calculate magnetic field
				if (q != 0)
					field->EFeld(y[0],y[1],y[2], V, E); // if particle has charge caculate electric field+potential
				if (q != 0){
					F[0] += q*(E[0] + y[4]*B[2][0] - y[5]*B[1][0]); // add Lorentz-force
					F[1] += q*(E[1] + y[5]*B[0][0] - y[3]*B[2][0]);
					F[2] += q*(E[2] + y[3]*B[1][0] - y[4]*B[0][0]);
				}
				if (mu != 0 && polarisation != 0){
					F[0] += polarisation*mu*B[3][1]; // add force on magnetic dipole moment
					F[1] += polarisation*mu*B[3][2];
					F[2] += polarisation*mu*B[3][3];
				}
			}
			long double rel = sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0))/m/ele_e; // relativstic factor 1/gamma/m
			dydx[3] = rel*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
			dydx[4] = rel*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
			dydx[5] = rel*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);
		};
		

		/**
		 * Check for reflection on surfaces or absorption in materials.
		 *
		 * Checks if a particle which flies from y1 to y2 in time x2-x1 hits a surface and is reflected or absorbed inside a material.
		 * If a surface is hit the routine splits the line segment y1->y2 on both sides of the collision point and calls itself recursively
		 * with the three new line segments as parameters. This is repeated until both split points are nearer than REFLECTION_TOLERANCE
		 * to the collision point (much like an bisection algorithm). For each line segment "Absorb" is called to check for absorption.
		 * For each short segment crossing a collision point "Reflect" is called to check for reflection.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment
		 * @param y2 End point of line segment
		 * @param iteration Iteration counter (incremented by recursive calls to avoid infinite loop)
		 * @return Returns true if particle was reflected/absorbed
		 */
		bool ReflectOrAbsorb(long double &x1, VecDoub_IO &y1, long double &x2, VecDoub_IO &y2, int iteration = 1){
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
					solid *sld = &geom->solids[coll.ID];
					if (colls.size() > 1 && coll.ID != (++colls.begin())->ID){
						// if particle hits two or more different surfaces at once
						printf("Reflection from two surfaces (%s & %s) at once!\n", sld->name.c_str(), geom->solids[(++colls.begin())->ID].name.c_str());
						printf("Check geometry for touching surfaces!");
						x2 = x1;
						y2 = y1;
						ID = ID_NRERROR;
						return true;
					}
//					cout << "Hit " << sld->ID << " - current solid " << currentsolids.rbegin()->ID << '\n';
					set<solid>::iterator currentsolid = currentsolids.find(*sld);
					Nhit++;
					bool resetintegration = false;
					if (distnormal < 0){ // particle is entering solid
						if (currentsolid != currentsolids.end()){ // if solid already is in currentsolids list something went wrong
							printf("Particle entering '%s' which it did enter before! Stopping it! Did you define overlapping solids with equal priorities?\n", geom->solids[coll.ID].name.c_str());
							x2 = x1;
							y2 = y1;
							ID = ID_NRERROR;
							return true;
						}
						if (sld->ID > currentsolids.rbegin()->ID){ // check for reflection if priority of new solid is higher than previously entered solids
							if (Reflect(x1, y1, x2, y2, coll.normal, &*currentsolids.rbegin(), sld, resetintegration))
								return true; // particle was reflected
						}
						currentsolids.insert(*sld); // when particle was not reflected insert new solid into currentsolids list
					}
					else if (distnormal > 0){ // particle is leaving solid
						if (currentsolid == currentsolids.end()){ // if solid is not in currentsolids list something went wrong
							printf("Particle inside '%s' which it did not enter before! Stopping it!\n", geom->solids[coll.ID].name.c_str());
							x2 = x1;
							y2 = y1;
							ID = ID_NRERROR;
							return true;
						}
						if (currentsolid->ID == currentsolids.rbegin()->ID){ // check only for reflection if hit wall belongs to solid with highest priority
							if (Reflect(x1, y1, x2, y2, coll.normal, sld, &*++currentsolids.rbegin(), resetintegration))
								return true; // particle was reflected
						}
						currentsolids.erase(*currentsolid); // if particle was not reflected, remove solid from currentsolids list
					}

					if (resetintegration) // if y2 was changed after reflection stop here and restart integration at x2
						return true;

					if (Absorb(x1, y1, x2, y2, &*currentsolids.rbegin())){ // check for absorption
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
						if (ReflectOrAbsorb(x1, y1, xbisect1, ybisect1, iteration+1)){ // recursive call for step before coll. point
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
						if (ReflectOrAbsorb(xbisect1, ybisect1, xbisect2, ybisect2, iteration+1)){ // recursive call for step over coll. point
							x2 = xbisect2;
							y2 = ybisect2;
							return true;
						}
					}

					if (ReflectOrAbsorb(xbisect2, ybisect2, x2, y2, iteration+1)) // recursive call for step after coll. point
						return true;
				}
			}
			else if (Absorb(x1, y1, x2, y2, &*currentsolids.rbegin())){ // if there was no collision: just check for absorption in solid with highest priority
				printf("Absorption!\n");
				return true;
			}
			return false;
		};
	

		/**
		 * Virtual routine to check for reflection on surfaces, has to be implemented separately for each derived particle.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment, may be altered
		 * @param y2 End point of line segment, may be altered
		 * @param normal Normal vector of hit surface
		 * @param leaving Solid that the particle is leaving
		 * @param entering Solid that the particle is entering
		 * @param resetintegration Tell integrator to restart integration at x2 because e.g. y2 was changed
		 * @return Returns true if particle was reflected
		 */
		virtual bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], const solid *leaving, const solid *entering, bool &resetintegration) = 0;


		/**
		 * Virtual routine to check for absorption in solids, has to be implemented separately for each derived particle.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment, may be altered
		 * @param y2 End point of line segment, may be altered
		 * @param currentsolid Solid through which the particle is moving
		 * @return Returns true if particle was absorbed
		 */
		virtual bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, const solid *currentsolid) = 0;


		/**
		 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
		 */
		virtual void Decay() = 0;

		/**
		 * Calculate potential energy of particle
		 *
		 * @param t Time
		 * @param y Coordinate vector
		 * @param field Pointer to TField structure for electric and magnetic potential
		 *
		 * @return Returns potential energy
		 */
		virtual long double Epot(const long double t, const long double y[3], TField *field = NULL){
			if ((q != 0 || mu != 0) &&field){
				long double B[4][4], E[3], V;
				field->BFeld(y[0],y[1],y[2],t,B);
				field->EFeld(y[0],y[1],y[2],V,E);
				return m*gravconst*y[2] + q/ele_e*V - polarisation*mu/ele_e*B[3][0];
			}
			return m*gravconst*y[2];
		};

		/**
		 * Calculates *end variables of particle at time x, position/velocity y and fields B,E,V.
		 *
		 * @param x Time
		 * @param y State vector (position+velocity)
		 * @param field Electric and magnetic fields (needed to determine total energy)
		 */
		void SetEndValues(long double x, VecDoub_I &y, TField *field){
			tend = x;
			for (int i = 0; i < 6; i++)
				yend[i] = y[i];
			long double vend = sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
			if (vend/c_0 < 0.01)
				Eend = 0.5*m*vend*vend;
			else{
				long double gammarel = 1/sqrt(1 - vend*vend/(c_0*c_0));
				Eend = c_0*c_0*m*(gammarel - 1);
			}
			Hend = Eend + Epot(x, &y[0], field);
			if (Hend > Hmax) Hmax = Hend;
		};


		/**
		 * Print *start and *end values to a file.
		 *
		 * @param file FILE-pointer to print into
		 */
		void Print(FILE *file){
			fprintf(file ,"%i %i %i %i "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG "
			               "%.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %i %i %LG "
			               "%i %.20LG %.20LG\n",
			               jobnumber, particlenumber, type, polarisation,
			               tstart, ystart[0], ystart[1], ystart[2],
			               ystart[3], ystart[4], ystart[5],
			               Hstart, Estart,
			               tend, yend[0], yend[1], yend[2],
			               yend[3], yend[4], yend[5],
			               Hend, Eend, ID, Nspinflip, comptime,
			               Nhit, trajlength, Hmax);
			               
			fflush(file);
		};
	

		/**
		 * Print current particle state into file to allow visualization of the particle's trajectory.
		 *
		 * @param trackfile FILE-pointer to print into
		 * @param x Time
		 * @param y State vector (position+velocity)
		 * @param h Last time step of ODE integrator
		 */
		void PrintIntegrationStep(FILE *trackfile, long double x, VecDoub_I &y, long double h){
			long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
			long double E[3] = {0,0,0};
			long double V = 0;
			long double logvlad = 0.0, logfrac = 0.0;
			printf("-");
			if (field){
				field->BFeld(y[0],y[1],y[2],x,B);
				field->EFeld(y[0],y[1],y[2],V,E);
			}
			SetEndValues(x,y,field);
			
			fprintf(trackfile,"%d %d %d "
							 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
							 "%.17LG %.17LG ",
							 particlenumber,type,polarisation,
							 x,y[0],y[1],y[2],y[3],y[4],y[5],
							 Hend,Eend);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					fprintf(trackfile,"%.17LG ",B[i][j]);
			fprintf(trackfile,"%.17LG %.17LG %.17LG %.17LG",E[0],E[1],E[2],V,h);

			fprintf(trackfile,"\n");
			fflush(trackfile);
		};
};



#endif // PARTICLE_H_
