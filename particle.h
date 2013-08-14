/**
 * \file
 * Particle class definitions.
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_


#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/time.h>

using namespace std;

#define MAX_SAMPLE_DIST 0.01 ///< max spatial distance of reflection checks, spin flip calculation, etc; longer integration steps will be interpolated
#define MIN_SAMPLE_DIST 0.005 ///< min spatial distance of track-entries

#include "nr/nr3.h"
#include "nr/stepper.h"
#include "nr/stepperdopr853.h"

#include "globals.h"
#include "fields.h"
#include "mc.h"
#include "geometry.h"
#include "bruteforce.h"
#include "ndist.h"
#include "adiabacity.h"

//#include "f.glueck/eH2.c"

/**
 * Basic particle class (virtual).
 *
 * Includes all functionality (integration, file output).
 * Special particles are derived from this class by declaring its own constructor, "Reflect", "Absorb" and "ForEachStep"
 */
struct TParticle{
	public:
		int protneut; ///< particle type (neutron: 1; proton: 2; electron: 6)
		int ID; ///< particle fate (0: not classified; 1: survived until end of simulation; 2: hit outer boundaries; 3: numerical error; 4: decayed; 5: did not find initial position; >=6: Absorption in solid with according ID)
		int particlenumber; ///< particle number
		const long double q; ///< charge [C]
		const long double m; ///< mass [eV/c^2]
		const long double mu; ///< magnetic moment [J/T]
		long double xstart; ///< start time
		long double xend; ///< max. simulation time (particle life time)
		long double dt; ///< actual simulation time
		long double ystart[6]; ///< state vector before integration
		long double yend[6]; ///< state vector after integration
		long double Hstart; ///< total start energy
		long double Hend; ///< total end energy
		long double Hmax; ///< max total energy
		long double Estart; ///< kinetic start energy
		long double Eend; ///< kinetic end energy
		long double rstart; ///< radial start coordinate
		long double rend; ///< radial end coordinate
		long double phistart; ///< start azimuth
		long double phiend; ///< end azimuth
		long double alphastart; ///< start angle between position vector and velocity vector projected onto xy-plane
		long double alphaend; ///< end angle between position vector and velocity vector projected onto xy-plane
		long double gammastart; ///< start angle between z axis and velocity vector
		long double gammaend; ///< end angle between z axis and velocity vector
		long double vstart; ///< absolute start velocity
		long double vend; ///< absolute end velocity
		long double trajlength; ///< trajectory length
		long double comptime; ///< computing time
		long double refltime; ///< reflection computing time
		long double vlad; ///< current spin flip probabilty by Vladimirskii
		long double vladmax; ///< max spin flip probability by Vladimirskii
		long double vladtotal; ///< product of all spin flip probabilities by Vladimirskii
		long double frac; ///< current adiabacity criterion (<<1 for adiabacity)
		long double thumbmax; ///< max adiabacity criterion
		long double BFsurvprob; ///< no-spinflip-probability determined by Bloch integration
		long double logBF; ///< log10(1-BFsurvprob)
		int polarisation; ///< polarisation of particle (-1,0,1)
		int nrefl; ///< number of reflection from wall
		int NSF; ///< number of spin flips
		int nsteps; ///< number of integration steps
		vector<TParticle*> secondaries; ///< list of secondary particles

		/**
		 * Generic constructor (empty), has to be implemented by each derived particle
		 */
		TParticle(int pp, long double qq, long double mm, long double mumu): protneut(pp), q(qq), m(mm), mu(mumu){ };

		/**
		 * Destructor, deletes secondaries
		 */
		~TParticle(){
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
		void Integrate(TGeometry &geometry, TMCGenerator &mcgen, TField *afield, FILE *ENDLOG, FILE *trackfile,
						FILE *SNAP = NULL, vector<float> *snapshots = NULL, int areflectlog = 0, FILE *aREFLECTLOG = NULL){
			geom = &geometry;
			mc = &mcgen;
			field = afield;
			reflectlog = areflectlog;
			REFLECTLOG = aREFLECTLOG;
			timeval clock_start, clock_end, refl_start, refl_end;
			gettimeofday(&clock_start, NULL); // start computing time measure
			unsigned int nextsnapshot = 0;
			while (snapshots && nextsnapshot < snapshots->size() && (*snapshots)[nextsnapshot] < xstart)
				nextsnapshot++;
			int perc = 0;
			printf("Teilchennummer: %i, Teilchensorte: %i\n",particlenumber,protneut);
			printf("r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG tau: %LG l: %LG\n",
				rstart, phistart/conv, ystart[2], vstart, alphastart/conv,gammastart/conv, Estart, xstart, xend, mc->MaxTrajLength(protneut));
		
			// set initial values for integrator
			long double x = xstart, x1, x2;
			VecDoub y(6,ystart), dydx(6), y1(6), y2(6);
			long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
			long double E[3] = {0,0,0};
			long double V = 0;
			long double h = 0.001/vstart; // first guess for stepsize
			derivs(x,y,dydx);
			if (trackfile)
				PrintIntegrationStep(trackfile, x, y, h); // print start into track file
			long double lastsave = x;
			// create integrator class (stepperdopr853 = 8th order Runge Kutta)
			stepper = new StepperDopr853<TParticle>(y, dydx, x, 1e-13, 0, true);// y, dydx, x, atol, rtol, dense output
			
			while (ID == ID_UNKNOWN){ // integrate as long as nothing happened to particle
				x1 = x; // save point before next step
				y1 = y;
				
				if (x + h > xstart + xend) 
					h = xstart + xend - x;	//If stepsize can overshoot, decrease.
				if (x + h > StorageTime)
					h = StorageTime - x;
				
				try{
					stepper->step(h, *this); // integrate one step
					nsteps++;
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
							printf("\n Snapshot at %LG s \n", xsnap);
							if (field){
								field->BFeld(ysnap[0],ysnap[1],ysnap[2],xsnap,B);
								field->EFeld(ysnap[0],ysnap[1],ysnap[2],V,E);
							}
							SetEndValues(xsnap,ysnap,B,E,V);
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
				
				if (trajlength/mc->MaxTrajLength(protneut) < (x - xstart)/xend)
					percent(x, xstart, xend, perc); // print percent of calculation
				else
					percent(trajlength, 0, mc->MaxTrajLength(protneut), perc); // print percent of calculation

				// x >= xstart + xend?
				if (ID == ID_UNKNOWN && x >= xstart + xend)
					ID = ID_DECAYED;
				else if (ID == ID_UNKNOWN && (x >= StorageTime || trajlength >= mc->MaxTrajLength(protneut)))
					ID = ID_NOT_FINISH;
			
				h = stepper->hnext; // set suggested stepsize for next step
			}
			
			gettimeofday(&clock_end, NULL);
			comptime = clock_end.tv_sec - clock_start.tv_sec + (long double)(clock_end.tv_usec - clock_start.tv_usec)/1e6;
			delete stepper;
			
			// Endwerte schreiben
			if (field){
				field->BFeld(y[0],y[1],y[2],x,B);
				field->EFeld(y[0],y[1],y[2],V,E);
			}
			SetEndValues(x,y,B,E,V);
			Print(ENDLOG);

			if (ID == ID_DECAYED)
				Decay();

			printf("Done!!\nBFFlipProb: %LG r: %LG z: %LG E: %LG Code: %i t: %LG tau: %LG l: %LG\n\n",
				(1 - BFsurvprob), rend, yend[2], Eend, ID, xstart + dt, dt, trajlength);
		};
	
	protected:
		list<solid*> currentsolids; ///< solids in which particle is currently inside
		TGeometry *geom; ///< TGeometry structure passed by "Integrate"
		TMCGenerator *mc; ///< TMCGenerator structure passed by "Integrate"
		TField *field; ///< TField structure passed by "Integrate"
		StepperDopr853<TParticle> *stepper; ///< ODE integrator
		int reflectlog; ///< Should reflections (1) or transmissions (2) or both (3) be logged? Passed by "Integrate"
		FILE *REFLECTLOG; ///< reflection log file passed by "Integrate"


		/**
		 * Set start values, should be called by every constructor.
		 *
		 * Set start values in cylindrical coordinates
		 *
		 * @param number Particle number
		 * @param t Start time
		 * @param tend Life time
		 * @param r Radial start coordinate
		 * @param phi Start azimuth
		 * @param z Start z coordinate
		 * @param Ekin Kinetic start energy
		 * @param alpha Start angle between position vector and velocity vector projected onto xy-plane
		 * @param gamma Start angle between z axis and velocity vector
		 * @param pol Start polarisation
		 * @param afield Electric and magnetic fields (needed to determine total energy)
		 */
		void Init(int number, long double t, long double tend, long double r, long double phi, long double z,
					long double Ekin, long double alpha, long double gamma, int pol, TField *afield){
			ID = ID_UNKNOWN;
			trajlength = comptime = refltime = 0;
			BFsurvprob = vladtotal = 1; vladmax = thumbmax = logBF = -9e99;
			nrefl = NSF = nsteps = 0;
			currentsolids.push_back(&defaultsolid);
			geom = NULL;
			mc = NULL;
			field = NULL;
			stepper = NULL;
			REFLECTLOG = NULL;

			particlenumber = number;
			polarisation = pol;
			xstart = t;
			xend = tend;
			dt = 0;
			rstart = rend = r;
			phistart = phiend = phi;
			Estart = Eend = Ekin;
			alphastart = alphaend = alpha;
			gammastart = gammaend = gamma;
			ystart[0] = yend[0] = rstart * cos(phistart);
			ystart[1] = yend[1] = rstart * sin(phistart);
			ystart[2] = yend[2] = z;

			Hstart = Estart + m*gravconst*ystart[2];
			if (afield){
				long double B[4][4], E[3], V;
				afield->BFeld(ystart[0],ystart[1],ystart[2],xstart,B);
				afield->EFeld(ystart[0],ystart[1],ystart[2],V,E);
				Hstart += q/ele_e*V - polarisation*mu/ele_e*B[3][0];
			}
			long double gammarel = Estart/m/c_0/c_0 + 1;
			if (gammarel < 1.0001)
				vstart = vend = sqrt(2*Estart/m);
			else
				vstart = vend = c_0 * sqrt(1-(1/(gammarel*gammarel)));
			Hmax = Hstart;

			ystart[3] = yend[3] = vstart * cos(phistart + alphastart) * sin(gammastart);
			ystart[4] = yend[4] = vstart * sin(phistart + alphastart) * sin(gammastart);
			ystart[5] = yend[5] = vstart * cos(gammastart);
		};


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
			list<TCollision> colls;
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
					else if (distnormal > 0 && find(currentsolids.rbegin(), currentsolids.rend(), sld) == currentsolids.rend()){
						// if particle is leaving solid and solid is not in currentsolids list something went wrong
						printf("Particle inside '%s' which it did not enter before! t = %LG - Trying to reflect back!\n", geom->solids[coll.ID].name.c_str(),x1-xstart);
						y1[3] = -y1[3];
						y1[4] = -y1[4];
						y1[5] = -y1[5];
						x2 = x1;
						y2 = y1;
//						ID = ID_NRERROR;
						return true;
					}
					else if (distnormal < 0 && find(currentsolids.rbegin(), currentsolids.rend(), sld) != currentsolids.rend()){
						// if particle is entering solid and solid is in currentsolids list something went wrong
						printf("Particle entering '%s' which it did enter before! Stopping it!\n", geom->solids[coll.ID].name.c_str());
						x2 = x1;
						y2 = y1;
						ID = ID_NRERROR;
						return true;
					}
					if (Reflect(x1, y1, x2, y2, coll.normal, coll.ID)) // check if particle should be reflected
						return true;
					else{ // if no reflection -> transmission through surface
						if (distnormal > 0)
							currentsolids.remove(sld); // if particle leaves solid: remove solid from currentsolids list
						else
							currentsolids.push_back(sld); // if particle enters solid: add solid to currentsolids list
						if (Absorb(x1, y1, x2, y2)){ // check for absorption
							printf("Absorption!\n");
							return true;
						}
					}
				}
				else{
					// else cut integration step right before and after first collision point
					// and call ReflectOrAbsorb again for each smaller step (quite similar to bisection algorithm)
					TCollision *c = &colls.front();
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
			else if (Absorb(x1, y1, x2, y2)){ // if there was no collision: just check for absorption
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
		 * @param solid Index of solid in the TGeometry::solids vector
		 * @return Returns true if particle was reflected
		 */
		virtual bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], unsigned solid) = 0;


		/**
		 * Virtual routine to check for absorption in solids, has to be implemented separately for each derived particle.
		 *
		 * @param x1 Start time of line segment
		 * @param y1 Start point of line segment
		 * @param x2 End time of line segment, may be altered
		 * @param y2 End point of line segment, may be altered
		 * @return Returns true if particle was absorbed
		 */
		virtual bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2) = 0;


		/**
		 * Virtual routine which is called when particles reaches its lifetime, has to be implemented for each derived particle.
		 */
		virtual void Decay() = 0;

		/**
		 * Calculates *end variables of particle at time x, position/velocity y and fields B,E,V.
		 *
		 * @param x Time
		 * @param y State vector (position+velocity)
		 * @param B Magnetic field vector, absolute value and its spatial derivatives (see TField)
		 * @param E Electric field vector
		 * @param V Electric potential
		 */
		void SetEndValues(long double x, VecDoub_I &y, long double B[4][4], long double E[3], long double V){
			dt = x - xstart;
			for (int i = 0; i < 6; i++)
				yend[i] = y[i];
			rend = sqrt(y[0]*y[0] + y[1]*y[1]);
			phiend = fmod(atan2(y[1],y[0]) + 2*pi, 2*pi);
			vend = sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
			if (vend/c_0 < 0.01)
				Eend = 0.5*m*vend*vend;
			else{
				long double gammarel = 1/sqrt(1 - vend*vend/(c_0*c_0));
				Eend = c_0*c_0*m*(gammarel - 1);
			}
			Hend = Eend + m*gravconst*yend[2] + q/ele_e*V - polarisation*mu/ele_e*B[3][0];
			if(vend > 0) gammaend = acos(y[5]/vend);
			else gammaend = 0;
			alphaend = fmod(atan2(y[4],y[3]) - phiend + 2*pi, 2*pi);	
			if (Hend > Hmax) Hmax = Hend;
		};


		/**
		 * Print *start and *end values to a file.
		 *
		 * @param file FILE-pointer to print into
		 */
		void Print(FILE *file){
			// calculate spin flip lifetime tauSF and influence on lifetime measurement 
			long double tauSF = (1 - BFsurvprob != 0) ? (-dt/log(BFsurvprob)) : 9e99;
			long double dtau = mc->tau_n - 1/(1/mc->tau_n + 1/tauSF) ;
			 
			int pol = 3;
			if (polarisation == -1) pol = POLARISATION_BAD;
			else if (polarisation == 1) pol = POLARISATION_GOOD;

			fprintf(file ,"%i %i %i %i "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %i %i %LG %.20LG "
			               "%i %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG\n",
			               jobnumber, particlenumber, protneut, pol,
			               xstart, ystart[0], ystart[1], ystart[2],
			               rstart, phistart, vstart, alphastart, gammastart, 
			               Hstart, Estart,
			               xstart + dt, yend[0], yend[1], yend[2],
			               rend, phiend, vend, alphaend, gammaend, dt,
			               Hend, Eend, ID, NSF, comptime, 1 - BFsurvprob,
			               nrefl, vladmax, vladtotal, thumbmax, trajlength,
			               Hmax, tauSF, dtau);
			               
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
			SetEndValues(x,y,B,E,V);
			
			if (vlad > 1e-99) 
				logvlad=log10(vlad);
			if (frac > 1e-99) 
				logfrac=log10(frac);
			int pol = 3;
			if (polarisation == -1) pol = POLARISATION_BAD;
			else if (polarisation == 1) pol = POLARISATION_GOOD;
	
			fprintf(trackfile,"%d %d %d "
							 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
							 "%.17LG %.17LG %.17LG ",
							 particlenumber,protneut,pol,
							 x,y[0],y[1],y[2],y[3],y[4],y[5],
							 vend,Hend,Eend);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					fprintf(trackfile,"%.17LG ",B[i][j]);
			fprintf(trackfile,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG",E[0],E[1],E[2],V,h,logvlad,logfrac,logBF);

			fprintf(trackfile,"\n");
			fflush(trackfile);
		};
};


/**
 * Proton particle class.
 *
 * Simulates a proton including gravitation and Lorentz-force
 */
struct TProton: TParticle{
public:
	/**
	 * Constructor, pass inital values directly.
	 *
		 * @param number Particle number
		 * @param t Start time
		 * @param tend Life time
		 * @param r Radial start coordinate
		 * @param phi Start azimuth
		 * @param z Start z coordinate
		 * @param Ekin Kinetic start energy
		 * @param alpha Start angle between position vector and velocity vector projected onto xy-plane
		 * @param gamma Start angle between z axis and velocity vector
		 * @param pol Start polarisation
		 * @param afield Electric and magnetic fields, needed to determine total energy, can be NULL
	 */
	TProton(int number, long double t, long double tend, long double r, long double phi, long double z,
				long double Ekin, long double alpha, long double gamma, int pol, TField *afield): TParticle(PROTON, ele_e, m_p, 0){
		Init(number, t, tend, r, phi, z, Ekin, alpha, gamma, pol, afield);
	};


	/**
	 * Constructor, create particle, set start values randomly according to all3inone.in.
	 *
	 * Sets start time ::xstart according to [SOURCE] in geometry.in.
	 * Sets kinetic start energy according to TMCGenerator::ProtonSpectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 *
	 * @param number Particle number
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TField used to calculate energies (can be NULL)
	 */
	TProton(int number, TSource &src, TMCGenerator &mcgen, TField *afield): TParticle(PROTON, ele_e, m_p, 0){
		xstart = mcgen.UniformDist(0,src.ActiveTime);

		Estart = mcgen.ProtonSpectrum();
		src.RandomPointInSourceVolume(mcgen, xstart, Estart, ystart[0], ystart[1], ystart[2], alphastart, gammastart);

		long double r = sqrt(ystart[0]*ystart[0] + ystart[1]*ystart[1]);
		long double phi = atan2(ystart[1],ystart[0]);

		Init(number, xstart, mcgen.LifeTime(PROTON), r, phi, ystart[2],
			Estart, alphastart, gammastart, polarisation, afield);
	};


	/**
	 * Constructor, pass initial values directly in cartesian coordinates
	 *
	 * @param number Particle number
	 * @param t Start time
	 * @param pos Start position in cartesian coordinates
	 * @param v Start velocity in cartesian coordinates
	 * @param pol Polarisation
	 * @param mcgen TMCGenerator used to get lifetime from
	 * @param afield TField used to calculate energies (can be NULL)
	 */
	TProton(int number, long double t, long double pos[3], long double v[3], int pol,
			TMCGenerator &mcgen, TField *afield): TParticle(PROTON,ele_e,m_p,0){
		vstart = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		phistart = atan2(pos[1],pos[0]);
		alphastart = fmod(atan2(v[1],v[0]) - phistart + 4*pi,2*pi);
		gammastart = acos(v[2]/vstart);
		Estart = 0.5*m*vstart*vstart;
		Init(number, t, mcgen.LifeTime(PROTON),
			sqrt(pos[0]*pos[0] + pos[1]*pos[1]), phistart, pos[2],
			Estart, alphastart, gammastart, pol, afield);
	};

protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Proton is never reflected
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param normal Normal vector of hit surface
	 * @param solid Index of solid in the TGeometry::solids vector
	 * @return Returns true if particle was reflected
	 */
	bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], unsigned solid){
		return false;
	};


	/**
	 * Checks for absorption
	 *
	 * Protons are immediately absorbed in solids other than ::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @return Returns true if particle was absorbed
	 */
	bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2){
		if (currentsolids.size() > 1){
			x2 = x1;
			y2 = y1;
			ID = currentsolids.back()->ID;
			return true;
		}
		return false;
	};


	/**
	 * Proton decay (not used)
	 */
	void Decay(){

	};
};


/**
 * Electron particle class.
 *
 * Simulates an electron including gravitation and Lorentz-force
 */
struct TElectron: TParticle{
public:
	/**
	 * Constructor, pass inital values directly.
	 *
	 * @param number Particle number
	 * @param t Start time
	 * @param tend Life time
	 * @param r Radial start coordinate
	 * @param phi Start azimuth
	 * @param z Start z coordinate
	 * @param Ekin Kinetic start energy
	 * @param alpha Start angle between position vector and velocity vector projected onto xy-plane
	 * @param gamma Start angle between z axis and velocity vector
	 * @param pol Start polarisation
	 * @param afield Electric and magnetic fields (needed to determine total energy)
	 */
	TElectron(int number, long double t, long double tend, long double r, long double phi, long double z,
				long double Ekin, long double alpha, long double gamma, int pol, TField *afield): TParticle(ELECTRON, -ele_e, m_e, 0){
		Init(number, t, tend, r, phi, z, Ekin, alpha, gamma, pol, afield);
	};


	/**
	 * Constructor, create particle, set start values randomly according to all3inone.in.
	 *
	 * Sets start time ::xstart according to [SOURCE] in geometry.in.
	 * Sets kinetic start energy according to TMCGenerator::ProtonSpectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 *
	 * @param number Particle number
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TField used to calculate energies
	 */
	TElectron(int number, TSource &src, TMCGenerator &mcgen, TField *afield): TParticle(ELECTRON, -ele_e, m_e, 0){
		xstart = mcgen.UniformDist(0,src.ActiveTime);

		Estart = mcgen.ElectronSpectrum();
		src.RandomPointInSourceVolume(mcgen, xstart, Estart, ystart[0], ystart[1], ystart[2], alphastart, gammastart);

		long double r = sqrt(ystart[0]*ystart[0] + ystart[1]*ystart[1]);
		long double phi = atan2(ystart[1],ystart[0]);

		Init(number, xstart, mcgen.LifeTime(ELECTRON), r, phi, ystart[2],
			Estart, alphastart, gammastart, polarisation, afield);
	};


	/**
	 * Constructor, pass initial values directly in cartesian coordinates
	 *
	 * @param number Particle number
	 * @param t Start time
	 * @param pos Start position in cartesian coordinates
	 * @param v Start velocity in cartesian coordinates
	 * @param pol Polarisation
	 * @param mcgen TMCGenerator used to get lifetime from
	 * @param afield TField used to calculate energies (can be NULL)
	 */
	TElectron(int number, long double t, long double pos[3], long double v[3], int pol,
			TMCGenerator &mcgen, TField *afield): TParticle(ELECTRON,-ele_e,m_e,0){
		vstart = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		phistart = atan2(pos[1],pos[0]);
		alphastart = fmod(atan2(v[1],v[0]) - phistart + 4*pi,2*pi);
		gammastart = acos(v[2]/vstart);
		Estart = (1/sqrt(1 - vstart*vstart/c_0/c_0) - 1)*m_e*c_0*c_0;
		Init(number, t, mcgen.LifeTime(ELECTRON),
			sqrt(pos[0]*pos[0] + pos[1]*pos[1]), phistart, pos[2],
			Estart, alphastart, gammastart, pol, afield);
	};

protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Electron is never reflected
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param normal Normal vector of hit surface
	 * @param solid Index of solid in the TGeometry::solids vector
	 * @return Returns true if particle was reflected
	 */
	bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], unsigned solid){
		return false;
	};


	/**
	 * Checks for absorption
	 *
	 * Electrons are immediately absorbed in solids other than ::defaultsolid
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @return Returns true if particle was absorbed
	 */
	bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2){
		if (currentsolids.size() > 1){
			x2 = x1;
			y2 = y1;
			ID = currentsolids.back()->ID;
			return true;
		}
/*		else{
			long double v = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			long double s = sqrt(pow(y2[0]- y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2));
			long double Ekin = m*c_0*c_0*(1/sqrt(1 - v*v/c_0/c_0) - 1);
			double Eloss, theta;
			int index;
			long double n = (1e-3*100)/(1.3806488e-23)/(77); // number density of rest gas (P [Pa])/(k [J/K])/(T [K])
			eH2(Ekin, s, n, &Eloss, &theta, &index); // did scattering happen?
			if (index != 0){
				printf("\nScattering %d, Eloss = %f\n",index,Eloss);
				double v_scat[3] = {y2[3], y2[4], y2[5]};
				eH2velocity(Eloss, theta, v_scat); // change velocity accordingly
				if (index == 3){ // on ionization create secondary electron
					TElectron *e = new TElectron(particlenumber, x2, &y2[0], &y2[3], polarisation, *mc, *field);
					secondaries.push_back(e);
				}
			}
		}
*/		return false;
	};


	/**
	 * Electron decay (not used)
	 */
	void Decay(){

	}
};


/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Constructor, pass inital values directly.
	 *
		 * @param number Particle number
		 * @param t Start time
		 * @param tend Life time
		 * @param r Radial start coordinate
		 * @param phi Start azimuth
		 * @param z Start z coordinate
		 * @param Ekin Kinetic start energy
		 * @param alpha Start angle between position vector and velocity vector projected onto xy-plane
		 * @param gamma Start angle between z axis and velocity vector
		 * @param pol Start polarisation
		 * @param afield Electric and magnetic fields (needed to determine total energy)
	 */
	TNeutron(int number, long double t, long double tend, long double r, long double phi, long double z,
				long double Ekin, long double alpha, long double gamma, int pol, TField *afield): TParticle(NEUTRON, 0, m_n, mu_nSI){
		Init(number, t, tend, r, phi, z, Ekin, alpha, gamma, pol, afield);
	};


	/**
	 * Constructor, create particle, set start values randomly according to all3inone.in.
	 *
	 * Sets start time ::xstart according to [SOURCE] in geometry.in.
	 * Sets start ::polarisation according to config.in.
	 * Sets start energy according to TMCGenerator::NeutronSpectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 * For volume sources the total energy ::Hstart is diced by TMCGenerator::NeutronSpectrum
	 * and the position search is repeated until the kinetic energy ::Estart is >0.
	 * Additionally the kinetic energy spectrum is weighted by sqrt(Estart) (see Golub/Richardson/Lamoreaux p. 82).
	 * For surface sources the kinetic energy ::Estart is diced and all positions are allowed.
	 *
	 * @param number Particle number
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TField used to calculate energies
	 */
	TNeutron(int number, TSource &src, TMCGenerator &mcgen, TField *afield): TParticle(NEUTRON, 0, m_n, mu_nSI){
		xstart = mcgen.UniformDist(0,src.ActiveTime);

		polarisation = mcgen.DicePolarisation();
		Hstart = mcgen.NeutronSpectrum();

		printf("Dice starting position for E_neutron = %LG neV ",Hstart*1e9);
		fflush(stdout);
		long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
		for (int nroll = 0; nroll <= MAX_DICE_ROLL; nroll++){ // try to create particle only MAX_DICE_ROLL times
			if (nroll % 1000000 == 0){
				printf("."); // print progress
				fflush(stdout);
			}

			src.RandomPointInSourceVolume(mcgen, xstart, Hstart, ystart[0], ystart[1], ystart[2], alphastart, gammastart);
			if (afield){
				afield->BFeld(ystart[0], ystart[1], ystart[2], xstart, B);
			}
			if (src.sourcemode == "volume" || src.sourcemode == "customvol"){ // dice H for volume source
				Estart = Hstart - m*gravconst*ystart[2] + polarisation*mu/ele_e*B[3][0];
				if (mcgen.UniformDist(0,1) > sqrt(Estart/Hstart))
					Estart = -1; // correct weighting of phase space (see Golub/Richardson/Lamoreaux p. 82)
			}
			else{ // dice E for surface source
				Estart = Hstart;
				Hstart = Estart + m*gravconst*ystart[2] - polarisation*mu/ele_e*B[3][0];
			}
			if (Estart >= 0){
				printf(" Found! %i dice(s) rolled\n", nroll+1);
				break;
			}
			else if (nroll >= MAX_DICE_ROLL){
				ID = ID_INITIAL_NOT_FOUND;
				printf(" ABORT: Failed %i times to find a compatible spot!! NO particle will be simulated!!\n\n", MAX_DICE_ROLL);
				return;
			}
		}
		long double r = sqrt(ystart[0]*ystart[0] + ystart[1]*ystart[1]);
		long double phi = atan2(ystart[1],ystart[0]);

		Init(number, xstart, mcgen.LifeTime(NEUTRON), r, phi, ystart[2],
			Estart, alphastart, gammastart, polarisation, afield);
	};

protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Each reflection has a certain chance of being diffuse.
	 * Diffuse reflection angles are cosine-distributed around the normal vector of the surface.
	 * Prints surface hits into ::REFLECTLOG if reflectlog is set in config.in.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, set to x1 if reflection happened
	 * @param y2 End point of line segment, returns reflected velocity
	 * @param normal Normal vector of hit surface
	 * @param ID Index of solid in the TGeometry::solids vector
	 * @return Returns true if particle was reflected
	 */
	bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], unsigned solid){
		material *mat = &geom->solids[solid].mat; // get material
		//long double r = sqrt(y1[0]*y1[0] + y1[1]*y1[1]);
		//printf("\nParticle hit %s (material %s) at t=%LG r=%LG z=%LG: ",geom->solids[ID].name.c_str(),mat->name.c_str(),x1,r,y1[2]);

		//************ handle different absorption characteristics of materials ****************

		long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		long double winkeben = 0, winksenkr = 0;
		bool log = false;
		long double prob = mc->UniformDist(0,1);
		long double Trans;
		if (vnormal < 0)
			Trans = Transmission(Enormal*1e9,mat->FermiReal,mat->FermiImag);
		else
			Trans = Transmission(Enormal*1e9 - mat->FermiReal, -mat->FermiReal, -mat->FermiImag);
		int refl = prob < Trans ? 0 : 1; // 0 for transmission, 1 for reflection

		//************** statistical absorption ***********
		if (!refl)
		{
			//printf("Transmission! Erefl=%LG neV Transprob=%LG\n",Enormal*1e9,Trans);
			if(reflectlog & 2)
				log = true;
		}
		else{
			//*************** specular reflexion ************
			prob = mc->UniformDist(0,1);
			nrefl++;
			if (prob >= mat->DiffProb)
			{
				//printf("Specular reflection! Erefl=%LG neV Transprop=%LG\n",Enormal*1e9,Trans);
				if(reflectlog & 1)
					log = true;
				x2 = x1;
				y2 = y1;
				y2[3] -= 2*vnormal*normal[0]; // reflect velocity
				y2[4] -= 2*vnormal*normal[1];
				y2[5] -= 2*vnormal*normal[2];
			}

			//************** diffuse reflection ************
			else
			{
				winkeben = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
				winksenkr = mc->SinCosDist(0, 0.5*pi);
				if (vnormal > 0) winksenkr += pi; // if normal points out of volume rotate by 180 degrees
				if(reflectlog & 1)
					log = true;
				x2 = x1;
				y2 = y1;
				y2[3] = vabs*cos(winkeben)*sin(winksenkr);	// new velocity with respect to z-axis
				y2[4] = vabs*sin(winkeben)*sin(winksenkr);
				y2[5] = vabs*cos(winksenkr);
				RotateVector(&y2[3],normal); // rotate coordinate system so, that new z-axis lies on normal
				//printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG Transprop=%LG\n",Enormal*1e9,winkeben/conv,winksenkr/conv,Trans);
			}
		}
		if (log){
			long double H = 0.5*m_n*vabs*vabs + m_n*gravconst*y1[2];
			if (field){
				long double B[4][4];
				field->BFeld(y1[0], y1[1], y1[2], x1, B);
				H += -polarisation*mu_nSI/ele_e*B[3][0];
			}
			int pol = 3;
			if (polarisation == -1) pol = POLARISATION_BAD;
			else if (polarisation == 1) pol = POLARISATION_GOOD;
			fprintf(REFLECTLOG, "%i %i %i %i %i %i "
								"%.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG "
								"%.10LG %.10LG %.10LG %.10LG "
								"%.10LG %.10LG %.10LG\n",
								jobnumber, particlenumber, protneut, pol, refl, geom->solids[solid].ID,
								x1,y1[0],y1[1],y1[2],y1[3],y1[4],y1[5],
								normal[0],normal[1],normal[2],H,
								winkeben,winksenkr,Trans);
			fflush(REFLECTLOG);
		}

		return refl;
	};


	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * If ::currentsolids contains more than one solid, it uses the one with highest absorption probability
	 * Caclulates adiabacity condition (::frac).
	 * Estimates spin flip probability according to Vladimirskii (::vlad)
	 * and by doing a bruteforce integration of the Bloch equation describing spin precession.
	 * Additionally it can print out a spatial neutron distribution matrix.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @return Returns true if particle was absorbed
	 */
	bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2){
		bool result = false;
/*		long double l = sqrt(pow(y2[0] - y1[0],2) + pow(y2[1] - y1[1],2) + pow(y2[2] - y1[2],2));
		long double v = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		long double E = 0.5*m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5])*1e9;
		long double absprob = 0;
		solid *sld = &defaultsolid;
		for (list<solid*>::iterator i = currentsolids.begin(); i != currentsolids.end(); i++){
			long double abs = Absorption(E, (*i)->mat.FermiReal, (*i)->mat.FermiImag, l);
			if (absprob < abs){
				absprob = abs;
				sld = *i;
			}
		}
		if (mc->UniformDist(0,1) < absprob){
			long double s = mc->UniformDist(0,l)/v;
			x2 = x1 + s*(x2 - x1);
			for (int i = 0; i < 6; i++)
				y2[i] = stepper->dense_out(i, x1 + s*(x2 - x1), stepper->hdid);
			ID = sld->ID;
			result = true;
		}*/
		for (list<solid*>::iterator i = currentsolids.begin(); i != currentsolids.end(); i++){
			if ((*i)->mat.FermiReal != 0 || (*i)->mat.FermiImag != 0){
				x2 = x1;
				y2 = y1;
				ID = (*i)->ID;
				result = true;
			}
		}

		// do special calculations for neutrons (spinflipcheck, snapshots, etc)
		if (neutdist == 1)
			fillndist(x1, y1, x2, y2); // write spatial neutron distribution

		if (field){
			long double B[4][4];
			field->BFeld(y1[0],y1[1],y1[2],x1,B);
			if (B[3][0] > 0){
				// spin flip properties according to Vladimirsky and thumbrule
				vlad = vladimirsky(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				frac = thumbrule(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				vladtotal *= 1-vlad;
				if (vlad > 1e-99)
					vladmax = max(vladmax,log10(vlad));
				if (frac > 1e-99)
					thumbmax = max(thumbmax,log10(frac));
			}

			long double B2[4][4];
			field->BFeld(y2[0],y2[1],y2[2],x2,B2);
			long double sp = BruteForceIntegration(x1,y1,B,x2,y2,B2,min(xstart+xend,StorageTime)); // integrate spinflip probability
			if (1-sp > 1e-30) logBF = log10(1-sp);
			else logBF = -99;
			BFsurvprob *= sp;
			// flip the spin with a probability of 1-BFsurvprob
			if (flipspin && mc->UniformDist(0,1) < 1-sp)
			{
				polarisation *= -1;
				NSF++;
				printf("\n The spin has flipped! Number of flips: %i\n",NSF);
				result = true;
			}
		}

		return result;
	};


	/**
	 * Neutron decay.
	 *
	 * Create decay proton and electron and add them to secondaries list
	 */
	void Decay(){
		long double v_p[3], v_e[3];
		TParticle *p;
		mc->MCZerfallsstartwerte(&yend[3], v_p, v_e);
		p = new TProton(particlenumber, xstart + dt, yend, v_p, polarisation, *mc, field);
		secondaries.push_back(p);
		p = new TElectron(particlenumber, xstart + dt, yend, v_e, polarisation, *mc, field);
		secondaries.push_back(p);
	};
};


#endif // PARTICLE_H_
