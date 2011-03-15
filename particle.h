#include <cmath>
#include <vector>
#include <sys/time.h>

using namespace std;

#define MAX_SAMPLE_DIST 0.01 // max spatial distance of track-entries, reflection checks, spin flip calculation, etc
#define MIN_SAMPLE_DIST 0.005 // min spatial distance of track-entries

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

struct TParticle{
	public:
		int kennz, particlenumber, protneut;
		long double xstart, xend, dt; // start time, max simulation time, actual simulation time
		long double ystart[6], yend[6];
		long double Hstart, Hend, Hmax, NeutEnergie; // start/end/max total energy, total energy of corresponding neutron
		long double Estart, Eend; //start/end kinetic energy
		long double rstart, rend, phistart, phiend, alphastart, alphaend, gammastart, gammaend, vstart, vend;
		long double trajlength, comptime, refltime;	// trajectory length, computing time, reflection computing time
		long double vlad, vladmax, vladtotal, frac, thumbmax, BFsurvprob; // spinflip probabilities, numerical error in decay calculation
		int polarisation; // polarisation of particle (-1,0,1)
		int nrefl, NSF, nsteps; // number of reflections/spin flips/"good" and "bad" integration steps 

		// constructor, pass inital values directly
		TParticle(int aprotneut, int number, long double t, long double tend, long double r, long double phi, long double z,
					long double Ekin, long double alpha, long double gamma, long double E_n, int pol, TField &afield){
			Init(aprotneut, number, t, tend, r, phi, z, Ekin, alpha, gamma, E_n, pol, afield);
		};

		// constructor, create particle, set start values randomly according to all3inone.in
		TParticle(int aprotneut, int number, TSource &src, TMCGenerator &mcgen, TField &afield){
			xstart = mcgen.StartTime(aprotneut);

			polarisation = mcgen.DicePolarisation();
			NeutEnergie = mcgen.NeutronSpectrum();
			long double Emax = NeutEnergie;
			if (polarisation > 0)
				NeutEnergie += src.Hmin_hfs;
			else if (polarisation < 0)
				NeutEnergie += src.Hmin_lfs;

			printf("Dice starting position for E_neutron = %LG neV ",NeutEnergie*1e9);
			fflush(stdout);
			long double B[4][4];
			for (int nroll = 0; nroll <= MAX_DICE_ROLL; nroll++){ // try to create particle only MAX_DICE_ROLL times
				if (nroll % 1000000 == 0){
					printf("."); // print progress
					fflush(stdout);
				}

				src.RandomPointInSourceVolume(mcgen, xstart, NeutEnergie, ystart[0], ystart[1], ystart[2], alphastart, gammastart);
				afield.BFeld(ystart[0], ystart[1], ystart[2], xstart, B);
				if (src.sourcemode == "volume" || src.sourcemode == "customvol"){ // dice H for volume source
					Estart = NeutEnergie - m_n*gravconst*ystart[2] + polarisation*mu_nSI/ele_e*B[3][0];
					if (mcgen.UniformDist(0,1) > sqrt(Estart/Emax))
						Estart = -1; // correct weighting of phase space (see Golub/Richardson/Lamoreaux p. 82)
				}
				else{ // dice E for surface source
					Estart = NeutEnergie;
					NeutEnergie = Estart + m_n*gravconst*ystart[2] - polarisation*mu_nSI/ele_e*B[3][0];
				}
				if (Estart >= 0){
					printf(" Found! %i dice(s) rolled\n", nroll+1);
					break;
				}
				else if (nroll >= MAX_DICE_ROLL){
					kennz = KENNZAHL_INITIAL_NOT_FOUND;
					printf(" ABORT: Failed %i times to find a compatible spot!! NO particle will be simulated!!\n\n", MAX_DICE_ROLL);
					return;
				}
			}
			switch (aprotneut){
				case PROTON: 
					Estart = mcgen.ProtonSpectrum();
					break;
				case ELECTRON:
					Estart = mcgen.ElectronSpectrum();
					break;
			}
			
			long double r = sqrt(ystart[0]*ystart[0] + ystart[1]*ystart[1]);
			long double phi = atan2(ystart[1],ystart[0]);
			
			Init(aprotneut, number, xstart, mcgen.LifeTime(aprotneut), r, phi, ystart[2], 
				Estart, alphastart, gammastart, NeutEnergie, polarisation, afield);
		};
	
		// constructor, create neutron-decay product with given velocity from monte-carlo-neutron-decay
		TParticle(int aprotneut, TParticle *neutron, long double v[3], TMCGenerator &mcgen, TField &afield){
			vstart = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			alphastart = fmod(atan2(v[1],v[0]) - neutron->phiend + 4*pi,2*pi);
			gammastart = acos(v[2]/vstart);
			switch(aprotneut){
				case PROTON:
					Estart = 0.5*m_p*vstart*vstart;
					break;
				case ELECTRON:
					Estart = (1/sqrt(1 - vstart*vstart/c_0/c_0) - 1)*m_e*c_0*c_0;
					break;
				default:
					printf("\nParticle type %d not allowed as neutron decay product!\n",protneut);
					exit(-1);
			}
			Init(aprotneut, neutron->particlenumber, neutron->xend, mcgen.LifeTime(aprotneut),
				neutron->rend, neutron->phiend, neutron->yend[2],
				Estart, alphastart, gammastart, neutron->Hstart, neutron->polarisation, afield);
		};

		// constructor, enter particle starting values interactively
		TParticle(int number, TField &afield){
			cout << "Particle type: ";
			cin >> protneut;
			particlenumber = number;
			cout << "Starting time [s]: ";
			cin >> xstart;
			cout << "Lifetime [s]: ";
			cin >> xend;
			cout << "r: ";
			cin >> rstart;
			cout << "phi [rad]: ";
			cin >> phistart;
			cout << "z: ";
			cin >> ystart[2];
			long double B[4][4], E[3], V;
			afield.BFeld(ystart[0],ystart[1],ystart[2],xstart,B);
			afield.EFeld(ystart[0],ystart[1],ystart[2],V,E);
			switch(protneut){
				case NEUTRON:
					cout << "Polarisation [-1 (lfs), 0, 1 (hfs)]: ";
					cin >> polarisation;
					do{
						if (Estart < 0) cout << "Neutron with energy " << Hstart << " not possible here! Try again..." << endl;
						cout << "Neutron energy [neV]: ";
						cin >> NeutEnergie;
						NeutEnergie *= 1e-9;
						Estart = NeutEnergie - m_n*gravconst*ystart[2] + polarisation*mu_nSI/ele_e*B[3][0];	
					}while (Estart < 0);
					break;
				case PROTON:
					cout << "Proton energy [eV]: ";
					cin >> Estart;
					break;
				case ELECTRON:
					cout << "Electron energy [keV]: ";
					cin >> Estart;
					Estart *= 1e3;
					break;
			}
			cout << "alpha [rad]: ";
			cin >> alphastart;
			cout << "gamma [rad]: ";
			cin >> gammastart;
			
			Init(protneut, number, xstart, xend, rstart, phistart, ystart[2], Estart, alphastart, gammastart, NeutEnergie, polarisation, afield);
		};
	
		// operator(), returns equations of motion, class TParticle is given to integrator which calls TParticle(x,y,dydx)
		void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx){
			derivs(x,y,dydx);
		};
		
		// integrate particle trajectory
		void Integrate(TGeometry &geometry, TMCGenerator &mcgen, TField &afield, FILE *ENDLOG, FILE *trackfile, 
						FILE *SNAP = NULL, set<int> *snapshots = NULL, int areflectlog = 0, FILE *aREFLECTLOG = NULL){
			geom = &geometry;
			mc = &mcgen;
			field = &afield;
			reflectlog = areflectlog;
			REFLECTLOG = aREFLECTLOG;
			timeval clock_start, clock_end, refl_start, refl_end;
			gettimeofday(&clock_start, NULL); // start computing time measure
			int snapshotsdone = -1, perc = 0;
			printf("Teilchennummer: %i, Teilchensorte: %i\n",particlenumber,protneut);
			printf("r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG l: %LG\n",
				rstart, phistart/conv, ystart[2], vstart, alphastart/conv,gammastart/conv, Estart, xend, mc->MaxTrajLength(protneut));
		
			// set initial values for integrator
			long double x = xstart, x1, x2;
			VecDoub y(6,ystart), dydx(6), y1(6), y2(6);
			long double B[4][4], E[3], V;	
			long double h = 0.001/vstart; // first guess for stepsize
			derivs(x,y,dydx);
			if (trackfile)
				PrintIntegrationStep(trackfile, x, y, h); // print start into track file
			long double lastsave = x;
			// create integrator class (stepperdopr853 = 8th order Runge Kutta)
			stepper = new StepperDopr853<TParticle>(y, dydx, x, 1e-13, 1e-13, true);// y, dydx, x, atol, rtol, dense output
			
			while (kennz == KENNZAHL_UNKNOWN){ // integrate as long as nothing happened to particle
				x1 = x; // save point before next step
				y1 = y;
				
				if (x + h > xstart + xend) 
					h = xstart + xend - x;	//If stepsize can overshoot, decrease.
				
				try{
					stepper->step(h, *this); // integrate one step
					nsteps++;
				}
				catch(...){ // catch Exceptions thrown by numerical recipes routines
					kennz = KENNZAHL_NRERROR;
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
						x = x2;
						y = y2;
					}
					gettimeofday(&refl_end, NULL);
					refltime += refl_end.tv_sec - refl_start.tv_sec + (long double)(refl_end.tv_usec - refl_start.tv_usec)/1e6;
				
					trajlength += sqrt(pow(y2[0] - y1[0],2) + pow(y2[1] - y1[1],2) + pow(y2[2] - y1[2],2)); // Trajectory length calculation		
		
					// do special calculations for neutrons (spinflipcheck, snapshots, etc)
					if (protneut == NEUTRON){
						field->BFeld(y1[0],y1[1],y1[2],x1,B);
						field->EFeld(y1[0],y1[1],y1[2],V,E);
						// spin flip properties according to Vladimirsky and thumbrule
/*						vlad = vladimirsky(B[0][0], B[1][0], B[2][0],
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
*/
						long double B2[4][4];
						field->BFeld(y2[0],y2[1],y2[2],x2,B2);
						long double sp = BruteForceIntegration(x1,y1,B,x2,y2,B2); // integrate spinflip pobability
						BFsurvprob *= sp;
						// flip the spin with a probability of 1-BFpol
						if (flipspin && mc->UniformDist(0,1) < 1-sp) 
						{
							polarisation *= -1;
							NSF++;
							printf("\n The spin has flipped! Number of flips: %i\n",NSF);
							x = x2;
							y = y2;
						}
					
						// take snapshots of neutrons at certain times
						if(SNAP && snapshots->count((int) x1) > 0 && (int) x1 != snapshotsdone)
						{
							printf("\n Snapshot at %LG s \n", x1);
							SetEndValues(x1,y1,B,E,V);
							Print(SNAP);
							snapshotsdone=(int) x1;
						}	
			
						if ((neutdist == 1) && (ExpPhase(x) == 4))
							fillndist(x1, y1, x2, y2); // write spatial neutron distribution
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
				if (x >= xstart + xend && mc->decay != 0 && protneut == NEUTRON)
					kennz = KENNZAHL_DECAYED;
				else if (x >= xstart + xend || x > StorageTime || trajlength >= mc->MaxTrajLength(protneut))
					kennz = KENNZAHL_NOT_FINISH;
			
				h = stepper->hnext; // set suggested stepsize for next step
			}
			
			gettimeofday(&clock_end, NULL);	
			comptime = clock_end.tv_sec - clock_start.tv_sec + (long double)(clock_end.tv_usec - clock_start.tv_usec)/1e6;
			delete stepper;
			
			// Endwerte schreiben
			field->BFeld(y[0],y[1],y[2],x,B);
			field->EFeld(y[0],y[1],y[2],V,E);
			SetEndValues(x,y,B,E,V);
			Print(ENDLOG);
			printf("Done!!\nBFFlipProb: %LG rend: %LG zend: %LG Eend: %LG Code: %i t: %LG l: %LG\n\n",
				(1 - BFsurvprob), rend, yend[2], Eend, kennz, dt, trajlength);
		};
	
	private:
		solid *currentsolid; // solid in which particle is currently moving
		TGeometry *geom;
		TMCGenerator *mc;	
		TField *field;	
		StepperDopr853<TParticle> *stepper;
		int reflectlog;
		FILE *REFLECTLOG;

		// set start values, called by every constructor
		void Init(int aprotneut, int number, long double t, long double tend, long double r, long double phi, long double z,
					long double Ekin, long double alpha, long double gamma, long double E_n, int pol, TField &afield){
			protneut = aprotneut;
			kennz = KENNZAHL_UNKNOWN; 
			xstart = xend = dt = 0;
			for (int i = 0; i < 6; i++)
				ystart[i] = yend[i] = 0;
			Hstart = Hend = Hmax = NeutEnergie = Estart = Eend = 0;
			rstart = rend = phistart = phiend = alphastart = alphaend = gammastart = gammaend = vstart = vend = 0;
			trajlength = comptime = refltime = 0;
			BFsurvprob = 1; vladtotal = vladmax = thumbmax = 0;
			polarisation = 0;
			nrefl = NSF = nsteps = 0;
			currentsolid = &defaultsolid;

			particlenumber = number;
			NeutEnergie = E_n;
			polarisation = pol;
			xstart = t;
			xend = tend;
			rstart = r;
			phistart = phi;
			Estart = Ekin;
			alphastart = alpha;
			gammastart = gamma;
			NeutEnergie = E_n;
			ystart[0] = rstart * cos(phistart);
			ystart[1] = rstart * sin(phistart);
			ystart[2] = z;
			long double B[4][4], E[3], V;
			afield.BFeld(ystart[0],ystart[1],ystart[2],xstart,B);
			afield.EFeld(ystart[0],ystart[1],ystart[2],V,E);
			switch(protneut){
				case NEUTRON:
					Hstart = NeutEnergie;
					Estart = Hstart - m_n*gravconst*ystart[2] + polarisation*mu_nSI/ele_e*B[3][0];	
					vstart = sqrt(2*Estart/m_n);
					break;
				case PROTON:
					Hstart = Estart + V + m_p*gravconst*ystart[2];	
					vstart = sqrt(2*Estart/m_p);
					break;
				case ELECTRON:
					Hstart = Estart - V;
					long double gammarel = Estart/m_e/c_0/c_0 + 1;
					vstart = c_0 * sqrt(1-(1/(gammarel*gammarel)));
					break;
			}
			ystart[3] = vstart * cos(phistart + alphastart) * sin(gammastart);
			ystart[4] = vstart * sin(phistart + alphastart) * sin(gammastart);
			ystart[5] = vstart * cos(gammastart);
		};

		// equations of motion dy/dx = f(x,y)
		void derivs(const Doub x, VecDoub_I &y, VecDoub_O &dydx){
			dydx[0] = y[3];
			dydx[1] = y[4];
			dydx[2] = y[5];
			long double B[4][4], E[3], V;
			switch (protneut){
				case NEUTRON:
					field->BFeld(y[0],y[1],y[2], x, B);
					dydx[3] = polarisation*mu_nSI/ele_e/m_n*B[3][1];
					dydx[4] = polarisation*mu_nSI/ele_e/m_n*B[3][2]; // [m/s^2]
					dydx[5] = polarisation*mu_nSI/ele_e/m_n*B[3][3] - gravconst;
					break;
				case PROTON:
					field->BFeld(y[0],y[1],y[2], x, B);
					field->EFeld(y[0],y[1],y[2], V, E);
					dydx[3] = (E[0] + y[4]*B[2][0] - y[5]*B[1][0])/m_p;
					dydx[4] = (E[1] + y[5]*B[0][0] - y[3]*B[2][0])/m_p; 
					dydx[5] = (E[2] + y[3]*B[1][0] - y[4]*B[0][0])/m_p - gravconst;
					break;
				case ELECTRON:
					field->BFeld(y[0],y[1],y[2], x, B);
					field->EFeld(y[0],y[1],y[2], V, E);
					long double gamma = 1/sqrt(1 - (y[3]*y[3] + y[4]*y[4] + y[5]*y[5])/(c_0*c_0));
					long double Qm0 = -1.0/(m_e*gamma);
					long double Fx = E[0] + y[4]*B[2][0] - y[5]*B[1][0];
					long double Fy = E[1] + y[5]*B[0][0] - y[3]*B[2][0];
					long double Fz = E[2] + y[3]*B[1][0] - y[4]*B[0][0]; 
					dydx[3] = Qm0*(Fx - (y[3]*y[3]*Fx + y[3]*y[4]*Fy + y[3]*y[5]*Fz)/c_0/c_0);
					dydx[4] = Qm0*(Fy - (y[4]*y[3]*Fx + y[4]*y[4]*Fy + y[4]*y[5]*Fz)/c_0/c_0);
					dydx[5] = Qm0*(Fz - (y[5]*y[3]*Fx + y[5]*y[4]*Fy + y[5]*y[5]*Fz)/c_0/c_0);
					break;
			}
		};
		
		// check for reflection or absorption
		bool ReflectOrAbsorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, int iteration = 1){
			if (!geom->CheckPoint(x1,&y1[0])){
				kennz = KENNZAHL_HIT_BOUNDARIES;
				return true;
			}
			list<TCollision> colls;
			if (geom->GetCollisions(x1, &y1[0], stepper->hdid, &y2[0], colls)){	// if there is a collision with a wall
				TCollision coll = *colls.begin();
				long double u[3] = {y2[0] - y1[0], y2[1] - y1[1], y2[2] - y1[2]};
				long double distnormal = u[0]*coll.normal[0] + u[1]*coll.normal[1] + u[2]*coll.normal[2];
				// if there is only one collision which is closer to y1 than REFLECT_TOLERANCE: reflect
				if ((colls.size() == 1 && coll.s*abs(distnormal) < REFLECT_TOLERANCE) || iteration > 99){
					if (colls.size() > 1) printf("Reflection from two surfaces (%s & %s) at once!", geom->solids[coll.ID].name.c_str(), geom->solids[(++colls.begin())->ID].name.c_str());
					if (Reflect(x1, y1, x2, y2, coll.normal, coll.ID)) return true;
					else if (Absorb(x1, y1, x2, y2)) return true;
				}
				else{ // else cut integration step right before collision points and call GetReflection again for each smaller step (quite similar to bisection algorithm)
					long double xbisect1 = x1;
					long double xbisect2 = x2;
					VecDoub ybisect1 = y1, ybisect2(6);
					for (list<TCollision>::iterator it = colls.begin(); it != colls.end(); it++){ // go through collisions
						xbisect2 = x1 + (x2 - x1)*it->s*(1 - 0.01*iteration); // cut integration step at this time
						for (int i = 0; i < 6; i++)
							ybisect2[i] = stepper->dense_out(i, xbisect2, stepper->hdid); // get point at cut time
						if (ReflectOrAbsorb(xbisect1, ybisect1, xbisect2, ybisect2, iteration+1)){ // recursive call for each smaller step
							x2 = xbisect2;
							y2 = ybisect2;
							return true;
						}
						xbisect1 = xbisect2;
						ybisect1 = ybisect2;
					}
					if (ReflectOrAbsorb(xbisect2, ybisect2, x2, y2, iteration+1)) // recursive call for last small step
						return true;
				}
			}
			else return (Absorb(x1, y1, x2, y2));
			return false;
		};
	
		// reflect velocity according to particle type and hit material
		bool Reflect(long double x1, VecDoub_I &y1, long double x2, VecDoub_IO &y2, long double normal[3], unsigned ID){
			material *mat = &geom->solids[ID].mat; // get material		
//			long double r = sqrt(y1[0]*y1[0] + y1[1]*y1[1]);
//			Log("\nParticle hit %s (material %s) at t=%LG r=%LG z=%LG: ",solids[ID].name.c_str(),mat->name.c_str(),x1,r,y1[2]);
				
			//************ absorption of electrons/protons ***********
			if (protneut != NEUTRON || !geom->reflekt[ExpPhase(x1)]){
				kennz = geom->solids[ID].kennz;
				printf("Absorption!");
				return true;
			}
			
			//************ handle different absorption characteristics of materials ****************
		
			long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
			if (vnormal > 0 && currentsolid != &geom->solids[ID]){
				printf("Particle inside '%s' which it did not enter before! Stopping it!\n", geom->solids[ID].name.c_str());
				kennz = KENNZAHL_NRERROR;
				return true;
			}
			long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
			long double winkeben = 0, winksenkr = 0;
			bool log = false;
			long double prob = mc->UniformDist(0,1);
			long double Trans=Transmission(Enormal*1e9,mat->FermiReal,mat->FermiImag);
			int refl = prob < Trans ? 0 : 1; // 0 for transmission, 1 for reflection
			
			//************** statistical absorption ***********
			if (!refl)
			{
				// if particle is transmitted through surface set the new medium as currentsolid
				//printf("Transmission! Erefl=%LG neV Transprob=%LG\n",Enormal*1e9,Trans);
				if (vnormal > 0)
					currentsolid = &defaultsolid;
				else
					currentsolid = &geom->solids[ID];
		       	if(reflectlog & 2)
		       		log = true;
			}	
			else{
				//*************** specular reflexion ************
				prob = mc->UniformDist(0,1);
				nrefl++;
				if ((geom->diffuse == 1) || ((geom->diffuse==3)&&(prob >= mat->DiffProb)))
				{
					//printf("Specular reflection! Erefl=%LG neV Transprop=%LG\n",Enormal*1e9,Trans);
					if(reflectlog & 1)
						log = true;
					y2 = y1;
					y2[3] -= 2*vnormal*normal[0]; // reflect velocity
					y2[4] -= 2*vnormal*normal[1];
					y2[5] -= 2*vnormal*normal[2];
				}

				//************** diffuse reflection ************
				else if ((geom->diffuse == 2) || ((geom->diffuse == 3)&&(prob < mat->DiffProb)))
				{
					winkeben = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
					winksenkr = mc->SinCosDist(0, 0.5*pi);
					if (vnormal > 0) winksenkr += pi; // if normal points out of volume rotate by 180 degrees
					if(reflectlog & 1)
						log = true;
					y2 = y1;
					y2[3] = vabs*cos(winkeben)*sin(winksenkr);	// new velocity with respect to z-axis
					y2[4] = vabs*sin(winkeben)*sin(winksenkr);
					y2[5] = vabs*cos(winksenkr);
					RotateVector(&y2[3],normal); // rotate coordinate system so, that new z-axis lies on normal
					//printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG Transprop=%LG\n",Enormal*1e9,winkeben/conv,winksenkr/conv,Trans);
				}
			}
			if (log){
				long double B[4][4];
				field->BFeld(y1[0], y1[1], y1[2], x1, B);
				long double H = 0.5*m_n*vabs*vabs + m_n*gravconst*y1[2] - polarisation*mu_nSI/ele_e*B[3][0];
				fprintf(REFLECTLOG, "%i %i %i %i %i "
									"%.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG "
									"%.10LG %.10LG %.10LG %.10LG "
									"%.10LG %.10LG %.10LG\n",
									jobnumber, particlenumber, protneut, refl, geom->solids[ID].kennz,
									x1,y1[0],y1[1],y1[2],y1[3],y1[4],y1[5],
									normal[0],normal[1],normal[2],H,
									winkeben,winksenkr,Trans);
			}

			return refl;
		};
		
		// check for absorption in current medium
		bool Absorb(long double x1, VecDoub_I &y1, long double x2, VecDoub_IO &y2){
			if (protneut == NEUTRON){
				long double l = sqrt(pow(y2[0] - y1[0],2) + pow(y2[1] - y1[1],2) + pow(y2[2] - y1[2],2));
				long double v = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
				long double E = 0.5*m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5])*1e9;
				if (mc->UniformDist(0,1) < Absorption(E, currentsolid->mat.FermiReal, currentsolid->mat.FermiImag, l)){
					long double s = mc->UniformDist(0,l)/v;
					x2 = x1 + s*(x2 - x1);
					for (int i = 0; i < 6; i++)
						y2[i] = stepper->dense_out(i, x1 + s*(x2 - x1), stepper->hdid);
					printf("Absorption in %s (material %s)!\n",currentsolid->name.c_str(), currentsolid->mat.name.c_str());
					kennz = currentsolid->kennz;
					return true;
				}
			}
			else if (currentsolid != &defaultsolid){
				x2 = x1;
				y2 = y1;
				kennz = currentsolid->kennz;
				return true;
			}
			return false;
		};
	
		// set *end-variables
		void SetEndValues(long double x, VecDoub_I &y, long double B[4][4], long double E[3], long double V){
			dt = x - xstart;
			for (int i = 0; i < 6; i++)
				yend[i] = y[i];
			rend = sqrt(y[0]*y[0] + y[1]*y[1]);
			phiend = fmod(atan2(y[1],y[0]) + 2*pi, 2*pi);
			vend = sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
			switch (protneut){
				case NEUTRON: 
					Eend = 0.5*m_n*vend*vend;
					Hend = Eend + m_n*gravconst*yend[2] - polarisation*mu_nSI/ele_e*B[3][0];
					break;
				case PROTON:
					Eend = 0.5*m_p*vend*vend;
					Hend = Eend + m_p*gravconst*yend[2] + V;
					break;
				case ELECTRON:
					Eend = c_0*c_0*m_e*(1/sqrt(1 - vend*vend/(c_0*c_0)) - 1);
					Hend = Eend - V;
					break;
			}
			if(vend > 0) gammaend = acos(y[5]/vend);
			else gammaend = 0;
			alphaend = fmod(atan2(y[4],y[3]) - phiend + 2*pi, 2*pi);	
			if (Hend > Hmax) Hmax = Hend;
		};

		// print end values to file
		void Print(FILE *file){
			// calculate spin flip lifetime tauSF and influence on lifetime measurement 
			long double tauSF = (1 - BFsurvprob != 0) ? (-dt/log(BFsurvprob)) : 9e99;
			long double dtau = mc->tau_n - 1/(1/mc->tau_n + 1/tauSF) ;
			 
			int pol = 3;
			if (polarisation == -1) pol = 1;
			else if (polarisation == 1) pol = 2;

			fprintf(file ,"%i %i %i %i "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %LG "
			               "%.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %i %i %LG %.20LG "
			               "%i %.20LG %.20LG %.20LG %.20LG "
			               "%.20LG %.20LG %.20LG\n",
			               jobnumber, particlenumber, protneut, pol,
			               xstart, ystart[0], ystart[1], ystart[2],
			               rstart, phistart, vstart, alphastart, gammastart, 
			               Hstart, Estart, NeutEnergie*1.0e9,
			               xstart + dt, yend[0], yend[1], yend[2],
			               rend, phiend, vend, alphaend, gammaend, dt,
			               Hend, Eend, kennz, NSF, comptime, 1 - BFsurvprob,
			               nrefl, vladmax, vladtotal, thumbmax, trajlength,
			               Hmax, tauSF, dtau);
			               
			fflush(file);
		};
	
		// print current integration step to file
		void PrintIntegrationStep(FILE *trackfile, long double x, VecDoub_I &y, long double h){
			long double B[4][4], E[3], V;
			long double logvlad = 0.0, logfrac = 0.0;
			printf("-");
			field->BFeld(y[0],y[1],y[2],x,B);
			field->EFeld(y[0],y[1],y[2],V,E);
			SetEndValues(x,y,B,E,V);
			
			if (vlad > 1e-99) 
				logvlad=log10(vlad);
			if (frac > 1e-99) 
				logfrac=log10(frac);
			int pol = 3;
			if (polarisation == -1) pol = 1;
			else if (polarisation == 1) pol = 2;
	
			fprintf(trackfile,"%d %d %d "
							 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
							 "%.17LG %.17LG %.17LG ",
							 particlenumber,protneut,pol,
							 x,y[0],y[1],y[2],y[3],y[4],y[5],
							 vend,Hend,Eend);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					fprintf(trackfile,"%.17LG ",B[i][j]);
			fprintf(trackfile,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %d",E[0],E[1],E[2],h,logvlad,logfrac,ExpPhase(x));
			fprintf(trackfile,"\n");
		};
};

