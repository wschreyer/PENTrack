#define KMDEF 1000
#define BFKMDEF 50000

// stdincludes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include <stddef.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <signal.h>
#include <list>
#include <vector>
#include <set>
#include <numeric>
#include <map>

// project includes
#include "mersenne/mt.h" // mersenne twister pseudo random number generator (rng) used for MonteCarlo alg.
#include "mc.h" // MonteCarlo related functions
#include "bruteforce.h"
#include "racetrack.h"
#include "adiabacity.h"
#include "fields.h"
#include "2dinterpolfeld.h"
#include "maplefield.h"
#include "BForbes.h"
#include "geometry.h"


// numeric recipes
#include "cel.h"
#include "nrutil.h"
#include "integrator.h" // the integrator functions

//defines
#define NEUTRON 1
#define PROTON 2
#define BF_ONLY 3
#define BF_CUT 4
#define ELECTRONS 6

#define POLARISATION_GOOD 1
#define POLARISATION_BAD 2
#define POLARISATION_NONE 3

#define OUTPUT_EVERYTHING 1
#define OUTPUT_ENDPOINTS 2
#define OUTPUT_NOTHING 5
#define OUTPUT_EVERYTHINGandSPIN 3
#define OUTPUT_ENDPOINTSandSPIN 4

#define KENNZAHL_UNKNOWN 0
#define KENNZAHL_NOT_FINISH 1
#define KENNZAHL_HIT_BOUNDARIES 2
#define KENNZAHL_LEFT_FIELD 3
#define KENNZAHL_DECAYED 4
#define KENNZAHL_INITIAL_NOT_FOUND 5

#define TRUE 1
#define FALSE 0
#define true 1
#define false 0

using namespace std;

//typedef unsigned short int bool;

//functions
void derivs(long double x, long double *y, long double *dydx);
void OpenFiles(int argc, char **argv); // Open in files
void PrepareParticle();
void IntegrateParticle(); // integrate particle trajectory
void PrintIntegrationStep(long double &timetemp);
void BruteForceIntegration(); // integrate spin flip probability
void initialStartbed();
void Startbed(int k);
void ausgabe(long double x2, long double *ystart, long double vend, long double H);
void prepndist(int k);
void fillndist(int k);
void outndist(int k);
void ConfigInit();
void PrintConfig();
void OutputState(long double *y, int l);
void PrintBFieldCut();	// print cut through BField to file
void PrintBField();	// print Bfield to file and investigate ramp heating
void CylKartCoord(long double Wr, long double Wphi, long double Wz, long double phi, long double *Wx0, long double *Wy0, long double *Wz0);
void KartCylCoord(long double Wx, long double Wy, long double Wz, long double phi, long double *Wr0, long double *Wphi0, long double *Wz0);
long double AbsValueCart(long double x, long double y, long double z);

// globals

// global file descriptors
extern FILE *LOGSCR, *OUTFILE1, *BFLOG, *ENDLOG, *STATEOUT, *STARTIN;

// files for in/output + paths
extern string inpath,outpath;
extern char mode_r[2],mode_rw[3],mode_w[2];

// physical constants
extern long double ele_e, Qm0;      //elementary charge in SI, charge/proton mass
extern long double gravconst, conv, mu0;      //g, Pi/180, permeability,
extern long double m_n, pi;  //neutron mass (eV/c^2), neutron magnetic moment (in eV/T), Pi
extern long double m_p;        //proton mass (eV/c^2), tempmass
extern long double m_e, c_0; //electron mass, lightspeed
extern long double hquer, mu_nSI;          // Neutron magn Mom (in J/T)
extern long double gamma_n;
extern long double tau;              // magn. moment of neutron/mass,  neutron lifetime
extern long double lengthconv, Bconv, Econv;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 												

//misc configuration
extern int MonteCarloAnzahl;   // user choice to use MC or not, number of particles for MC simulation
extern int reflekt, bfeldwahl, protneut, Racetracks;       //user choice for reflecting walls, B-field, prot oder neutrons, experiment mode
extern int SaveIntermediate;               // 1: reflections shall be logged, save intermediate steps of Runge Kutta?
extern int polarisation, polarisationsave, ausgabewunsch, ausgabewunschsave, Feldcount; // user choice for polarisation of neutrons und Ausgabewunsch
extern int runge;                            // Runge-Kutta or Bulirsch-Stoer?  set to Runge right now!!!
extern long double Ibar;                // current through rod
extern int diffuse; // diffuse reflection switch
extern int neutdist;

// Fields
extern list<string> fieldvaltab;	// filenames of field table-files
extern long double dBrdr, dBrdz, dBzdr, dBzdz, Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
extern long double dBphidr, dBphidz, dBrdphi, dBzdphi, dBphidphi;
extern long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field
extern long double BFeldSkal, EFeldSkal, BFeldSkalGlobal;                        // parameter to scale the magnetic field for ramping, scale electric field, Global: also scale ramping etc...
extern long double DiceRodField, RodFieldMultiplicator;
extern long double BCutPlanePoint[3], BCutPlaneNormalAlpha, BCutPlaneNormalGamma, BCutPlaneSampleDist; // plane for cut through BField
extern int BCutPlaneSampleCount;
extern long double Emin_n;

// particles
extern long double M, H, mu_n, mumB;                               // total energy of particle
extern long double ystart[7], xstart;       //z-component of velocity, initial and intermediate values of y[9]
extern int iMC;                             //  counter for MonteCarloSim
extern bool noparticle; // true if there was NO particle simulated
extern long double  x1, x2;                         // start and endtime handed over to integrator
extern long int nrefl; // reflection counter
extern long double trajlengthsum;
extern long double Hstart, Hend, Hmax;     //maximum energy
extern long double gammarel, NeutEnergie; // relativistic gamma factor, neutron energy

struct decayinfo				// containing all data from neutron decay for the emerging proton and electron
{	unsigned short int on;		// (!=0) if neutron decay is allowed // "do the neutrons decay?"
								// (==2) if neutron decay in a proton and an electron 
	unsigned short int ed;		// (!=0) if the previous neutron decayed // "did the neutron decay?"
	long int counter;
	int Npolarisation;			// = 'polarisation' of decayed neutron
	long double Nr;				// = 'ystart[1]' (rend) of decayed neutron
	long double Nphi;			// = 'phiend' of decayed neutron
	long double Nz;				// = 'ystart[3]' (zend) of decayed neutron
	long double Nv;				// = 'vend' of decayed neutron
	long double Nalpha;			// = 'alphaend' of decayed neutron
	long double Ngamma;			// = 'gammaend' of decayed neutron
	long double Nx;				// = 'x2' (t) of decayed neutron
	long double NH;				// = 'H' (total energie) of decayed neutron (gravitational + kinetic + B-potential) [neV]
	long double PEnergie;		// = kinetic energy of decay proton
	long double Palpha;			// = 'alphastart' of decay proton
	long double Pgamma;			// = 'gammastart' of decay proton
	long double EEnergie;		// = kinetic energy of decay electron
	long double Ealpha;			// = 'alphastart' of decay electron
	long double Egamma;			// = 'gammastart' of decay electron
	long double error;			// decay error (square root of delta s)
};
extern struct decayinfo decay;	// containing all data from neutron decay for the emerging proton and electron

// initial values of particle
extern long double EnergieS, EnergieE, Energie;    //initial energy range
extern long double r_n, phi_n, z_n, v_n;                //initial particle coordinates
extern long double alpha, gammaa, hmin;                  //initial angle to x-Achse, zo z-Achse, Schrittweite
extern long double alphas, gammas;   //initial values from
extern long double alphae, gammae;   //initial values to
extern long double dalpha, dgamma;   //initial values step
extern long double delx;                           // initial timestep for the integrator

extern long double delx_n;      // shorter timestep for BruteForce integration
extern int stopall;                            //  if stopall=1: stop particle

struct initial												// record for initial values of all 3 particle types (for more details see "all3inone.in")
{	long double EnergieS, EnergieE;				// E
	long double alphas, alphae;						// alpha
	long double gammas, gammae;						// gamma
	long double delx, xend;									// delx | xend
};
extern struct initial nini, pini, eini;						// one 'initial' for each particle type

// final values of particle
extern int kennz;                                  // ending code
extern long double vend, gammaend, alphaend, phiend, xend;    //endvalues for particle

//integrator params
extern long double  eps, h1;  // desired accuracy in ODEINT: normal, for polarisation and phase, trial time step
extern int nvar, nok, nbad;                            // in ODEINT: number of variables in Derivs, good und bad steps

//set the duration of the experiment
extern long double FillingTime;										//Filling in UCN
extern long double CleaningTime;                        // cleaning without field
extern long double RampUpTime;                          // ramping up coils
extern long double FullFieldTime;                       // storing in full field
extern long double RampDownTime;                        // ramping down coils
extern long double EmptyingTime;                        // emptying without field
extern long double StorageTime;                     // time when ramping down shall start, if xend > storage time, let neutron decay
extern int ffBruteForce,ffreflekt,ffspinflipcheck;  // fullfield
extern int ruBruteForce,rureflekt,ruspinflipcheck; // rampup
extern int rdBruteForce,rdreflekt,rdspinflipcheck;  // rampdown
extern int fiBruteForce,fireflekt,fispinflipcheck;  // filling
extern int coBruteForce,coreflekt,cospinflipcheck;  // counting UCN
extern int clBruteForce,clreflekt,clspinflipcheck;  // cleaning time

// Spintracking
extern int spinflipcheck;                          // user choice to check adiabacity
extern long double vlad, vladtotal, frac;                   // adiabacity after Vladimirsky
extern long double vladmax; // maximum values of spinflip prob
extern long double thumbmax;

// file output
extern int jobnumber;
extern unsigned long int monthinmilliseconds; // RandomSeed
extern long Zeilencount;
extern int Filecount, p;                                   // counts the output files, counter for last line written in outs, random generator temp value
extern long double BahnPointSaveTime;               // default 2e-7; 0=1e-19 not changed at run time, time between two lines written in outs

// variables for BruteForce Bloch integration BEGIN
extern long double *BFtime, **BFField;    // time, Bx, By, Bz, r, z array
extern int offset, BFkount, BFindex;			// counter in BFarray, offset off BFarray, maximum index of intermediate values , index in BFarray;
extern long double *BFBws;                    // BFpolarisation
extern long double BFBmin, BFBminmem,BFTargetB;     // smallest value of Babs during step, Babs < BFTargetB => integrate,
extern unsigned short int BruteForce, firstint, flipspin;  // enable BruteForce?, flipthespin?
extern long double I_n[4], **BFypFields;        // Spinvector, intermediate field values in BFodeint
extern long BFZeilencount; extern int BFFilecount;                  // to control output filename of BF
extern long double BFflipprob, BFsurvprob;                         // spinflip pobability

//for the output of intermediate steps
extern int kmax, BFkmax;                                         // number of steps for intermediate output
extern long double nintcalls, ntotalsteps;					// counters to determine average steps per integrator call
extern int kount, hfs, NSF;                               // counter for intermediate output steps, highfieldseeker, numberofspinflips
extern long double *xp,**yp, *BFxp, **BFyp, dxsav, BFdxsav;          // Arrays for intermediate output
extern long double **Bp,**Ep;                            // Arrays for intermediate output

// incorporate B-fieldoszillations into the code
extern int FieldOscillation;        // turn field oscillation on if 1
extern long double OscillationFraction, OscillationFrequency;


extern mt_state_t *v_mt_state; //mersenne twister state var
