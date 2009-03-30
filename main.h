#define KMDEF 1000
#define BFKMDEF 3000

// stdincludes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <stddef.h>
#include <string.h>
#include <iostream.h>
#include <sstream>


// root stuff
#include "TROOT.h"/*
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"*/
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"


// project includes
#include "mersenne/mt.h" // mersenne twister pseudo random number generator (rng) used for MonteCarlo alg.
#include "mc.h" // MonteCarlo related functions
#include "neutron.h"
#include "bruteforce.h"
#include "2dinterpolfeld.h"
#include "adiabacity.h"
#include "reflect.h"
#include "maplefield.h"
#include "BForbes.h"
#include "geomcheck.h"


// numeric recipes
#include "cel.h"
#include "nrutil.h"
#include "integrator.h" // the integrator functions

//defines
#define NEUTRON 1
#define PROTON 2
#define BF_ONLY 3
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
#define KENNZAHL_HIT_BOTTOM 2
#define KENNZAHL_HIT_WALL 3
#define KENNZAHL_ESCAPE_TOP 4
#define KENNZAHL_ESCAPE_SLIT 5
#define KENNZAHL_DETECTOR_BOTTOM 6
#define KENNZAHL_ABSORBED 7
#define KENNZAHL_STILL_THERE 8
#define KENNZAHL_EATEN 9
#define KENNZAHL_DETECTOR_TOP 10

#define TRUE 1
#define FALSE 0
#define true 1
#define false 0

using namespace std;

//typedef unsigned short int bool;

//functions
extern void derivs(long double x, long double *y, long double *dydx);
extern void OpenFiles(int argc, char **argv); // Open in files
extern void PrepareBField(); // Read fieldval tab and preinterpolated or read coils.cond
extern void IntegrateParticle(); // integrate particle trajectory
extern void BruteForceIntegration(); // integrate spin flip probability
extern void initialStartbed();
extern void Startbed(int k);
extern void ausgabe(long double x2, long double *ystart, long double vend, long double H);
extern void prepndist(int k);
extern void fillndist(int k);
extern void outndist(int k);
extern void ConfigInit(void);
extern void OutputState(long double *y, int l);
extern void csleep(int sec);

// globals

// global file descriptors
extern FILE *LOGSCR, *OUTFILE1, *REFLECTLOG, *BFLOG, *TESTLOG, *ENDLOG, *FIN, *STATEOUT, *STARTIN;

// files for in/output + paths
extern string inpath,outpath;
extern char mode_r[2],mode_rw[3],mode_w[2];

// physical constants
extern long double ele_e, Qm0;      //elementary charge in SI, charge/proton mass
extern long double gravconst, conv, mu0;      //g, Pi/180, permeability,
extern long double m_n, mu_n, pi;  //neutron mass (eV/c^2), neutron magnetic moment (in eV/T), Pi
extern long double M,m_p;        //proton mass (eV/c^2), tempmass
extern long double m_e, c_0, gammarel, NeutEnergie; //electron mass, lightspeed, relativistic gamma factor
extern long double hquer, mu_nSI;          // Neutron magn Mom (in J/T)
extern long double gamma_n;
extern long double mumB, tau;              // magn. moment of neutron/mass,  neutron lifetime
extern long double lengthconv, Bconv, Econv;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 												

//misc configuration
extern int MonteCarlo, MonteCarloAnzahl;   // user choice to use MC or not, number of particles for MC simulation
extern int reflekt, Efeldwahl, bfeldwahl, protneut, expmode, Racetracks;       //user choice for reflecting walls, B-field, prot oder neutrons, experiment mode
extern int reflektlog, SaveIntermediate;               // 1: reflections shall be logged, save intermediate steps of Runge Kutta?
extern int polarisation, polarisationsave, ausgabewunsch, ausgabewunschsave, Feldcount; // user choice for polarisation of neutrons und Ausgabewunsch
extern int runge;                            // Runge-Kutta or Bulirsch-Stoer?  set to Runge right now!!!
extern long double Ibar;                // current through rod
extern int diffuse; // diffuse reflection switch
extern long double DiffProb, diffuprob; // property of diffuse reflection 0.125
extern unsigned short int nodelay, slit, DetOpen, decay;                // delays for monte carlo, is there an entrance slit?, do the neutrons decay?;
extern long double lossprob, epsi, AbsProb;   // Lossprobability per Wallbounce, Distance from the Wall within which Reflekt is called, probability of absorption at absorber

// Fields
extern long double dBrdr, dBrdz, dBzdr, dBzdz, Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
extern long double dBphidr, dBphidz, dBrdphi, dBzdphi, dBphidphi;
extern long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field
extern long double BFeldSkal, EFeldSkal, BFeldSkalGlobal;                        // parameter to scale the magnetic field for ramping, scale electric field, Global: also scale ramping etc...
extern long double EFeldSkalSave, BFeldSkalGlobalSave;    // temperorary variables to save initial values
extern int n, m;                               // number of colums and rows in the input B-field matrices
extern long double Vflux, Bre0, Bphie0, Bze0, Be0, Bemax, FluxStep, CritAngle, ElecAngleB, IncidentAngle, DetEnergy, RodFieldMultiplicator;
extern long double DiceRodField;
extern long double VolumeB[200];   // Volume[E] accessible to neutrons at energy E without and with B-field

// particles
extern long double H;                               // total energy of particle
extern long double projz, ystart[7], ysave[7], xstart;       //z-component of velocity, initial and intermediate values of y[9]
extern int iMC;                             //  counter for MonteCarloSim
extern long double  x1, x2;                         // start and endtime handed over to integrator
extern long int kennz0,kennz1,kennz2,kennz3,kennz4,kennz5,kennz6,kennz7,kennz8,kennz9,kennz10,kennz11,kennz12, kennz99,nrefl; // Counter for the particle codes
extern long double trajlength, trajlengthsum, ytemp1, ytemp3, ytemp5;
extern unsigned short int TrajectoryLength;
extern long double Hstart, Hend, Hmax;     //maximum energy

// initial values of particle
extern long double EnergieS,dEnergie, EnergieE, Energie;    //initial energy range
extern long double r_n, phi_n, z_n, v_n;                //initial particle coordinates
extern long double alpha, gammaa,phia, hmin;                  //initial angle to x-Achse, zo z-Achse, Schrittweite
extern long double phis,r_ns, z_ns, v_ns, alphas, gammas;   //initial values from
extern long double phie,r_ne, z_ne, v_ne, alphae, gammae;   //initial values to
extern long double dphi,dr_n, dz_n, dv_n, dalpha, dgamma;   //initial values step
extern long double delx;                           // initial timestep for the integrator

extern long double vr_n, vphi_n, vz_n, vtemp;          //velocity, vtemp: Geschw.komp in xy-Ebene
extern long double delx_n;      // shorter timestep for BruteForce integration
extern int stopall;                            //  if stopall=1: stop particle

struct initial												// record for initial values of all 3 particle types (for more details see "all3inone.in")
{	long double EnergieS, dEnergie, EnergieE;				// E
	long double zs, dz, ze;									// z
	long double rs, dr, re;									// r
	long double phis, dphi, phie;							// phi
	long double alphas, dalpha, alphae;						// alpha
	long double gammas, dgamma, gammae;						// gamma
	long double delx, xend;									// delx | xend
};
extern struct initial nini, pini, eini;						// one 'initial' for each particle type

// final values of particle
extern int kennz;                                  // ending code
extern long double vend, vtest, gammaend, alphaend, phiend, xend;    //endvalues for particle

// geometry of bottle
extern long double ri, ra;                      //dimensions of torus for Maple Field
extern long double zmax, zmin;    //dimensions of torus
extern long double hlid;    // height of a possible lid to be put on the storage bottle [m]
extern long double R, rmax, rmin;
extern long double LueckeR, LueckeZ, Luecke;      // size of the gab in the outer left corner (m)
extern long double wanddicke, wandinnen;                        // Dicke des Bereichs innerhalb der Spulen, der ben�tigt wird
extern long double detz, detrmin, detrmax;          // where is the detector?
extern long double StorVolrmin,StorVolrmax,StorVolzmin, StorVolzmax;
extern long double FillChannelrmin,FillChannelrmax,FillChannelzmin, FillChannelzmax,FillChannelBlockageAngle;
extern long double Bufferrmin,Bufferrmax,Bufferzmin, Bufferzmax;
extern long double DetVolrmin,DetVolrmax,DetVolzmin,DetVolzmax;
extern long double DetConertop, DetConerbot,DetConezbot,DetConeztop;
extern long double FillConertop,FillConerbot,FillConeztop,FillConezbot;
extern long double UCNdetradius,UCNdetr,UCNdetphi;
extern long double UCNentrancermax;
extern long double RoundBottomCornerCenterr,RoundBottomCornerCenterz,RoundBottomCornerradius;

//particle distribution
extern long double rdist, zdist;
extern long double *ndistr, *ndistz, **ndistW;                                            // matrix for probability of finding particle
extern int v,w, neutdist;                                                            // dimension of matrix above, user choice for ndist
extern long double *yyy, *yyy1,*yyy2,*yyy12;        //rectangle with values at the 4 corners for interpolation
extern long double dr, dz;                          // distance between field values in r and z-direction

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
extern long double storagetime;                     // time when ramping down shall start, if xend > storage time, let neutron decay
extern int ffslit,ffBruteForce,ffreflekt,ffspinflipcheck,ffDetOpen;  // fullfield
extern int ruslit,ruBruteForce,rureflekt,ruspinflipcheck,ruDetOpen; // rampup
extern int rdslit,rdBruteForce,rdreflekt,rdspinflipcheck,rdDetOpen;  // rampdown
extern int fislit,fiBruteForce,fireflekt,fispinflipcheck,fiDetOpen;  // filling
extern int coslit,coBruteForce,coreflekt,cospinflipcheck,coDetOpen;  // counting UCN
extern int clslit,clBruteForce,clreflekt,clspinflipcheck,clDetOpen;  // cleaning time

// Spintracking
extern int spinflipcheck;                          // user choice to check adiabacity
extern long double vlad, vladtotal, frac;                   // adiabacity after Vladimirsky
extern long double vladmax; // maximum values of spinflip prob
extern long double thumbmax;
extern long double Bxcoor, Bycoor, Bzcoor;    // B-field in cart Labor coord, cart coord of vector for spin coor sys
//long double matoranorm, matoratime, matoratemptimeb, matoratemptimee,
//            matoratempprob = 1.0;                     // spin flip probabiltiy per second
//long double directprob = 1.0, directtime = 0.0, writeprob = 1.0,

// file output
extern int jobnumber;
extern long Zeilencount;
extern int Filecount, p;                                   // counts the output files, counter for last line written in outs, random generator temp value
extern long double BahnPointSaveTime;               // default 2e-7; 0=1e-19 not changed at run time, time between two lines written in outs
extern char msg[500];

// data for material storage
extern long double FPrealNocado, FPimNocado;     // real and imaginary part of fermi potential for milk tubes
extern long double FPrealPE, FPimPE;  			// for polyethylene
extern long double FPrealTi , FPimTi;			// for titanium
extern long double FPrealCu, FPimCu;			// for Copper
extern long double FPrealCsI, FPimCsI;      // CsI on detector
extern long double FPrealDLC, FPimDLC;      // diamond like carbon
extern int AbsorberChoice;    // 1: PE, 2: Ti

// variables for BruteForce Bloch integration BEGIN
extern long double *BFtime, **BFField;    // time, Bx, By, Bz, r, z array
extern int offset, BFkount, BFindex;			// counter in BFarray, offset off BFarray, maximum index of intermediate values , index in BFarray;
extern long double *BFBws;                    // BFpolarisation
extern long double BFBmin, BFTargetB;     // smallest value of Babs during step, Babs < BFTargetB => integrate,
extern long double BFBxcoor, BFBycoor, BFBzcoor;        // cartesian coord of B field
extern unsigned short int BruteForce, firstint, alwayscont, flipspin;  // enable BruteForce?, flipthespin?
extern long double I_n[4], **BFypFields;        // Spinvector, intermediate field values in BFodeint
extern long BFZeilencount; extern int BFFilecount;                  // to control output filename of BF
extern long double BFflipprob, BFsurvprob;                         // spinflip pobability

 // Absorber integrieren
extern long double abszmin;
extern long double abszmax;
extern long double absrmin; // 48.5, if absorber shall be effectiv, 50.0 if not
extern long double absrmax;
extern long double absphimin;
extern long double absphimax;
extern long double Mf, Pf;  // real and imaginary part if neutron fermi potential of absorber, absorption cross section [1e-28 m^2]
extern int NoAbsorption, AbsorberHits;

//for the output of intermediate steps
extern int kmax, BFkmax;                                         // number of steps for intermediate output
extern long double nintcalls, ntotalsteps;					// counters to determine average steps per integrator call
extern int kount, hfs, NSF;                                            // counter for intermediate output steps, highfieldseeker, numberofspinflips
extern long double *xp,**yp, *BFxp, **BFyp, dxsav;          // Arrays for intermediate output
extern long double **Bp,**Ep;                            // Arrays for intermediate output

// incorporate B-fieldoszillations into the code
extern int FieldOscillation;        // turn field oscillation on if 1
extern long double OscillationFraction, OscillationFrequency;

// coil data for Forbes method
extern long double rFo, phiFo, zFo, aFo,  bFo,  R_0Fo,  J_0Fo, zoffsetFo;
extern int sign1, sign2;
extern long double C1a, C1b, C1R_0, C1J_0, C1zoffset;
extern long double aF[100], bF[100], R_0[100], zoffset[100], J_0[100];
extern int CoilNr;   // number of coils read in

// blank variables
extern long int blankint;
extern long double time_temp;

// define racetrack current bars
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the straight wire 
// current from outside in
extern long double Bars_1r[14],Bars_1phi[14],Bars_1z[14],Bars_2r[14],Bars_2phi[14],Bars_2z[14]; 




extern mt_state_t *v_mt_state; //mersenne twister state var
