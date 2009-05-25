
/*********************************************
pnTracker

**********************************************/
#include "main.h"
 

// files for in/output + paths
FILE *LOGSCR = NULL, *OUTFILE1 = NULL, *REFLECTLOG = NULL, *BFLOG = NULL, *TESTLOG = NULL, *ENDLOG = NULL, *FIN = NULL, *STATEOUT = NULL, *STARTIN = NULL;
string inpath, outpath;
char mode_r[2] = "r",mode_rw[3] = "rw",mode_w[2] = "w"; // modes for fopen()

// physical constants
long double ele_e=1.602176487E-19, Qm0=1.602176487E-19/1.672621637E-27;      //elementary charge in SI, charge/proton mass 																					
long double gravconst=9.80665, conv=0.01745329251, mu0= 1.25663706144e-6;      //g, Pi/180, permeability,
long double m_n=1.674927211E-27/1.602176487E-19, mu_n, pi=3.141592655359;  //neutron mass (eV/c^2), neutron magnetic moment (in eV/T), Pi
long double M,m_p=1.672621637E-27/1.602176487E-19;        //proton mass (eV/c^2), tempmass
long double m_e = 9.10938215e-31/1.602176487E-19, c_0 = 299792458; //electron mass (eV/c^2), lightspeed
long double hquer=1.05457266e-34, mu_nSI=0.96623641e-26;          // Neutron magn Mom (in J/T)
long double gamma_n = 1.83247185e8;            
long double mumB, tau=885.7;              // magn. moment of neutron/mass,  neutron lifetime [s]
long double lengthconv = 0.01 , Bconv = 1e-4, Econv = 1e2;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 
												// cm => m,  Gauss => Tesla,   V/cm => V/m    Bconv temporr von 1e-4 auf 1 gesetzt

// misc configurations
int MonteCarlo=0, MonteCarloAnzahl=1;   // user choice to use MC or not, number of particles for MC simulation
int reflekt=0, Efeldwahl, bfeldwahl, protneut, expmode=1,Racetracks=2;       //user choice for reflecting walls, B-field, prot or neutrons, experiment mode
int reflektlog = 0, SaveIntermediate=0;                // 1: reflections shall be logged, save intermediate steps of Runge Kutta?
int polarisation=0, polarisationsave=0, ausgabewunsch=5, ausgabewunschsave; // user choice for polarisation of neutrons und Ausgabewunsch
long double Ibar= 2250.;                // B-field strength, current through rod
int runge;                            // Runge-Kutta or Bulirsch-Stoer?  set to Runge right now!!!
int diffuse; //diffuse reflection switch
long double DiffProb = 0.16, diffuprob; // property of diffuse reflection 0.125
unsigned short int nodelay, slit = 0, DetOpen=0; // delays for monte carlo, is there an entrance slit?
long double lossprob = 5.0e-4, epsi= 0, AbsProb = 0;   // Lossprobability per Wallbounce, Distance from the Wall within which Reflekt is called, probability of absorption at absorber

// Fields
long double dBrdr, dBrdz, dBrdphi=0.0, dBphidr, dBphidz, dBphidphi=0.0, dBzdr, dBzdz, dBzdphi=0.0;
long double Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field
long double BFeldSkal = 1.0, EFeldSkal = 1.0, BFeldSkalGlobal = 1.0;          // parameter to scale the magnetic field for ramping, scale electric field, Global: also scale ramping etc...
long double EFeldSkalSave, BFeldSkalGlobalSave;    // temperorary variables to save initial values
int n, m;                               // number of colums and rows in the input B-field matrices
long double Vflux, Bre0, Bphie0, Bze0, Be0, Bemax, FluxStep=0.001, CritAngle, ElecAngleB, IncidentAngle, DetEnergy, RodFieldMultiplicator = 0.0;
long double DiceRodField=0;
long double VolumeB[200] = {0.0};   // Volume[E] accessible to neutrons at energy E without and with B-field
long double BCutPlanePoint[3], BCutPlaneNormalAlpha, BCutPlaneNormalGamma, BCutPlaneSampleDist;	// plane for B field cut
int BCutPlaneSampleCount;

// particles
long double H;                               // total energy of particle
long double projz, ystart[7], ysave[7], xstart = 0;       //z-component of velocity, initial and intermediate values of y[9]
long double  x1, x2;                         // start and endtime handed over to integrator
int iMC;                             //  counter for MonteCarloSim
long double trajlength, trajlengthsum, ytemp1, ytemp3, ytemp5;
unsigned short int TrajectoryLength=1;
long double Hstart, Hend,Hmax;     //maximum energy
long double gammarel, NeutEnergie;	//relativistic gamma factor, Neutron energy

struct decayinfo decay; //containing all data from neutron decay for the emerging proton and electron

// inital values of particle
long double EnergieS,dEnergie, EnergieE, Energie;    //initial energy range
long double r_n, phi_n, z_n, v_n;                //initial particle coordinates
long double alpha, gammaa, hmin;                  //initial angle to x-Achse, zo z-Achse, Schrittweite
long double phis,r_ns, z_ns, v_ns, alphas, gammas;   //initial values from
long double phie,r_ne, z_ne, v_ne, alphae, gammae;   //initial values to
long double dphi,dr_n, dz_n, dv_n, dalpha, dgamma;   //initial values step
long double delx;                            // initial timestep for the integrator

long double vr_n, vphi_n, vz_n, vtemp;          //velocity, vtemp: Geschw.komp in xy-Ebene
long double delx_n=0.0;
int stopall=0, Feldcount=0;                            //  if stopall=1: stop particle

struct initial nini, pini, eini; 	// one 'initial' for each particle type

// final values of particle
int kennz;                                  // ending code
long double vend, vtest, gammaend, alphaend, phiend, xend;    //endvalues for particle
long int kennz0[3]={0},kennz1[3]={0},kennz2[3]={0},kennz3[3]={0},kennz4[3]={0},kennz5[3]={0},kennz6[3]={0},kennz7[3]={0},kennz8[3]={0},kennz9[3]={0},kennz10[3]={0},kennz11[3]={0},kennz12[3]={0},kennz99[3]={0},nrefl; // Counter for the particle codes

// geometry of bottle
long double ri= 0.12, ra= 0.48;                      //dimensions of torus for Maple field
long double zmax=1.2, zmin=0.0;    //dimensions of torus
long double hlid;    // height of a possible lid to be put on the storage bottle [m]
long double rmax, rmin;
long double LueckeR=0.001, LueckeZ=0.05, Luecke=0.05;      // size of the gap in the outer left corner (m)
long double wanddicke, wandinnen;                        // Dicke des Bereichs innerhalb der Spulen, der bentigt wird
long double detz, detrmin, detrmax;          // where is the detector?
long double StorVolrmin=0.129,StorVolrmax=0.488,StorVolzmin=0.01, StorVolzmax=1.345;   // main storage hollow cylinder
long double FillChannelrmin=0.4935,FillChannelrmax=0.54131,FillChannelzmin=-0.188, FillChannelzmax=0.02543,FillChannelBlockageAngle=6.1;     // filling channel at at outer bottom edge, Angle: 4 times on every corner there is the racetrack coil crossing
long double Bufferrmin=0,Bufferrmax=0.54131,Bufferzmin=-0.45, Bufferzmax=-0.188;   // buffervolume below, which has filling opening and detector embedded
long double DetVolrmin=0.0889,DetVolrmax=0.488,DetVolzmin=1.345,DetVolzmax=1.565;  // proton detector volume
long double DetConerbot=0.488, DetConertop=0.3559,DetConezbot=1.345,DetConeztop=1.565; // cone for focussing coils in this space
long double FillConerbot=0.54131,FillConertop=0.488,FillConezbot=-0.022,FillConeztop=0.02543; // cone to reflect UCN in filling channel
long double UCNdetradius=0.06,UCNdetr=0.45,UCNdetphi=0;  // UCN detector in buffer volume - phi in Bogenmass
long double UCNentrancermax=0.1;  // UCN entrance tube at bottom center of buffer volume
long double RoundBottomCornerCenterr,RoundBottomCornerCenterz,RoundBottomCornerradius=0.02; // round entrance corner

// particle distribution
long double *ndistr = NULL, *ndistz = NULL, **ndistW = NULL;    // matrix for probability of finding particle
int v=300,w=1200, neutdist;                                       // dimension of matrix above, user choice for ndist
long double *yyy = NULL, *yyy1 = NULL, *yyy2 = NULL, *yyy12 = NULL;        //rectangle with values at the 4 corners for interpolation

// integrator params
long double  eps, h1;  // desired accuracy in ODEINT: normal, for polarisation and phase, trial time step
int nvar, nok, nbad;                            // in ODEINT: number of variables in Derivs, good und bad steps

//set the duration of the experiment
long double FillingTime = 0;										// filling time, entrance open
long double CleaningTime = 0;                        // cleaning without field
long double RampUpTime = 0;                          // ramping up coils
long double FullFieldTime = 1000;                       // storing in full field
long double RampDownTime = 5;                        // ramping down coils
long double EmptyingTime = 100;                        // emptying without field
long double storagetime = 1500.0;                     // time when ramping down shall start, if xend > storage time, let neutron decay
int ffslit,ffBruteForce,ffreflekt,ffspinflipcheck,ffDetOpen;  // fullfield
int ruslit,ruBruteForce,rureflekt,ruspinflipcheck,ruDetOpen; // rampup
int rdslit,rdBruteForce,rdreflekt,rdspinflipcheck,rdDetOpen;  // rampdown
int fislit,fiBruteForce,fireflekt,fispinflipcheck,fiDetOpen;  // filling
int coslit,coBruteForce,coreflekt,cospinflipcheck,coDetOpen;  // counting UCN
int clslit,clBruteForce,clreflekt,clspinflipcheck,clDetOpen;  // cleaning time
int jobnumber;
unsigned long int monthinmilliseconds; // RandomSeed

// Spintracking
int spinflipcheck = 0;                          // user choice to check adiabacity
long double vlad=0.0, vladtotal = 1.0, frac;                   // adiabacity after Vladimirsky
long double vladmax=0.0; // maximum values of spinflip prob
long double thumbmax=0.0;
long double Bxcoor, Bycoor, Bzcoor;    // B-field in cart Labor coord, cart coord of vector for spin coor sys
//long double matoranorm, matoratime, matoratemptimeb, matoratemptimee,
//            matoratempprob = 1.0;                     // spin flip probabiltiy per second
//long double directprob = 1.0, directtime = 0.0, writeprob = 1.0,

// file output
long Zeilencount;
int Filecount=1, p;                                   // counts the output files, counter for last line written in outs, random generator temp value
long double BahnPointSaveTime = 5.0e-7;               // default 2e-7; 0=1e-19 not changed at run time, time between two lines written in outs
char msg[500];

// data for material storage
long double FPrealNocado = 183.04, FPimNocado = 0.018985481;     // real and imaginary part of fermi potential for milk tubes
long double FPrealPE = -8.56, FPimPE = 0.001912531;  			// for polyethylene
long double FPrealTi = -50.76, FPimTi = 0.024983971;			// for titanium
long double FPrealCu = 169.98, FPimCu = 0.023134523;			// for Copper
long double FPrealCsI = 29.51, FPimCsI = 0.03;     // CsI on proton detector
long double FPrealDLC=256, FPimDLC=0.00182;      // diamond like carbon: mu = 1e-4 values from PSI measurement[J. Res. Natl. Inst. Stand. Technol. 110, 279-281 (2005)]
int AbsorberChoice = 1;    // 1: PE, 2: Ti

// variables for BruteForce Bloch integration BEGIN
long double *BFtime=NULL, **BFField=NULL;   // time, Bx, By, Bz, r, z array
int offset=0, BFkount, BFindex = 3;			// counter in BFarray, offset of BFarray, maximum index of intermediate values , index in BFarray
long double BFpol, *BFBws=NULL;                    // BFpolarisation
long double BFBmin = 10.0, BFBminmem=10.0, BFTargetB=0.1;     // smallest value of Babs during step, memorize value during value accumulation, Babs < BFTargetB => integrate,
long double BFBxcoor, BFBycoor, BFBzcoor;        // cartesian coord of B field
unsigned short int BruteForce = 0, firstint = 1, flipspin=1;  // enable BruteForce?,
long double I_n[4], **BFypFields=NULL;        // Spinvector, intermediate field values in BFodeint
long BFZeilencount; int BFFilecount=1;                  // to control output filename of BF
long double BFflipprob = 0.0, BFsurvprob=1.0;                // spinflip probability, survival (non-flip) probability

 // Absorber integrieren
long double abszmin = 0.6;
long double abszmax = 0.8;
long double absrmin = 0.285; // 48.5, if absorber shall be effectiv, 50.0 if not
long double absrmax = 0.29;
long double absphimin = 0.0;
long double absphimax = 360.0;
long double Mf = -8.56, Pf = 0.00191; // real and imaginary part if neutron fermi potential of absorber, absorption cross section [1e-28 m^2]
int NoAbsorption = 0, AbsorberHits = 0;

//for the output of intermediate steps
//#define KMDEF 1000
int BFNrIntermediate=BFKMDEF;    // number of steps for intermediate output in the case of 
int kmax=0, BFkmax = BFKMDEF;                                         // number of steps for intermediate output
long double nintcalls=0, ntotalsteps=0;     // counters to determine average steps per integrator call
int kount, hfs, NSF = 0;                                            // counter for intermediate output steps, highfieldseeker: +1 yes, -1 lowfieldseeker, number of spin flips
long double *xp=NULL,**yp=NULL, *BFxp=NULL, **BFyp=NULL, dxsav=0, BFdxsav=0;          // Arrays for intermediate output
long double **Bp=NULL,**Ep=NULL;                            // Arrays for intermediate output
// END output of intermediate steps 

// incorporate B-fieldoszillations into the code
int FieldOscillation = 0;        // turn field oscillation on if 1
long double OscillationFraction = 1e-4, OscillationFrequency = 1;    // Frequency in Hz

// coil data for Forbes method
long double rFo, phiFo, zFo, aFo,  bFo,  R_0Fo,  J_0Fo, zoffsetFo;
int sign1, sign2;
long double C1a = 0.02, C1b=0.015, C1R_0 = 0.515, C1J_0=3e8, C1zoffset = 0;
long double aF[100], bF[100], R_0[100], zoffset[100], J_0[100];
int CoilNr=0;   // number of coils read in

// blank variables
long int blankint;
long double time_temp;

// define racetrack current bars
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the straight wire 
long double Bars_1r[14],Bars_1phi[14],Bars_1z[14],Bars_2r[14],Bars_2phi[14],Bars_2z[14];  // lower horizontal 0 deg





mt_state_t *v_mt_state = NULL; //mersenne twister state var

// uebergabe: jobnumber inpath outpath                 paths without last slash
int main(int argc, char **argv){
	time_t mytime;
	tm *monthday;
	timeval daysec;
	ldiv_t divresult; 

	// for random numbers we need a statevar + we need to set an initial seed
	mt_state_t mtstate;
	v_mt_state = &mtstate;
	mytime = time(NULL);
	
	monthday = localtime(&mytime);
	monthinmilliseconds = (unsigned long int)(monthday->tm_mday)*24*60*60*1000; // add day in ms
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_hour)*60*60*1000; // add hour in ms
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_min)*60*1000; // add minute in ms 
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_sec)*1000;  // add second in ms	
	gettimeofday(&daysec, 0);
	divresult = div((long int)(daysec.tv_usec), (long int)(1000));
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(divresult.quot); // add milliseconds
	
	mt_set (v_mt_state, monthinmilliseconds);
	
	// setting some default values
	nvar=6;           // number of variables
	eps=1.0e-13;      // desired relative precision for tracking of n and p 10^-13
	dxsav=1e-5;  //Kleinster ausgegebener zeitschritt Neutronen
	hmin= 0;       // minimum stepsize for runge kutta
	
	// compute values for RoundBottomCorner
	RoundBottomCornerCenterr=FillChannelrmin-RoundBottomCornerradius;
	RoundBottomCornerCenterz=StorVolzmin-RoundBottomCornerradius;	
	// END compute values for RoundBottomCorner
	
	
	// current bar test
	// current from outside in
	Bars_1r[1]=0.60;Bars_1phi[1]=0.0;	Bars_1z[1]=-0.15;	Bars_2r[1]=0.0;	Bars_2phi[1]=0.0;	Bars_2z[1]=-0.15;   // lower horizontal 0 deg
	Bars_1r[2]=0.60;Bars_1phi[2]=pi/2.0;Bars_1z[2]=-0.15;	Bars_2r[2]=0.0;	Bars_2phi[2]=pi/2.0;Bars_2z[2]=-0.15; // lower horizontal 90 deg
	Bars_1r[3]=0.60;Bars_1phi[3]=pi;	Bars_1z[3]=-0.15;	Bars_2r[3]=0.0;	Bars_2phi[3]=pi;	Bars_2z[3]=-0.15; // lower horizontal 180 deg
	Bars_1r[4]=0.60;Bars_1phi[4]=pi*1.5;Bars_1z[4]=-0.15;	Bars_2r[4]=0.0;	Bars_2phi[4]=pi*1.5;Bars_2z[4]=-0.15; // lower horizontal 270 deg
	// current from inside out
	Bars_1r[5]=0.0;Bars_1phi[5]=0.0;	Bars_1z[5]=1.35;	Bars_2r[5]=0.6;	Bars_2phi[5]=0.0;	Bars_2z[5]=1.35;   // upper horizontal 0 deg
	Bars_1r[6]=0.0;Bars_1phi[6]=pi/2.0;	Bars_1z[6]=1.35;	Bars_2r[6]=0.6;	Bars_2phi[6]=pi/2.0;Bars_2z[6]=1.35;  // upper horizontal 90 deg
	Bars_1r[7]=0.0;Bars_1phi[7]=pi;		Bars_1z[7]=1.35;	Bars_2r[7]=0.6;	Bars_2phi[7]=pi;	Bars_2z[7]=1.35;  // upper horizontal 180 deg
	Bars_1r[8]=0.0;Bars_1phi[8]=pi*1.5;	Bars_1z[8]=1.35;	Bars_2r[8]=0.6;	Bars_2phi[8]=pi*1.5;Bars_2z[8]=1.35;  // upper horizontal 270 deg
	// current from high to low
	Bars_1r[9]=0.60;Bars_1phi[9]=0;		Bars_1z[9]=1.35;	Bars_2r[9]=0.6; Bars_2phi[9]=0;		Bars_2z[9]=-0.15;  //outer current 0 deg
	Bars_1r[10]=0.6;Bars_1phi[10]=pi/2;	Bars_1z[10]=1.35;	Bars_2r[10]=0.6;Bars_2phi[10]=pi/2;	Bars_2z[10]=-0.15; //outer current 90 deg
	Bars_1r[11]=0.6;Bars_1phi[11]=pi;	Bars_1z[11]=1.35;	Bars_2r[11]=0.6;Bars_2phi[11]=pi;	Bars_2z[11]=-0.15; //outer current 180 deg
	Bars_1r[12]=0.6;Bars_1phi[12]=pi*1.5;Bars_1z[12]=1.35;	Bars_2r[12]=0.6;Bars_2phi[12]=pi*1.5;Bars_2z[12]=-0.15; //outer current 270 deg
	// current from low to high
	Bars_1r[13]=0.0;Bars_1phi[13]=0.0;	Bars_1z[13]=-0.15;	Bars_2r[13]=0;	Bars_2phi[13]=0.0;	Bars_2z[13]=1.35;   // center current 4 TIMES THE CURRENT OF OTHERS!!!
  
  /*
	long double err, testr = 0.25, result;
	ystart[5]=pi/2;
	ystart[3]=-0.1;
	//protneut=2;
	fprintf(TESTLOG,"coord Br Bphi Bz \n");
	for(int p=-100;p<200;p++)
	{
	Br=0; Bphi=0;Bz=0;dBrdr =0;dBrdphi=0;dBrdz=0;dBphidr=0;dBphidphi =0;dBphidz=0;dBzdr=0;dBzdphi=0;dBzdz=0;
	BarRaceTrack(0.1, -pi/3, p/100.0, 2250);	
	cout << "z " << p/100.0 << " (" << Br << ", " << Bphi << ", " << Bz << ")" << endl;
	//cout << "Ableitungen" << dBrdr << " " << dBrdphi << " " << dBrdz << " " << dBphidr << " " << dBphidphi << " " << dBphidz << " " << dBzdr << " " << dBzdphi << " " << dBzdz << endl;
	fprintf(TESTLOG,"%.17LG %.17LG %.17LG %.17LG\n",(long double) p/100.0, Br,Bphi,Bz);	
	
		//result = dfridr(LowHor0_r, 0.25, 0.01, &err);
	}
	
	cout << "All B (" << Br << ", " << Bphi << ", " << Bz << ")" << endl;
	cout << "Componentwise B (" << BarRaceTrack_Br(0.2, pi/5, 0.2, 2250) << ", " << BarRaceTrack_Bphi(0.2, pi/5, 0.2, 2250) << ", " << BarRaceTrack_Bz(0.2, pi/5, 0.2, 2250) << ")"<< endl;
	cout << "Ableitungen" << dBrdr << " " << dBrdphi << " " << dBrdz << " " << dBphidr << " " << dBphidphi << " " << dBphidz << " " << dBzdr << " " << dBzdphi << " " << dBzdz << endl;
	// cin >> blankint;
	
	//return 0;
	// current bar test end
	*/
	
	// globals init end
	
	if(argc>3) // if user supplied 3 args (outputfilestamp, inpath, outpath)
	{
		outpath = argv[3]; // set the output path pointer
		inpath = argv[2]; // same with input path pointer
		jobnumber = atoi(argv[1]); // stamp for output filenames
	}
	else if(argc>2) // if user supplied 2 args (outputfilestamp, inpath)
	{
		inpath = argv[2]; // input path pointer set
		jobnumber = atoi(argv[1]); 
		outpath = "./out"; // setting outpath to default
	}
	else if(argc==2) // if user supplied 1 arg (outputfilestamp)
	{
		jobnumber=atoi(argv[1]);
		outpath = "./out";
		inpath = "./in";
	}
	else // no args supplied
	{
		jobnumber=0;
		outpath = "./out";
		inpath = "./in";
	}
	
	// initial step ... reading userinput, inputfiles etc ...
	printf(
	" ################################################################\n"
	" ###                 Welcome to PNTracker,                    ###\n"
	" ###     the tracking program for neutrons and protons        ###\n"
	" ################################################################\n");
	
	ConfigInit();
	OpenFiles(argc, argv);	// Open .in and .out files and write headers
	
	//printf("\nMonteCarlo: %i\n MonteCarloAnzahl %i \n", MonteCarlo, MonteCarloAnzahl);
	
	// allocate vectors and matrizes for BruteForce only if necessary
	if(BruteForce||clBruteForce||coBruteForce||fiBruteForce||ruBruteForce||ffBruteForce||rdBruteForce)
	{
		BFtime=dvector(0,BFNrIntermediate);		
		BFField=dmatrix(1,5,0,BFNrIntermediate);
		BFBws=dvector(0,BFNrIntermediate);
		BFypFields=dmatrix(1,3,0,BFNrIntermediate);
		BFxp=dvector(0,BFNrIntermediate);
		BFyp=dmatrix(1,3,0,BFNrIntermediate);		
	}
	
	// allocate vector intermediate values when desired by the user
	if((ausgabewunsch==1)||(ausgabewunsch==3)||BruteForce||ruBruteForce||ffBruteForce||rdBruteForce||neutdist)
	{
		xp=dvector(1,kmax);
		yp=dmatrix(1,6,1,kmax);         
		Bp=dmatrix(1,13,1,kmax);
		Ep=dmatrix(1,2,1,kmax);
	}
	
	//printconfig();
		
	if (bfeldwahl != 1){
		PrepareBField();	// read fieldval.tab or coils.cond
	}

	Startbed(1); // read in starting values of particle
	
	PrintConfig();
	
	PrepareParticle(); // setting values depending on particle type

	if((protneut == NEUTRON)&&((ausgabewunsch == OUTPUT_EVERYTHINGandSPIN)||(ausgabewunsch == OUTPUT_ENDPOINTSandSPIN)))
	{
		//fprintf(TESTLOG,"t log(pol) spinflipprob\n");
		fprintf(BFLOG,"t Babs Polar logPolar Ixnorm Iynorm Iznorm Bxnorm Bynorm Bznorm dx dy dz\n");
		// for field interpolation tests fprintf(BFLOG,"r z Br Bphi Bz BrTab BphiTab BzTab deltabr_abs deltabphi_abs deltabz_abs deltabr_rel deltabphi_rel deltabz_rel\n"); 
	}
	
	//end initial step
	
	/*printf("\nm=%i \nn=%i\n",m,n);
	for(int i=3;i<=m-1;i=i+2){
		for(int j=3;j<=n-1;j=j+2){
			//if(rind[i]>0.275 && rind[i]<0.305 && zind[j]>0.045 && zind[j]<0.075){
			if(rind[i]>0.13 && rind[i]<0.49 && zind[j]>0.01 && zind[j]<1.1){
				BFeld(rind[i],0.0,zind[j],500.0);
				fprintf(BFLOG,"%.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG %.16LG\n",rind[i],zind[j],Br,Bphi,Bz,BrTab[i][j],BphiTab[i][j],BzTab[i][j],fabsl(BrTab[i][j]-Br),fabsl(BphiTab[i][j]-Bphi),fabsl(BzTab[i][j]-Bz),fabsl((BrTab[i][j]-Br)/BrTab[i][j]),fabsl((BphiTab[i][j]-Bphi)/BphiTab[i][j]),fabsl((BzTab[i][j]-Bz)/BzTab[i][j]));
			}
		}
		
	}
	//return;*/
	
	decay.counter = 0; // no neutron had decayed yet (no chance)

	iMC = 1;
	do
	{      
		for(Energie=EnergieS; Energie<=EnergieE; Energie+=dEnergie)
		{
			for(r_n=r_ns; r_n<=r_ne; r_n+=dr_n)
			{
				for(z_n=z_ns; z_n<=z_ne; z_n+=dz_n)
				{
					for(alpha=alphas; alpha<=alphae; alpha+=dalpha)
					{
						for(gammaa=gammas; gammaa<=gammae; gammaa+=dgamma)
						{
							for(phi_n=phis; phi_n<=phie; phi_n+=dphi)
							{
								IntegrateParticle();
							}							
						}
					}
				}
			}
		}
		if((decay.on==2) && decay.ed)
		{	switch(protneut) // simulation sequenz in decay case: neutron, proton, electron
			{	case NEUTRON:	protneut = PROTON;
								decay.counter++;
								break;	
				case PROTON:	protneut = ELECTRONS;
								break;
				case ELECTRONS:	protneut = NEUTRON;
								iMC++; // !!
								break;
			}
			initialStartbed(); // initalizes initial values from record to used variables
			PrepareParticle(); // setting values depending on particle type
		}
		else
		{	iMC++; // !!
		}
	} while(iMC <= MonteCarloAnzahl);
	
	if (neutdist == 1) outndist(1);   // Neutronenverteilung in der Flasche ausgeben
	OutputCodes(iMC);

	// for investigating ramp heating of neutrons, volume accessible to neutrons with and
	// without B-field is calculated and the heating approximated by thermodynamical means
	if (protneut == BF_ONLY){
		fprintf(LOGSCR,"\nEnergie [neV], Volumen ohne B-Feld, mit B-Feld, 'Erwaermung'");
		int i;
		long double Volume;
		for (i = 0; i <= EnergieE; i++) 
		{
			Volume = ((i * 1.0e-9 / (M * gravconst))-wanddicke) * pi * (r_ne*r_ne-r_ns*r_ns);
			// isentropische zustandsnderung, kappa=5/3
			fprintf(LOGSCR,"\n%i %.17LG %.17LG %.17LG",i,Volume,VolumeB[i],i * powl((Volume/VolumeB[i]),(2.0/3.0)) - i);
		}
	}
	//printf("The B field time is:%.17LG\n",timer1);
	//printf("The integration time is:%.17LG\n",timer2);
	printf("Integrator used (1 Bulirsch Stoer, 2 Runge Kutta): %d \n", runge);
	fprintf(LOGSCR,"Integrator used (1 Bulirsch Stoer, 2 Runge Kutta): %d \n ", runge);
	// printf("We spent %.17LG seconds for BF-field interpolation\n",timer3);
	printf("The integrator was called: %LF times with %LF internal steps on average. \n", nintcalls,ntotalsteps/nintcalls);
	fprintf(LOGSCR,"The integrator was called: %LF times with %LF internal steps on average. \n", nintcalls,ntotalsteps/nintcalls);
	printf("That's it... Have a nice day!\n");
	fprintf(LOGSCR,"That's it... Have a nice day!\n");
	
	// cleanup ... lassen wir bleiben macht linux fuer uns *hoff*	
	/*if(LOGSCR != NULL)
		fclose(LOGSCR);
	if(ausgabewunsch == OUTPUT_EVERYTHING)
		fclose(OUTFILE1);
	if(REFLECTLOG != NULL)
		fclose(REFLECTLOG);
	if(ENDLOG != NULL)
		fclose(ENDLOG);
	if(BFLOG != NULL)
		fclose(BFLOG);
	if(TESTLOG != NULL)
		fclose(TESTLOG);
	if(BFLOG != NULL)
		fclose(BFLOG);*/
		
	return 0;
}

//======== setting values depending on particle type ==================================================================
void PrepareParticle()
{	
	switch(protneut)
	{	case NEUTRON:	M = m_n;			// [eV/c^2]
						h1 = 5e-5;			// guess for the first step size of runge kutta
						dxsav = 1e-5;		// kleinster ausgegebener zeitschritt Neutronen
						EnergieS = EnergieS*1.e-9;
						EnergieE = EnergieE*1.e-9;
						dEnergie = dEnergie*1.e-9;
						if(polarisationsave == POLARISATION_GOOD) // in ev/T
						{	hfs = -1;
							mu_n = hfs * mu_nSI / ele_e; // => mu_n is negative
							mumB = mu_n/M;	// [c^2/T]
						}
						else if(polarisationsave == POLARISATION_BAD) // in ev/T
						{	hfs = 1;
							mu_n = hfs * mu_nSI / ele_e; // => mu_n is positive
							mumB = mu_n/M;	// [c^2/T]
						}
						else if(polarisationsave == POLARISATION_NONE) // in ev/T
						{	hfs = 0;
							mu_n = 0;
							mumB = mu_n/M;	// [c^2/T]
						}
						decay.ed = 0;
						break;		
		case PROTON:	M = m_p;			// [eV/c^2]
						Qm0 = 1.0/M;
						h1 = 1e-8;		// guess for the first step size of runge kutta
						dxsav = 1e-10;	// kleinster ausgabeschritt der zwischenwerte im integrator
						BahnPointSaveTime = 1e-8;
						reflekt = 0;
						break;		
		case BF_ONLY:	fprintf(ENDLOG,"r phi z Br Bphi Bz 0 0 0 \n");
						break;		
		case BF_CUT:	PrintBFieldCut();
						exit(0);
						break;		
		case ELECTRONS:	M=m_e;			// [eV/c^2]
						Qm0 = -1.0/M;
						h1 = 2e-10;		// guess for the first step size of runge kutta
						dxsav = 2e-12;	// kleinster ausgabeschritt der zwischenwerte im integrator
						BahnPointSaveTime = 5e-12;
						reflekt = 0;
						//fprintf(ENDLOG,"rstart zstart vr vphi vz ElecAngleB Dethit? CritAngle Ekin Br0 Bz0 Babs0 Babsm rend zend Babsend Vdiff EnergyonDet IncAngle\n");
						break;
	}
}
//======== end of PrepareParticle =====================================================================================

//int derivsaufrufe=0;
void derivs(long double x, long double *y, long double *dydx){
	
	Entkommen(y, x, H);
	// call the B-field
	BFeld(y[1],y[5],y[3], x);
	//if((bfeldwahl==0)&&(Bphi==0))
	//			OutputState(y,1);
	
	if (protneut != NEUTRON)
		EFeld(y[1],y[5],y[3]);
	
	//fprintf(LOGSCR,"fc: %i\n",++derivsaufrufe);
	
	if (protneut == NEUTRON)
	{ // neutron equations of motion
		// Bahnverfolgung ausgeschaltet => gerade Bahn der Teilchen
		/*              
		dydx[1]= y[2];
		dydx[2]= 0;
		dydx[3]= y[4];
		dydx[4]= 0;
		dydx[5]= y[6];
		dydx[6]= 0;
		*/
		
		dydx[1]= y[2];
		dydx[2]= y[1]*(y[6]*y[6])+mumB*dBdr;
		dydx[3]= y[4];
		dydx[4]= mumB*dBdz-gravconst;
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+mumB/y[1]*dBdphi/y[1];
		
	}
	else if(protneut == PROTON)
	{ // equations for the proton
		dydx[1]= y[2];
		dydx[2]= y[1]*(y[6]*y[6])+Qm0*(Bz*y[1]*y[6]-Bphi*y[4]+Er);
		dydx[3]= y[4];
		dydx[4]= Qm0*(Ez+Bphi*y[2]-y[1]*Br*y[6]);
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+Qm0*(Ephi+Br*y[4]-Bz*y[2])/y[1];
	}
	else if(protneut == ELECTRONS)
	{ // equations for the electron (relativistic mass)
		Qm0 = -1.0/M*sqrtl(1-(y[2]*y[2]+y[1]*y[1]*y[6]*y[6]+y[4]*y[4])/(c_0*c_0));
		dydx[1]= y[2];
		dydx[2]= y[1]*(y[6]*y[6])+Qm0*(Bz*y[1]*y[6]-Bphi*y[4]+Er);
		dydx[3]= y[4];
		dydx[4]= Qm0*(Ez+Bphi*y[2]-y[1]*Br*y[6]);
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+Qm0*(Ephi+Br*y[4]-Bz*y[2])/y[1];
	}
	return;
}

void OpenFiles(int argc, char **argv){
	// printing the path into the vars
	ostringstream logscrfile;
	logscrfile << outpath << "/" << jobnumber << "log.out";
	LOGSCR = fopen(logscrfile.str().c_str(),mode_w);

	fprintf(LOGSCR,
	" ################################################################\n"
	" ###                 Welcome to PNTracker,                    ###\n"
	" ###     the tracking program for neutrons and protons        ###\n"
	" ################################################################\n");


	if(reflektlog == 1){
		ostringstream reflectlogfile;
		reflectlogfile << outpath << "/" << jobnumber << "reflect.out";
		REFLECTLOG = fopen(reflectlogfile.str().c_str(),mode_w);
		fprintf(REFLECTLOG,"t r z phi x y diffuse vabs Eges Erefl winkeben winksenkr vr vz vtang phidot dvabs\n"); // Header for Reflection File
	}
	if((ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)||(ausgabewunsch==OUTPUT_ENDPOINTSandSPIN))
	{
		ostringstream BFoutfile1;
		BFoutfile1 << outpath << "/" << jobnumber << "BF1.out";
		BFLOG = fopen(BFoutfile1.str().c_str(),mode_w);
	}
	
	// Endpunkte
	ostringstream endlogfile;
	endlogfile << outpath << "/" << jobnumber << "end.out";
	ENDLOG = fopen(endlogfile.str().c_str(),mode_w);
	if (protneut != BF_ONLY) 
	{
        fprintf(ENDLOG,"jobnumber RandomSeed protneut polarisation tstart rstart phistart zstart NeutEnergie vstart alphastart "
        			   "gammastart rend phiend zend vend alphaend gammaend t H kennz "
        			   "NSF RodFieldMult BFflipprob AnzahlRefl vladmax vladtotal thumbmax trajlength Hdiff Hmax "
        			   "AbsorberHits BFeldSkal EFeldSkal lossprob tauSF dtau\n");
	}		
	
	// Print track to file
	if ((ausgabewunsch == OUTPUT_EVERYTHING)||(ausgabewunsch == OUTPUT_EVERYTHINGandSPIN))
	{ 
		SaveIntermediate=1; // turn on saving of intermediate values in integrator
		kmax=KMDEF;
		ostringstream wholetrackfile;
		wholetrackfile << outpath << "/" << jobnumber << "track1.out";
		OUTFILE1 = fopen(wholetrackfile.str().c_str(),mode_w);       // open outfile neut001.out
		Zeilencount=0;
		fprintf(OUTFILE1,"Teilchen protneut t r drdt z dzdt phi dphidt x y "
						 "v H Br dBrdr dBrdphi dBrdz Bphi dBphidr dBphidphi dBphidz "
						 "Bz dBzdr dBzdphi dBzdz Babs Er Ez timestep logvlad logthumb\n");
	}

	// read starting values from file start.in   NOT IMPLEMENTED YET
	// open this file
	if (MonteCarlo==2)
	{
		string path;
		path = inpath + "/start.in";
		FILE *STARTIN = fopen (path.c_str(),mode_r);
		if (STARTIN == NULL) exit(-1);        // Fehlerbehandlung
		fgets(msg,500,STARTIN);
	}
}

void PrepareBField(){
	if ((bfeldwahl == 0)||(bfeldwahl == 2))
	{
        printf("\nPreparing the electromagnetic fields... \n");
        PrepIntpol(1);          // read input table file with E and B fields
		printf("allocating space for preinterpolation ... (about %.4LG MB)\n",(long double)n*m*12*16*3/1024/1024);

		// doing the preinterpolation
		Preinterpol(1);						
		
		/*
		// interpolation test
		wholetrackfile = (char*)malloc((outpathlength+20)*sizeof(char));
		sprintf(wholetrackfile, "%s/%06dinterpol.out", outpath, jobnumber);
		OUTFILE1 = fopen(wholetrackfile,mode_w);       // open outfile neut001.out
		Zeilencount=0;
		fprintf(OUTFILE1,"r z Br dBrdr dBrdphi dBrdz Bz dBzdr dBzdphi dBzdz Babs\n");
		//fprintf(OUTFILE1,"indr indz c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34 c41 c42 c43 c44 \n");
		BFeldSkal=1;
		long double rtst = 0.45;
		for(long double ztst = 0.137; ztst<=0.15; ztst=ztst+0.0001)
		{
			//BInterpol(rtst,0,ztst);
			BInterpol(rtst,0,ztst);
			fprintf(OUTFILE1,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG \n",
													rtst, ztst,Br,dBrdr,dBrdphi,dBrdz,Bz,dBzdr,dBzdphi,dBzdz,Bws);
			
		}
		return 0;
		// ENDE interpolation test*/
		
	}
	else if(bfeldwahl==4)
	{
		ReadMagnets();
	
		printf("\n \n Test of integration\n");
		//	long double TestInt;
		BFeldSkal=1.0;
		sign1 = 1, sign2 = 1;
		for (int a = 0;a<1;a++)
		{
		
			BFeld(0.3,0,0.1, 500.0);
			cout << "T" << endl;
		}
		printf("Br = %.17LG \n",Br);
		printf("dBrdr = %.17LG \n",dBrdr);
		printf("dBrdz = %.17LG \n",dBrdz);
		printf("Bz = %.17LG \n",Bz);
		printf("dBzdr = %.17LG \n",dBzdr);
		printf("dBzdz = %.17LG \n",dBzdz);	
	}	
}

void IntegrateParticle(){
	int schritte=0,iii=0, perc=0;   // ?? , ??, percentage of particle done counter
	unsigned short int DetHit=0;
		// reset some values for new particle
	stopall=0;
	kennz=0; // not categorized yet									
	// initial values for Brute-Force Spinintegration 
		BFpol = 0.5;
		I_n[3]=0.5; I_n[2]=I_n[1]=0;
		offset = 0;
		BFsurvprob = 1.0;
		BFflipprob = 0.0;
	Hmax=0.0;
	NoAbsorption = 0;
	AbsorberHits = 0;									
	x2=x1= 0.;     //set time to zero
	
	// randomly chosen start parameters
	if(MonteCarlo==1)
	{			
		int BadParticle=1;
		while(BadParticle==1)
		{
			MCStartwerte(delx);   // MonteCarlo Startwert und Lebensdauer fr Teilchen festlegen
			if(noparticle) break;
			long double ytemp[7]={r_n,0,z_n,0,phi_n,0};
			BadParticle = GeomCheck(ytemp,r_n,0,z_n,0,phi_n,0,0.0);
			if (BadParticle==1)
				cout << "Not within boundaries at start... redicing...!" << endl;
		}
		x2=xstart;    // set time starting value										
	}	
	
	// remember starting energy of particle
	Hstart = Energie;
	
	if(BruteForce)
	{
		NSF=0;
		firstint = 1;                 // after BF step, it will not be the first BF step any more...
		if(protneut == NEUTRON){
	
			if(polarisation==POLARISATION_GOOD){
				hfs = -1;
				mu_n=hfs * mu_nSI / ele_e;
				mumB= mu_n/M;  // [c^2/T]
			}else if(polarisation==POLARISATION_BAD){  // in ev/T
				hfs = 1;
				mu_n=hfs * mu_nSI / ele_e;
				mumB= mu_n/M;  // [c^2/T]
			} else if (polarisation==POLARISATION_NONE){
				hfs = 0;
				mu_n=0;
				mumB= mu_n/M;  // [c^2/T]
			}  // in ev/T
		}
	}
	
	// start parameters from a file with format of end.out file
	if (MonteCarlo==2)
	{			
		do
		{
			fgets(msg,250,STARTIN);
		}while(!feof(STARTIN));
		if(protneut == PROTON)
		{
			long double rtmp=r_n, phitmp=phi_n,ztmp=z_n;  // store position values for later 
			MCStartwerte(delx);  // determine energy and velocity vector of proton randomly...
			r_n=rtmp; phi_n=phitmp; ztmp=z_n;   // put position values read in back in place
			
		}
	}
	
	//-------- If there is no particle to be simulated, make a note in the ENDLOG and go on with the next particle. --------
	if(noparticle)
	{	xstart = 0.0;
		r_n = 0.0;
		phi_n = 0.0;
		z_n = 0.0;
		v_n = 0.0;
		alpha  = 0.0;
		gammaa = 0.0;
		ystart[1] = 0.0;
		ystart[2] = 0.0;
		phiend = 0.0;
		ystart[3] = 0.0;
		vend = 0.0;
		alphaend = 0.0;
		gammaend = 0.0;
		x2 = 0.0;
		H = 0.0;
		NSF = 0;
		nrefl = 0;
		vladmax = 0.0;
		thumbmax = 0.0;
		trajlengthsum = 0.0;
		Hstart = 0.0;
		Hmax = 0.0;
		AbsorberHits = 0;
		// output in the ENDLOG
		ausgabe(0.0, ystart, 0.0, 0.0);
		// to ensure only ONE loop
		if (MonteCarlo)
		{	Energie = EnergieE + 1;
			z_n = z_ne + 1;
			r_n = r_ne + 1;
			alpha = alphae + 1;
			gammaa = gammae + 1;
			phi_n = phie + 1;
		}
		// reset 'noparticle' to false 		
		noparticle = false;
		return;
	}
	//-------- Finished ----------------------------------------------------------------------------------------------------
	
	if (protneut != BF_ONLY)
	{
		printf("\nRodFieldMultiplicator: %.17LG\n",RodFieldMultiplicator);
		fprintf(LOGSCR,"\nRodFieldMultiplicator: %.17LG\n",RodFieldMultiplicator);
		printf("Feldcount = %i\n\n",Feldcount);
		fprintf(LOGSCR,"Feldcount = %i\n\n",Feldcount);
	}

	trajlengthsum = 0;
	nrefl=0;
	kennz=KENNZAHL_UNKNOWN;  // 0 
	stopall=0;
	
	Feldcount=0;

	if (protneut == NEUTRON)
	{// Spinflipwahrscheinlichkeiten zurcksetzen
		thumbmax=0.0;
		vladtotal = 1.0;
		vladmax=0.0;									
		BFeld(r_n,phi_n*conv,z_n, 0.0);							
		long double Ekin=Energie-M*gravconst*z_n+mu_n*Bws;      // kin Energie = Anfangsen. - Pot Energie + Energie im B-Feld
		if(Ekin>=0.0)
		{
			v_n=powl(2.0/M*(Ekin),0.5);
			if (bfeldwahl == 3) v_n=1.0;    // for analytic input field set v manually
			ausgabewunsch=ausgabewunschsave;
		}else
		{
			v_n=0.0;
			stopall=1;
			ausgabewunsch=5;
			printf("\nEkin: %.17LG  smaller than Zero!!! \n",Ekin);
			fprintf(LOGSCR,"\nEkin: %.17LG  smaller than Zero!!! \n",Ekin);
			//iMC--;
			 if (!nodelay)
				csleep(1);
			
			if (MonteCarlo)
			{
				Energie=EnergieE+1;
				z_n=z_ne+1;
				r_n=r_ne+1;
				alpha=alphae+1;
				gammaa=gammae+1;												
			}
				
			return;
		}
	}
	else if(protneut == PROTON)
	{                // Proton
		v_n=powl(2*Energie/m_p,0.5);     // Gravitiationsenergie vernachlssigbar, auch B-Feld Energie
		cout << "Proton: Energy: " << Energie << " v= " << v_n << " m/s ";
	}
	else if(protneut == ELECTRONS)
	{
		gammarel = Energie/m_e/c_0/c_0 + 1;
		v_n = c_0 * sqrtl(1-(1/(gammarel*gammarel)));     // relatistic velocity from E_kin in eV
		//cout << "Energie: " << Energie << " m_e: " << m_e << " v_n: " << v_n << endl;
		//sleep(5);										
	}
	
	if(protneut != BF_ONLY)
	{
		projz= cosl(conv*gammaa);  // projection of velocity on z-axis
		vz_n= v_n*projz;						// multiplied by the velocity
		vtemp= v_n*sinl(conv*gammaa);  // projection of velocity on x-y plane
		vr_n=  vtemp*cosl(conv*(alpha-phi_n));  // 
		vphi_n=vtemp*sinl(conv*(alpha-phi_n));
		ystart[1]= r_n; ystart[2]= vr_n;           // fill array for ODEint integrator
		ystart[3]= z_n; ystart[4]= vz_n;
		ystart[5]= conv*phi_n;
		if (TrajectoryLength)
		{ // Trajectory length
			ytemp1=ystart[1]; 
			ytemp3=ystart[3]; 
			ytemp5=ystart[5];      
		}
		if (r_n!=0.) 
			ystart[6]= vphi_n/r_n;
		else if (r_n==0.)
			ystart[6]= 0.;
		
										
		//x2=x1= 12.600989800684586; // only temporary
		Entkommen(ystart,  x2, H);                 // Test ob das Teilchen am Anfang schon in einem falschen Bereich ist
	}
	
	//Hier wird nur der erste Punkt geschrieben
	if(protneut == NEUTRON)                // n
		H= (m_n*gravconst*ystart[3]+0.5*m_n*v_n*v_n-mu_n*Bws)*1E9 ;       // Energie in neV
	else if(protneut == PROTON)           // p
	{
		H= (0.5*m_p*v_n*v_n);                                           // Energie in eV
		cout << " Energy: " << H << " eV" << endl;
	}
	else if(protneut == ELECTRONS)           // e-
	{
		H= c_0*c_0  * m_e * ( 1/sqrtl(1-v_n*v_n/(c_0*c_0)) - 1) ;                                        // rel Energie in eV
		cout << " Energy: " << H << " eV" << endl;
	}
	Hmax = H;
	
	// do integration for neutrons, protons or electrons 
	if(protneut == NEUTRON || protneut == PROTON || protneut == ELECTRONS)
	{
/*		if(BruteForce)
		{
			// set initial value of polarisation vector: parallel to BField, length 0.5
			I_n[1]= Bxcoor/Bws*0.5;
			I_n[2]= Bycoor/Bws*0.5;
			I_n[3]= Bzcoor/Bws*0.5;
			
			//I_n[1]=.002326711314962561; I_n[2]=.007738274230667618; I_n[3]=.9999673522302555; // only temporary
		}
*/		
		
		printf("Teilchennummer: %i\n",iMC);
		fprintf(LOGSCR,"Teilchennummer: %i\n",iMC);
		if(decay.on == 2)
		{	printf("Teilchensorte : %i\n", protneut);
			fprintf(LOGSCR,"Teilchensorte : %i\n", protneut);
		}
		printf("r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG\n",r_n,phi_n,z_n,v_n,alpha,gammaa,H,xend);
		fprintf(LOGSCR,"r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG\n",r_n,phi_n,z_n,v_n,alpha,gammaa,H,xend);

		schritte = 0;   // zhlt die integrationsschritte mit
		

		//-----------------------------------------------------
		// Schleife fr ein Teilchen, bis die Zeit aus ist oder das Teilchen entkommt
		long double timetemp = 0;                                    // temporre Variable, Zeit wann letzter Schritt in outs geschrieben wurde
		do
		{
			 //if(x2 >= 35.5)  OutputState(ystart,1);
			
			BFeld(ystart[1],ystart[5],ystart[3], x2);

			if( (Bws<(BFTargetB+0.1) ) && BruteForce)
				delx_n = delx/10;
			else if((Bws<BFTargetB) && BruteForce)
				delx_n = delx/100;
			else
				delx_n = delx;

			x1= x2; x2=x1+delx_n;                 // determine start and endtime of step
			schritte++;
			if (TrajectoryLength) ytemp1 = ystart[1]; ytemp3 = ystart[3]; ytemp5 = ystart[5];    // for trajectory length calculation
														
			// put phi (ystart[5]) back to [-2Pi,2Pi]
			if(ystart[5]>(2.0*pi))
				ystart[5]=ystart[5]-2*pi;
			if(ystart[5]<(-2.0*pi))
				ystart[5]=ystart[5]+2*pi;
			
			//mytime1 = clock();
			//###################### Integrationsroutine #####################
			if (runge==2)  (*odeint) (ystart,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);           // runge kutta step
			if (runge==1)  (*odeint) (ystart,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);        // bulirsch stoer step
														// (ystart: input vector | nvar: number of variables | x1, x2: start and end time | eps: precision to be achieved | h1: guess for first stepsize | hmin: mininum stepsize | nok,nbad: number of good and bad steps taken | derivs: function for differential equation to be integrated, rkqs, bsstep: integrator to be used (runge kutta, bulirsch stoer) )
			//###################### Integrationsroutine #####################
			nintcalls++;
			ntotalsteps=ntotalsteps+kount;
			
			
			// write status of particle to console
			// one point per 2 percent of endtime
			if((x2-xstart)*100/(xend-xstart)>perc)
			{
				perc++;
				if(perc%10==0)
				{
					printf("%i%%",perc);
					fflush(stdout);
				}
				if((perc%10!=0)&&(perc%2==0)) 
				{
					printf(".");
					fflush(stdout);
				}
			}
			
			vend = sqrtl(fabsl(ystart[2]*ystart[2]+ystart[1]*ystart[1]*ystart[6]*ystart[6]+ystart[4]*ystart[4]));
			if(protneut == NEUTRON)                // n
				H = (M*gravconst*ystart[3]+0.5*M*vend*vend-mu_n*Bws)*1E9 ;       // Energie in neV
			else if(protneut == PROTON)           // p			
				H= (0.5*m_p*vend*vend);           // Energie in eV for p						
			else if(protneut == ELECTRONS)           // e
				H= c_0*c_0  * M *  (1/sqrtl(1-vend*vend/(c_0*c_0))-1);                                        // rel Energie in eV
		
			if (H>Hmax) Hmax=H;

			if ((neutdist == 1)&&(protneut == NEUTRON))
				fillndist(1);

			
			if ((!BruteForce) && (protneut != PROTON) && (protneut != ELECTRONS))
				BahnPointSaveTime = 1e-3;
			else if (BruteForce)
				BahnPointSaveTime = 1e-4;
			
			if(spinflipcheck==3)
				BahnPointSaveTime = 1e-4;
			

			if (BruteForce)
			{
				BruteForceIntegration();
			}
			
			//Ausgabe der Zwischenwerte aus odeint
			int klauf;
			long double logvlad = 0.0, logfrac = 0.0;
			if (((ausgabewunsch==OUTPUT_EVERYTHING)||(ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)) && ((x2-x1)>=BahnPointSaveTime)){
				for (klauf=2;klauf<=kount;klauf++){
					if ((xp[klauf]-time_temp)>=BahnPointSaveTime)
					{
						
						printf("-");
						fflush(stdout);
						// Ausgabewerte berechnen
						vend = sqrtl(fabsl(yp[2][klauf]*yp[2][klauf]+yp[1][klauf]*yp[1][klauf]*yp[6][klauf]*yp[6][klauf]+yp[4][klauf]*yp[4][klauf]));
						if(protneut == NEUTRON)
							H = (M*gravconst*yp[3][klauf]+0.5*M*vend*vend-mu_n*Bp[13][klauf])*1E9 ;    // mu_n negative for low-field seekers
						else if(protneut == PROTON)           // p			
							H= (0.5*m_p*vend*vend);           // Energie in eV for p		
						else if(protneut == ELECTRONS)  
							H= c_0*c_0  * M * (1/sqrtl(1-vend*vend/(c_0*c_0))-1);                                        // rel Energie in eV
						
						if (spinflipcheck==2)
						{
							if (vlad>1e-99) 
								logvlad=log10l(vlad);
							if (frac>1e-99) 
								logfrac=log10l(frac);
						}
						
						//cout << "Br " << Bp[1][klauf] << endl;
						fprintf(OUTFILE1,"%d %d %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
										 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
										 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG \n",
										 iMC,protneut,xp[klauf],yp[1][klauf],yp[2][klauf],yp[3][klauf],yp[4][klauf],yp[5][klauf],yp[6][klauf],yp[1][klauf]*cosl(yp[5][klauf]),yp[1][klauf]*sinl(yp[5][klauf]),
										 vend,H,Bp[1][klauf],Bp[2][klauf],Bp[3][klauf],Bp[4][klauf],Bp[5][klauf],Bp[6][klauf], Bp[7][klauf],Bp[8][klauf],
										 Bp[9][klauf],Bp[10][klauf],Bp[11][klauf],Bp[12][klauf],Bp[13][klauf],Ep[1][klauf],Ep[2][klauf],x2-x1,logvlad,logfrac);
						//fprintf(OUTFILE1,"%LG\n",xp[klauf]);
						fflush(OUTFILE1);
						Zeilencount++;
						time_temp = xp[klauf];
					}
				}
			}
			else if (((ausgabewunsch==OUTPUT_EVERYTHING)||(ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)) && ((x2-x1)<BahnPointSaveTime)){
				if((x2-timetemp)>=BahnPointSaveTime){
					printf(":");
					fflush(stdout);
					// Ausgabewerte berechnen
					BFeld(ystart[1],ystart[5],ystart[3], x2);
					vend    = sqrtl(fabsl(ystart[2]*ystart[2]+ystart[1]*ystart[1]*ystart[6]*ystart[6]+ystart[4]*ystart[4]));
					if(protneut == NEUTRON) 
						H = (M*gravconst*ystart[3]+0.5*M*vend*vend-mu_n*Bws)*1E9 ;    // mu_n negative for low-field seekers
					else if(protneut == PROTON)           // p			
						H= (0.5*m_p*vend*vend);           // Energie in eV for p		
					else if(protneut == ELECTRONS)  
						H= c_0*c_0  * M * (1/sqrtl(1-vend*vend/(c_0*c_0))-1);                                        // rel Energie in eV
					if (spinflipcheck == 2){
						if (vlad>1e-99) 
							logvlad=log10l(vlad);
						if (frac>1e-99) 
							logfrac=log10l(frac);
					}
					fprintf(OUTFILE1,"%d %d %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
									 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
									 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG \n",
									 iMC,protneut,x2,ystart[1],ystart[2],ystart[3],ystart[4],ystart[5],ystart[6],ystart[1]*cosl(ystart[5]),ystart[1]*sinl(ystart[5]),
									 sqrtl(ystart[2]*ystart[2]+ystart[6]*ystart[6]*ystart[1]*ystart[1]+ystart[4]*ystart[4]),H,Br,dBrdr,dBrdphi,dBrdz,Bphi,dBphidr,dBphidphi,dBphidz,
									 Bz,dBzdr,dBzdphi,dBzdz,Bws,Er,Ez,x2-x1,logvlad,logfrac);
					fflush(OUTFILE1);
					Zeilencount++;
					timetemp = x2;

				}
			}//Ende Ausgabe
	
			if (Zeilencount>40000)
			{
				fclose(OUTFILE1);
				Filecount++;
				ostringstream wholetrackfile;
				wholetrackfile << outpath << "/" << jobnumber << "track" << Filecount << ".out";
				OUTFILE1=fopen(wholetrackfile.str().c_str(),mode_w);
				fprintf(OUTFILE1,"Teilchen protneut t r drdt z dzdt phi dphidt x y "
								 "v H Matora Br dBrdr dBrdphi dBrdz Bphi dBphidr "
								 "dBphidphi dBphidz Bz dBzdr dBzdphi dBzdz Babs Polar Er Ez "
								 "timestep Bcheck logvlad logthumb\n");
				printf(" ##");
				printf(wholetrackfile.str().c_str());
				printf("## \n");
				fprintf(LOGSCR," ##");
				fprintf(LOGSCR,wholetrackfile.str().c_str());
				fprintf(LOGSCR,"## \n");
				Zeilencount=1;
			}

			if ((BFZeilencount>50000) && BruteForce)
			{
				fclose(BFLOG);
				BFFilecount++;
				ostringstream BFoutfile1;
				BFoutfile1 << outpath << "/" << jobnumber << "BF" << BFFilecount << ".out";
				if(!(BFLOG = fopen(BFoutfile1.str().c_str(),mode_w))) 
				{
					perror("fopen");
					exit(1);
				}
				fprintf(BFLOG,"t Babs Polar logPolar Ix Iy Iz Bx By Bz\n");
				printf(" ##");
				printf(BFoutfile1.str().c_str());
				printf("## \n");
				fprintf(LOGSCR," ##");
				fprintf(LOGSCR,BFoutfile1.str().c_str());
				fprintf(LOGSCR,"## \n");
				BFZeilencount=1;
			}
		
			Entkommen(ystart,  x2, H);  // check if particle should end!
			fflush(LOGSCR);
			
			//if(!((x2<=xend)&&((ystart[3]>=zmin)||(reflekt==1))&&(ystart[3]<=zmax)&& ((ystart[1]>=rmin)||(reflekt==1)) &&((ystart[1]<=rmax)||reflekt==1)&&(stopall==0))) 
			//	OutputState(ystart,1);
			
		}while (((x2-xstart)<=xend)&&(!stopall)); // end integration do - loop
		// END of loop for one partice
		
		
		timetemp = 0.0;         // set times for points on track to zero for new particle
		time_temp = 0.0;

		vend    = sqrtl(fabsl(ystart[2]*ystart[2]+ystart[1]*ystart[1]*ystart[6]*ystart[6]+ystart[4]*ystart[4]));
		long double phitemp = ((ystart[5])/conv);     // calculate end angle
		phiend  = fmodl(phitemp, 360.);   // in degree
		if (phiend<0)                    // from 0 to 360
			phiend=360.0 + phiend;

		if(protneut == NEUTRON)                // n
			H = (M*gravconst*ystart[3]+0.5*M*vend*vend-mu_n*Bws)*1E9 ;       // Energie in neV
		else if(protneut == PROTON)           // p			
				H= (0.5*m_p*vend*vend);           // Energie in eV for p		
		else if(protneut == ELECTRONS)           // p,e
			H= c_0*c_0  * M * (1/sqrtl(1-v_n*v_n/(c_0*c_0))-1);                                        // rel Energie in eV

		ausgabe(x2,ystart, vend, H);// Endwerte schreiben

		printf("Done!!\nBFFlipProb: %.17LG rend: %.17LG zend: %.17LG Eend: %.17LG Code: %i t: %.17LG\n",(BFflipprob),ystart[1],ystart[3],H,kennz,x2);
		fprintf(LOGSCR,"Done!!\nBFFlipProb: %.17LG rend: %.17LG zend: %.17LG Eend: %.17LG Code: %i t: %.17LG\n",(BFflipprob),ystart[1],ystart[3],H,kennz,x2);
		
		
		IncrementCodes(kennz);

		 
	} // end proton neutron calc
	
	 // calculation of electrons only through following field lines and magnetic mirror effect formula
	// TURNED OFF!!!
	if (protneut == 99)
	{
		BFeld(r_n, phi_n*conv, z_n, 500.0);
		Bre0=Br; Bphie0 = Bphi; Bze0= Bz;
		DetHit = CalcFluxLine(r_n,phi_n*conv,z_n, FluxStep);
		if (DetHit) {
			CritAngle = CalcCritAngle(r_n,phi_n*conv,z_n,Energie);
			IncidentAngle = CalcIncidentAngle (r_n, phi_n*conv,z_n, vr_n, vphi_n, vz_n, Bre0, Bphie0,  Bze0, Bemax);
			DetEnergy = Energie - Vflux;
		}else if (!DetHit) {
			CritAngle = 0;
			IncidentAngle = 0;
			DetEnergy = 0;
			ElecAngleB = 0;
		}
		fprintf(ENDLOG,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %i %.17LG %.17LG %.17LG "
					   "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
					   "%.17LG\n",
					   r_n,z_n,vr_n,vphi_n,vz_n,ElecAngleB,DetHit,CritAngle,Energie,Bre0,
					   Bze0,Be0,Bemax,ystart[1],ystart[3],Bws,Vflux,DetEnergy,IncidentAngle,BFeldSkal,
					   EFeldSkal);
	}

	
	if(protneut == BF_ONLY) // B-Feld Ausgabe
	{
		//Entkommen(ystart,  x2, H);
		BFeldSkal = 1.0;
		BFeld(r_n, phi_n, z_n, 500.0);
		Bws=sqrtl(Br*Br+Bz*Bz+Bphi*Bphi);
		//fprintf(ENDLOG,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG\n",
		// x:=r y=0 z:=z
		fprintf(ENDLOG,"%LG %G %LG %LG %LG %LG %G %G %LG \n",r_n*100,0.0,z_n*100,Br*1e4,Bphi*1e4,Bz*1e4,0.0,0.0,Bws);
		//printf("B");
		//printf("%LG %G %LG %LG %LG %LG %G %G %G\n", r_n*100.0,0.0,z_n*100.0,Br*1e4,Bphi*1e4,Bz*1e4,0.0,0.0,0.0);
		cout << "r= " << r_n << " z= " << z_n << " Br= " << Br << " T, Bz= " << Bz << " T"  << endl;
		
		// Ramp Heating Analysis
		long double EnTest;
		for (Energie = 0; Energie <= EnergieE; Energie++){
			EnTest = Energie*1.0e-9 - M*gravconst*z_n + mu_n * Bws;
			if (EnTest >= 0){
				iii = (int) Energie;
				// add the volume segment to the volume that is accessible to a neutron with energy Energie
				VolumeB[iii] = VolumeB[iii] + pi * dz_n * ((r_n+0.5*dr_n)*(r_n+0.5*dr_n) - (r_n-0.5*dr_n)*(r_n-0.5*dr_n));
			}
		}

	}// Ende B-Feld Ausgabe


	// sorgt dafr, dass bei MonteCarlo die for Schleifen nur 1mal durchlaufen werden
	if (MonteCarlo){
		Energie=EnergieE+1;
		z_n=z_ne+1;
		r_n=r_ne+1;
		alpha=alphae+1;
		gammaa=gammae+1;			
		phi_n=phie+1;
	}
}

void BruteForceIntegration(){
	// Array fr BruteForce Integration wird gebildet
	int klauf, klaufstart;
	BFdxsav=5e-7;  // timestep for intermediate output
	bool BFPolmin;
	
	BFPolmin = (BFBmin < BFTargetB);
	// if at last step there was BFintegration => true
	// if there was no integration => false
	
	BFBmin = 10; // set to a value higher than all real bws values
	
	for (klauf=1;klauf<=kount;klauf++)
	{    // go through intermediate values
		BFBmin = min(BFBmin, Bp[13][klauf]); // write out smallest value of Bws
	}
	
	if(BFBmin<BFBminmem)    // accumulate the smallest value Babs for which BF integration is done
		BFBminmem=BFBmin;
	
	if ((BFBmin>BFTargetB)&&(BFPolmin))
	{    // output of polarisation after BF int completed
		BFsurvprob = (BFpol+0.5) * BFsurvprob;
		BFflipprob = 1-BFsurvprob;		// update spinflip probability after passing low field region
		// flip the spin with the probability BFflipprob
		if (flipspin){
			// if (rando < 0.5) { // for testing, remove 
			if (mt_get_double(v_mt_state) < (1-(BFpol+0.5))) 
			{
				hfs *= -1;
				mu_n=hfs * mu_nSI / ele_e;
				mumB= mu_n/M;
				NSF++; 
				printf("\n BFpol: %LG  The spin has flipped! Number of flips: %i\n",BFpol,NSF);			
			}
		}
		//fprintf(TESTLOG,"%.17LG %.17LG %.17LG\n",BFtime[offset],BFlogpol,(1-BFflipprob));
	}

	if (BFBmin>BFTargetB)
	{   // => no BF integration will take place, so for the next one, the polarisation vector will be set parallel to magnetic field
		firstint = 1;
	}

	if (BFBmin<BFTargetB)
	{   // check if this value is worth for Bloch integration 
		klaufstart = 1;// start with index 1												
		
		if (offset>0)
		{// start with 2 if there are already values in the arrays, because otherwise two values could be the same
			klaufstart = 2;
		}			
		/*
		// output start position and velocities of BruteForce
		if(klaufstart == 1)
		{
			gammaend=atan2l(sqrtl(powl(yp[2][2],2)+powl(yp[1][2]*yp[6][2],2)),yp[4][2])/conv;
			fprintf(LOGSCR,"\n r:%.17LG phi:%.17LG z:%.17LG H:%.17LG alpha:%.17LG gamma:%.17LG \n ",
							yp[1][2],yp[5][2]/conv,yp[3][2],(M*gravconst*yp[3][2]+0.5*M*fabsl(yp[2][2]*yp[2][2]+yp[1][2]*yp[1][2]*yp[6][2]*yp[6][2]+yp[4][2]*yp[4][2])-mu_n*Bp[13][2])*1E9, atanl(yp[6][2]*yp[1][2]/yp[2][2])/conv, gammaend);
		}*/
		
													
		for (klauf=klaufstart;klauf<=kount;klauf++)
		{    // build array for Bloch integration
			//if(xp[klauf]<=BFtime[offset+klauf-(klaufstart-1)])
			//	klauf++;
			BFtime[offset+klauf-(klaufstart-1)]=xp[klauf];
			// transform cylindrical into kartesian lokal coordinates
			CylKartCoord(Bp[1][klauf], Bp[5][klauf], Bp[9][klauf],  yp[5][klauf], &BFBxcoor, &BFBycoor, &BFBzcoor);
			BFField[1][offset+klauf-(klaufstart-1)]=BFBxcoor;
			BFField[2][offset+klauf-(klaufstart-1)]=BFBycoor;
			BFField[3][offset+klauf-(klaufstart-1)]=BFBzcoor;
			BFField[4][offset+klauf-(klaufstart-1)]=yp[1][klauf];    // r
			BFField[5][offset+klauf-(klaufstart-1)]=yp[3][klauf];    // z
			
			// exemplary spinflip investigation
			//BFField[1][offset+klauf-(klaufstart-1)]=10787.1388395727981340375148505*xp[klauf]*xp[klauf]*xp[klauf]-379485.201099151545103425695947*xp[klauf]*xp[klauf]+4425265.48802544343656564311135*xp[klauf]-17089608.2309508308489746319934;
			//BFField[2][offset+klauf-(klaufstart-1)]=16511.3272584589837067791527764*xp[klauf]*xp[klauf]*xp[klauf]-594543.289553430286823139813885*xp[klauf]*xp[klauf]+7118405.94609221980769652287442*xp[klauf]-28331085.0061183671391854839742;
			//BFField[3][offset+klauf-(klaufstart-1)]=1944281.06055634049312257407394*xp[klauf]*xp[klauf]*xp[klauf]-69404618.8278242196709266829602*xp[klauf]*xp[klauf]+822964791.430415938909499801820*xp[klauf]-3239972252.28819208686600696119;
			
																	// Spinflipper test
			//BFField[1][offset+klauf-(klaufstart-1)]= B1 * cosl(xp[klauf]*gamma_n * 1.0e-3);
			//BFField[2][offset+klauf-(klaufstart-1)]= B1 * sinl(xp[klauf]*gamma_n * 1.0e-3);
			//BFField[3][offset+klauf-(klaufstart-1)]= 1.0e-3;
		}
		offset = offset + kount-(klaufstart-1);
	}


	
	// Perform integration
	if ((((BFBmin>=BFTargetB)||((x2-xstart)>=xend))&&(offset>=10))||(offset>=2000))
	{                      
		if (firstint){
			BFBws[1] = sqrtl(powl(BFField[1][1],2)+powl(BFField[2][1],2)+powl(BFField[3][1],2));
			I_n[1]= (BFField[1][1]/BFBws[1])*0.5;
			I_n[2]= (BFField[2][1]/BFBws[1])*0.5;
			I_n[3]= (BFField[3][1]/BFBws[1])*0.5;
			
			//printf("Bvector before %LG %LG %LG Babs %LG \n",BFField[1][1],BFField[2][1],BFField[3][1],BFBws[1]);
			//printf("Spinvector before %LG %LG %LG  \n",I_n[1],I_n[2],I_n[3]);													
			printf(" BF starttime %.6LG, BFflipprob before %LG, offset %i, delx_n %LG, Bmin %LG, |I| before %LG ",BFtime[1], BFflipprob, offset, delx_n, BFBmin,sqrtl(powl(I_n[1],2)+powl(I_n[2],2)+powl(I_n[3],2)));
		}													
		
		/*
		for(int itmp=1;itmp<=offset;itmp++)
		{
			printf("%d t %LG Bx %LG By %LG Bz %LG \n",itmp,BFtime[itmp],BFField[1][itmp],BFField[2][itmp],BFField[3][itmp],BFBws[1]);
		}
		*/
		
		//    (eingangsvektor,nvar,xbeg,xend,rel.accur,begstepsize,hmin,&nok,&nbad,derivs,bsstep)
		(*BFodeintrk)(I_n,3,BFtime[1],BFtime[offset],1e-13,1e-5,0,&nok,&nbad,BFderivs,BFrkqs);
		//(*BFodeintrk) (I_n,3,BFtime[1],2e-7,1e-13,1e-5,0,&nok,&nbad,BFderivs,BFrkqs);
		printf("|I| after %LG, BF endtime %LG, intsteps taken %d, Babsmin %LG\n",sqrtl(powl(I_n[1],2)+powl(I_n[2],2)+powl(I_n[3],2)), BFtime[offset], nok, BFBminmem);
		firstint = 0;                 // after BF step, it will not be the first BF step any more...	
		BFBminmem=10;
		// calculate polarisation at end of step BFpol = (I_n*B/|B|) in [-1/2,1/2]
		BFpol = 	(BFyp[1][BFkount]* BFypFields[1][BFkount] + 
						BFyp[2][BFkount]* BFypFields[2][BFkount] + 
						BFyp[3][BFkount]* BFypFields[3][BFkount])
							/sqrtl(powl(BFypFields[1][BFkount],2) + 
										powl(BFypFields[2][BFkount],2) + 
										powl(BFypFields[3][BFkount],2));
		//printf("Bvector after %LG %LG %LG Babs %LG \n",BFypFields[1][BFkount],BFypFields[2][BFkount],BFypFields[3][BFkount],sqrtl(powl(BFypFields[1][BFkount],2) + 
		//								powl(BFypFields[2][BFkount],2) + 
		//								powl(BFypFields[3][BFkount],2)));
		//printf("Spinvector after %LG %LG %LG  BFpol %LG \n",I_n[1],I_n[2],I_n[3],BFpol);		
		
		
		
		// only print interesting values where deviation is > 1e-6 
		//if ((1-BFpol)>1.0e-6){
		if((ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)||(ausgabewunsch==OUTPUT_ENDPOINTSandSPIN))
		{
			for (int BFcount=2; BFcount<=(BFkount);BFcount++)
			{
				BFBws[BFcount] =  sqrtl(BFypFields[1][BFcount]*BFypFields[1][BFcount] + BFypFields[2][BFcount]*BFypFields[2][BFcount] + BFypFields[3][BFcount]*BFypFields[3][BFcount]);
				BFpol = (BFyp[1][BFcount]* BFypFields[1][BFcount] + BFyp[2][BFcount]* BFypFields[2][BFcount] + BFyp[3][BFcount]* BFypFields[3][BFcount])/sqrtl(BFypFields[1][BFcount]*BFypFields[1][BFcount] + BFypFields[2][BFcount]*BFypFields[2][BFcount] + BFypFields[3][BFcount]*BFypFields[3][BFcount]);
				long double BFlogpol;
				if (BFpol<0.5) 
					BFlogpol = log10l(0.5-BFpol);
				else if (BFpol==0.5) 
					BFlogpol = 0.0;
				fprintf(BFLOG,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG \n",
								BFxp[BFcount],BFBws[BFcount],BFpol,BFlogpol,BFyp[1][BFcount]*2,BFyp[2][BFcount]*2,BFyp[3][BFcount]*2,BFypFields[1][BFcount]/BFBws[BFcount],BFypFields[2][BFcount]/BFBws[BFcount],BFypFields[3][BFcount]/BFBws[BFcount], BFyp[1][BFcount]*2-BFypFields[1][BFcount]/BFBws[BFcount],
				BFyp[2][BFcount]*2-BFypFields[2][BFcount]/BFBws[BFcount],
				BFyp[3][BFcount]*2-BFypFields[3][BFcount]/BFBws[BFcount]);
				
				BFZeilencount++;
			}  
		}
		//}
		
		offset = 0;                       // after BF step start with new array
													
	}			
// END brute force integration
}
