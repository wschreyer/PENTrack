
/*********************************************
pnTracker

**********************************************/
#include "main.h"
 

// files for in/output + paths
FILE *OUTFILE1 = NULL, *BFLOG = NULL, *SNAP = NULL;
TFile *treefile;
TNtupleD *endtree, *tracktree;
string inpath, outpath;	// "in" and "out" directories
char mode_r[2] = "r",mode_rw[3] = "rw",mode_w[2] = "w"; // modes for fopen()

// physical constants
long double ele_e=1.602176487E-19, Qm0=1.602176487E-19/1.672621637E-27;      //elementary charge in SI, charge/proton mass 																					
long double gravconst=9.80665, conv=0.01745329251, mu0= 1.25663706144e-6;      //g [m/s^2], Pi/180, permeability,
long double m_n=1.674927211E-27/1.602176487E-19, pi=3.141592655359;  //neutron mass (eV/c^2), neutron magnetic moment (in eV/T), Pi
long double m_p=1.672621637E-27/1.602176487E-19;        //proton mass (eV/c^2), tempmass
long double m_e = 9.10938215e-31/1.602176487E-19, c_0 = 299792458; //electron mass (eV/c^2), lightspeed
long double hquer=1.05457266e-34, mu_nSI=0.96623641e-26;          // Neutron magn Mom (in J/T)
long double gamma_n = 1.83247185e8;            
long double tau=885.7;              // magn. moment of neutron/mass,  neutron lifetime [s]
long double lengthconv = 0.01 , Bconv = 1e-4, Econv = 1e2;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 
												// cm => m,  Gauss => Tesla,   V/cm => V/m    Bconv temporr von 1e-4 auf 1 gesetzt

// misc configurations
int MonteCarloAnzahl=1;   // user choice to use MC or not, number of particles for MC simulation
int reflekt=0, bfeldwahl, protneut, Racetracks=2;       //user choice for reflecting walls, B-field, prot or neutrons, experiment mode
int polarisation=0, polarisationsave=0, ausgabewunsch=5, ausgabewunschsave; // user choice for polarisation of neutrons und Ausgabewunsch
int WriteTree = 0;
int snapshot=0, snapshotsdone=0;  // make snapshots of the neutron at specified times
//vector<double> snapshots;
set<int> snapshots;
long double Ibar= 2250.;                // B-field strength, current through rod
int runge;                            // Runge-Kutta or Bulirsch-Stoer?  set to Runge right now!!!
int diffuse; //diffuse reflection switch
int neutdist;

// Fields
list<string> fieldvaltab; // filenames of field table files
long double dBrdr, dBrdz, dBrdphi=0.0, dBphidr, dBphidz, dBphidphi=0.0, dBzdr, dBzdz, dBzdphi=0.0;
long double Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field
long double BFeldSkal = 1.0, EFeldSkal = 1.0, BFeldSkalGlobal = 1.0;          // parameter to scale the magnetic field for ramping, scale electric field, Global: also scale ramping etc...
long double DiceRodField=0, RodFieldMultiplicator;
long double BCutPlanePoint[3], BCutPlaneNormalAlpha, BCutPlaneNormalGamma, BCutPlaneSampleDist;	// plane for B field cut
int BCutPlaneSampleCount;
long double Emin_n = 1e30;  // minimum energy of neutrons in the B-field

// particles
long double M, H, mu_n, mumB;                               // mass, total energy, magnetic moment, moment per mass
long double ystart[7], xstart = 0;       //initial and intermediate values of y[7], start time
long double  x1, x2;                         // start and endtime handed over to integrator
int iMC;                             //  counter for MonteCarloSim
long double trajlengthsum;			// trajectory length
long double Hstart, Hend,Hmax;     //maximum energy
long double NeutEnergie;	//relativistic gamma factor, Neutron energy

struct decayinfo decay; //containing all data from neutron decay for the emerging proton and electron

// inital values of particle
long double EnergieS, EnergieE, Energie;    //initial energy range
long double r_n, phi_n, z_n, v_n;                //initial particle coordinates
long double alpha, gammaa, hmin;                  //initial angle to x-Achse, zo z-Achse, Schrittweite
long double alphas, gammas;   //initial values from
long double alphae, gammae;   //initial values to
long double delx;                            // initial timestep for the integrator

long double delx_n=0.0;
int stopall=0, Feldcount=0;                            //  if stopall=1: stop particle

initial nini, pini, eini; 	// one 'initial' for each particle type

// final values of particle
int kennz;                                  // ending code
long double vend, gammaend, alphaend, phiend, xend;    //endvalues for particle
long int kennz0[3]={0},kennz1[3]={0},kennz2[3]={0},kennz3[3]={0},kennz4[3]={0},kennz5[3]={0},kennz6[3]={0},kennz7[3]={0},kennz8[3]={0},kennz9[3]={0},kennz10[3]={0},kennz11[3]={0},kennz12[3]={0},kennz99[3]={0},nrefl; // Counter for the particle codes

// integrator params
long double  eps, h1;  // desired accuracy in ODEINT: normal, for polarisation and phase, trial time step
int nvar, nok, nbad;                            // in ODEINT: number of variables in Derivs, good und bad steps

//set the duration of the experiment
long double FillingTime = 0;										// filling time, entrance open
long double CleaningTime = 0;                        // cleaning without field
long double RampUpTime = 0;                          // ramping up coils
long double FullFieldTime = 1000;                       // storing in full field
long double RampDownTime = 5;                        // ramping down coils
long double EmptyingTime = 0;                        // emptying without field
long double StorageTime = 1500.0;                     // time when ramping down shall start, if xend > storage time, let neutron decay
int ExpPhase = 0;  															// current experiment phase
int ffBruteForce,ffreflekt,ffspinflipcheck;  // fullfield
int ruBruteForce,rureflekt,ruspinflipcheck; // rampup
int rdBruteForce,rdreflekt,rdspinflipcheck;  // rampdown
int fiBruteForce,fireflekt,fispinflipcheck;  // filling
int coBruteForce,coreflekt,cospinflipcheck;  // counting UCN
int clBruteForce,clreflekt,clspinflipcheck;  // cleaning time
int jobnumber;
unsigned long int monthinmilliseconds; // RandomSeed

// Spintracking
int spinflipcheck = 0;                          // user choice to check adiabacity
long double vlad=0.0, vladtotal = 1.0, frac;                   // adiabacity after Vladimirsky
long double vladmax=0.0; // maximum values of spinflip prob
long double thumbmax=0.0;

// file output
long Zeilencount;
int Filecount=1, p;                                   // counts the output files, counter for last line written in outs, random generator temp value
long double BahnPointSaveTime = 5.0e-7;               // default 2e-7; 0=1e-19 not changed at run time, time between two lines written in outs

// variables for BruteForce Bloch integration BEGIN
long double *BFtime=NULL, **BFField=NULL;   // time, Bx, By, Bz, r, z array
int offset=0, BFkount, BFindex = 3;			// counter in BFarray, offset of BFarray, maximum index of intermediate values , index in BFarray
long double BFpol, *BFBws=NULL;                    // BFpolarisation
long double BFBmin = 10.0, BFBminmem=10.0, BFTargetB=0.1;     // smallest value of Babs during step, memorize value during value accumulation, Babs < BFTargetB => integrate,
unsigned short int BruteForce = 0, firstint = 1, flipspin=1;  // enable BruteForce?,
long double I_n[4], **BFypFields=NULL;        // Spinvector, intermediate field values in BFodeint
long BFZeilencount; int BFFilecount=1;                  // to control output filename of BF
long double BFflipprob = 0.0, BFsurvprob=1.0;                // spinflip probability, survival (non-flip) probability

//for the output of intermediate steps
//#define KMDEF 1000
int BFNrIntermediate=BFKMDEF;    // number of steps for intermediate output in the case of 
int BFkmax = BFKMDEF;                                         // number of steps for intermediate output
long double nintcalls=0, ntotalsteps=0;     // counters to determine average steps per integrator call
int kount, hfs, NSF = 0;                                            // counter for intermediate output steps, highfieldseeker: +1 yes, -1 lowfieldseeker, number of spin flips
long double *xp=NULL,**yp=NULL, *BFxp=NULL, **BFyp=NULL, BFdxsav=0;          // Arrays for intermediate output
long double **Bp=NULL,**Ep=NULL;                            // Arrays for intermediate output
// END output of intermediate steps 

// incorporate B-fieldoszillations into the code
int FieldOscillation = 0;        // turn field oscillation on if 1
long double OscillationFraction = 1e-4, OscillationFrequency = 1;    // Frequency in Hz

// time statistics
float InitTime = 0, ReflectionTime = 0, IntegratorTime = 0, DiceTime = 0;

mt_state_t *v_mt_state = NULL; //mersenne twister state var


/*
 * handler:		catch_alarm(int sig)
 * 				terminates a program if a specific signal occurs
 * type:		void
 * var sig:		signalnumber which called the handler; to get the right number
 * 				for corresponding signals have a look <man signal.h>.
 * 				e.g: <SIGFPE> is connected to number 8
 * return: 		noreturn
 */
void catch_alarm (int sig){
	Log("Program was terminated, because Signal %i occured", sig);
	exit(1);
}

/*	
 *	function:	int file_rights(char path)
 *	It checks the rights of a file to be sure one is allowed to
 *	write on it. If it fails, the program will be terminated.
 *
 *	ATTENTION: The whole path is required
 */
int file_rights(char *path){
	
	struct stat attribut;
	
	stat(path,&attribut);
	//if(!(attribut.st_mode & S_IWUSR){ //Windows-Bedingung
	if(!(attribut.st_mode & S_IWUSR) || !(attribut.st_mode & S_IWGRP) || !(attribut.st_mode & S_IWOTH)){ //Linux-Bedingung
		Log("Die Datei <%s> besitzt keine Schreibrechte", path);
		exit(1);
	}
	return 0;
}


// uebergabe: jobnumber inpath outpath                 paths without last slash
int main(int argc, char **argv){
	//Initialize signal-analizing
	signal (SIGUSR1, catch_alarm);
	signal (SIGUSR2, catch_alarm);
	signal (SIGXCPU, catch_alarm);   
	
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
	hmin= 0;       // minimum stepsize for runge kutta
	
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
	ConfigInit();
	OpenFiles(argc, argv);	// Open .in and .out files and write headers
		
	Log(
	" ################################################################\n"
	" ###                 Welcome to PNTracker,                    ###\n"
	" ###     the tracking program for neutrons and protons        ###\n"
	" ################################################################\n");

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
	xp=dvector(1,kmax);
	yp=dmatrix(1,6,1,kmax);         
	Bp=dmatrix(1,13,1,kmax);
	Ep=dmatrix(1,2,1,kmax);
	
	LoadGeometry();	// read STL-files		
	
	PrepareFields();	// read fieldval.tab or coils.cond


	Startbed(1); // read in starting values of particle
	
	PrintConfig();
	
	PrepareParticle(); // setting values depending on particle type
	
	InitTime = (1.*clock())/CLOCKS_PER_SEC;

	decay.counter = 0; // no neutron had decayed yet (no chance)

	iMC = 1;
	do
	{      
		IntegrateParticle();
		
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
	
	Log("Integrator used (1 Bulirsch Stoer, 2 Runge Kutta): %d \n", runge);
	Log("The integrator was called: %LF times with %LF internal steps on average. \n", nintcalls,ntotalsteps/nintcalls);
	float SimulationTime = (1.*clock())/CLOCKS_PER_SEC - InitTime;
	IntegratorTime -= ReflectionTime;
	Log("Init: %.2fs, Simulation: %.2fs, Integrator: %.2fs (%.2f%%), Reflection: %.2fs (%.2f%%), Dicing: %.4fs\n",
			InitTime, SimulationTime, IntegratorTime, IntegratorTime*100/SimulationTime, ReflectionTime, ReflectionTime*100/SimulationTime, DiceTime);
	Log("That's it... Have a nice day!\n");
	
	FreeFields();
	
	if (WriteTree){
	//	endtree->Print();
	//	tracktree->Print();
		treefile->Write();
		treefile->Close();
	}

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
						decay.error = 0;
						break;		
		case PROTON:	M = m_p;			// [eV/c^2]
						Qm0 = 1.0/M;
						h1 = 1e-8;		// guess for the first step size of runge kutta
						BahnPointSaveTime = 1e-8;
						reflekt = 0;
						break;		
		case BF_ONLY:	PrintBField();
						exit(0);
						break;		
		case BF_CUT:	PrintBFieldCut();
						exit(0);
						break;		
		case ELECTRONS:	M=m_e;			// [eV/c^2]
						Qm0 = -1.0/M;
						h1 = 2e-10;		// guess for the first step size of runge kutta
						BahnPointSaveTime = 1e-11;
						reflekt = 0;
						break;
	}
}
//======== end of PrepareParticle =====================================================================================

void derivs(long double x, long double *y, long double *dydx){
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
		dydx[4]= mumB*dBdz-gravconst; // [m/s^2]
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+mumB/y[1]*dBdphi/y[1];
		
	}
	else if(protneut == PROTON)
	{ // equations for the proton
		dydx[1]= y[2];
		dydx[2]= y[1]*(y[6]*y[6])+Qm0*(Bz*y[1]*y[6]-Bphi*y[4]+Er);
		dydx[3]= y[4];
		dydx[4]= Qm0*(Ez+Bphi*y[2]-y[1]*Br*y[6])-gravconst; 
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+Qm0*(Ephi+Br*y[4]-Bz*y[2])/y[1];
	}
	else if(protneut == ELECTRONS)
	{ // equations for the electron (relativistic mass)
		Qm0 = -1.0/M*sqrtl(1-(y[2]*y[2]+y[1]*y[1]*y[6]*y[6]+y[4]*y[4])/(c_0*c_0));
		dydx[1]= y[2];
		dydx[2]= y[1]*(y[6]*y[6])+Qm0*(Bz*y[1]*y[6]-Bphi*y[4]+Er);
		dydx[3]= y[4];
		dydx[4]= Qm0*(Ez+Bphi*y[2]-y[1]*Br*y[6])-gravconst;
		dydx[5]= y[6];
		dydx[6]= -2*y[2]*y[6]/y[1]+Qm0*(Ephi+Br*y[4]-Bz*y[2])/y[1];
	}
	return;
}


void IntegrateParticle(){
	int perc=0;   // percentage of particle done counter
		// reset some values for new particle
	stopall=0;
	kennz=KENNZAHL_UNKNOWN; // not categorized yet									
	// initial values for Brute-Force Spinintegration 
		BFpol = 0.5;
		I_n[3]=0.5; I_n[2]=I_n[1]=0;
		offset = 0;
		BFsurvprob = 1.0;
		BFflipprob = 0.0;
	Hmax=0.0;
	x2=x1= 0.;     //set time to zero
	snapshotsdone=0;
	
	timeval dicestart, diceend;
	gettimeofday(&dicestart, NULL);
	MCStartwerte(delx);   // MonteCarlo Startwert und Lebensdauer fr Teilchen festlegen
	gettimeofday(&diceend, NULL);
	DiceTime += diceend.tv_sec - dicestart.tv_sec + float(diceend.tv_usec - dicestart.tv_usec)/1e6;
	x2=xstart;    // set time starting value										
	
	
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
		// output in the ENDLOG
		ausgabe(0.0, ystart, 0.0, 0.0);
		// reset 'noparticle' to false 		
		noparticle = false;
		return;
	}
	//-------- Finished ----------------------------------------------------------------------------------------------------
	
	Log("\nRodFieldMultiplicator: %.17LG\n",RodFieldMultiplicator);
	Log("Feldcount = %i\n\n",Feldcount);

	trajlengthsum = 0;
	nrefl=0;
	kennz=KENNZAHL_UNKNOWN;
	stopall=0;
	
	Feldcount=0;

	if (protneut == NEUTRON)
	{// Spinflipwahrscheinlichkeiten zurcksetzen
		thumbmax=0.0;
		vladtotal = 1.0;
		vladmax=0.0;									
		BFeld(r_n,phi_n,z_n, 0.0);							
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
			Log("\nEkin: %.17LG  smaller than Zero!!! \n",Ekin);
			//iMC--;
			
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
		long double gammarel = Energie/m_e/c_0/c_0 + 1;
		v_n = c_0 * sqrtl(1-(1/(gammarel*gammarel)));     // relatistic velocity from E_kin in eV
		//cout << "Energie: " << Energie << " m_e: " << m_e << " v_n: " << v_n << endl;
		//sleep(5);										
	}
	
	long double projz= cosl(gammaa);  // projection of velocity on z-axis
	long double vtemp= v_n*sinl(gammaa);  // projection of velocity on x-y plane
	ystart[1]= r_n;           // fill array for ODEint integrator
	ystart[2]= vtemp*cosl(alpha-phi_n);
	ystart[3]= z_n;
	ystart[4]= v_n*projz;
	ystart[5]= phi_n;
	if (r_n!=0.) 
		ystart[6]= vtemp*sinl(alpha-phi_n)/r_n;
	else if (r_n==0.)
		ystart[6]= 0.;
	
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
		
		Log("Teilchennummer: %i\n",iMC);
		if(decay.on == 2)
			Log("Teilchensorte : %i\n", protneut);
		Log("r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG\n",
			r_n, phi_n/conv, z_n, v_n, alpha/conv, gammaa/conv, H, xend);

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
														
			// put phi (ystart[5]) back to [-2Pi,2Pi]
			if(ystart[5]>(2.0*pi))
				ystart[5]=ystart[5]-2*pi;
			if(ystart[5]<(-2.0*pi))
				ystart[5]=ystart[5]+2*pi;
			
			timeval intstart, intend;
			gettimeofday(&intstart, NULL);	
			//###################### Integrationsroutine #####################
			if (runge==2)  odeint(ystart,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,&derivs,&rkqs);           // runge kutta step
			else if (runge==1)  odeint(ystart,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,&derivs,&bsstep);        // bulirsch stoer step
			// (ystart: input vector | nvar: number of variables | x1, x2: start and end time | eps: precision to be achieved | h1: guess for first stepsize | hmin: mininum stepsize | nok,nbad: number of good and bad steps taken | derivs: function for differential equation to be integrated, rkqs, bsstep: integrator to be used (runge kutta, bulirsch stoer) )
			//###################### Integrationsroutine #####################
			gettimeofday(&intend, NULL);
			IntegratorTime += intend.tv_sec - intstart.tv_sec + float(intend.tv_usec - intstart.tv_usec)/1e6;
			nintcalls++;
		
			// spin flip properties according to Vladimirsky and thumbrule
			for (int i = 1; i < kount; i++){
				if (spinflipcheck == 2){   
					// y[6]: phidot
					vlad = vladimirsky(yp[1][i], Bp[1][i], Bp[5][i], Bp[9][i], 
									   Bp[2][i], Bp[3][i], Bp[4][i], Bp[6][i], Bp[7][i], Bp[8][i], Bp[10][i], Bp[11][i], Bp[12][i],
									   yp[2][i], yp[6][i], yp[4][i]);
					frac = thumbrule(Bp[1][i], Bp[5][i], Bp[9][i], 
									 Bp[2][i], Bp[3][i], Bp[4][i], Bp[6][i], Bp[7][i], Bp[8][i], Bp[10][i], Bp[11][i], Bp[12][i],
									 yp[2][i], yp[6][i], yp[4][i]);
					if (vlad > 1e-99){
						vladtotal = vladtotal * (1-vlad);
						if (vladtotal < 0.9999)
							Log(" VladShit (%LG) at t= %.17LG\n",vladtotal,xp[i]);
					}
					if ((vlad > vladmax)&&(vlad > 1e-99))
						vladmax=log10l(vlad);
					if ((frac > thumbmax)&&(frac > 1e-99))
						thumbmax=log10l(frac);
				}
			}
								
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
			
			PrintIntegrationStep(timetemp);
			
			// take snapshots of neutrons at certain times
			if(snapshot==1  && protneut==1)
				Snapshooter(x2,ystart, H);
			
		}while (((x2-xstart)<=xend) && (!stopall) && (x2 <= StorageTime)); // end integration do - loop
		// END of loop for one partice
		
		vend    = sqrtl(fabsl(ystart[2]*ystart[2]+ystart[1]*ystart[1]*ystart[6]*ystart[6]+ystart[4]*ystart[4]));
		phiend = fmod(ystart[5], 2*pi);

		if(protneut == NEUTRON)                // n
			H = (M*gravconst*ystart[3]+0.5*M*vend*vend-mu_n*Bws)*1E9 ;       // Energie in neV
		else if(protneut == PROTON)           // p			
				H= (0.5*m_p*vend*vend);           // Energie in eV for p		
		else if(protneut == ELECTRONS)           // p,e
			H= c_0*c_0  * M * (1/sqrtl(1-v_n*v_n/(c_0*c_0))-1);                                        // rel Energie in eV

		ausgabe(x2,ystart, vend, H);// Endwerte schreiben

		Log("Done!!\nBFFlipProb: %.17LG rend: %.17LG zend: %.17LG Eend: %.17LG Code: %i t: %.17LG\n",(BFflipprob),ystart[1],ystart[3],H,kennz,x2);		
		
		IncrementCodes(kennz);
		 
	} // end proton neutron calc
}

void BruteForceIntegration(){
	// Array fr BruteForce Integration wird gebildet
	int klauf, klaufstart;
	BFdxsav=5e-7;  // timestep for intermediate output
	bool BFPolmin = (BFBmin < BFTargetB);
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
			gammaend=atan2l(sqrtl(powl(yp[2][2],2)+powl(yp[1][2]*yp[6][2],2)),yp[4][2]);
			fprintf(LOGSCR,"\n r:%.17LG phi:%.17LG z:%.17LG H:%.17LG alpha:%.17LG gamma:%.17LG \n ",
							yp[1][2],yp[5][2]/conv,yp[3][2],(M*gravconst*yp[3][2]+0.5*M*fabsl(yp[2][2]*yp[2][2]+yp[1][2]*yp[1][2]*yp[6][2]*yp[6][2]+yp[4][2]*yp[4][2])-mu_n*Bp[13][2])*1E9, atanl(yp[6][2]*yp[1][2]/yp[2][2])/conv, gammaend/conv);
		}*/
		
													
		for (klauf=klaufstart;klauf<=kount;klauf++)
		{    // build array for Bloch integration
			//if(xp[klauf]<=BFtime[offset+klauf-(klaufstart-1)])
			//	klauf++;
			BFtime[offset+klauf-(klaufstart-1)]=xp[klauf];
			// transform cylindrical into kartesian lokal coordinates
			long double BFBxcoor, BFBycoor, BFBzcoor;
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

//Ausgabe der Zwischenwerte aus odeint
void PrintIntegrationStep(long double &timetemp){
	long double logvlad = 0.0, logfrac = 0.0;
	if ((ausgabewunsch==OUTPUT_EVERYTHING)||(ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)){
		for (int klauf=1;klauf<kount;klauf++){
			if ((xp[klauf]-timetemp)>=BahnPointSaveTime)
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
				fprintf(OUTFILE1,"%d %d %d %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
								 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
								 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %d \n",
								 iMC,protneut,polarisation,xp[klauf],yp[1][klauf],yp[2][klauf],yp[3][klauf],yp[4][klauf],yp[5][klauf],yp[6][klauf],yp[1][klauf]*cosl(yp[5][klauf]),yp[1][klauf]*sinl(yp[5][klauf]),
								 vend,H,Bp[1][klauf],Bp[2][klauf],Bp[3][klauf],Bp[4][klauf],Bp[5][klauf],Bp[6][klauf], Bp[7][klauf],Bp[8][klauf],
								 Bp[9][klauf],Bp[10][klauf],Bp[11][klauf],Bp[12][klauf],Bp[13][klauf],Ep[1][klauf],Ep[2][klauf],x2-x1,logvlad,logfrac,ExpPhase);
				//fprintf(OUTFILE1,"%LG\n",xp[klauf]);
				
				if (WriteTree){
					double outvars[] = {iMC, protneut, polarisation, xp[klauf], yp[1][klauf], yp[2][klauf], yp[3][klauf], yp[4][klauf], yp[5][klauf], yp[6][klauf], yp[1][klauf]*cosl(yp[5][klauf]), yp[1][klauf]*sinl(yp[5][klauf]),
										vend, H, Bp[1][klauf], Bp[2][klauf], Bp[3][klauf], Bp[4][klauf], Bp[5][klauf], Bp[6][klauf],  Bp[7][klauf], Bp[8][klauf],
										Bp[9][klauf], Bp[10][klauf], Bp[11][klauf], Bp[12][klauf], Bp[13][klauf], Ep[1][klauf], Ep[2][klauf], x2-x1, logvlad, logfrac, ExpPhase};
					tracktree->Fill(outvars);
				}
				
				fflush(OUTFILE1);
				Zeilencount++;
				timetemp = xp[klauf];
			}
		}
	}

	if (Zeilencount>40000)
	{
		fclose(OUTFILE1);
		Filecount++;
		ostringstream wholetrackfile;
		wholetrackfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "track" << setw(3) << Filecount << setw(0) << ".out";
		OUTFILE1=fopen(wholetrackfile.str().c_str(),mode_w);
		fprintf(OUTFILE1,"Teilchen protneut polarisation t r drdt z dzdt phi dphidt x y "
						 "v H Matora Br dBrdr dBrdphi dBrdz Bphi dBphidr "
						 "dBphidphi dBphidz Bz dBzdr dBzdphi dBzdz Babs Polar Er Ez "
						 "timestep Bcheck logvlad logthumb ExpPhase \n");
		Log(" ##");
		Log(wholetrackfile.str().c_str());
		Log("## \n");
		Zeilencount=1;
	}

	if ((BFZeilencount>50000) && BruteForce)
	{
		fclose(BFLOG);
		BFFilecount++;
		ostringstream BFoutfile1;
		BFoutfile1 << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "BF" << setw(3) << BFFilecount << setw(0) << ".out";
		if(!(BFLOG = fopen(BFoutfile1.str().c_str(),mode_w))) 
		{
			perror("fopen");
			exit(1);
		}
		fprintf(BFLOG,"t Babs Polar logPolar Ix Iy Iz Bx By Bz\n");
		Log(" ##");
		Log(BFoutfile1.str().c_str());
		Log("## \n");
		BFZeilencount=1;
	}
}
