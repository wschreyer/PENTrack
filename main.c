
/*********************************************
pnTracker

**********************************************/
#include <iostream>
#include <cmath>
#include <csignal>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sys/time.h>

#include "particle.h"
#include "globals.h"
#include "fields.h"
#include "geometry.h"
#include "mc.h" 
#include "integrator.h"
#include "bruteforce.h"
#include "adiabacity.h"
#include "ndist.h"
#include "nrutil.h"
#include "racetrack.h"

using namespace std;

#define TINY 1.0e-30

void ConfigInit(); // read config.in
void IntegrateParticle(TParticle *particle); // integrate particle trajectory
void SetEndValues(TParticle *particle, long double x, long double *y, long double H);
void derivs(long double x, long double *y, long double *dydx); // calculates derivatives of state vector
void SetExpPhase(long double t);
void SpinFlipCheck(TParticle *particle, long double x, long double *y);
void OpenFiles(int snapshot); // open out-files and write file headers
void PrintConfig(); // print info on config
void PrintParticle(FILE *file, TParticle *particle);
void OutputCodes(vector<int> kennz_counter[3]);
void PrintIntegrationStep(TParticle *particle, long double x, long double h, long double *y, long double H); // print step to track-file


FILE *SNAP = NULL, *ENDLOG = NULL;
long Zeilencount;
int Filecount=0;                                   // counts the output files, counter for last line written in outs, random generator temp value
FILE *OUTFILE1 = NULL;


// integrator params
const int nvar = 6; // number of variables in derivs
const long double eps = 1.0e-13, hmin = 0;      // desired relative precision for tracking of n and p 10^-13, minimum step size

int MonteCarloAnzahl=1;   // number of particles for MC simulation
int protneut;       //user choice for reflecting walls, B-field, prot or neutrons, experiment mode
int runge;                            // Runge-Kutta or Bulirsch-Stoer?  set to Runge right now!!!
unsigned short int flipspin=1;  // MonteCarlo for splinflip

//parameters for different experiment phases
int ffBruteForce,ffreflekt,ffspinflipcheck;	// fullfield
int ruBruteForce,rureflekt,ruspinflipcheck;	// rampup
int rdBruteForce,rdreflekt,rdspinflipcheck;	// rampdown
int fiBruteForce,fireflekt,fispinflipcheck;	// filling
int coBruteForce,coreflekt,cospinflipcheck;	// counting UCN
int clBruteForce,clreflekt,clspinflipcheck;	// cleaning time

unsigned short int BruteForce = 0;
int neutdist=0;
int spinflipcheck = 0;                          // user choice to check adiabacity
int snapshot=0;  // make snapshots of the neutron at specified times
set<int> snapshots;
vector<int> kennz_counter[3];
float InitTime = 0, ReflectionTime = 0, IntegratorTime = 0, DiceTime = 0; // time statistics

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


// uebergabe: jobnumber inpath outpath                 paths without last slash
int main(int argc, char **argv){
	//Initialize signal-analizing
	signal (SIGUSR1, catch_alarm);
	signal (SIGUSR2, catch_alarm);
	signal (SIGXCPU, catch_alarm);   
	
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
	OpenFiles(snapshot);	// Open .in and .out files and write headers
		
	Log(
	" ################################################################\n"
	" ###                 Welcome to PNTracker,                    ###\n"
	" ###     the tracking program for neutrons and protons        ###\n"
	" ################################################################\n");

	//printf("\nMonteCarlo: %i\n MonteCarloAnzahl %i \n", MonteCarlo, MonteCarloAnzahl);

	LoadGeometry((inpath + "/geometry.in").c_str());	// read STL-files		
	
	LoadMCGenerator((inpath + "/all3inone.in").c_str());	// set random seed and load inital value limits
	
	PrintConfig();
	
	InitTime = (1.*clock())/CLOCKS_PER_SEC;

	switch(protneut)
	{	
		case BF_ONLY:	PrintBField((outpath+"/BF.out").c_str()); // estimate ramp heating
						exit(0);
						break;		
		case BF_CUT:	PrintBFieldCut((outpath+"/BFCut.out").c_str()); // print cut through B field
						exit(0);
						break;		
	}

	int ntotalsteps = 0;     // counters to determine average steps per integrator call
	TParticle *particle;
	for (int iMC = 1; iMC <= MonteCarloAnzahl; iMC++)
	{      
		switch(protneut)
		{	
			case NEUTRON:	particle = new TParticle(NEUTRON, iMC);
							break;		
			case PROTON:	particle = new TParticle(PROTON, iMC);
							break;		
			case ELECTRON:	particle = new TParticle(ELECTRON, iMC);	
							break;
		}

		timeval dicestart, diceend;
		gettimeofday(&dicestart, NULL);
		MCStartwerte(particle); // dice initial values and write them into particle
		gettimeofday(&diceend, NULL);
		DiceTime += diceend.tv_sec - dicestart.tv_sec + float(diceend.tv_usec - dicestart.tv_usec)/1e6;

		IntegrateParticle(particle); // integrate particle trajectory
		ntotalsteps += particle->nok + particle->nbad;
		
		if(particle->kennz == KENNZAHL_DECAYED && decay == 2)
		{	
			TParticle p(PROTON, iMC);
			TParticle e(ELECTRON, iMC);
			MCZerfallsstartwerte(particle, &p, &e); // dice decay products
			IntegrateParticle(&p); // integrate decay products
			IntegrateParticle(&e);
			ntotalsteps += p.nok + p.nbad + e.nok + e.nbad;
		}
		delete particle;
	}

	ostringstream path;
	path << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "ndist.out";	
	if (neutdist == 1) outndist(path.str().c_str());   // Neutronenverteilung in der Flasche ausgeben

	OutputCodes(kennz_counter); // print particle kennzahlen
	
	Log("Integrator used (1 Bulirsch Stoer, 2 Runge Kutta): %d \n", runge);
	Log("The integrator made %d steps. \n", ntotalsteps);
	float SimulationTime = (1.*clock())/CLOCKS_PER_SEC - InitTime;
	IntegratorTime -= ReflectionTime;
	Log("Init: %.2fs, Simulation: %.2fs, Integrator: %.2fs (%.2f%%), Reflection: %.2fs (%.2f%%), Dicing: %.4fs\n",
			InitTime, SimulationTime, IntegratorTime, IntegratorTime*100/SimulationTime, ReflectionTime, ReflectionTime*100/SimulationTime, DiceTime);
	Log("That's it... Have a nice day!\n");
	
	FreeFields();
	
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

// read config.in
void ConfigInit(void){
	/* setting default values */
	protneut = NEUTRON;
	runge = 2;
	polarisation = 1;
	spinflipcheck = 0;
	BruteForce = 0;
	neutdist = 0;
	reflekt = 1;
	diffuse = 3;
	bfeldwahl = 0;
	ausgabewunsch = 2;
	MonteCarloAnzahl = 1;
	/*end default values*/
	
	ifstream infile((inpath+"/config.in").c_str());
	map<string, map<string, string> > config;
	char c;
	string rest,section,key;
	while (infile.good() && (infile >> ws) && (c = infile.peek())){
		if (c == '[' && infile.ignore()){
			if (infile.peek() == '/'){
				section = "";
//				cout << endl;
			}
			else{
				getline(infile, section, ']');
//				cout << "section: " << section << endl;
			}
			getline(infile,rest);
		}
		else if (c == '#')
			getline(infile,rest);
		else if (section != ""){
			infile >> key;
			getline(infile,config[section][key]);
//			cout << key << ":" << config[section][key];
		}
		else
			getline(infile,rest);
	}
	
	istringstream(config["global"]["protneut"])				>> protneut;
	istringstream(config["global"]["BruteForce"])			>> BruteForce;
	istringstream(config["global"]["reflekt"])				>> reflekt;
	istringstream(config["global"]["spinflipcheck"])		>> spinflipcheck;
	istringstream(config["global"]["runge"])				>> runge;
	istringstream(config["global"]["polarisation"])			>> polarisation;
	istringstream(config["global"]["neutdist"])				>> neutdist;
	istringstream(config["global"]["diffuse"])				>> diffuse;
	istringstream(config["global"]["bfeldwahl"])			>> bfeldwahl;
	istringstream(config["global"]["BFeldSkalGlobal"])		>> BFeldSkalGlobal;
	istringstream(config["global"]["EFeldSkal"])			>> EFeldSkal;
	istringstream(config["global"]["Racetracks"])			>> Racetracks;
	istringstream(config["global"]["Ibar"])					>> Ibar;	
	istringstream(config["global"]["ausgabewunsch"])		>> ausgabewunsch;
	
	istringstream(config["global"]["snapshot"])				>> snapshot;
	int tmp_snapshot;
	istringstream isnapshots(config["global"]["snapshots"]);	
	do{
		isnapshots >> tmp_snapshot;
		if (isnapshots.good()) snapshots.insert(tmp_snapshot);
	}while(isnapshots.good());
	
	istringstream(config["global"]["reflektlog"])			>> reflektlog;
	istringstream(config["global"]["MonteCarloAnzahl"])		>> MonteCarloAnzahl;
	istringstream(config["global"]["FillingTime"])			>> FillingTime;
	istringstream(config["global"]["CleaningTime"])			>> CleaningTime;
	istringstream(config["global"]["RampUpTime"])			>> RampUpTime;
	istringstream(config["global"]["FullFieldTime"])		>> FullFieldTime;
	istringstream(config["global"]["RampDownTime"])			>> RampDownTime;
	istringstream(config["global"]["EmptyingTime"])			>> EmptyingTime;
	istringstream(config["global"]["storagetime"])			>> StorageTime;
	istringstream(config["global"]["BFTargetB"])			>> BFTargetB;
	istringstream(config["global"]["decay"])				>> decay;
	istringstream(config["global"]["tau"])					>> tau;	
	istringstream(config["global"]["decayoffset"])			>> decayoffset;	
	istringstream(config["global"]["BCutPlane"])			>> BCutPlanePoint[0] >> BCutPlanePoint[1] >> BCutPlanePoint[2] 
															>> BCutPlaneNormalAlpha >> BCutPlaneNormalGamma >> BCutPlaneSampleDist >> BCutPlaneSampleCount;

	
	istringstream(config["filling"]["BruteForce"])			>> fiBruteForce;
	istringstream(config["filling"]["reflekt"])				>> fireflekt;
	istringstream(config["filling"]["spinflipcheck"])		>> fispinflipcheck;
	                                    
	istringstream(config["cleaning"]["BruteForce"])			>> clBruteForce;
	istringstream(config["cleaning"]["reflekt"])			>> clreflekt;
	istringstream(config["cleaning"]["spinflipcheck"])		>> clspinflipcheck;

	istringstream(config["rampup"]["BruteForce"])			>> ruBruteForce;
	istringstream(config["rampup"]["reflekt"])				>> rureflekt;
	istringstream(config["rampup"]["spinflipcheck"])		>> ruspinflipcheck;
	                                    
	istringstream(config["fullfield"]["BruteForce"])		>> ffBruteForce;
	istringstream(config["fullfield"]["reflekt"])			>> ffreflekt;
	istringstream(config["fullfield"]["spinflipcheck"])		>> ffspinflipcheck;
	                                    
	istringstream(config["rampdown"]["BruteForce"])			>> rdBruteForce	;
	istringstream(config["rampdown"]["reflekt"])			>> rdreflekt;
	istringstream(config["rampdown"]["spinflipcheck"])		>> rdspinflipcheck;
	                                    
	istringstream(config["counting"]["BruteForce"])			>> coBruteForce;
	istringstream(config["counting"]["reflekt"])			>> coreflekt;
	istringstream(config["counting"]["spinflipcheck"])		>> cospinflipcheck;


	if(protneut == NEUTRON)
	{
		if (neutdist == 1) prepndist();
	}
	
	return;
}

// calculates derivatives, passed to integration stepper
void derivs(long double x, long double *y, long double *dydx, void *params){
	((TParticle*)params)->derivs(x,y,dydx);
	return;
}

// integrate particle trajectory
void IntegrateParticle(TParticle *particle){	
	Log("Feldcount = %i\n\n",Feldcount);

	int snapshotsdone=0;
	Feldcount=0;
	
	long double y[nvar+1];
	particle->Fillystart(y);

	Log("Teilchennummer: %i\n",particle->particlenumber);
	if(decay == 2)
		Log("Teilchensorte : %i\n", particle->protneut);
	Log("r: %LG phi: %LG z: %LG v: %LG alpha: %LG gamma: %LG E: %LG t: %LG\n",
		particle->rstart, particle->phistart/conv, particle->zstart, particle->vstart, particle->alphastart/conv,
		particle->gammastart/conv, particle->Hstart, particle->xend);

	timeval intstart, intend;
	gettimeofday(&intstart, NULL);	
	//###################### Integrationsroutine #####################
	int i;
	long double x = 0, hnext = 0, hdid = 0, h = particle->h1; // current time, hnext: suggestion for next time step, hdid: actual timestep done by stepper 
	long double yscal[nvar+1], dydx[nvar+1]; // yscal needed for error estimation, dydx holds derivatives of y
	int perc=0;   // percentage of particle done counter
	long double timetemp = 0; // temporre Variable, Zeit wann letzter Schritt in outs geschrieben wurde
	long double xprev, yprev[nvar+1], H = particle->Hstart; // save last position to check for reflection, calculate trajectory length etc., energy at current position
	int itercount = 0; // iteration counter for reflection check

	x=particle->xstart;
	SetExpPhase(x); // set experiment phase parameters (rampup, fullfield, etc)
	
	while (particle->kennz == KENNZAHL_UNKNOWN){ // while nothing happened to particle 
		derivs(x,y,dydx,particle);
		for (i=1;i<=nvar;i++){
			yscal[i]=fabsl(y[i])+fabsl(dydx[i]*h)+TINY; //Scaling used to monitor accuracy.
		}
		
		// spin flip properties according to Vladimirsky and thumbrule
		if ((spinflipcheck == 2) && (particle->protneut == NEUTRON)){
			SpinFlipCheck(particle, x, y);
		}

		xprev = x; // save current point before step
		for (i=1;i<=nvar;i++) {
			yprev[i] = y[i];
		}
	
		if (particle->xstart + particle->xend < x + h) 
			h = particle->xstart + particle->xend - x;	//If stepsize can overshoot, decrease.
		if (h > particle->hmax) h = particle->hmax;
		if((Bws<BFTargetB) && BruteForce && h > 1e-5)
			h = 1e-5;		
		else if( (Bws<(BFTargetB+0.1) ) && BruteForce && h > 1e-4)
			h = 1e-4;
			
		try{
			if (runge == 1) bsstep(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,particle,derivs);        // runge kutta or bulirsch stoer step
			else if (runge == 2) rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,particle,derivs);        // runge kutta or bulirsch stoer step
		}
		catch(...){ // catch exception thrown by numerical recipes routines
			particle->kennz = KENNZAHL_NRERROR;
		}
		
		timeval reflectstart, reflectend;
		gettimeofday(&reflectstart, NULL);	
		if (ReflectCheck(particle, xprev, hdid, yprev, y, itercount)){ // check if particle should be reflected
			// if ReflectCheck returns true, reset last step and repeat with smaller timestep
			for (i = 1; i <= nvar; i++) y[i] = yprev[i];
			x = xprev;
			h = hdid;
			itercount++;
			gettimeofday(&reflectend, NULL);
			ReflectionTime += reflectend.tv_sec - reflectstart.tv_sec + float(reflectend.tv_usec - reflectstart.tv_usec)/1e6;
			continue;
		}
		itercount = 0;		
		gettimeofday(&reflectend, NULL);
		ReflectionTime += reflectend.tv_sec - reflectstart.tv_sec + float(reflectend.tv_usec - reflectstart.tv_usec)/1e6;
		
		SetExpPhase(x);
		
		// Trajectory length calculation
//		long double trajlength = sqrtl(pow(yprev[1] - y[1],2) + pow(yprev[3] - y[3],2) + pow(yprev[5]*yprev[1] - y[5]*y[1],2));
//		particle->trajlength += trajlength;
		particle->trajlength += sqrt(pow(yprev[1]*cos(yprev[5]) - y[1]*cos(y[5]),2) + 
									pow(yprev[1]*sin(yprev[5]) - y[1]*sin(y[5]),2) + 
									pow(yprev[3] - y[3],2));
		
		H = particle->Energy(x,y);
		if (H > particle->Hmax) particle->Hmax = H;
	
		AbsorbCheck(particle,yprev,y,H); // check if particle is absorbed in current medium
	
		if (hdid == h) 
			++particle->nok; 
		else 
			++particle->nbad;
	
		percent(x, particle->xstart, particle->xend, perc); // print percent of calculation
	
		long double BahnPointSaveTime = particle->BahnPointSaveTime;
		if ((!BruteForce) && (particle->protneut == NEUTRON))
			BahnPointSaveTime = 1e-3;
		else if (BruteForce)
			BahnPointSaveTime = 1e-4;		
		if(spinflipcheck==3)
			BahnPointSaveTime = 1e-4;	
		if ((x-timetemp) >= BahnPointSaveTime)
		{
			PrintIntegrationStep(particle, x, hdid, y, H); // print integration step in track file
			timetemp = x;
		}

		// take snapshots of neutrons at certain times
		if(snapshot==1  && particle->protneut==1 && snapshots.count((int) x) > 0 && (int) x != snapshotsdone)
		{
			Log("\n Snapshot at %LG s \n", x);
			SetEndValues(particle, x, y, H);
			PrintParticle(SNAP, particle);
			snapshotsdone=(int) x;	
		}	

		if (BruteForce)
		{
			long double BFsurvprob = BruteForceIntegration(x,y,Br,Bphi,Bz,Bws); // integrate spinflip pobability
			particle->BFsurvprob *= BFsurvprob;
			// flip the spin with the probability 1-BFpol
			if (flipspin && UniformDist(0,1) < 1-BFsurvprob) 
			{
				particle->polarisation *= -1;
				particle->NSF++; 
				printf("\n The spin has flipped! Number of flips: %i\n",particle->NSF);			
			}
		}
		
		if ((neutdist == 1)&&(particle->protneut == NEUTRON) /*&&(ExpPhase==4)*/ )
			fillndist(xprev, yprev, x, y); // write spatial neutron distribution
		
		if (x >= particle->xstart + particle->xend && decay != 0 && particle->protneut == NEUTRON)
			particle->kennz = KENNZAHL_DECAYED;
		else if (x >= particle->xstart + particle->xend || x > StorageTime)
			particle->kennz = KENNZAHL_NOT_FINISH;
			
	
		if (fabsl(hnext) <= hmin) 
			nrerror("Step size too small in ODEINT");
		h=hnext;	

	}
	SetEndValues(particle, x, y, H);
	PrintParticle(ENDLOG, particle);// Endwerte schreiben
	Log("Done!!\nBFFlipProb: %.17LG rend: %.17LG zend: %.17LG Eend: %.17LG Code: %i t: %.17LG l: %LG\n",
		(1 - particle->BFsurvprob),y[1],y[3],H,particle->kennz,x, particle->trajlength);		
	//###################### Integrationsroutine #####################
	gettimeofday(&intend, NULL);
	IntegratorTime += intend.tv_sec - intstart.tv_sec + float(intend.tv_usec - intstart.tv_usec)/1e6;
			
	kennz_counter[particle->protneut % 3][particle->kennz]++; // increase kennz-counter
}

// spin flip properties according to Vladimirsky and thumbrule
void SpinFlipCheck(TParticle *particle, long double x, long double *y){
	// y[6]: phidot
	particle->vlad = vladimirsky(y[1], Br, Bphi, Bz, 
					   dBrdr, dBrdphi, dBrdz, dBphidr, dBphidphi, dBphidz, dBzdr, dBzdphi, dBzdz, Bws,
					   y[2], y[6], y[4]);
	particle->frac = thumbrule(Br, Bphi, Bz, 
					   dBrdr, dBrdphi, dBrdz, dBphidr, dBphidphi, dBphidz, dBzdr, dBzdphi, dBzdz, Bws,
					   y[2], y[6], y[4]);
	if (particle->vlad > 1e-99){
		particle->vladtotal *= 1-particle->vlad;
		if (particle->vladtotal < 0.9999)
			Log(" VladShit (%LG) at t= %.17LG\n",particle->vladtotal,x);
	}
	if ((particle->vlad > particle->vladmax)&&(particle->vlad > 1e-99))
		particle->vladmax = log10(particle->vlad);
	if ((particle->frac > particle->thumbmax)&&(particle->frac > 1e-99))
		particle->thumbmax = log10(particle->frac);	
}

// write values into end*-variables of particle
void SetEndValues(TParticle *particle, long double x, long double *y, long double H){
	particle->dt = x - particle->xstart;
	particle->rend = y[1];
	particle->phiend = fmod(y[5],2*pi);
	particle->zend = y[3];
	particle->Hend = H;
	particle->vend = sqrt(y[2]*y[2] + y[4]*y[4] + y[6]*y[1]*y[6]*y[1]);
	if(particle->vend > 0) particle->gammaend = acos(y[4]/particle->vend);
	else particle->gammaend = 0;
	particle->alphaend = atan2(y[6]*y[1],y[2]);	
}

// set experiment paramters (rampup, fullfield, etc)
void SetExpPhase(long double t){
	if ((t < FillingTime)&&(FillingTime>0))
	{      // filling in neutrons
		BruteForce = fiBruteForce; // no spin tracking
		reflekt = fireflekt;
		spinflipcheck = fispinflipcheck;		
		ExpPhase=1;
	}
	
	else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
	{      // spectrum cleaning
		BruteForce = clBruteForce; // no spin tracking
		reflekt = clreflekt;
		spinflipcheck = clspinflipcheck;
		ExpPhase=2;
	}
	else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
	{   // ramping up field
		BruteForce = ruBruteForce; 
		reflekt = rureflekt;
		spinflipcheck = ruspinflipcheck;
		ExpPhase=3;
	}
	else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
	{  // storage time
		BruteForce = ffBruteForce;  // start spin tracking
		reflekt = ffreflekt;
		spinflipcheck = ffspinflipcheck;
		ExpPhase=4;
	}
	else if ((t >= (RampUpTime+CleaningTime+FillingTime+FullFieldTime)) && (t < (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime)) && (RampDownTime != 0))
	{   // ramping down field
		BruteForce = rdBruteForce; // no spin tracking, neutrons shall stay in trap through reflection
		reflekt = rdreflekt;
		spinflipcheck=rdspinflipcheck;
		ExpPhase=5;
	}
	else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
	{      // emptying of neutrons into detector
		BruteForce = coBruteForce;
		spinflipcheck= cospinflipcheck;
		ExpPhase=6;
	}
}


// create out-files and write file headers
void OpenFiles(int snapshot){
	ostringstream logscrfile;
	logscrfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "log.out";
	LOGSCR = fopen(logscrfile.str().c_str(),"w");

	// Endpunkte
	ostringstream endlogfile;
	endlogfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "end.out";
	ENDLOG = fopen(endlogfile.str().c_str(), "w");
	
    fprintf(ENDLOG,"jobnumber RandomSeed teilchen protneut polarisation "
                   "tstart rstart phistart zstart NeutEnergie "
                   "vstart alphastart gammastart decayerror "
                   "rend phiend zend "
                   "vend alphaend gammaend tend dt "
                   "H kennz NSF RodFieldMult BFflipprob "
                   "AnzahlRefl vladmax vladtotal thumbmax trajlength "
                   "Hstart Hmax BFeldSkal EFeldSkal tauSF dtau\n");
	
	// make snapshots as specified time into snapshot.out
	if (snapshot==1)
	{ 
		ostringstream snapshotfile;
		snapshotfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "snapshot.out";
		SNAP = fopen(snapshotfile.str().c_str(),"w");       // open outfile neut001.out
        fprintf(SNAP,"jobnumber RandomSeed teilchen protneut polarisation "
                       "tstart rstart phistart zstart NeutEnergie "
                       "vstart alphastart gammastart decayerror "
                       "rend phiend zend "
                       "vend alphaend gammaend tend dt "
                       "H kennz NSF RodFieldMult BFflipprob "
                       "AnzahlRefl vladmax vladtotal thumbmax trajlength "
                       "Hstart Hmax BFeldSkal EFeldSkal tauSF dtau\n");
	}
	
}


// write the current configuration of the programm to screen and to the logfile
void PrintConfig(void)
{
	// ONSCREEN
	Log("The following parameters were set for program execution: \n");
	Log("protneut: %i => (1) neutron, (2) proton tracking, (3) magnetic field evaluation \n",protneut);
	if (protneut == NEUTRON)
	{
		if (polarisation == 1)
		{
			Log("Low field seekers are tracked! \n");		
		}
		else if (polarisation == 2)
		{
			Log("High field seekers are tracked! \n");		
		}
		else if (polarisation == 3)
		{
			Log("No polarisation! \n");		
		}
		else if (polarisation == 4)
		{
			Log("Polarisation is diced! \n");		
		}
		
		if (fireflekt == 1)
			Log("Filling Time %LG \n Reflection is on, ", FillingTime );
		else if (fireflekt == 0)
			Log("Filling Time %LG \n Reflection is off, ",FillingTime);
		if (fiBruteForce == 1)
			Log("Brute Force on, ");
		else if (fiBruteForce == 0)
			Log("Brute Force off, ");		
		
		if (clreflekt == 1)
			Log("Cleaning Time %LG \n Reflection is on, ", CleaningTime );
		else if (clreflekt == 0)
			Log("Cleaning  Time %LG \n Reflection is off, ",CleaningTime);
		if (clBruteForce == 1)
			Log("Brute Force on, ");
		else if (clBruteForce == 0)
			Log("Brute Force off, ");
		
		if (rureflekt == 1)
			Log("Ramp Up Time %LG \n Reflection is on, ", RampUpTime );
		else if (rureflekt == 0)
			Log("Ramp Up  Time %LG \n Reflection is off, ",RampUpTime);
		if (ruBruteForce == 1)
			Log("Brute Force on, ");
		else if (ruBruteForce == 0)
			Log("Brute Force off, ");
		
		if (ffreflekt == 1)
			Log("Full Field Time %LG \n Reflection is on, ", FullFieldTime );
		else if (ffreflekt == 0)
			Log("Full Field Time %LG \n Reflection is off, ",FullFieldTime);
		if (ffBruteForce == 1)
			Log("Brute Force on, ");
		else if (ffBruteForce == 0)
			Log("Brute Force off, ");
		
		if (rdreflekt == 1)
			Log("Ramp Down Time  %LG \n Reflection is on, ", RampDownTime );
		else if (rdreflekt == 0)
			Log("Ramp Down Time  %LG \n Reflection is off, ",RampDownTime);
		if (rdBruteForce == 1)
			Log("Brute Force on, ");
		else if (rdBruteForce == 0)
			Log("Brute Force off, ");
		
		if (coreflekt == 1)
			Log("Counting Time  \n Reflection is on, ");
		else if (coreflekt == 0)
			Log("Counting Time  \n Reflection is off, ");
		if (coBruteForce == 1)
			Log("Brute Force on, ");
		else if (coBruteForce == 0)
			Log("Brute Force off, ");
		
			
	}
	Log("The choice of Bfield is: %i  (0) interpolated field, (1) no field, (2) check interpolation routine \n",bfeldwahl);
	Log("The choice of reflection is: %i  (1) specular, (2) diffuse, (3) statistically specular or diffuse \n",diffuse);
	Log("Choice of output: %i  (1) Endpoints and track (2) only endpoints (3) Endpoints, track and spin \n (4) Endpoints and spin (5) nothing \n ", ausgabewunsch);
}


//Ausgabe der Zwischenwerte aus odeint
void PrintIntegrationStep(TParticle *particle, long double x, long double h, long double *y, long double H){
	if ((ausgabewunsch==OUTPUT_EVERYTHING)||(ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)){
		if (Zeilencount>40000){
			fclose(OUTFILE1);
			OUTFILE1 = NULL;
		}
			
		if (!OUTFILE1){
			Filecount++;
			ostringstream wholetrackfile;
			wholetrackfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "track" << setw(3) << Filecount << setw(0) << ".out";
			OUTFILE1=fopen(wholetrackfile.str().c_str(),"w");
			fprintf(OUTFILE1,"Teilchen protneut polarisation t r drdt z dzdt phi dphidt x y "
							 "v H Br dBrdr dBrdphi dBrdz Bphi dBphidr "
							 "dBphidphi dBphidz Bz dBzdr dBzdphi dBzdz Babs Er Ez "
							 "timestep logvlad logthumb ExpPhase \n");
			Log(" ##");
			Log(wholetrackfile.str().c_str());
			Log("## \n");
			Zeilencount=1;
		}
	
			
		printf("-");
		fflush(stdout);
		
		long double v = sqrt(y[2]*y[2] + y[4]*y[4] + y[6]*y[1]*y[6]*y[1]);
		
		long double logvlad = 0.0, logfrac = 0.0;
		if (particle->vlad > 1e-99) 
			logvlad=log10(particle->vlad);
		if (particle->frac > 1e-99) 
			logfrac=log10(particle->frac);
		
		int pol;
		if (particle->polarisation == -1) pol = 1;
		else if (particle->polarisation == 0) pol = 3;
		else if (particle->polarisation == 1) pol = 2;
		
		//cout << "Br " << Bp[1][klauf] << endl;
		fprintf(OUTFILE1,"%d %d %d %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
						 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG "
						 "%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %d \n",
						 particle->particlenumber,particle->protneut,pol,x,y[1],y[2],y[3],y[4],y[5],y[6],y[1]*cos(y[5]),y[1]*sin(y[5]),
						 v,H,Br,dBrdr,dBrdphi,dBrdz,Bphi,dBphidr, dBphidphi,dBphidz,
						 Bz,dBzdr,dBzdphi,dBzdz,Bws,Er,Ez,h,logvlad,logfrac,ExpPhase);
		//fprintf(OUTFILE1,"%LG\n",xp[klauf]);
		
		fflush(OUTFILE1);
		Zeilencount++;
	}

}

// print particle start- and end-values into file
void PrintParticle(FILE *file, TParticle *particle){
	// calculate spin flip lifetime tauSF and influence on lifetime measurement 
	long double tauSF = (1 - particle->BFsurvprob != 0) ? (-(particle->dt)/log(particle->BFsurvprob)) : 9e99;
	long double dtau=tau-1/(1/tau+1/tauSF) ;
	 
	int pol;
	if (particle->polarisation == -1) pol = 1;
	else if (particle->polarisation == 0) pol = 3;
	else if (particle->polarisation == 1) pol = 2;

	fprintf(file ,"%i %li %i %i %i "
	               "%LG %LG %LG %LG %LG "
	               "%LG %LG %LG %LG "
	               "%LG %LG %LG "
	               "%LG %LG %LG %LG %LG "
	               "%LG %i %i %f %LG "
	               "%i %LG %LG %LG %LG "
	               "%LG %LG %LG %LG %LG %LG\n",
	               jobnumber, monthinmilliseconds, particle->particlenumber, particle->protneut, pol,
	               particle->xstart, particle->rstart, particle->phistart, particle->zstart, particle->NeutEnergie*1.0e9,
	               particle->vstart, particle->alphastart, particle->gammastart, particle->decayerror,
	               particle->rend, particle->phiend, particle->zend,
	               particle->vend, particle->alphaend, particle->gammaend, particle->xstart + particle->xend, particle->dt,
	               particle->Hend, particle->kennz, particle->NSF, 1.0, 1 - particle->BFsurvprob,
	               particle->nrefl, particle->vladmax, particle->vladtotal, particle->thumbmax, particle->trajlength,
	               particle->Hstart, particle->Hmax, BFeldSkal, EFeldSkal, tauSF, dtau);
	               
	fflush(file);
}


//======== Output of the endcodes ==========================================================================================
void OutputCodes(vector<int> kennz_counter[3]){
	int ncount = accumulate(kennz_counter[1].begin(),kennz_counter[1].end(),0);
	int pcount = accumulate(kennz_counter[2].begin(),kennz_counter[2].end(),0);
	int ecount = accumulate(kennz_counter[0].begin(),kennz_counter[0].end(),0);
	Log("\nThe calculations of %li particle(s) yielded:\n"
	       "endcode:  of %4i neutron(s) ; of %4i proton(s) ; of %4i electron(s)\n"
	       "   0 %12i %20i %19i 		(were not categorized)\n"
	       "   1 %12i %20i %19i 		(did not finish)\n"
	       "   2 %12i %20i %19i 		(hit outer boundaries)\n"
	       "   3 %12i %20i %19i 		(produced a numerical error (most likely step size underflow))\n"
	       "   4 %12i %20i %19i 		(decayed)\n"
	       "   5 %12i %20i %19i 		(found no initial position)\n",
	       ncount + pcount + ecount,
	       ncount, pcount, ecount,
	       kennz_counter[1][0], kennz_counter[2][0], kennz_counter[0][0],
	       kennz_counter[1][1], kennz_counter[2][1], kennz_counter[0][1],
	       kennz_counter[1][2], kennz_counter[2][2], kennz_counter[0][2],
	       kennz_counter[1][3], kennz_counter[2][3], kennz_counter[0][3],
	       kennz_counter[1][4], kennz_counter[2][4], kennz_counter[0][4],
	       kennz_counter[1][5], kennz_counter[2][5], kennz_counter[0][5]);
	for (unsigned i = 6; i < kennz_counter[0].size(); i++){
		string solidnames;
		for (vector<solid>::iterator it = solids.begin(); it != solids.end(); it++)
			if (it->kennz == i)
				solidnames += '/' + it->name;
		Log("  %2i %12i %20i %19i		(were statistically absorbed by %s)\n",
				i,kennz_counter[1][i],kennz_counter[2][i],kennz_counter[0][i],solidnames.c_str()+1);
	}
}
//======== end of OutputCodes ==============================================================================================


/*
// print debug info
void OutputState(TParticle *particle, long double x, long double h, long double *y){
	ostringstream stateoutfile;
	stateoutfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "state.out";
	FILE *STATEOUT = fopen(stateoutfile.str().c_str(),"w");	
	printf("\n \n Something happened! \n \n ");
	fprintf(STATEOUT,"In this file the state of the programm will be written, when something happens for debug porposes... \n");
	fprintf(STATEOUT,"Jobnumber: %i \n", jobnumber);
	fprintf(STATEOUT,"Particle number: %i \n", particle->particlenumber);
	fprintf(STATEOUT,"Starting energy: %.10LG neV\n", particle->NeutEnergie*1e9);
	fprintf(STATEOUT,"Time: %.10LG s\n", x);
	fprintf(STATEOUT,"Endtime: %.10LG s\n", particle->xend);
	fprintf(STATEOUT,"Particle ID: %i \n", particle->kennz);
	fprintf(STATEOUT,"Current energy: %.10LG neV\n", particle->Energy(x,y));
	fprintf(STATEOUT,"Spatial coordinates (y[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", y[1], y[5], y[3]);	
	fprintf(STATEOUT,"Spatial coordinates (ystart[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", particle->rstart, particle->phistart, particle->zstart);	
	fprintf(STATEOUT,"Velocities: rdot: %.10LG, phidot: %.10LG, zdot: %.10LG, vabs: %.10LG \n", y[2], y[6], y[4], sqrtl(fabsl(y[2]*y[2]+y[1]*y[1]*y[6]*y[6]+y[4]*y[4])));
	fprintf(STATEOUT,"Magnetic Field: Br: %.10LG, Bphi: %.10LG, Bz: %.10LG, Babs: %.10LG \n", Br,Bphi,Bz,Bws);
	fprintf(STATEOUT,"Magnetic Field Derivatives: \n  %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG \n", dBrdr,dBrdphi,dBrdz,dBphidr,dBphidphi,dBphidz,dBzdr,dBzdphi,dBzdz);
	fprintf(STATEOUT,"Reflection state (1 is on): %i \n",reflekt);
	fprintf(STATEOUT,"Diffuse reflection state (1 is specular, 2 diffuse, 3 statistical): %i \n",diffuse);
	fprintf(STATEOUT,"timestep: %.10LG\n", h);
	fprintf(STATEOUT,"Reflekt while Rampdown: %i\n", rdreflekt);
	fclose(STATEOUT);
	particle->kennz = 88;
//	exit(-1);
}
*/

