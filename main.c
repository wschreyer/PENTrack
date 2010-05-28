#include "main.h"
//#include "lines.h"
/*
b = 0  kein bereich alle zeilen werden ignoriert bis auf [ startende
b = 1  globale settings
b = 2  ramp up
b = 3  full field
b = 4  ramp down
*/

long double *ndistr = NULL, *ndistz = NULL, **ndistW = NULL;    // matrix for probability of finding particle
int v=300,w=1200;                                       // dimension of matrix above, user choice for ndist

FILE *LOGSCR = NULL, *ENDLOG = NULL;


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
	ausgabewunschsave = 2;
	MonteCarloAnzahl = 1;
	M = m_n;
	mu_n = mu_nSI / ele_e * -1;
	mumB = mu_n/M;
	decay.on = 0;
	decay.ed = 0;
	decay.error = 0;
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
	istringstream(config["global"]["ausgabewunsch"])		>> ausgabewunsch;
	istringstream(config["global"]["snapshot"])		>> snapshot;
	int tmp_snapshot;
	istringstream isnapshots(config["global"]["snapshots"]);	
	do{
		isnapshots >> tmp_snapshot;
		if (isnapshots.good()){
//			cout << "tmp_snapshot " << tmp_snapshot << endl;
//			cout << "snapshot " << *
			snapshots.insert(tmp_snapshot);
//			.first << endl;					
		}			
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
	istringstream(config["global"]["decay"])				>> decay.on;
	istringstream(config["global"]["tau"])					>> tau;	
	istringstream(config["global"]["decayoffset"])		>> decayoffset;	

	
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


	if((decay.on == 2) && ((protneut == PROTON) || (protneut == ELECTRONS))) protneut = NEUTRON;
	
	if(decay.on)
	{	decay.ed = 0;
		decay.error = 0;
	}
	
	polarisationsave = polarisation; // save polarisation choise
	
	if(protneut == NEUTRON)
	{
		if (neutdist == 1) prepndist(1);
	}
	
	ausgabewunschsave=ausgabewunsch;
	return;
}


// create out-files and write file headers
void OpenFiles(int argc, char **argv){
	ostringstream logscrfile;
	logscrfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "log.out";
	LOGSCR = fopen(logscrfile.str().c_str(),mode_w);

	if((ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)||(ausgabewunsch==OUTPUT_ENDPOINTSandSPIN))
	{
		ostringstream BFoutfile1;
		BFoutfile1 << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "BF1.out";
		BFLOG = fopen(BFoutfile1.str().c_str(),mode_w);
		fprintf(BFLOG,"t Babs Polar logPolar Ixnorm Iynorm Iznorm Bxnorm Bynorm Bznorm dx dy dz\n");
	}
	
	// Endpunkte
	ostringstream endlogfile;
	endlogfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "end.out";
	ENDLOG = fopen(endlogfile.str().c_str(),mode_w);
	
	if (protneut != BF_ONLY) 
	{
        fprintf(ENDLOG,"jobnumber RandomSeed teilchen protneut polarisation "
                       "tstart rstart phistart zstart NeutEnergie "
                       "vstart alphastart gammastart decayerror "
                       "rend phiend zend "
                       "vend alphaend gammaend tend dt "
                       "H kennz NSF RodFieldMult BFflipprob "
                       "AnzahlRefl vladmax vladtotal thumbmax trajlength "
                       "Hstart Hmax BFeldSkal EFeldSkal tauSF dtau\n");
	}		
	
	// Print track to file
	if ((ausgabewunsch == OUTPUT_EVERYTHING)||(ausgabewunsch == OUTPUT_EVERYTHINGandSPIN))
	{ 
		ostringstream wholetrackfile;
		wholetrackfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "track001.out";
		OUTFILE1 = fopen(wholetrackfile.str().c_str(),mode_w);       // open outfile neut001.out

		Zeilencount=0;
		fprintf(OUTFILE1,"Teilchen protneut polarisation t r drdt z dzdt phi dphidt x y "
						 "v H Br dBrdr dBrdphi dBrdz Bphi dBphidr dBphidphi dBphidz "
						 "Bz dBzdr dBzdphi dBzdz Babs Er Ez timestep logvlad logthumb ExpPhase\n");
						 
	}
	
	// make snapshots as specified time into snapshot.out
	if (snapshot==1)
	{ 
		ostringstream snapshotfile;
		snapshotfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "snapshot.out";
		SNAP = fopen(snapshotfile.str().c_str(),mode_w);       // open outfile neut001.out
		fprintf(SNAP,"jobnumber Teilchennummer protneut polarisation "
                       "tstart rstart phistart zstart NeutEnergie "
                       "vstart alphastart gammastart "
                       "r phi z x y "
                       "v alphaend gammaend tend dt "
                       "H kennz NSF RodFieldMult BFflipprob "
                       "AnzahlRefl trajlength "
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


//======== initalizes initial values from record to used variables ==============================================
void initialStartbed()
{	initial particleini;
	long double edm; // energie dimension multiplikator
	switch(protneut)
	{	case NEUTRON:	particleini = nini;
						edm = 1e-9; 
						break;	
		case PROTON:	particleini = pini;
						edm = 1;
						break;
		case ELECTRONS:	particleini = eini;
						edm = 1e+3;
						break;
	}
	EnergieS = particleini.EnergieS * edm; // [eV]
	EnergieE = particleini.EnergieE * edm; // [eV]
	alphas = particleini.alphas;
	alphae = particleini.alphae;
	gammas = particleini.gammas;
	gammae = particleini.gammae;
	delx = particleini.delx;
	xend = particleini.xend;
	//return;
}
//======== end of innitialStartbed ==============================================================================

//======== read in of the initial values (out of all3inone.in) and the dimensions (out of dimensions.in) ========
void Startbed(int k)
{	int i = 0, ncont;
    char cline[200];    
	string path;

	// setting path for all3inone.in
	path = inpath + "/all3inone.in";

	// creating 'inistream' to all3inone.in
	FILE *inistream = fopen(path.c_str(), mode_r);
	if(inistream == NULL)
	{	cout << "Can't open " << path << "\n";
		exit(-1);
	}
	
	// looping  over the lines of all3inone.in and reading them in
	for(i = 1; i < 31; i++)
	{	fgets(cline, 200, inistream);
		ncont = 42;
		switch (i)
		{	case  2:	ncont = sscanf(cline, "%LG %LG %LG ", &DiceRodField, &RodFieldMultiplicator, &Ibar);
						break;
			case  4:	ncont = sscanf(cline, "%LG %LG ", &nini.EnergieS, &nini.EnergieE);
						break;
			case  5:	ncont = sscanf(cline, "%LG %LG ", &nini.alphas, &nini.alphae);
						nini.alphas *= conv;
						nini.alphae *= conv;
						break;
			case  6:	ncont = sscanf(cline, "%LG %LG ", &nini.gammas, &nini.gammae);
						nini.gammas *= conv;
						nini.gammae *= conv;
						break;
			case 7:		ncont = sscanf(cline, "%LG %LG ", &nini.delx, &nini.xend);
						break;
			case 9:		ncont = sscanf(cline, "%LG %LG ", &pini.EnergieS, &pini.EnergieE);
						break;
			case 10:	ncont = sscanf(cline, "%LG %LG ", &pini.alphas, &pini.alphae);
						pini.alphas *= conv;
						pini.alphae *= conv;
						break;
			case 11:	ncont = sscanf(cline, "%LG %LG ", &pini.gammas, &pini.gammae);
						pini.gammas *= conv;
						pini.gammae *= conv;
						break;
			case 12:	ncont = sscanf(cline, "%LG %LG ", &pini.delx, &pini.xend);
						break;
			case 14:	ncont = sscanf(cline, "%LG %LG ", &eini.EnergieS, &eini.EnergieE);
						break;
			case 15:	ncont = sscanf(cline, "%LG %LG ", &eini.alphas, &eini.alphae);
						eini.alphas *= conv;
						eini.alphae *= conv;
						break;
			case 16:	ncont = sscanf(cline, "%LG %LG ", &eini.gammas, &eini.gammae);
						eini.gammas *= conv;
						eini.gammae *= conv;
						break;
			case 17:	ncont = sscanf(cline, "%LG %LG ", &eini.delx, &eini.xend);
						break;
			case 19:	ncont = sscanf(cline, "%LG %LG %LG ", &BCutPlanePoint[0], &BCutPlanePoint[1], &BCutPlanePoint[2]);
						break;
			case 20:	ncont = sscanf(cline, "%LG %LG ", &BCutPlaneNormalAlpha, &BCutPlaneNormalGamma);
						BCutPlaneNormalAlpha *= conv;
						BCutPlaneNormalGamma *= conv;
						break;
			case 21:	ncont = sscanf(cline, "%LG %u ", &BCutPlaneSampleDist, &BCutPlaneSampleCount);
						break;
		}
		if(ncont < 1) Log("an error occourd while reading the %i. item in line %i of all3inone.in", ncont, i);
    }
    
    // closing 'inistream' to all3inone.in
 	fclose(inistream);
	
	// initalizes initial values from record to used variables
	initialStartbed();

	// minimal necessary energy (Emin_n) > maximum input energy (nini.EnergieE)?? EXIT !!
	if(Emin_n*1e9 > nini.EnergieE)
	{	Log("\n\n\nERROR: Emin_n (= %.5LG neV) > nini.EnergieE (= %.5LG neV)\n"
		             "       Check the energy range of the neutrons and / or the B-field configuration!\n\n"
		             "EXIT!!\n", (Emin_n*1e9), nini.EnergieE);
		exit(-1);
	}
	
	// check if lower neutron energy boundary is possible, if not set to possible value
	if((FillingTime == 0) && (CleaningTime == 0) && (RampUpTime == 0) && ((FullFieldTime != 0) || (RampDownTime != 0)))
	{	nini.EnergieS = max(Emin_n*1e9, nini.EnergieS); // higher energy of the neutron
		if(protneut == NEUTRON) EnergieS = nini.EnergieS; // reinitializes initial value
	}

	// starting value (S) > final value (E)?? EXIT!!
	for(int i = 1; i < 4; i++)
	{	struct initial particleini;
		switch(i)
		{	case 1:	particleini = nini;
					break;	
			case 2:	particleini = pini;
					break;
			case 3:	particleini = eini;
					break;
		}
		if((particleini.EnergieS > particleini.EnergieE) || (particleini.alphas > particleini.alphae) || (particleini.gammas > particleini.gammae))
		{	Log("\n\n\nERROR: Two or more initial values are inconsistent.\n"
			             "       Check ALL starting and final values in all3inone.in!\n\n"
			             "EXIT!!\n");
			exit(-1);
		}
	}

	// writing parameters to screen
	Log("\nStart parameters:\n"
	       "  Energy (min, max): %.17LG neV/eV/keV, %.17LG neV/eV/keV\n"
	       "  Maximum runtime: %.17LG s\n"
	       "  Alpha (min, max): %.17LG degree, %.17LG degree\n"
	       "  Gamma (min, max): %.17LG degree, %.17LG degree\n"
	       "  Filling time: %LG s\n"
	       "  Cleaning time: %.17LG s\n"
	       "  Ramp up time: %.17LG s\n"
	       "  Full field time: %.17LG s\n"
	       "  RampDownTime: %.17LG s\n"
	       "  B field scaling factor: %.17LG\n"
	       "  E field scaling factor: %.17LG\n",
	       EnergieS, EnergieE,
	       xend,
	       alphas/conv, alphae/conv,
	       gammas/conv,  gammae/conv,
	       FillingTime,
	       CleaningTime,
	       RampUpTime,
	       FullFieldTime,
	       RampDownTime,
	       BFeldSkalGlobal,
	       EFeldSkal);
}
//======== end of Startbed ======================================================================================


// works like printf, prints to logfile and, if jobnumber == 0, to stdout 
int Log(const char* format, ...){
	va_list args;
	va_start(args, format);
	if (jobnumber == 0){
		vprintf(format, args);
		fflush(stdout);
	} 
	int result = vfprintf(LOGSCR, format, args);
	fflush(LOGSCR);
	va_end(args);
	return result;
}


// write finished particle to endlog
void ausgabe(long double x2, long double *ystart, long double vend, long double H)
{	
	long double dt = x2 - xstart; // simulation time dt

	if(vend>0) gammaend= acosl(ystart[4]/vend);
	else gammaend=0;
	alphaend= atan2l(ystart[6]*ystart[1],ystart[2]);

	if(dt>=xend) 
	{                                                           //Zeit abgelaufen
		if (decay.on && (protneut == NEUTRON))
		{	kennz = KENNZAHL_DECAYED; // neutron decayed
			decay.ed = 1;
			decay.Npolarisation = polarisation;
			decay.Nr = ystart[1];
			decay.Nphi = phiend;
			decay.Nz = ystart[3];
			decay.Nv = vend;
			decay.Nalpha = alphaend;
			decay.Ngamma = gammaend;
			decay.Nx = x2;
			decay.NH = H; // [neV]
		}
		else
		{	kennz = KENNZAHL_NOT_FINISH; // particle survived until xend
		}			
	}
	else if (x2 >= StorageTime)
		kennz = KENNZAHL_NOT_FINISH;

	// calculate spin flip lifetime tauSF and influence on lifetime measurement 
	long double tauSF = (BFflipprob != 0) ? (-x2/logl(1-BFflipprob)) : 9e99;
	long double dtau=tau-1/(1/tau+1/tauSF) ;
	 
	// output of end values
	fprintf(ENDLOG,"%i %li %i %i %i "
	               "%LG %LG %LG %LG %LG "
	               "%LG %LG %LG %LG "
	               "%LG %LG %LG "
	               "%LG %LG %LG %LG %LG "
	               "%LG %i %i %LG %LG "
	               "%li %LG %LG %LG %LG "
	               "%LG %LG %LG %LG %LG %LG\n",
	               jobnumber, monthinmilliseconds, iMC, protneut, polarisation,
	               xstart, r_n, phi_n, z_n, NeutEnergie*1.0e9,
	               v_n, alpha, gammaa, decay.error,
	               ystart[1], phiend, ystart[3],
	               vend, alphaend, gammaend, x2, dt,
	               H, kennz, NSF, RodFieldMultiplicator, BFflipprob,
	               nrefl, vladmax, vladtotal, thumbmax, trajlengthsum,
	               Hstart, Hmax, BFeldSkal, EFeldSkal, tauSF, dtau);
	               
	fflush(ENDLOG);
    return;
}


// Diese Funktionen werten die in yp[][] gespeichterten Werte aus und berechnen daraus die
// Neutronenverteilung in der r-z-Ebene. Dies wird �ber die Auswertung der Aufenthalts-
// wahrscheinlichkeit des Neutrons in Quadraten erreicht.
// Wie bei der Interpolation gibt es einen Vektor f�r die r-Werte der Ecken der Quadrate,
// eine f�r die z-Werte und eine Matrix in der die Verteilung abgespeichert wird

// Vorbereitung der r und z-Vektoren
void prepndist(int k)
{
long double dooo, increment;
//v=300;       //  2mm steps in r-direction 0-60cm
//w=1200;      //  z: -0.5 bis1.7
increment = 0.6 / v;  // Schrittweite des Gitters (0.002m)

ndistr=dvector(0,v);       // r-Werte
ndistz=dvector(0,w);       // z-Werte
ndistW=dmatrix(0,v,0,w);   // Verteilung

        // r-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
        dooo=0.0;
        for(int i=0;i<=v;i++)
        {
				ndistr[i] = dooo;
			  dooo = dooo + increment;
        }

        // z-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
        dooo=-0.5;
        for(int i=0;i<=w;i++)
        {
				ndistz[i] = dooo;
				dooo = dooo + increment;
        }

        // Matrix der Verteilung auf 0 setzen
        for(int i=0;i<=v;i++)
        {
                for(int j=0;j<=w;j++)
                {
                ndistW[i][j]=0.0;
                }
        }

 }

// Bef�llung der Verteilungsmatrix
void fillndist(int k)
{
	
	int ir1, iz1, ir2, iz2;

	hunt(ndistr, v, yp[1][1], &ir1);     // Index des ersten Punktes suchen
	hunt(ndistz, w, yp[3][1], &iz1);
	//ir1 = (yp[1][1] - conv_rA) / (conv_rB);
	//iz1 = (yp[3][1] - conv_zA) / (conv_zB);
	if((ir1>v)||(iz1>w)) printf("Ndist Error %d %d \n",ir1,iz1);
	
	ir2=ir1; iz2=iz1;

	for(int klauf=2;klauf<=kount;klauf++)
	{
		ir1=ir2; iz1=iz2;
		hunt(ndistr, v, yp[1][klauf], &ir2);   // Index des n�chsten Punktes suchen
		hunt(ndistz, w, yp[3][klauf], &iz2);
		//ir2 = (yp[1][klauf] - conv_rA) / (conv_rB);
		//iz2 = (yp[3][klauf] - conv_zA) / (conv_zB);
//		printf("%i,%i,",ir2,iz2);
		if((ir1>v)||(iz2>w)) printf("Ndist Error %d %d \n",ir2,iz2);
		
		if ((ir2==ir1) && (iz2==iz1))                     // sind die Pkte im gleichen Quadrat?
		{
			ndistW[ir1][iz1]=ndistW[ir1][iz1] + xp[klauf]-xp[klauf-1];            // Zeit zum Quadrat dazuaddieren
		}
		else
		{
	
			ndistW[ir1][iz1]=ndistW[ir1][iz1] + (xp[klauf]-xp[klauf-1])/2;        // H�lfte der Zeit zum ersten Quadrat
			ndistW[ir2][iz2]=ndistW[ir2][iz2] + (xp[klauf]-xp[klauf-1])/2;        // H�lfte der Zeit zum anderen
		}
	}

}

// Ausgabe in ndist.out
void outndist(int k)
{
	//printf("\nOutputting the particle spacial distribution... \n");
		
	ostringstream path;
	path << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "ndist.out";
	FILE *NDIST=fopen(path.str().c_str(),mode_w);
	int Treffer = 0;
	
	//FILE *NDIST=fopen("ndist.out",mode_w);
	
	fprintf(NDIST,"Rindex Zindex Rmtlpkt Zmtlpkt Whk Treffer\n");
	long double increment = ndistr[1] - ndistr[0];
	int i,j;
	for(i=0;i<=v;i++){
		for(j=0;j<=w;j++){
			if (ndistW[i][j] != 0) Treffer =1;
			if (ndistW[i][j] == 0) Treffer =0;	
			fprintf(NDIST,"%i %i %.5LG %.5LG %.17LG %i\n",i,j,(ndistr[i]+(increment/2.0)),(ndistz[j]+(increment/2.0)), ndistW[i][j], Treffer);
		}
	}
	fclose(NDIST);
	//printf("\nDone outputting the particle spacial distribution! \n");
	
}


// print debug info
void OutputState(long double *y, int l){
	ostringstream stateoutfile;
	stateoutfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "state.out";
	FILE *STATEOUT = fopen(stateoutfile.str().c_str(),mode_w);	
	printf("\n \n Something happened! \n \n ");
	fprintf(STATEOUT,"In this file the state of the programm will be written, when something happens for debug porposes... \n");
	fprintf(STATEOUT,"Jobnumber: %i \n", jobnumber);
	fprintf(STATEOUT,"Particle number: %i \n", iMC);
	fprintf(STATEOUT,"Starting energy: %.10LG neV\n", NeutEnergie*1e9);
	fprintf(STATEOUT,"Time: %.10LG s\n", x2);
	fprintf(STATEOUT,"Endtime: %.10LG s\n", xend);
	fprintf(STATEOUT,"Particle ID: %i \n", kennz);
	fprintf(STATEOUT,"Current energy: %.10LG neV\n", H);
	fprintf(STATEOUT,"Spatial coordinates (y[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", y[1], y[5], y[3]);	
	fprintf(STATEOUT,"Spatial coordinates (ystart[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", ystart[1], ystart[5], ystart[3]);	
	fprintf(STATEOUT,"Velocities: rdot: %.10LG, phidot: %.10LG, zdot: %.10LG, vabs: %.10LG \n", y[2], y[6], y[4], sqrtl(fabsl(y[2]*y[2]+y[1]*y[1]*y[6]*y[6]+y[4]*y[4])));
	fprintf(STATEOUT,"Magnetic Field: Br: %.10LG, Bphi: %.10LG, Bz: %.10LG, Babs: %.10LG \n", Br,Bphi,Bz,Bws);
	fprintf(STATEOUT,"Magnetic Field Derivatives: \n  %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG \n", dBrdr,dBrdphi,dBrdz,dBphidr,dBphidphi,dBphidz,dBzdr,dBzdphi,dBzdz);
	fprintf(STATEOUT,"Reflection state (1 is on): %i \n",reflekt);
	fprintf(STATEOUT,"Diffuse reflection state (1 is specular, 2 diffuse, 3 statistical): %i \n",diffuse);
	fprintf(STATEOUT,"timestep: %.10LG\n", (x2-x1));
	fprintf(STATEOUT,"Reflekt while Rampdown: %i\n", rdreflekt);
	fclose(STATEOUT);
	kennz=88;
	stopall =1;
	exit(-1);
	iMC = MonteCarloAnzahl +1;
}


// print cut through BField to file
void PrintBFieldCut(){
	// transform plane parameters to cartesian coords
	long double P[3] = {BCutPlanePoint[0]*cosl(BCutPlanePoint[1]), BCutPlanePoint[0]*sinl(BCutPlanePoint[1]), BCutPlanePoint[2]};
	long double n[3] = {cosl(BCutPlanePoint[1]+BCutPlaneNormalAlpha)*sinl(BCutPlaneNormalGamma), 
						sinl(BCutPlanePoint[1]+BCutPlaneNormalAlpha)*sinl(BCutPlaneNormalGamma), 
						cosl(BCutPlaneNormalGamma)};

	// project x/y/z-axes on plane for direction vectors u,v
	long double u[3], v[3];
	if (n[0] < 0.9){		// n not parallel to x-axis
		u[0] = 1-n[0]*n[0];
		u[1] = -n[0]*n[1];
		u[2] = -n[0]*n[2];
	}
	else{					// if n parallel to x-axis use z-axis for u
		u[0] = -n[2]*n[0];
		u[1] = -n[2]*n[1];
		u[2] = 1-n[2]*n[1];
	}
	if (n[1] < 0.9){		// n not parallel to y-axis
		v[0] = -n[1]*n[0];
		v[1] = 1-n[1]*n[1];
		v[2] = -n[1]*n[2];
	}
	else{					// if n parallel to y-axis use z-axis for v
		v[0] = -n[2]*n[0];
		v[1] = -n[2]*n[1];
		v[2] = 1-n[2]*n[1];
	}
	
	int i,j,k;
	// normalize u,v
	long double uabs = sqrtl(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
	long double vabs = sqrtl(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	for (i = 0; i < 3; i++){
		u[i] /= uabs;
		v[i] /= vabs;
	}
	
	// print BField to file
	string path = outpath + "/BFCut.out";
	FILE *cutfile = fopen(path.c_str(), mode_w);
	if (!cutfile){
		Log("Could not open %s!",path.c_str());
		exit(-1);
	}
	fprintf(cutfile, "r phi z x y Br Bphi Bz Bx By dBrdr dBrdphi dBrdz dBphidr dBphidphi dBphidz dBzdr dBzdphi dBzdz Babs dBdr dBdphi dBdz Er Ephi Ez Ex Ey dErdr dErdz dEphidr dEphidz dEzdr dEzdz\n");
	
	long double Pp[3],r,phi,Bx,By,Ex,Ey;
	float start = clock();
	for (i = -BCutPlaneSampleCount/2; i <= BCutPlaneSampleCount/2; i++) {
		for (j = -BCutPlaneSampleCount/2; j <= BCutPlaneSampleCount/2; j++){
			for (k = 0; k < 3; k++)
				Pp[k] = P[k] + i*BCutPlaneSampleDist*u[k] + j*BCutPlaneSampleDist*v[k];
			r = sqrtl(Pp[0]*Pp[0]+Pp[1]*Pp[1]);
			phi = atan2(Pp[1],Pp[0]);
			BFeld(r, phi, Pp[2], 0);
			Bx = Br*cosl(phi) - Bphi*sinl(phi);
			By = Br*sinl(phi) + Bphi*cosl(phi);
			fprintf(cutfile, "%LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG ",
							  r,phi,Pp[2],Pp[0],Pp[1],Br,Bphi,Bz,Bx,By,dBrdr,dBrdphi,dBrdz,dBphidr,dBphidphi,dBphidz,dBzdr,dBzdphi,dBzdz,Bws,dBdr,dBdphi,dBdz);
			EFeld(r,phi,Pp[2]);
			Ex = Er*cos(phi) - Ephi*sin(phi);
			Ey = Er*sin(phi) + Ephi*cos(phi);
			fprintf(cutfile, "%LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
							  Er,Ephi,Ez,Ex,Ey,dErdr,dErdz,dEphidr,dEphidz,dEzdr,dEzdz);
		}
	}
	start = (clock() - start)/CLOCKS_PER_SEC;
	fclose(cutfile);
	printf("Called BFeld and EFeld %u times in %fs (%fms per call)",BCutPlaneSampleCount*BCutPlaneSampleCount, start, start/BCutPlaneSampleCount/BCutPlaneSampleCount);
}


// calculate heating of the neutron gas by field rampup
void PrintBField(){
	fprintf(ENDLOG,"r phi z Br Bphi Bz 0 0 Babs\n");
	BFeldSkal = 1.0;
	int E;
	long double rmin = 0.12, rmax = 0.5, zmin = 0, zmax = 1.2;
	long double dr = 0.1, dz = 0.1;
	long double VolumeB[(int)(nini.EnergieE) + 1];
	for (E = 0; E <= nini.EnergieE; E++) VolumeB[E] = 0;
	
	long double EnTest;
	for (long double r = rmin; r <= rmax; r += dr){
		for (long double z = zmin; z <= zmax; z += dz){
			BFeld(r, 0, z, 500.0);
			fprintf(ENDLOG,"%LG %G %LG %LG %LG %LG %G %G %LG \n",r/lengthconv,0.0,z/lengthconv,Br/Bconv,Bphi/Bconv,Bz/Bconv,0.0,0.0,Bws/Bconv);
			cout << "r= " << r << " z= " << z << " Br= " << Br << " T, Bz= " << Bz << " T"  << endl;
			
			// Ramp Heating Analysis
			for (E = 0; E <= nini.EnergieE; E++){
				EnTest = E*1.0e-9 - m_n*gravconst*z + mu_n * Bws;
				if (EnTest >= 0){
					// add the volume segment to the volume that is accessible to a neutron with energy Energie
					VolumeB[E] = VolumeB[E] + pi * dz * ((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
				}
			}
		}
	}

	// for investigating ramp heating of neutrons, volume accessible to neutrons with and
	// without B-field is calculated and the heating approximated by thermodynamical means
	Log("\nEnergie [neV], Volumen ohne B-Feld, mit B-Feld, 'Erwaermung'");
	long double Volume;
	for (E = 0; E <= nini.EnergieE; E++) 
	{
		Volume = ((E * 1.0e-9 / (m_n * gravconst))) * pi * (rmax*rmax-rmin*rmin);
		// isentropische zustandsnderung, kappa=5/3
		Log("\n%i %.17LG %.17LG %.17LG",E,Volume,VolumeB[E],E * powl((Volume/VolumeB[E]),(2.0/3.0)) - E);
	}
}


//-----------------------------------------------------------------------------
// Vektor W (kein Ortsvektor, sondern B-Feld o.�hnliches) von zyl in kart koord xyz umrechnen ## Wr, Wphi, Wz: Eingabewerte ## Wx0, Wy0, Wz0 werden ver�ndert von CylKartCoord
// phi ist dabei die winkelkoordinate des aktuellen ortsvektors
void CylKartCoord(long double Wr, long double Wphi, long double Wz, long double phi, long double *Wx0, long double *Wy0, long double *Wz0){
	(*Wx0) = cosl(phi) * Wr - sinl(phi) * Wphi;
	(*Wy0) = sinl(phi) * Wr + cosl(phi) * Wphi;
	(*Wz0) = Wz;
	return;
}


//-----------------------------------------------------------------------------
// Vektor W (kein Ortsvektor, sondern B-Feld o.�hnliches) von kart in zyl koord xyz umrechnen ## Wx, Wy, Wz: Eingabewerte ## Wr0, Wphi0, Wz0 werden ver�ndert von KartCylCoord
// phi ist dabei die winkelkoordinate des aktuellen ortsvektors
void KartCylCoord(long double Wx, long double Wy, long double Wz, long double phi, long double *Wr0, long double *Wphi0, long double *Wz0){
	(*Wr0) = cosl(phi) * Wx + sinl(phi) * Wy;
	(*Wphi0) = (-1) * sinl(phi) * Wx + cosl(phi) * Wy;
	(*Wz0) = Wz;
	return;
}


// return absolute value of a 3 vector in cartesian coordinates
long double AbsValueCart(long double x, long double y, long double z)
{
	return sqrt(x*x + y*y + z*z);
}
