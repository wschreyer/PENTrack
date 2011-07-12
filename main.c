
/*********************************************
pnTracker

**********************************************/
#include <iostream>
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
#include "bruteforce.h"
#include "ndist.h"

using namespace std;

void ConfigInit(); // read config.in
void OpenFiles(FILE *&endlog, FILE *&tracklog, FILE *&snap, FILE *&reflectlog);
void OutputCodes(vector<int> kennz_counter[3]); // print simulation summary at program exit
void PrintBFieldCut(const char *outfile, TField &field); // evaluate fields on given plane and write to outfile
void PrintBField(const char *outfile, TField &field);
void PrintGeometry(const char *outfile, TGeometry &geom); // do many random collisionchecks and write all collisions to outfile


int MonteCarloAnzahl=1; ///< number of particles for MC simulation (read from config)
int protneut; ///< type of particle which shall be simulated (read from config)
int reflektlog = 0; ///< write reflections (1), transmissions (2) or both (3) to file? (read from config)
vector<float> snapshots; ///< times when to take snapshots (read from config)
vector<int> kennz_counter[3]; ///< Array of three vectors (for n, p and e) to count particle kennzahlen
int polarisation = POLARISATION_GOOD; ///< polarisation of neutrons (read from config)
int decay = 2; ///< should neutrons decay? (no: 0; yes: 1; yes, with simulation of decay particles: 2) (read from config)
long double decayoffset = 0; ///< start neutron decay timer after decayoffset seconds (read from config)
long double tau_n = 885.7; ///< life time of neutrons (read from config)
long double BCutPlanePoint[9]; ///< 3 points on plane for field slice (read from config)
int BCutPlaneSampleCount1; ///< number of field samples in BCutPlanePoint[3..5]-BCutPlanePoint[0..2] direction (read from config)
int BCutPlaneSampleCount2; ///< number of field samples in BCutPlanePoint[6..8]-BCutPlanePoint[0..2] direction (read from config)

/**
 * Catch signals.
 *
 * terminates a program if a specific signal occurs
 *
 * @param sig signalnumber which called the handler; to get the right number
 * 				for corresponding signals have a look "man signal.h".
 * 				e.g: "SIGFPE" is connected to number 8
 */
void catch_alarm (int sig){
	printf("Program was terminated, because Signal %i occured\n", sig);
	exit(1);
}


/**
 * main function.
 *
 * @param argc Number of parameters passed via the command line
 * @param argv Array of parameters passed via the command line (./Track [jobnumber [configpath [outputpath]]])
 * @return Return 0 on success, value !=0 on failure
 *
 */
int main(int argc, char **argv){
	//Initialize signal-analizing
	signal (SIGUSR1, catch_alarm);
	signal (SIGUSR2, catch_alarm);
	signal (SIGXCPU, catch_alarm);   
	
	jobnumber=0;
	outpath = "./out";
	inpath = "./in";
	if(argc>1) // if user supplied at least 1 arg (outputfilestamp)
		jobnumber=atoi(argv[1]);
	if(argc>2) // if user supplied 2 or more args (outputfilestamp, inpath)
		inpath = argv[2]; // input path pointer set
	if(argc>3) // if user supplied all 3 args (outputfilestamp, inpath, outpath)
		outpath = argv[3]; // set the output path pointer
	
	// read config.in
	ConfigInit();

	if(protneut == NEUTRON){
		if (neutdist == 1) prepndist(); // prepare for neutron distribution-calculation
	}
	else
		neutdist = 0;
	
	cout << "Loading fields..." << endl;
	// load field configuration from geometry.in
	TField field((inpath + "/geometry.in").c_str());

	switch(protneut)
	{
		case BF_ONLY:	PrintBField((outpath+"/BF.out").c_str(), field); // estimate ramp heating
						return 0;
		case BF_CUT:	PrintBFieldCut((outpath+"/BFCut.out").c_str(), field); // print cut through B field
						return 0;
	}


	cout << "Loading geometry..." << endl;
	//load geometry configuration from geometry.in
	TGeometry geom((inpath + "/geometry.in").c_str(), kennz_counter);
	
	if (protneut == GEOMETRY){
		// print random points on walls in file to visualize geometry
		PrintGeometry((outpath+"/geometry.out").c_str(), geom);
		return 0;
	}
	
	cout << "Loading source..." << endl;
	// load source configuration from geometry.in
	TSource source((inpath + "/geometry.in").c_str(), geom, field);
	
	cout << "Loading random number generator..." << endl;
	// load random number generator from all3inone.in
	TMCGenerator mc((inpath + "/all3inone.in").c_str(), polarisation, decay, decayoffset, tau_n);
	
	FILE *endlog = NULL, *tracklog = NULL, *snap = NULL, *reflectlog = NULL;
	// open output files according to config
	OpenFiles(endlog, tracklog, snap, reflectlog);	

	int ntotalsteps = 0;     // counters to determine average steps per integrator call
	float InitTime = (1.*clock())/CLOCKS_PER_SEC, DiceTime = 0, IntegratorTime = 0, ReflTime = 0; // time statistics

	// simulation time counter
	timeval simstart, simend;
	gettimeofday(&simstart, NULL);	

	printf(
	" ################################################################\n"
	" ###                 Welcome to PNTracker,                    ###\n"
	" ###     the tracking program for neutrons and protons        ###\n"
	" ################################################################\n");

/*
	stringstream filename;
	filename << "in/42_0063eout2000m_" << jobnumber << ".out";
	ifstream infile(filename.str().c_str());
	if (!infile.is_open()){
		printf("\ninfile %s not found!\n",filename.str().c_str());
		exit(-1);
	}
	infile.ignore(1024*1024, '\n');
	int i = 0;
	long double r,phi,z,phieuler,thetaeuler,E_n,Ekin,dt,dummy;
	while (infile.good()){
		i++;
		infile >> r >> phi >> z >> phieuler >> thetaeuler >> E_n >> Ekin >> dummy >> dummy >> dummy >>  dummy >> dummy >> dummy >> dummy >> dt;
		infile.ignore(1024*1024, '\n');
		TParticle particle(ELECTRON, i, 0, dt, r, phi*conv, z, Ekin, (phieuler-phi)*conv, thetaeuler*conv, E_n, 0, field);
		particle.Integrate(geom, mc, field, endlog, tracklog, snap, &snapshots, reflectlog);
		kennz_counter[particle.protneut % 3][particle.kennz]++; // increase kennz-counter
		ntotalsteps += particle.nsteps;
		IntegratorTime += particle.comptime;
		ReflTime += particle.refltime;
		infile >> ws;
	}
*/
	TNeutron *neutron = NULL;
	for (int iMC = 1; iMC <= MonteCarloAnzahl; iMC++)
	{      
		timeval dicestart, diceend;
		gettimeofday(&dicestart, NULL); // measure time to find initial point
		neutron = new TNeutron(iMC, source, mc, field); // always create a neutron first, p/e will be simulated from its decay
		gettimeofday(&diceend, NULL);
		DiceTime += diceend.tv_sec - dicestart.tv_sec + (float)(diceend.tv_usec - dicestart.tv_usec)/1e6;

		if (protneut == NEUTRON){ // if neutron should be simulated
			neutron->Integrate(geom, mc, field, endlog, tracklog, snap, &snapshots, reflektlog, reflectlog);
			kennz_counter[NEUTRON][neutron->kennz]++; // increment counters
			ntotalsteps += neutron->nsteps;
			IntegratorTime += neutron->comptime;
			ReflTime += neutron->refltime;
		}
		
		if(protneut == PROTON || protneut == ELECTRON ||
			(neutron->kennz == KENNZAHL_DECAYED && decay == 2)) // if only p/e should be simulated or neutron decayed
		{	
			long double v_p[3], v_e[3];
			mc.MCZerfallsstartwerte(&neutron->yend[3], v_p, v_e); // get velocity spectrum from random number generator

			if (protneut != ELECTRON){ // if not only electrons should be simulated
				TProton p(neutron, v_p, mc, field); // create proton
				p.Integrate(geom, mc, field, endlog, tracklog); // integrate proton
				kennz_counter[2][p.kennz]++; // increment counters
				ntotalsteps += p.nsteps;
				IntegratorTime += p.comptime;
				ReflTime += p.refltime;
			}

			if (protneut != PROTON){ // if not only protons should be simulated
				TElectron e(neutron, v_e, mc, field); // create electron
				e.Integrate(geom, mc, field, endlog, tracklog); // integrate electron
				kennz_counter[0][e.kennz]++; // increment counters
				ntotalsteps += e.nsteps;
				IntegratorTime += e.comptime;
				ReflTime += e.refltime;
			}
			
		}
		delete neutron;
	}


	OutputCodes(kennz_counter); // print particle kennzahlen
	
	// print statistics
	printf("The integrator made %d steps. \n", ntotalsteps);
	gettimeofday(&simend, NULL);
	float SimulationTime = simend.tv_sec - simstart.tv_sec + (float)(simend.tv_usec - simstart.tv_usec)/1e6;
	printf("Init: %.2fs, Simulation: %.2fs, Integrator: %.2fs (%.2f%%) including Reflection: %.2fs (%.f%%), Dicing: %.4fs\n",
			InitTime, SimulationTime, IntegratorTime, IntegratorTime*100/SimulationTime, ReflTime, ReflTime*100/SimulationTime, DiceTime);
	printf("That's it... Have a nice day!\n");
	

	ostringstream fileprefix;
	fileprefix << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0);
	if (neutdist == 1) outndist((fileprefix.str() + "ndist.out").c_str());   // print neutron distribution into file

	// close files
	if (tracklog)
		fclose(tracklog);
	if (endlog)
		fclose(endlog);
	if (snap)
		fclose(snap);
	if (reflectlog)
		fclose(reflectlog);
		
	return 0;
}

/**
 * Read config file.
 */
void ConfigInit(void){
	/* setting default values */
	protneut = NEUTRON;
	polarisation = POLARISATION_GOOD;
	neutdist = 0;
	ausgabewunsch = 2;
	MonteCarloAnzahl = 1;
	/*end default values*/
	
	/* read lines in config.in into map */
	ifstream infile((inpath+"/config.in").c_str());
	map<string, map<string, string> > config;
	char c;
	string rest,section,key;
	while (infile.good() && (infile >> ws) && (c = infile.peek())){
		if (c == '[' && infile.ignore()){
			if (infile.peek() == '/'){
				section = "";
				printf("\n");
			}
			else{
				getline(infile, section, ']');
				printf("section : %s\n",section.c_str());
			}
			getline(infile,rest);
		}
		else if (c == '#')
			getline(infile,rest);
		else if (section != ""){
			infile >> key;
			getline(infile,config[section][key]);
			printf("%s:%s\n", key.c_str(), config[section][key].c_str());
		}
		else
			getline(infile,rest);
	}
	
	/* read variables from map by casting strings in map into istringstreams and extracting value with ">>"-operator */
	istringstream(config["global"]["protneut"])				>> protneut;
	istringstream(config["global"]["polarisation"])			>> polarisation;
	istringstream(config["global"]["neutdist"])				>> neutdist;
	istringstream(config["global"]["ausgabewunsch"])		>> ausgabewunsch;
	
	int snapshot;
	istringstream(config["global"]["snapshot"])				>> snapshot;
	if (snapshot){
		float tmp_snapshot;
		istringstream isnapshots(config["global"]["snapshots"]);	
		do{
			isnapshots >> tmp_snapshot;
			if (isnapshots.good()) snapshots.push_back(tmp_snapshot);
		}while(isnapshots.good());
	}
	
	int BFstart, BFend;
	istringstream iBruteForce(config["global"]["BruteForce"]);
	do{
		iBruteForce >> BFstart >> BFend;
		if (iBruteForce.good()){
			BFtimes.push_back(BFstart);
			BFtimes.push_back(BFend);
		}
	}while(iBruteForce.good());

	istringstream(config["global"]["reflektlog"])			>> reflektlog;
	istringstream(config["global"]["MonteCarloAnzahl"])		>> MonteCarloAnzahl;
	istringstream(config["global"]["storagetime"])			>> StorageTime;
	istringstream(config["global"]["BFTargetB"])			>> BFTargetB;
	istringstream(config["global"]["decay"])				>> decay;
	istringstream(config["global"]["tau"])					>> tau_n;	
	istringstream(config["global"]["decayoffset"])			>> decayoffset;	
	istringstream(config["global"]["BCutPlane"])			>> BCutPlanePoint[0] >> BCutPlanePoint[1] >> BCutPlanePoint[2] 
															>> BCutPlanePoint[3] >> BCutPlanePoint[4] >> BCutPlanePoint[5]
															>> BCutPlanePoint[6] >> BCutPlanePoint[7] >> BCutPlanePoint[8]
															>> BCutPlaneSampleCount1 >> BCutPlaneSampleCount2;
}

/**
 * Open output files.
 *
 * @param endlog Returns FILE-pointer to particles' end point log
 * @param tracklog Returns FILE-pointer to particles' trajectory log (if specified in config)
 * @param snap Return FILE-pointer to print snapshots at specified times (if specified in config)
 * @param reflectlog Return FILE-pointer to reflection and transmission log (if specified in config)
 */
void OpenFiles(FILE *&endlog, FILE *&tracklog, FILE *&snap, FILE *&reflectlog){
	// **************** create log files ****************
	ostringstream fileprefix;
	fileprefix << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0);

	// open endlog
	endlog = fopen((fileprefix.str() + "end.out").c_str(),"w");
    if (!endlog){
    	printf("\n%s not found!\n",(fileprefix.str() + "end.out").c_str());
    	exit(-1);
    }
    // print endlog header
	fprintf(endlog,"jobnumber teilchen protneut polarisation "
                   "tstart xstart ystart zstart "
                   "rstart phistart vstart alphastart gammastart "
                   "Hstart Estart "
                   "tend xend yend zend "
                   "rend phiend vend alphaend gammaend dt "
                   "Hend Eend kennz NSF ComputingTime BFflipprob "
                   "AnzahlRefl vladmax vladtotal thumbmax trajlength "
                   "Hmax tauSF dtau\n");

	if (protneut == NEUTRON && snapshots.size() > 0){
		// open snapshot log
		snap = fopen((fileprefix.str() + "snapshot.out").c_str(),"w");
	    if (!snap){
	    	printf("\n%s not found!\n",(fileprefix.str() + "snapshot.out").c_str());
	    	exit(-1);
	    }
	    // print snapshot log header
		fprintf(snap,"jobnumber teilchen protneut polarisation "
	                   "tstart xstart ystart zstart "
	                   "rstart phistart vstart alphastart gammastart "
	                   "Hstart Estart "
	                   "tend xend yend zend "
	                   "rend phiend vend alphaend gammaend dt "
	                   "Hend Eend kennz NSF ComputingTime BFflipprob "
	                   "AnzahlRefl vladmax vladtotal thumbmax trajlength "
	                   "Hmax tauSF dtau\n");
	}

	if ((ausgabewunsch==OUTPUT_EVERYTHING)||(ausgabewunsch==OUTPUT_EVERYTHINGandSPIN)){
		// open track log
		tracklog = fopen((fileprefix.str() + "track.out").c_str(),"w");
	    if (!tracklog){
	    	printf("\n%s not found!\n",(fileprefix.str() + "track.out").c_str());
	    	exit(-1);
	    }
	    // print track log header
		fprintf(tracklog,"Teilchen protneut polarisation t x y z dxdt dydt dzdt "
							"v H E Bx dBxdx dBxdy dBxdz By dBydx "
							 "dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V "
							 "timestep logvlad logthumb logBF\n");
	}

	if (protneut == NEUTRON && reflektlog){
		// open reflect log
		reflectlog = fopen((fileprefix.str() + "reflect.out").c_str(),"w");
	    if (!reflectlog){
	    	printf("\n%s not found!\n",(fileprefix.str() + "reflect.out").c_str());
	    	exit(-1);
	    }
	    // print reflect log header
		fprintf(reflectlog,"jobnumber teilchen protneut polarisation reflection solid t x y z vx vy vz nx ny nz H winkeben winksenkr Transprob\n"); // Header for Reflection File
	}
	// **************** end create log files ****************	
}


/**
 * Print final particles statistics.
 *
 * @param kennz_counter Array of three vectors (for n, p and e) containing the count of particles with each kennzahl
 */
//======== Output of the endcodes ==========================================================================================
void OutputCodes(vector<int> kennz_counter[3]){
	int ncount = accumulate(kennz_counter[1].begin(),kennz_counter[1].end(),0);
	int pcount = accumulate(kennz_counter[2].begin(),kennz_counter[2].end(),0);
	int ecount = accumulate(kennz_counter[0].begin(),kennz_counter[0].end(),0);
	printf("\nThe calculations of %i particle(s) yielded:\n"
	       "endcode:  of %4i neutron(s) ; of %4i proton(s) ; of %4i electron(s)\n"
	       "   0 %12i %20i %19i 		(were not categorized)\n"
	       "   1 %12i %20i %19i 		(did not finish)\n"
	       "   2 %12i %20i %19i 		(hit outer boundaries)\n"
	       "   3 %12i %20i %19i 		(produced a numerical error)\n"
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
		printf("  %2i %12i %20i %19i		(were statistically absorbed)\n",
				i,kennz_counter[1][i],kennz_counter[2][i],kennz_counter[0][i]);
	}
}


/**
 * Print planar slice of fields into a file.
 *
 * The slice plane is given by three points BCutPlayPoint[0..8] on the plane
 *
 * @param outfile filename of result file
 * @param field TField structure which should be evaluated
 */
void PrintBFieldCut(const char *outfile, TField &field){
	// get directional vectors from points on plane by u = p2-p1, v = p3-p1
	long double u[3] = {BCutPlanePoint[3] - BCutPlanePoint[0], BCutPlanePoint[4] - BCutPlanePoint[1], BCutPlanePoint[5] - BCutPlanePoint[2]};
	long double v[3] = {BCutPlanePoint[6] - BCutPlanePoint[0], BCutPlanePoint[7] - BCutPlanePoint[1], BCutPlanePoint[8] - BCutPlanePoint[2]};
	
	// open output file
	FILE *cutfile = fopen(outfile, "w");
	if (!cutfile){
		printf("Could not open %s!",outfile);
		exit(-1);
	}
	// print file header
	fprintf(cutfile, "x y z Bx dBxdx dBxdy dBxdz By dBydx dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V\n");
	
	long double Pp[3],B[4][4],Ei[3],V;
	float start = clock(); // do some time statistics
	// sample field BCutPlaneSmapleCount1 times in u-direction and BCutPlaneSampleCount2 time in v-direction
	for (int i = 0; i < BCutPlaneSampleCount1; i++) {
		for (int j = 0; j < BCutPlaneSampleCount2; j++){
			for (int k = 0; k < 3; k++)
				Pp[k] = BCutPlanePoint[k] + i*u[k]/BCutPlaneSampleCount1 + j*v[k]/BCutPlaneSampleCount2;
			// print B-/E-Field to file
			fprintf(cutfile, "%LG %LG %LG ", Pp[0],Pp[1],Pp[2]);
			
			field.BFeld(Pp[0], Pp[1], Pp[2], 0, B);
			for (int k = 0; k < 4; k++)
				for (int l = 0; l < 4; l++)
					fprintf(cutfile, "%LG ",B[k][l]);

			field.EFeld(Pp[0], Pp[1], Pp[2], V, Ei);
			fprintf(cutfile, "%LG %LG %LG %LG\n",
							  Ei[0],Ei[1],Ei[2],V);
		}
	}
	start = (clock() - start)/CLOCKS_PER_SEC;
	//close file
	fclose(cutfile);
	// print time statistics
	printf("Called BFeld and EFeld %u times in %fs (%fms per call)\n",BCutPlaneSampleCount1*BCutPlaneSampleCount2, start, start/BCutPlaneSampleCount1/BCutPlaneSampleCount2*1000);
}


/**
 * Ramp Heating Analysis.
 *
 * "Count" phase space for each energy bin and calculate "heating" of the neutrons due to
 * phase space compression by magnetic field ramping
 *
 * @param outfile Filename of output file
 * @param field TField structure which should be evaluated
 */
void PrintBField(const char *outfile, TField &field){
	// print BField to file
	FILE *bfile = fopen(outfile, "w");
	if (!bfile){
		printf("Could not open %s!",outfile);
		exit(-1);
	}

	fprintf(bfile,"r phi z Bx By Bz 0 0 Babs\n");
	long double rmin = 0.12, rmax = 0.5, zmin = 0, zmax = 1.2;
	int E;
	const int Emax = 108;
	long double dr = 0.1, dz = 0.1;
	long double VolumeB[Emax];
	for (E = 0; E <= Emax; E++) VolumeB[E] = 0;
	
	long double EnTest, B[4][4];
	// sample space in cylindrical pattern
	for (long double r = rmin; r <= rmax; r += dr){
		for (long double z = zmin; z <= zmax; z += dz){
			field.BFeld(r, 0, z, 500.0, B); // evaluate field
			// print field values
			fprintf(bfile,"%LG %G %LG %LG %LG %LG %G %G %LG \n",r,0.0,z,B[0][0],B[1][0],B[2][0],0.0,0.0,B[3][0]);
			printf("r=%LG, z=%LG, Br=%LG T, Bz=%LG T\n",r,z,B[0][0],B[2][0]);
			
			// Ramp Heating Analysis
			for (E = 0; E <= Emax; E++){
				EnTest = E*1.0e-9 - m_n*gravconst*z - mu_nSI/ele_e * B[3][0];
				if (EnTest >= 0){
					// add the volume segment to the volume that is accessible to a neutron with energy Energie
					VolumeB[E] = VolumeB[E] + pi * dz * ((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
				}
			}
		}
	}

	// for investigating ramp heating of neutrons, volume accessible to neutrons with and
	// without B-field is calculated and the heating approximated by thermodynamical means
	printf("\nEnergie [neV], Volumen ohne B-Feld, mit B-Feld, 'Erwaermung'");
	long double Volume;
	for (E = 0; E <= Emax; E++) 
	{
		Volume = ((E * 1.0e-9 / (m_n * gravconst))) * pi * (rmax*rmax-rmin*rmin);
		// isentropische zustandsnderung, kappa=5/3
		printf("\n%i %.17LG %.17LG %.17LG",E,Volume,VolumeB[E],E * powl((Volume/VolumeB[E]),(2.0/3.0)) - E);
	}
}


/**
 * Sample geometry randomly to visualize it.
 *
 * Creates random line segments and prints every intersection point with a surface
 * into outfile
 *
 * @param outfile File name of output file
 * @param geom TGeometry structure which shall be sampled
 */
void PrintGeometry(const char *outfile, TGeometry &geom){
    long double p1[3], p2[3];
    long double theta, phi;
    // create count line segments with length raylength
    unsigned count = 1000000, collcount = 0, raylength = 1;

    ofstream f(outfile);
    f << "x y z" << endl; // print file header

    list<TCollision> c;
    srand(time(NULL));
    float timer = clock(), colltimer = 0;
    for (unsigned i = 0; i < count; i++){
    	// random segment start point
        p1[0] = (((long double)rand())/RAND_MAX) * (geom.kdtree->hi[0] - geom.kdtree->lo[0]) + geom.kdtree->lo[0];
        p1[1] = (((long double)rand())/RAND_MAX) * (geom.kdtree->hi[1] - geom.kdtree->lo[1]) + geom.kdtree->lo[1];
        p1[2] = (((long double)rand())/RAND_MAX) * (geom.kdtree->hi[2] - geom.kdtree->lo[2]) + geom.kdtree->lo[2];
		// random segment direction
        theta = (long double)rand()/RAND_MAX*pi;
		phi = (long double)rand()/RAND_MAX*2*pi;
		// translate direction and length into segment end point
		p2[0] = p1[0] + raylength*sin(theta)*cos(phi);
		p2[1] = p1[1] + raylength*sin(theta)*sin(phi);
		p2[2] = p1[2] + raylength*cos(theta);

		timeval collstart,collend;
		gettimeofday(&collstart,NULL);
		if (geom.kdtree->Collision(p1,p2,c)){ // check if segment intersected with surfaces
			collcount++;
			for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){ // print all intersection points into file
				f << p1[0] + i->s*(p2[0]-p1[0]) << " " << p1[1] + i->s*(p2[1] - p1[1]) << " " << p1[2] + i->s*(p2[2] - p1[2]) << endl;
			}
		}
		gettimeofday(&collend,NULL);
		colltimer += (collend.tv_sec - collstart.tv_sec)*1e6+collend.tv_usec - collstart.tv_usec;
    }
    timer = (clock() - timer)*1000/CLOCKS_PER_SEC;
    // print some time statistics
    printf("%u tests, %u collisions in %fms (%fms per Test, %fms per Collision)\n",count,collcount,timer,timer/count,colltimer/collcount/1000);
    f.close();	
}
