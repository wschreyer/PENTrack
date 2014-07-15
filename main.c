/**
 * \file
 * Main program.
 *
 * Create particles to your liking...
 */

#include <cstdlib>
#include <cstdio>
#include <csignal>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sys/time.h>

#include "nr/nr3.h"
#include "nr/interp_1d.h"
#include "nr/interp_linear.h"
#include "nr/interp_2d.h"
#include "nr/odeint.h"
#include "nr/stepper.h"
#include "nr/stepperdopr853.h"

#include "particle.h"
#include "neutron.h"
#include "proton.h"
#include "electron.h"
#include "globals.h"
#include "fields.h"
#include "geometry.h"
#include "source.h"
#include "mc.h" 
#include "bruteforce.h"
#include "ndist.h"

void ConfigInit(TConfig &config); // read config.in
void OutputCodes(map<string, map<int, int> > &ID_counter); // print simulation summary at program exit
void PrintBFieldCut(const char *outfile, TFieldManager &field); // evaluate fields on given plane and write to outfile
void PrintBField(const char *outfile, TFieldManager &field);
void PrintGeometry(const char *outfile, TGeometry &geom); // do many random collisionchecks and write all collisions to outfile


long double SimTime = 1500.; ///< max. simulation time
int simcount = 1; ///< number of particles for MC simulation (read from config)
int simtype = PARTICLE; ///< type of particle which shall be simulated (read from config)
int secondaries = 1; ///< should secondary particles be simulated? (read from config)
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
	if ((argc > 1) && (strcmp(argv[1], "-h") == 0)){
		cout << "Usage:\nPENTrack [jobnumber [path/to/in/files [path/to/out/files]]]" << endl;
		exit(0);
	}

	//Initialize signal-analizing
	signal (SIGUSR1, catch_alarm);
	signal (SIGUSR2, catch_alarm);
	signal (SIGXCPU, catch_alarm);   
	
	jobnumber = 0;
	outpath = "./out";
	string inpath = "./in";
	if(argc>1) // if user supplied at least 1 arg (outputfilestamp)
		istringstream(argv[1]) >> jobnumber;
	if(argc>2) // if user supplied 2 or more args (outputfilestamp, inpath)
		inpath = argv[2]; // input path pointer set
	if(argc>3) // if user supplied all 3 args (outputfilestamp, inpath, outpath)
		outpath = argv[3]; // set the output path pointer
	
	TConfig configin;
	ReadInFile(string(inpath + "/config.in").c_str(), configin);
	TConfig geometryin;
	ReadInFile(string(inpath + "/geometry.in").c_str(), geometryin);
	TConfig particlein;
	ReadInFile(string(inpath + "/particle.in").c_str(), particlein); // read particle specific log configuration from particle.in
	for (TConfig::iterator i = particlein.begin(); i != particlein.end(); i++){
		if (i->first != "all"){
			i->second = particlein["all"]; // set all particle specific settings to the "all" settings
		}
	}
	ReadInFile(string(inpath+"/particle.in").c_str(), particlein); // read again to overwrite "all" settings with particle specific settings

	// read config.in
	ConfigInit(configin);

	if(simtype == PARTICLE){
		if (neutdist == 1) prepndist(); // prepare for neutron distribution-calculation
	}
	
	cout << "Loading fields...\n";
	// load field configuration from geometry.in
	TFieldManager field(geometryin);

	switch(simtype)
	{
		case BF_ONLY:	PrintBField(string(outpath+"/BF.out").c_str(), field); // estimate ramp heating
						return 0;
		case BF_CUT:	PrintBFieldCut(string(outpath+"/BFCut.out").c_str(), field); // print cut through B field
						return 0;
	}


	cout << "Loading geometry...\n";
	//load geometry configuration from geometry.in
	TGeometry geom(geometryin);
	
	if (simtype == GEOMETRY){
		// print random points on walls in file to visualize geometry
		PrintGeometry(string(outpath+"/geometry.out").c_str(), geom);
		return 0;
	}
	
	cout << "Loading source...\n";
	// load source configuration from geometry.in
	TSource source(geometryin, geom, field);
	
	cout << "Loading random number generator...\n";
	// load random number generator from all3inone.in
	TMCGenerator mc(string(inpath + "/particle.in").c_str());
	
	int ntotalsteps = 0;     // counters to determine average steps per integrator call
	float InitTime = (1.*clock())/CLOCKS_PER_SEC, DiceTime = 0, IntegratorTime = 0, ReflTime = 0; // time statistics

	// simulation time counter
	timespec simstart, simend;
	clock_gettime(CLOCK_REALTIME, &simstart);

	printf(
	" ########################################################################\n"
	" ###                      Welcome to PENTrack,                        ###\n"
	" ### a simulation tool for ultra-cold neutrons, protons and electrons ###\n"
	" ########################################################################\n");

	map<string, map<int, int> > ID_counter; // 2D map to store number of each ID for each particle type

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
		ID_counter[particle.protneut % 3][particle.ID]++; // increase ID-counter
		ntotalsteps += particle.nsteps;
		IntegratorTime += particle.comptime;
		ReflTime += particle.refltime;
		infile >> ws;
	}
*/
	if (simtype == PARTICLE){ // if proton or neutron shall be simulated
		for (int iMC = 1; iMC <= simcount; iMC++)
		{
			TParticle *p = source.CreateParticle(mc, geom, &field);
			p->Integrate(SimTime, geom, mc, &field, particlein[p->name]); // integrate particle
			ID_counter[p->name][p->ID]++; // increment counters
			ntotalsteps += p->Nstep;
			IntegratorTime += p->inttime;
			ReflTime += p->refltime;

			if (secondaries == 1){
				for (vector<TParticle*>::iterator i = p->secondaries.begin(); i != p->secondaries.end(); i++){
					(*i)->Integrate(SimTime, geom, mc, &field, particlein[(*i)->name]); // integrate secondary particles
					ID_counter[(*i)->name][(*i)->ID]++;
					ntotalsteps += (*i)->Nstep;
					IntegratorTime += (*i)->inttime;
					ReflTime += (*i)->refltime;
				}
			}

			delete p;
		}
	}
	else{
		printf("\nDon't know simtype %i! Exiting...\n",simtype);
		exit(-1);
	}



	OutputCodes(ID_counter); // print particle IDs
	
	// print statistics
	printf("The integrator made %d steps. \n", ntotalsteps);
	clock_gettime(CLOCK_REALTIME, &simend);
	float SimulationTime = simend.tv_sec - simstart.tv_sec + (float)(simend.tv_nsec - simstart.tv_nsec)/1e9;
	printf("Init: %.2fs, Simulation: %.2fs, Integrator: %.2fs (%.2f%%), Reflection: %.2fs (%.f%%), Dicing: %.4fs\n",
			InitTime, SimulationTime, IntegratorTime, IntegratorTime*100/SimulationTime, ReflTime, ReflTime*100/SimulationTime, DiceTime);
	printf("That's it... Have a nice day!\n");
	

	ostringstream fileprefix;
	fileprefix << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0);
	if (neutdist == 1) outndist((fileprefix.str() + "ndist.out").c_str());   // print neutron distribution into file

	return 0;
}


/**
 * Read config file.
 *
 * @param config TConfig struct containing [global] options map
 */
void ConfigInit(TConfig &config){
	/* setting default values */
	simtype = PARTICLE;
	neutdist = 0;
	simcount = 1;
	/*end default values*/

	/* read variables from map by casting strings in map into istringstreams and extracting value with ">>"-operator */
	istringstream(config["global"]["simtype"])		>> simtype;
	istringstream(config["global"]["neutdist"])		>> neutdist;
	

	istringstream(config["global"]["simcount"])		>> simcount;
	istringstream(config["global"]["simtime"])		>> SimTime;
	istringstream(config["global"]["secondaries"])	>> secondaries;
	istringstream(config["global"]["BCutPlane"])	>> BCutPlanePoint[0] >> BCutPlanePoint[1] >> BCutPlanePoint[2]
													>> BCutPlanePoint[3] >> BCutPlanePoint[4] >> BCutPlanePoint[5]
													>> BCutPlanePoint[6] >> BCutPlanePoint[7] >> BCutPlanePoint[8]
													>> BCutPlaneSampleCount1 >> BCutPlaneSampleCount2;
}


/**
 * Print final particles statistics.
 */
void OutputCodes(map<string, map<int, int> > &ID_counter){
	cout << "\nThe simulated particles suffered following fates:\n";
	for (map<string, map<int, int> >::iterator i = ID_counter.begin(); i != ID_counter.end(); i++){
		map<int, int> counts = i->second;
		const char *name = i->first.c_str();
		printf("%4i: %6i %10s(s) were not categorized\n",		 0, counts[ 0], name);
		printf("%4i: %6i %10s(s) did not finish\n",				-1, counts[-1], name);
		printf("%4i: %6i %10s(s) hit outer boundaries\n",		-2, counts[-2], name);
		printf("%4i: %6i %10s(s) produced a numerical error\n", -3, counts[-3], name);
		printf("%4i: %6i %10s(s) decayed\n",					-4, counts[-4], name);
		printf("%4i: %6i %10s(s) found no initial position\n",	-5, counts[-5], name);
		for (map<int, int>::iterator j = counts.begin(); j != counts.end(); j++){
			if (j->first > 0)
				printf("%4i: %6i %10s(s) were absorbed\n", j->first, j->second, name);
		}
		printf("\n");
	}
}


/**
 * Print planar slice of fields into a file.
 *
 * The slice plane is given by three points BCutPlayPoint[0..8] on the plane
 *
 * @param outfile filename of result file
 * @param field TFieldManager structure which should be evaluated
 */
void PrintBFieldCut(const char *outfile, TFieldManager &field){
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
			
			field.BField(Pp[0], Pp[1], Pp[2], 0, B);
			for (int k = 0; k < 4; k++)
				for (int l = 0; l < 4; l++)
					fprintf(cutfile, "%LG ",B[k][l]);

			field.EField(Pp[0], Pp[1], Pp[2], 0, V, Ei);
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
void PrintBField(const char *outfile, TFieldManager &field){
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
	long double VolumeB[Emax + 1];
	for (E = 0; E <= Emax; E++) VolumeB[E] = 0;
	
	long double EnTest, B[4][4];
	// sample space in cylindrical pattern
	for (long double r = rmin; r <= rmax; r += dr){
		for (long double z = zmin; z <= zmax; z += dz){
			field.BField(r, 0, z, 500.0, B); // evaluate field
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
		printf("\n%i %.17LG %.17LG %.17LG",E,Volume,VolumeB[E],E * pow((Volume/VolumeB[E]),(2.0L/3.0)) - E);
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
    f << "x y z" << '\n'; // print file header

    set<TCollision> c;
    srand(time(NULL));
	timespec collstart,collend;
	clock_gettime(CLOCK_REALTIME, &collstart);
	for (unsigned i = 0; i < count; i++){
    	// random segment start point
        for (int j = 0; j < 3; j++)
        	p1[j] = (long double)rand()/RAND_MAX * (geom.mesh.tree.bbox().max(j) - geom.mesh.tree.bbox().min(j)) + geom.mesh.tree.bbox().min(j);
		// random segment direction
        theta = (long double)rand()/RAND_MAX*pi;
		phi = (long double)rand()/RAND_MAX*2*pi;
		// translate direction and length into segment end point
		p2[0] = p1[0] + raylength*sin(theta)*cos(phi);
		p2[1] = p1[1] + raylength*sin(theta)*sin(phi);
		p2[2] = p1[2] + raylength*cos(theta);

		if (geom.mesh.Collision(p1,p2,c)){ // check if segment intersected with surfaces
			collcount++;
			for (set<TCollision>::iterator i = c.begin(); i != c.end(); i++){ // print all intersection points into file
				f << p1[0] + i->s*(p2[0]-p1[0]) << " " << p1[1] + i->s*(p2[1] - p1[1]) << " " << p1[2] + i->s*(p2[2] - p1[2]) << '\n';
			}
		}
    }
	clock_gettime(CLOCK_REALTIME, &collend);
	float colltimer = (collend.tv_sec - collstart.tv_sec)*1e9 + collend.tv_nsec - collstart.tv_nsec;
    // print some time statistics
    printf("%u tests, %u collisions in %fms (%fms per Test, %fms per Collision)\n",count,collcount,colltimer,colltimer/count/1e6,colltimer/collcount/1e6);
    f.close();	
}
