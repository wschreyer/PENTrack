/**
 * \file
 * Main program.
 *
 * Create particles to your liking...
 */

#include <csignal>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <boost/format.hpp>

using namespace std;

#include "particle.h"
#include "config.h"
#include "fields.h"
#include "geometry.h"
#include "source.h"
#include "mc.h" 
#include "microroughness.h"

TConfig ConfigInit(int argc, char **argv); // read config.in
void OutputCodes(const map<string, map<int, int> > &ID_counter); // print simulation summary at program exit
void PrintBFieldCut(TConfig &config, const boost::filesystem::path &outfile, const TFieldManager &field); // evaluate fields on given plane and write to outfile
void PrintBField(const boost::filesystem::path &outfile, const TFieldManager &field);
void PrintGeometry(const boost::filesystem::path &outfile, TGeometry &geom); // do many random collisionchecks and write all collisions to outfile
void PrintMROutAngle(TConfig &config, const boost::filesystem::path &outpath); // produce a 3d table of the MR-DRP for each outgoing solid angle
void PrintMRThetaIEnergy(TConfig &config, const boost::filesystem::path &outpath); // produce a 3d table of the total (integrated) MR-DRP for a given incident angle and energy


double SimTime = 1500.; ///< max. simulation time
int simcount = 1; ///< number of particles for MC simulation (read from config)
simType simtype = PARTICLE; ///< type of particle which shall be simulated (read from config)
int secondaries = 1; ///< should secondary particles be simulated? (read from config)
uint64_t seed = 0; ///< random seed used for random-number generator (generated from high-resolution clock)

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
 * @param argv Array of parameters passed via the command line (./Track [jobnumber [configpath [outputpath [seed]]]])
 * @return Return 0 on success, value !=0 on failure
 *
 */
int main(int argc, char **argv){
	if ((argc > 1) && (strcmp(argv[1], "-h") == 0)){
		cout << "Usage:\nPENTrack [jobnumber [location/of/config.in [path/to/out/files [seed]]]]" << endl;
		return 0;
	}

	cout <<
	" #######################################################################\n"
	" ###                      Welcome to PENTrack,                       ###\n"
	" ### a simulation tool for ultracold neutrons, protons and electrons ###\n"
	" #######################################################################\n";

	//Initialize signal-analizing
	signal (SIGINT, catch_alarm);
	signal (SIGUSR1, catch_alarm);
	signal (SIGUSR2, catch_alarm);
	signal (SIGXCPU, catch_alarm);
	
	// read config
	TConfig configin = ConfigInit(argc, argv);

	if (simtype == MR_THETA_OUT_ANGLE){
		PrintMROutAngle(configin, outpath);
		return 0;
	}
	else if (simtype == MR_THETA_I_ENERGY){
		PrintMRThetaIEnergy(configin, outpath);
		return 0;
	}


	cout << "Loading fields...\n";
	// load field configuration from geometry.in
	TFieldManager field(configin);

	if (simtype == BF_ONLY){
		PrintBField(outpath / "BF.out", field); // estimate ramp heating
		return 0;
	}
	else if (simtype == BF_CUT){
		PrintBFieldCut(configin, outpath / "BFCut.out", field); // print cut through B field
		return 0;
	}


	cout << "Loading geometry...\n";
	//load geometry configuration from geometry.in
	TGeometry geom(configin);
	
	if (simtype == GEOMETRY){
		// print random points on walls in file to visualize geometry
		PrintGeometry(outpath / "geometry.out", geom);
		return 0;
	}
	
	cout << "Loading random number generator...\n";
	// load random number generator
	TMCGenerator mc;
	if (seed == 0){
		// get high-resolution timestamp to generate seed
		using namespace std::chrono;
		seed = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
	}
	std::cout << "Random Seed: " << seed << "\n\n";
	mc.seed(seed);


	
	cout << "Loading source...\n";
	// load source configuration from geometry.in
	unique_ptr<TParticleSource> source(CreateParticleSource(configin, geom));

	int ntotalsteps = 0;     // counters to determine average steps per integrator call
	float InitTime = (1.*clock())/CLOCKS_PER_SEC; // time statistics

	// simulation time counter
	chrono::time_point<chrono::steady_clock> simstart = chrono::steady_clock::now();

	map<string, map<int, int> > ID_counter; // 2D map to store number of each ID for each particle type

	if (simtype == PARTICLE){ // if proton or neutron shall be simulated
		for (int iMC = 1; iMC <= simcount; iMC++)
		{
			TParticle *p = source->CreateParticle(mc, geom, field);
			p->Integrate(SimTime, configin[p->GetName()], mc, geom, field); // integrate particle
			ID_counter[p->GetName()][p->GetStopID()]++; // increment counters
			ntotalsteps += p->GetNumberOfSteps();

			if (secondaries == 1){
				std::vector<TParticle*> secondaries = p->GetSecondaryParticles();
				for (auto i = secondaries.begin(); i != secondaries.end(); i++){
					(*i)->Integrate(SimTime, configin[(*i)->GetName()], mc, geom, field); // integrate secondary particles
					ID_counter[(*i)->GetName()][(*i)->GetStopID()]++;
					ntotalsteps += (*i)->GetNumberOfSteps();
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
	chrono::time_point<chrono::steady_clock> simend = chrono::steady_clock::now();
	float SimulationTime = chrono::duration_cast<chrono::milliseconds>(simend - simstart).count()/1000.;
	printf("Init: %.2fs, Simulation: %.2fs\n",
			InitTime, SimulationTime);
	printf("That's it... Have a nice day!\n");
	

	ostringstream fileprefix;
	fileprefix << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0);

	return 0;
}


/**
 * Read config file.
 *
 * @param argc Number of command line parameters
 * @param argv Array command line parameters
 * 
 * @return Returns TConfig struct containing options map
 */
TConfig ConfigInit(int argc, char **argv){
	/* setting default values */
	jobnumber = 0;
	outpath = boost::filesystem::absolute("out/");
	configpath = boost::filesystem::absolute("in/config.in");
	seed = 0;
	simtype = PARTICLE;
	simcount = 1;
	/*end default values*/

	if(argc>1) // if user supplied at least 1 arg (jobnumber)
		istringstream(argv[1]) >> jobnumber;
	if(argc>2){ // if user supplied 2 or more args (jobnumber, configpath)
		configpath = boost::filesystem::absolute(argv[2]); // input path pointer set
		if (boost::filesystem::is_directory(configpath))
			configpath /= "config.in";
	}
	if(argc>3) // if user supplied 3 or more args (jobnumber, configpath, outpath)
		outpath = boost::filesystem::absolute(argv[3]); // set the output path pointer
	if (argc>4) // if user supplied 4 or more args (jobnumber, configpath, outpath, seed)
		istringstream(argv[4]) >> seed;
	
	TConfig config(configpath.native());
	config.convert(configpath.native());

	/* read variables from map by casting strings in map into istringstreams and extracting value with ">>"-operator */
	int stype;
	istringstream(config["GLOBAL"]["simtype"])		>> stype;
	simtype = static_cast<simType>(stype);
	istringstream(config["GLOBAL"]["simcount"])		>> simcount;
	istringstream(config["GLOBAL"]["simtime"])		>> SimTime;
	istringstream(config["GLOBAL"]["secondaries"])	>> secondaries;
	
	// add default parameters from PARTICLES section to each individual particle's parameters
	for (auto i = config["PARTICLES"].begin(); i != config["PARTICLES"].end(); ++i){
		config["neutron"].insert(*i);
		config["proton"].insert(*i);
		config["electron"].insert(*i);
		config["mercury"].insert(*i);
		config["xenon"].insert(*i);
	}
	
	return config;
}

/**
 * 
 * Output a table containing the MR diffuse reflection probability for the specified range of solid angles from the config.in file
 *
 * @param config TConfig class containing parameters
 * @param outpath The file name of the file to which results will be printed
 *  
 * Other params are read in from the config.in file
*/
void PrintMROutAngle(TConfig &config, const boost::filesystem::path &outpath) {
	vector<double> MRSolidAngleDRPParams; ///< params to output the  MR-DRP values for given theta_inc and phi_inc [Fermi potential, incident neutron energy, b (nm), w (nm), theta_i] (read from config)
	istringstream ss(config["GLOBAL"]["MRSolidAngleDRP"]);
	copy(istream_iterator<double>(ss), istream_iterator<double>(), back_inserter(MRSolidAngleDRPParams));
	if (MRSolidAngleDRPParams.size() != 5)
		throw std::runtime_error("Incorrect number of parameters to print micro-roughness distribution!");
	
	/** Create a material struct that defines vacuum (the leaving material) and the material the neutron is being reflected from (the entering material **/
	material matEnter = { "reflection surface material" , MRSolidAngleDRPParams[0], 0, 0, 0, MRSolidAngleDRPParams[2], MRSolidAngleDRPParams[3] };
	material matLeav = { "vacuum material", 0, 0, 0, 0, 0, 0 };
	
	/** Create a solid object that the neutron is leaving and entering based on the materials created in the previous step **/
	solid solEnter = { "path", "reflection solid", matEnter, 2 }; // no ignore times (priority = 2) 
	solid solLeav = { "path", "vacuum solid", matLeav, 1 }; //no ignore times (priority = 1 )
	
	ostringstream oss;
	oss << "MR-SldAngDRP" << "-F" << MRSolidAngleDRPParams[0] << "-En" << MRSolidAngleDRPParams[1] << "-b" << MRSolidAngleDRPParams[2] << "-w" << MRSolidAngleDRPParams[3] << "-th" << MRSolidAngleDRPParams[4] << ".out"; 
 	boost::filesystem::path fileName = outpath / oss.str();
	
	cout << "\nGenerating table of MR diffuse reflection probability for all solid angles in " << fileName << "...\n";	
	
	ofstream mrproboutfile(fileName.c_str());
	if (!mrproboutfile)
		throw ((boost::format("Could not open %1%!") % fileName.c_str()).str());

	//print file header
	mrproboutfile << "phi_out theta_out mrdrp\n";
	double theta_inc = MRSolidAngleDRPParams[4];
	
	//determine neutron velocity corresponding to the energy and create a state_type vector from it
	double vabs = sqrt(2*MRSolidAngleDRPParams[1]*1e-9/m_n);
	double v[3] = {0, vabs*sin(theta_inc), -vabs*cos(theta_inc)};
	double norm[] = { 0, 0, 1 };

	//write the mrprob values to the output file
	for (double phi=-pi; phi<pi; phi+=(2*pi)/100) {
		
		for (double theta=0; theta<pi/2; theta+=(pi/2)/100) {
			//the sin(theta) factor is needed to normalize for different size of surface elements in spherical coordinates
			double mrprob = MR::MRDist(false, false, v, norm, solLeav, solEnter, theta, phi)*sin(theta);
			mrproboutfile << phi << ' ' << theta << ' ' << mrprob << '\n';
		}
		for (double theta=0; theta<pi/2; theta+=(pi/2)/100) {
			//the sin(theta) factor is needed to normalize for different size of surface elements in spherical coordinates
			double mrprob = MR::MRDist(true, false, v, norm, solLeav, solEnter, theta, phi)*sin(theta);
			mrproboutfile << phi << ' ' << pi - theta << ' ' << mrprob << '\n';
		}
	}
} // end PrintMRThetaIEnergy

/**
 * 
 * Output a table in root format giving the total MR DRP for a set of incident theta angles and neutron energy. 
 *
 * @param config TConfig class containing parameters
 * @param outpath The file name to which the results will be printed
*/
void PrintMRThetaIEnergy(TConfig &config, const boost::filesystem::path &outpath) {
	vector<double> MRThetaIEnergyParams; ///< params for which to output the integrated MR-DRP values [Fermi potential, b (nm), w (nm), theta_i_start, theta_i_end, neut_energy_start, neut_energy_end] (read from config.in)
	istringstream ss(config["GLOBAL"]["MRThetaIEnergy"]);
	copy(istream_iterator<double>(ss), istream_iterator<double>(), back_inserter(MRThetaIEnergyParams));
	if (MRThetaIEnergyParams.size() != 7)
		throw std::runtime_error("Incorrect number of parameters to print total micro-roughness-scattering probability!");


	/** Create a material struct that defines vacuum (the leaving material) and the material the neutron is being reflected from (the entering material **/
	material matEnter = { "reflection surface material", MRThetaIEnergyParams[0], 0, 0, 0, MRThetaIEnergyParams[1], MRThetaIEnergyParams[2] };
	material matLeav = { "reflection surface material", 0, 0, 0, 0, 0, 0  };

	/** Create a solid object that the neutron is leaving and entering based on the materials created in the previous step **/
	solid solEnter = { "path", "reflection solid", matEnter, 2 }; // no ignore times (priority = 2) 
	solid solLeav = { "path", "vacuum solid", matLeav, 1 }; //no ignore times (priority = 1 )
	
	ostringstream oss;
	oss << "MR-Tot-DRP" << "-F" << MRThetaIEnergyParams[0] << "-b" << MRThetaIEnergyParams[1] << "-w" << MRThetaIEnergyParams[2] << ".out"; 
 	boost::filesystem::path fileName = outpath / oss.str();	
	
	cout << "\nGenerating table of integrated MR diffuse reflection probability for different incident angle and energy in " << fileName << "...\n";	

	ofstream mroutfile(fileName.c_str());
	if (!mroutfile)
		throw ((boost::format("Could not open %1%!") % fileName).str());

	mroutfile << "theta_i neut_en totmrdrp\n"; //print the header
	
	//define the min and max values of the following for loop
	double theta_start = MRThetaIEnergyParams[3];
	double theta_end = MRThetaIEnergyParams[4];
	double neute_start = MRThetaIEnergyParams[5];
	double neute_end = MRThetaIEnergyParams[6];
	int prevProg=0;
	double norm[] = { 0, 0, 1 };
 	//write the integrated mrprob values to the output file 
	for (double theta = theta_start; theta<theta_end; theta += (theta_end-theta_start)/100) {
		//since this part can be slow it is helpful to monitor the progress
		PrintPercent((theta - theta_start)/(theta_end - theta_start), prevProg);

		for (double energy = neute_start; energy<neute_end; energy += (neute_end-neute_start)/100) {
			//determine neutron velocity corresponding to the energy and create a state_type vector from it
			double vabs = sqrt(2*energy*1e-9/m_n);
			double v[3] = {0, vabs*sin(theta), -vabs*cos(theta)};
			//the sin(theta) factor is needed to noramlize for different sizes of surface elements in spherical coordinates
			double totmrprob = MR::MRProb(false, v, norm, solLeav, solEnter);
			mroutfile << theta << ' ' << energy << ' ' << totmrprob << '\n';
		}
	}
	cout << '\n';
} // end PrintMRThetaIEnergy


/**
 * Print final particles statistics.
 *
 * @param ID_counter A list of counters indicating the numbers of particles with each stopID.
 */
void OutputCodes(const map<string, map<int, int> > &ID_counter){
	cout << "\nThe simulated particles suffered following fates:\n";
	for (auto i = ID_counter.begin(); i != ID_counter.end(); i++){
		map<int, int> counts = i->second;
		const char *name = i->first.c_str();
		printf("%4i: %6i %10s(s) were absorbed on a surface\n",	 2, counts[ 2], name);
		printf("%4i: %6i %10s(s) were absorbed in a material\n", 1, counts[ 1], name);
		printf("%4i: %6i %10s(s) were not categorized\n",		 0, counts[ 0], name);
		printf("%4i: %6i %10s(s) did not finish\n",				-1, counts[-1], name);
		printf("%4i: %6i %10s(s) hit outer boundaries\n",		-2, counts[-2], name);
		printf("%4i: %6i %10s(s) produced integration error\n", -3, counts[-3], name);
		printf("%4i: %6i %10s(s) decayed\n",					-4, counts[-4], name);
		printf("%4i: %6i %10s(s) found no initial position\n",	-5, counts[-5], name);
		printf("%4i: %6i %10s(s) encountered CGAL error\n",		-6, counts[-6], name);
		printf("%4i: %6i %10s(s) encountered geometry error\n",	-7, counts[-7], name);
		printf("\n");
	}
}


/**
 * Print planar slice of fields into a file.
 *
 * The slice plane is given by three points BCutPlayPoint[0..8] on the plane
 *
 * @param config TConfig class containing cut parameters
 * @param outfile filename of result file
 * @param field TFieldManager structure which should be evaluated
 */
void PrintBFieldCut(TConfig &config, const boost::filesystem::path &outfile, const TFieldManager &field){
	double BCutPlanePoint[9]; ///< 3 points on plane for field slice (read from config)
	int BCutPlaneSampleCount1; ///< number of field samples in BCutPlanePoint[3..5]-BCutPlanePoint[0..2] direction (read from config)
	int BCutPlaneSampleCount2; ///< number of field samples in BCutPlanePoint[6..8]-BCutPlanePoint[0..2] direction (read from config)
	istringstream(config["GLOBAL"]["BCutPlane"])	>> BCutPlanePoint[0] >> BCutPlanePoint[1] >> BCutPlanePoint[2]
													>> BCutPlanePoint[3] >> BCutPlanePoint[4] >> BCutPlanePoint[5]
													>> BCutPlanePoint[6] >> BCutPlanePoint[7] >> BCutPlanePoint[8]
													>> BCutPlaneSampleCount1 >> BCutPlaneSampleCount2;

	// get directional vectors from points on plane by u = p2-p1, v = p3-p1
	double u[3] = {BCutPlanePoint[3] - BCutPlanePoint[0], BCutPlanePoint[4] - BCutPlanePoint[1], BCutPlanePoint[5] - BCutPlanePoint[2]};
	double v[3] = {BCutPlanePoint[6] - BCutPlanePoint[0], BCutPlanePoint[7] - BCutPlanePoint[1], BCutPlanePoint[8] - BCutPlanePoint[2]};
	
	// open output file
	ofstream cutfile(outfile.c_str());
	if (!cutfile){
		std::cout << "Could not open " << outfile << "!\n";
		exit(-1);
	}
	// print file header
	cutfile << "x y z Bx dBxdx dBxdy dBxdz By dBydx dBydy dBydz Bz dBzdx dBzdy dBzdz Ex Ey Ez V\n";
	
	double Pp[3];
	double B[3], dBidxj[3][3],Ei[3],V;
	float start = clock(); // do some time statistics
	// sample field BCutPlaneSmapleCount1 times in u-direction and BCutPlaneSampleCount2 time in v-direction
	for (int i = 0; i < BCutPlaneSampleCount1; i++) {
		for (int j = 0; j < BCutPlaneSampleCount2; j++){
			for (int k = 0; k < 3; k++)
				Pp[k] = BCutPlanePoint[k] + i*u[k]/BCutPlaneSampleCount1 + j*v[k]/BCutPlaneSampleCount2;
			// print B-/E-Field to file
			cutfile << Pp[0] << " " << Pp[1] << " " << Pp[2] << "\n";
			
			field.BField(Pp[0], Pp[1], Pp[2], 0, B, dBidxj);
			for (int k = 0; k < 3; k++){
				cutfile << B[k] << " ";
				for (int l = 0; l < 3; l++)
					cutfile << dBidxj[k][l] << " ";
			}

			field.EField(Pp[0], Pp[1], Pp[2], 0, V, Ei);
			cutfile << Ei[0] << " " << Ei[1] << " " << Ei[2] << " " << V << "\n";
		}
	}
	start = (clock() - start)/CLOCKS_PER_SEC;
	//close file
	cutfile.close();
	// print time statistics
	printf("\nWrote magnetic and electric fields %u times into %s in %fs (%fms per call)\n",BCutPlaneSampleCount1*BCutPlaneSampleCount2, outfile.c_str(), start, start/BCutPlaneSampleCount1/BCutPlaneSampleCount2*1000);
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
void PrintBField(const boost::filesystem::path &outfile, const TFieldManager &field){
	// print BField to file
	ofstream bfile(outfile.c_str());
	if (!bfile){
		std::cout << "Could not open " << outfile << "!\n";
		exit(-1);
	}

	bfile << "r phi z Bx By Bz 0 0 Babs\n";
	double rmin = 0.12, rmax = 0.5, zmin = 0, zmax = 1.2;
	int E;
	const int Emax = 108;
	double dr = 0.1, dz = 0.1;
	double VolumeB[Emax + 1];
	for (E = 0; E <= Emax; E++) VolumeB[E] = 0;
	
	double EnTest;
	double B[3];
	// sample space in cylindrical pattern
	for (double r = rmin; r <= rmax; r += dr){
		for (double z = zmin; z <= zmax; z += dz){
			field.BField(r, 0, z, 500.0, B); // evaluate field
			double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
			// print field values
			bfile << r << " " << 0.0 << " " << z << " " << B[0] << " " << B[1] << " " << B[2] << " " << 0.0 << " " << 0.0 << " " << Babs << "\n";
			std::cout << "r=" << r << ", z=" << z << ", Br=" << B[0] << " T, Bz=" << B[2] << " T\n";
			
			// Ramp Heating Analysis
			for (E = 0; E <= Emax; E++){
				EnTest = E*1.0e-9 - m_n*gravconst*z - mu_nSI/ele_e * Babs;
				if (EnTest >= 0){
					// add the volume segment to the volume that is accessible to a neutron with energy Energie
					VolumeB[E] = VolumeB[E] + pi * dz * ((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
				}
			}
		}
	}
	bfile.close();

	// for investigating ramp heating of neutrons, volume accessible to neutrons with and
	// without B-field is calculated and the heating approximated by thermodynamical means
	printf("\nEnergie [neV], Volumen ohne B-Feld, mit B-Feld, 'Erwaermung'");
	double Volume;
	for (E = 0; E <= Emax; E++) 
	{
		Volume = ((E * 1.0e-9 / (m_n * gravconst))) * pi * (rmax*rmax-rmin*rmin);
		// isentropische zustandsnderung, kappa=5/3
		printf("\n%i %.17g %.17g %.17g",E,Volume,VolumeB[E],E * pow((Volume/VolumeB[E]),(2.0/3.0)) - E);
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
void PrintGeometry(const boost::filesystem::path &outfile, TGeometry &geom){
    double p1[3], p2[3];
    double theta, phi;
    // create count line segments with length raylength
    unsigned count = 1000000, collcount = 0, raylength = 1;

    ofstream f(outfile.c_str());
    f << "x y z ID" << '\n'; // print file header

    CBox bbox = geom.mesh.GetBoundingBox();
    srand(time(NULL));
	chrono::time_point<chrono::steady_clock> collstart = chrono::steady_clock::now();
	for (unsigned i = 0; i < count; i++){
    	// random segment start point
        for (int j = 0; j < 3; j++)
        	p1[j] = (double)rand()/RAND_MAX * (bbox.max(j) - bbox.min(j)) + bbox.min(j);
		// random segment direction
        theta = (double)rand()/RAND_MAX*pi;
		phi = (double)rand()/RAND_MAX*2*pi;
		// translate direction and length into segment end point
		p2[0] = p1[0] + raylength*sin(theta)*cos(phi);
		p2[1] = p1[1] + raylength*sin(theta)*sin(phi);
		p2[2] = p1[2] + raylength*cos(theta);

	    multimap<TCollision, bool> c;
	    geom.GetCollisions(0,p1,0,p2,c);
		collcount += c.size();

		for (auto i = c.begin(); i != c.end(); i++)
			f << p1[0] + i->first.s*(p2[0]-p1[0]) << " " << p1[1] + i->first.s*(p2[1] - p1[1]) << " " << p1[2] + i->first.s*(p2[2] - p1[2]) << " " << i->first.ID << '\n'; // print all intersection points into file
    }
	chrono::time_point<chrono::steady_clock> collend = chrono::steady_clock::now();
	float colltimer = chrono::duration_cast<chrono::nanoseconds>(collend - collstart).count();
    // print some time statistics
    printf("Wrote %u tests, %u collisions in %fms (%fms per Test, %fms per Collision) into %s\n",count,collcount,colltimer/1e6,colltimer/count/1e6,colltimer/collcount/1e6, outfile.c_str());
    f.close();	
}
