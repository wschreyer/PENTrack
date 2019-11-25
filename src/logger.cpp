#include "logger.h"

#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

unique_ptr<TLogger> logger;

void CreateLogger(TConfig& config){
	logger = std::unique_ptr<TLogger>(new TLogger(config));
}

void TLogger::Print(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y, const state_type &spin, const TGeometry &geom, const TFieldManager &field, const std::string& suffix){
	bool log;
	istringstream(config[p->GetName()][suffix + "log"]) >> log;
	if (not log)
		return;

	ofstream& file = logstreams[p->GetName() + suffix];
	if (!file.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << p->GetName() << suffix << ".out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		file.open(outfile.c_str());
		if (!file.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		file <<	"jobnumber particle "
					"tstart xstart ystart zstart "
					"vxstart vystart vzstart polstart "
					"Sxstart Systart Szstart "
					"Hstart Estart Bstart Ustart solidstart "
					"tend xend yend zend "
					"vxend vyend vzend polend "
					"Sxend Syend Szend "
					"Hend Eend Bend Uend solidend "
					"stopID Nspinflip spinflipprob "
					"Nhit Nstep trajlength Hmax wL\n";
		file << std::setprecision(std::numeric_limits<double>::digits10); // need maximum precision for wL and delwL 
	}
//	cout << "Printing status\n";

	value_type E = p->GetKineticEnergy(&y[3]);
	double B[3], Ei[3], V;

	value_type tstart = p->GetInitialTime();
	state_type ystart = p->GetInitialState();
	state_type spinstart = p->GetInitialSpin();
	field.BField(ystart[0], ystart[1], ystart[2], tstart, B);
	field.EField(ystart[0], ystart[1], ystart[2], tstart, V, Ei);

	file	<< jobnumber << " " << p->GetParticleNumber() << " "
			<< tstart << " " << ystart[0] << " " << ystart[1] << " " << ystart[2] << " "
			<< ystart[3] << " " << ystart[4] << " " << ystart[5] << " " << ystart[7] << " "
			<< spinstart[0] << " " << spinstart[1] << " " << spinstart[2] << " " << p->GetInitialTotalEnergy(geom, field) << " " << p->GetInitialKineticEnergy() << " "
			<< sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) << " " << V << " " << p->GetInitialSolid().ID << " ";

	double H;
	solid sld = geom.GetSolid(x, &y[0]);
	H = E + p->GetPotentialEnergy(x, y, field, sld);

	field.BField(y[0], y[1], y[2], x, B);
	field.EField(y[0], y[1], y[2], x, V, Ei);

	double wL = 0;
	if (spin[3] > 0)
		wL = spin[4]/spin[3];

	file	<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< y[3] << " " << y[4] << " " << y[5] << " " << y[7] << " "
			<< spin[0] << " " << spin[1] << " " << spin[2] << " " << H << " " << E << " "
			<< sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) << " " << V << " " << sld.ID << " "
			<< p->GetStopID() << " " << p->GetNumberOfSpinflips() << " " << 1 - p->GetNoSpinFlipProbability() << " "
			<< p->GetNumberOfHits() << " " << p->GetNumberOfSteps() << " " << y[8] << " " << p->GetMaxTotalEnergy() << " " << wL << '\n';
}

void TLogger::PrintSnapshot(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x2, const state_type &y2,
		const state_type &spin, const dense_stepper_type& stepper, const TGeometry &geom, const TFieldManager &field){
	bool log;
	istringstream(config[p->GetName()]["snapshotlog"]) >> log;
	if (not log)
		return;
	istringstream snapshottimes(config[p->GetName()]["snapshots"]);
	auto tsnap = find_if(istream_iterator<double>(snapshottimes), istream_iterator<double>(), [&](const double& tsnapshot){ return x1 <= tsnapshot and tsnapshot < x2; });
	if (tsnap != istream_iterator<double>()){
		state_type ysnap(STATE_VARIABLES);
		stepper.calc_state(*tsnap, ysnap);
		Print(p, *tsnap, ysnap, spin, geom, field, "snapshot");
	}
}

void TLogger::PrintTrack(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x, const state_type& y,
		const state_type &spin, const solid &sld, const TFieldManager &field){
	bool log;
	double interval;
	istringstream(config[p->GetName()]["tracklog"]) >> log;
	istringstream(config[p->GetName()]["trackloginterval"]) >> interval;
	if (not log or interval <= 0)
		return;

	if (y[8] > 0 and int(y1[8]/interval) == int(y[8]/interval)) // if this is not the first point and tracklength did not cross an integer multiple of trackloginterval
		return;

	ofstream &trackfile = logstreams[p->GetName() + "track"];
	if (!trackfile.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << p->GetName() << "track.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		trackfile.open(outfile.c_str());
		if (!trackfile.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		trackfile << 	"jobnumber particle polarisation "
						"t x y z vx vy vz "
						"H E Bx dBxdx dBxdy dBxdz By dBydx "
						"dBydy dBydz Bz dBzdx dBzdy dBzdz Ex Ey Ez V\n";
		trackfile.precision(10);
	}

//	cout << "-";
	double B[3] = {0,0,0};
	double dBidxj[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double E[3] = {0,0,0};
	double V = 0;
	field.BField(y[0],y[1],y[2],x,B, dBidxj);
	field.EField(y[0],y[1],y[2],x,V,E);
	value_type Ek = p->GetKineticEnergy(&y[3]);
	value_type H = Ek + p->GetPotentialEnergy(x, y, field, sld);

	trackfile << jobnumber << " " << p->GetParticleNumber() << " " << y[7] << " "
				<< x << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " "
				<< H << " " << Ek << " ";
	for (int i = 0; i < 3; i++){
		trackfile << B[i] << " ";
		for (int j = 0; j < 3; j++)
			trackfile << dBidxj[i][j] << " ";
	}
	trackfile << E[0] << " " << E[1] << " " << E[2] << " " << V << '\n';
}


void TLogger::PrintHit(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering){
	bool log;
	istringstream(config[p->GetName()]["hitlog"]) >> log;
	if (not log)
		return;

	ofstream &hitfile = logstreams[p->GetName() + "hit"];
	if (!hitfile.is_open()){
		ostringstream filename;
		filename << setw(12) << setfill('0') << jobnumber << p->GetName() << "hit.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		cout << "Creating " << outfile << '\n';
		hitfile.open(outfile.c_str());
		if (!hitfile.is_open()){
			throw std::runtime_error("Could not create" + outfile.native());
		}
		hitfile << "jobnumber particle "
					"t x y z v1x v1y v1z pol1 "
					"v2x v2y v2z pol2 "
					"nx ny nz solid1 solid2\n";
		hitfile.precision(10);
	}

//	cout << ":";
	hitfile << jobnumber << " " << p->GetParticleNumber() << " "
			<< x << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << y1[5] << " " << y1[7] << " "
			<< y2[3] << " " << y2[4] << " " << y2[5] << " " << y2[7] << " "
			<< normal[0] << " " << normal[1] << " " << normal[2] << " " << leaving.ID << " " << entering.ID << '\n';
}


void TLogger::PrintSpin(const std::unique_ptr<TParticle>& p, const value_type x, const dense_stepper_type& spinstepper,
		const dense_stepper_type &trajectory_stepper, const TFieldManager &field){
	bool log;
	double interval;
	istringstream(config[p->GetName()]["spinlog"]) >> log;
	istringstream(config[p->GetName()]["spinloginterval"]) >> interval;
	if (not log or interval <= 0)
		return;

	double x1 = spinstepper.previous_time();
	if (x > x1 and int(x1/interval) == int(x/interval)) // if time did not cross an integer multiple of spinloginterval
		return;

	ofstream &spinfile = logstreams[p->GetName() + "spin"];
	if (!spinfile.is_open()){
		std::ostringstream filename;
		filename << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << p->GetName() << "spin.out";
		boost::filesystem::path outfile = outpath / filename.str();
//		std::cout << "Creating " << outfile << '\n';
		spinfile.open(outfile.c_str());
		if(!spinfile.is_open())
		{
			throw std::runtime_error("Could not open " + outfile.native());
		}

		//need the maximum accuracy in spinoutlog for the larmor frequency to see any difference
		spinfile << std::setprecision(std::numeric_limits<double>::digits10);
		// spinfile << "jobnumber particle t Sx Sy Sz Wx Wy Wz\n";
		spinfile << "jobnumber particle t x y z Sx Sy Sz Wx Wy Wz Bx By Bz\n";
	}

	double B[3] = {0,0,0};
	state_type y(STATE_VARIABLES);
	trajectory_stepper.calc_state(x, y);
	field.BField(y[0],y[1],y[2],x,B);
	double Omega[3];
	p->SpinPrecessionAxis(x, trajectory_stepper, field, Omega[0], Omega[1], Omega[2]);

	state_type spin(SPIN_STATE_VARIABLES);
	spinstepper.calc_state(x, spin);

	// spinfile << jobnumber << " " << particlenumber << " "
	// 		<< x << " " << spin[0] << " " << spin[1] << " " << spin[2] << " "
	// 		<< Omega[0] << " " << Omega[1] << " " << Omega[2] << "\n";
	spinfile << jobnumber << " " << p->GetParticleNumber() << " "
			<< x << " " << y[0] << " " << y[1] << " " << y[2] << " "
			<< spin[0] << " " << spin[1] << " " << spin[2] << " "
			<< Omega[0] << " " << Omega[1] << " " << Omega[2] << " "
			<< B[0] << " " << B[1] << " " << B[2] << "\n";
}
