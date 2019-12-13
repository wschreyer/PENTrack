#include "logger.h"

#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

std::unique_ptr<TLogger> CreateLogger(TConfig& config){
    bool ROOTlog;
    istringstream(config["GLOBAL"]["ROOTlog"]) >> ROOTlog;
    if (ROOTlog){
        #ifdef USEROOT
            return std::unique_ptr<TLogger>(new TROOTLogger(config));
        #else
            throw runtime_error("ROOTlog is set but PENTrack was compiled without ROOT support!");
        #endif
    }
    else
	    return std::unique_ptr<TLogger>(new TTextLogger(config));
}


void TLogger::Print(const std::unique_ptr<TParticle>& p, const double &t, const TStep &stepper,
        const TGeometry &geom, const TFieldManager &field, const std::string suffix){
    bool log;
    istringstream(config[p->GetName()][suffix + "log"]) >> log;
    if (not log)
        return;

    vector<string> titles = {"jobnumber", "particle",
                             "tstart", "xstart", "ystart", "zstart",
                             "vxstart", "vystart", "vzstart", "polstart",
                             "Sxstart", "Systart", "Szstart",
                             "Hstart", "Estart", "Bstart",
                             "Ustart", "solidstart",
                             "tend", "xend", "yend", "zend",
                             "vxend", "vyend", "vzend", "polend",
                             "Sxend", "Syend", "Szend",
                             "Hend", "Eend", "Bend", "Uend", "solidend",
                             "stopID", "Nspinflip", "spinflipprob",
                             "Nhit", "Nstep", "trajlength", "Hmax", "wL"};

    auto pos = stepper.GetPosition(t);
    auto v = stepper.GetVelocity(t);
    double E = p->GetKineticEnergy(v);
    double Bstart[3], Eistart[3], Vstart;

    double tstart = p->GetInitialTime();
    auto posstart = p->GetInitialPosition();
    auto vstart = p->GetInitialVelocity();
    auto spinstart = p->GetInitialSpin();
    field.BField(posstart[0], posstart[1], posstart[2], tstart, Bstart);
    field.EField(posstart[0], posstart[1], posstart[2], tstart, Vstart, Eistart);

    double H;
    solid sld = geom.GetSolid(t, &pos[0]);
    H = E + p->GetPotentialEnergy(t, stepper.GetPosition(t), stepper.GetVelocity(t), stepper.GetPolarization(t), field, sld);

    double B[3], Ei[3], V;
    field.BField(pos[0], pos[1], pos[2], t, B);
    field.EField(pos[0], pos[1], pos[2], t, V, Ei);

    auto spin = stepper.GetSpin();
    double wL = stepper.GetSpinPhase();

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           tstart, posstart[0], posstart[1], posstart[2],
                           vstart[0], vstart[1], vstart[2], static_cast<double>(p->GetInitialPolarization()),
                           spinstart[0], spinstart[1], spinstart[2],
                           p->GetInitialTotalEnergy(geom, field), p->GetInitialKineticEnergy(),
                           sqrt(Bstart[0]*Bstart[0] + Bstart[1]*Bstart[1] + Bstart[2]*Bstart[2]), Vstart, static_cast<double>(p->GetInitialSolid().ID),
                           t, pos[0], pos[1], pos[2],
                           v[0], v[1], v[2], static_cast<double>(stepper.GetPolarization(t)),
                           spin[0], spin[1], spin[2], H, E,
                           sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]), V, static_cast<double>(sld.ID),
                           static_cast<double>(p->GetStopID()), static_cast<double>(p->GetNumberOfSpinflips()), 1 - p->GetNoSpinFlipProbability(),
                           static_cast<double>(p->GetNumberOfHits()), static_cast<double>(p->GetNumberOfSteps()), stepper.GetPathLength(), p->GetMaxTotalEnergy(), wL};

    Log(p->GetName(), suffix, titles, vars);
}

void TLogger::PrintSnapshot(const std::unique_ptr<TParticle>& p, const TStep &stepper,
                            const TGeometry &geom, const TFieldManager &field){
    bool log;
    istringstream(config[p->GetName()]["snapshotlog"]) >> log;
    if (not log)
        return;
    istringstream snapshottimes(config[p->GetName()]["snapshots"]);
    auto tsnap = find_if(istream_iterator<double>(snapshottimes), istream_iterator<double>(), [&](const double& tsnapshot){ return stepper.GetStartTime() <= tsnapshot and tsnapshot < stepper.GetTime(); });
    if (tsnap != istream_iterator<double>()){
        Print(p, *tsnap, stepper, geom, field, "snapshot");
    }
}

void TLogger::PrintTrack(const std::unique_ptr<TParticle>& p, const TStep &stepper,
                        const solid &sld, const TFieldManager &field){
    bool log;
    double interval;
    istringstream(config[p->GetName()]["tracklog"]) >> log;
    istringstream(config[p->GetName()]["trackloginterval"]) >> interval;
    if (not log or interval <= 0)
        return;

    
    if (stepper.GetPathLength() > 0 and int(stepper.GetPathLength(stepper.GetStartTime())/interval) == int(stepper.GetPathLength()/interval)) // if this is the first point or tracklength did cross an integer multiple of trackloginterval
        return;

    vector<string> titles = {"jobnumber", "particle", "polarisation",
                             "t", "x", "y", "z", "vx", "vy", "vz",
                             "H", "E", "Bx", "dBxdx", "dBxdy", "dBxdz", "By", "dBydx",
                             "dBydy", "dBydz", "Bz", "dBzdx", "dBzdy", "dBzdz", "Ex", "Ey", "Ez", "V"};

    double t = stepper.GetTime();
    auto pos = stepper.GetPosition();
    auto v = stepper.GetVelocity();
    double B[3] = {0,0,0};
    double dBidxj[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double E[3] = {0,0,0};
    double V = 0;
    field.BField(pos[0], pos[1], pos[2], t, B, dBidxj);
    field.EField(pos[0], pos[1], pos[2], t, V, E);
    double Ek = p->GetKineticEnergy(stepper.GetVelocity());
    double H = Ek + p->GetPotentialEnergy(stepper.GetTime(), stepper.GetPosition(), stepper.GetVelocity(), static_cast<double>(stepper.GetPolarization()), field, sld);

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()), static_cast<double>(stepper.GetPolarization()),
                           t, pos[0], pos[1], pos[2], v[0], v[1], v[2],
                           H, Ek};
    for (int i = 0; i < 3; i++){
        vars.push_back(B[i]);
        for (int j = 0; j < 3; j++)
            vars.push_back(dBidxj[i][j]);
    }
    for (int i = 0; i < 3; ++i)
        vars.push_back(E[i]);
    vars.push_back(V);

    Log(p->GetName(), "track", titles, vars);
}

void TLogger::PrintHit(const std::unique_ptr<TParticle>& p, const TStep &stepper, const double *normal, const solid &leaving, const solid &entering){
    bool log;
    istringstream(config[p->GetName()]["hitlog"]) >> log;
    if (not log)
        return;

    vector<string> titles = {"jobnumber", "particle",
                             "t", "x", "y", "z", "v1x", "v1y", "v1z", "pol1",
                             "v2x", "v2y", "v2z", "pol2",
                             "nx", "ny", "nz", "solid1", "solid2"};

    double t = stepper.GetStartTime();
    auto pos1 = stepper.GetPosition(t);
    auto v1 = stepper.GetVelocity(t);
    auto v2 = stepper.GetVelocity();
    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           t, pos1[0], pos1[1], pos1[2], v1[0], v1[1], v1[2], static_cast<double>(stepper.GetPolarization(t)),
                           v2[0], v2[1], v2[2], static_cast<double>(stepper.GetPolarization()),
                           normal[0], normal[1], normal[2], static_cast<double>(leaving.ID), static_cast<double>(entering.ID)};

    Log(p->GetName(), "hit", titles, vars);
}

void TLogger::PrintSpin(const std::unique_ptr<TParticle>& p, const double& t, const std::array<double, 3>& spin,
               const TStep &trajectory_stepper, const TFieldManager &field) {
    bool log;
    istringstream(config[p->GetName()]["spinlog"]) >> log;
    if (not log)
        return;

    double interval;
    istringstream(config[p->GetName()]["spinloginterval"]) >> interval;
    if (t > lastspinlog and t < lastspinlog + interval)
        return;
    lastspinlog = t;

    vector<string> titles = {"jobnumber", "particle",
                             "t", "x", "y", "z",
                             "Sx", "Sy", "Sz",
                             "Wx", "Wy", "Wz",
                             "Bx", "By", "Bz"};

    double B[3] = {0,0,0};
    auto pos = trajectory_stepper.GetPosition(t);
    field.BField(pos[0], pos[1], pos[2], t, B);
    auto Omega = p->SpinPrecessionAxis(t, trajectory_stepper, field);

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           t, pos[0], pos[1], pos[2],
                           spin[0], spin[1], spin[2],
                           Omega[0], Omega[1], Omega[2],
                           B[0], B[1], B[2]};

    Log(p->GetName(), "spin", titles, vars);
}




void TTextLogger::Log(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){
    ofstream &file = logstreams[particlename + suffix];
    if (!file.is_open()){
        std::ostringstream filename;
        filename << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << particlename << suffix << ".out";
        boost::filesystem::path outfile = outpath / filename.str();
//		std::cout << "Creating " << outfile << '\n';
        file.open(outfile.c_str());
        if(!file.is_open())
        {
            throw std::runtime_error("Could not open " + outfile.native());
        }

        file << std::setprecision(std::numeric_limits<double>::digits10);
        copy(titles.begin(), titles.end(), ostream_iterator<string>(file, " "));
        file << '\n';
    }

    copy(vars.begin(), vars.end(), ostream_iterator<double>(file, " "));
    file << '\n';
}


TROOTLogger::TROOTLogger(TConfig& aconfig){
    config = aconfig;
    lastspinlog = numeric_limits<double>::lowest();
    ostringstream filename;
    filename << setw(12) << std::setfill('0') << jobnumber << ".root";
    boost::filesystem::path outfile = outpath / filename.str();
    ROOTfile = new TFile(outfile.c_str(), "RECREATE");
    if (not ROOTfile->IsOpen())
        throw std::runtime_error("Could not open " + outfile.native());
}

void TROOTLogger::Log(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){
    string name = particlename + suffix;
    TNtupleD* tree = static_cast<TNtupleD*>(ROOTfile->Get(name.c_str()));
    if (not tree){
        ostringstream varliststr;
        copy(titles.begin(), titles.end(), ostream_iterator<string>(varliststr, ":"));
        string varlist = varliststr.str();
        varlist.pop_back();
        tree = new TNtupleD(name.c_str(), name.c_str(), varlist.c_str());
    }
    tree->Fill(&vars[0]);
}

TROOTLogger::~TROOTLogger(){
    ROOTfile->Write();
    delete ROOTfile;
}