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


void TLogger::Print(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y, const state_type &spin,
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

    value_type E = p->GetKineticEnergy(&y[3]);
    double Bstart[3], Eistart[3], Vstart;

    value_type tstart = p->GetInitialTime();
    state_type ystart = p->GetInitialState();
    state_type spinstart = p->GetInitialSpin();
    field.BField(ystart[0], ystart[1], ystart[2], tstart, Bstart);
    field.EField(ystart[0], ystart[1], ystart[2], tstart, Vstart, Eistart);

    double H;
    solid sld = geom.GetSolid(x, &y[0]);
    H = E + p->GetPotentialEnergy(x, y, field, sld);

    double B[3], Ei[3], V;
    field.BField(y[0], y[1], y[2], x, B);
    field.EField(y[0], y[1], y[2], x, V, Ei);

    double wL = 0;
    if (spin[3] > 0)
        wL = spin[4]/spin[3];

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           tstart, ystart[0], ystart[1], ystart[2],
                           ystart[3], ystart[4], ystart[5], ystart[7],
                           spinstart[0], spinstart[1], spinstart[2],
                           p->GetInitialTotalEnergy(geom, field), p->GetInitialKineticEnergy(),
                           sqrt(Bstart[0]*Bstart[0] + Bstart[1]*Bstart[1] + Bstart[2]*Bstart[2]), Vstart, static_cast<double>(p->GetInitialSolid().ID),
                           x, y[0], y[1], y[2],
                           y[3], y[4], y[5], y[7],
                           spin[0], spin[1], spin[2], H, E,
                           sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]), V, static_cast<double>(sld.ID),
                           static_cast<double>(p->GetStopID()), static_cast<double>(p->GetNumberOfSpinflips()), 1 - p->GetNoSpinFlipProbability(),
                           static_cast<double>(p->GetNumberOfHits()), static_cast<double>(p->GetNumberOfSteps()), y[8], p->GetMaxTotalEnergy(), wL};

    Log(p->GetName(), suffix, titles, vars);
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

    if (y[8] > 0 and int(y1[8]/interval) == int(y[8]/interval)) // if this is the first point or tracklength did cross an integer multiple of trackloginterval
        return;

    vector<string> titles = {"jobnumber", "particle", "polarisation",
                             "t", "x", "y", "z", "vx", "vy", "vz",
                             "H", "E", "Bx", "dBxdx", "dBxdy", "dBxdz", "By", "dBydx",
                             "dBydy", "dBydz", "Bz", "dBzdx", "dBzdy", "dBzdz", "Ex", "Ey", "Ez", "V"};

    double B[3] = {0,0,0};
    double dBidxj[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double E[3] = {0,0,0};
    double V = 0;
    field.BField(y[0],y[1],y[2],x,B, dBidxj);
    field.EField(y[0],y[1],y[2],x,V,E);
    value_type Ek = p->GetKineticEnergy(&y[3]);
    value_type H = Ek + p->GetPotentialEnergy(x, y, field, sld);

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()), y[7],
                           x, y[0], y[1], y[2], y[3], y[4], y[5],
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

void TLogger::PrintHit(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering){
    bool log;
    istringstream(config[p->GetName()]["hitlog"]) >> log;
    if (not log)
        return;

    vector<string> titles = {"jobnumber", "particle",
                             "t", "x", "y", "z", "v1x", "v1y", "v1z", "pol1",
                             "v2x", "v2y", "v2z", "pol2",
                             "nx", "ny", "nz", "solid1", "solid2"};

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           x, y1[0], y1[1], y1[2], y1[3], y1[4], y1[5], y1[7],
                           y2[3], y2[4], y2[5], y2[7],
                           normal[0], normal[1], normal[2], static_cast<double>(leaving.ID), static_cast<double>(entering.ID)};

    Log(p->GetName(), "hit", titles, vars);
}

void TLogger::PrintSpin(const std::unique_ptr<TParticle>& p, const value_type x, const dense_stepper_type& spinstepper,
               const dense_stepper_type &trajectory_stepper, const TFieldManager &field) {
    bool log;
    double interval;
    istringstream(config[p->GetName()]["spinlog"]) >> log;
    istringstream(config[p->GetName()]["spinloginterval"]) >> interval;
    if (not log or interval <= 0)
        return;

    double x1 = spinstepper.previous_time();
    if (x > x1 and int(x1 / interval) == int(x / interval)) // if time crossed an integer multiple of spinloginterval
        return;

    vector<string> titles = {"jobnumber", "particle",
                             "t", "x", "y", "z",
                             "Sx", "Sy", "Sz",
                             "Wx", "Wy", "Wz",
                             "Bx", "By", "Bz"};

    double B[3] = {0,0,0};
    state_type y(STATE_VARIABLES);
    trajectory_stepper.calc_state(x, y);
    field.BField(y[0],y[1],y[2],x,B);
    double Omega[3];
    p->SpinPrecessionAxis(x, trajectory_stepper, field, Omega[0], Omega[1], Omega[2]);

    state_type spin(SPIN_STATE_VARIABLES);
    spinstepper.calc_state(x, spin);

    vector<double> vars = {static_cast<double>(jobnumber), static_cast<double>(p->GetParticleNumber()),
                           x, y[0], y[1], y[2],
                           spin[0], spin[1], spin[2],
                           Omega[0], Omega[1], Omega[2],
                           B[0], B[1], B[2]};

    Log(p->GetName(), "spin", titles, vars);
}




void TTextLogger::Log(std::string particlename, std::string suffix, std::vector<std::string> titles, std::vector<double> vars){
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
    ostringstream filename;
    filename << setw(12) << std::setfill('0') << jobnumber << ".root";
    boost::filesystem::path outfile = outpath / filename.str();
    ROOTfile = new TFile(outfile.c_str(), "RECREATE");
    if (not ROOTfile->IsOpen())
        throw std::runtime_error("Could not open " + outfile.native());
}

void TROOTLogger::Log(std::string particlename, std::string suffix, std::vector<std::string> titles, std::vector<double> vars){
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