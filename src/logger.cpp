#include "logger.h"

#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

std::unique_ptr<TLogger> CreateLogger(TConfig& config){
    bool ROOTlog = false;
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

    map<string, double> variables = {{"jobnumber", static_cast<double>(jobnumber)},
                                     {"particle", static_cast<double>(p->GetParticleNumber())},
                                     {"m", p->GetMass()},
                                     {"q", p->GetCharge()},
                                     {"mu", p->GetMagneticMoment()},
                                     {"tstart", tstart},
                                     {"xstart", ystart[0]},
                                     {"ystart", ystart[1]},
                                     {"zstart", ystart[2]},
                                     {"vxstart", ystart[3]},
                                     {"vystart", ystart[4]},
                                     {"vzstart", ystart[5]},
                                     {"polstart", ystart[7]},
                                     {"Sxstart", spinstart[0]},
                                     {"Systart", spinstart[1]},
                                     {"Szstart", spinstart[2]},
                                     {"Hstart", p->GetInitialTotalEnergy(geom, field)},
                                     {"Estart", p->GetInitialKineticEnergy()},
                                     {"Bstart", sqrt(Bstart[0]*Bstart[0] + Bstart[1]*Bstart[1] + Bstart[2]*Bstart[2])},
                                     {"Ustart", Vstart},
                                     {"solidstart", static_cast<double>(p->GetInitialSolid().ID)},
                                     {"tend", x},
                                     {"xend", y[0]},
                                     {"yend", y[1]},
                                     {"zend", y[2]},
                                     {"vxend", y[3]},
                                     {"vyend", y[4]},
                                     {"vzend", y[5]},
                                     {"polend", y[7]},
                                     {"Sxend", spin[0]},
                                     {"Syend", spin[1]},
                                     {"Szend", spin[2]},
                                     {"Hend", H},
                                     {"Eend", E},
                                     {"Bend", sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])},
                                     {"Uend", V},
                                     {"solidend", static_cast<double>(sld.ID)},
                                     {"stopID", static_cast<double>(p->GetStopID())},
                                     {"Nspinflip", static_cast<double>(p->GetNumberOfSpinflips())},
                                     {"spinflipprob", 1 - p->GetNoSpinFlipProbability()},
                                     {"Nhit", static_cast<double>(p->GetNumberOfHits())},
                                     {"Nstep", static_cast<double>(p->GetNumberOfSteps())},
                                     {"propert", y[6]},
                                     {"trajlength", y[8]},
                                     {"Hmax", p->GetMaxTotalEnergy()},
                                     {"wL", wL}
                                    };

    vector<string> default_titles = {"jobnumber", "particle",
                                     "tstart", "xstart", "ystart", "zstart", "vxstart", "vystart", "vzstart", "polstart",
                                     "Sxstart", "Systart", "Szstart", "Hstart", "Estart", "Bstart", "Ustart", "solidstart",
                                     "tend", "xend", "yend", "zend", "vxend", "vyend", "vzend", "polend",
                                     "Sxend", "Syend", "Szend", "Hend", "Eend", "Bend", "Uend", 
                                     "solidend", "stopID", "Nspinflip", "spinflipprob", "Nhit", "Nstep", "propert", "trajlength", "Hmax", "wL"};

    Log(p->GetName(), suffix, variables, default_titles);
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

    double B[3] = {0,0,0};
    double dBidxj[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double E[3] = {0,0,0};
    double V = 0;
    field.BField(y[0],y[1],y[2],x,B, dBidxj);
    field.EField(y[0],y[1],y[2],x,V,E);
    value_type Ek = p->GetKineticEnergy(&y[3]);
    value_type H = Ek + p->GetPotentialEnergy(x, y, field, sld);

    map<string, double> variables = {{"jobnumber", static_cast<double>(jobnumber)},
                                     {"particle", static_cast<double>(p->GetParticleNumber())},
                                     {"polarisation", y[7]},
                                     {"t", x},
                                     {"x", y[0]},
                                     {"y", y[1]},
                                     {"z", y[2]},
                                     {"vx", y[3]},
                                     {"vy", y[4]},
                                     {"vz", y[5]},
                                     {"H", H},
                                     {"E", Ek},
                                     {"Bx", B[0]},
                                     {"dBxdx", dBidxj[0][0]},
                                     {"dBxdy", dBidxj[0][1]},
                                     {"dBxdz", dBidxj[0][2]},
                                     {"By", B[1]},
                                     {"dBydx", dBidxj[1][0]},
                                     {"dBydy", dBidxj[1][1]},
                                     {"dBydz", dBidxj[1][2]},
                                     {"Bz", B[2]},
                                     {"dBzdx", dBidxj[2][0]},
                                     {"dBzdy", dBidxj[2][1]},
                                     {"dBzdz", dBidxj[2][2]},
                                     {"Ex", E[0]},
                                     {"Ey", E[1]},
                                     {"Ez", E[2]},
                                     {"V", V}
                                    };

    vector<string> default_titles = {"jobnumber", "particle",
                                     "polarisation", "t", "x", "y", "z", "vx", "vy", "vz", "H", "E",
                                     "Bx", "dBxdx", "dBxdy", "dBxdz", "By", "dBydx", "dBydy", "dBydz", "Bz", "dBzdx", "dBzdy", "dBzdz",
                                     "Ex", "Ey", "Ez", "V"};

    Log(p->GetName(), "track", variables, default_titles);
}

void TLogger::PrintHit(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y1, const state_type &y2, const double *normal, const solid &leaving, const solid &entering){
    bool log;
    istringstream(config[p->GetName()]["hitlog"]) >> log;
    if (not log)
        return;

    map<string, double> variables = {{"jobnumber", static_cast<double>(jobnumber)},
                                     {"particle", static_cast<double>(p->GetParticleNumber())},
                                     {"t", x},
                                     {"x", y1[0]},
                                     {"y", y1[1]},
                                     {"z", y1[2]},
                                     {"v1x", y1[3]},
                                     {"v1y", y1[4]},
                                     {"v1z", y1[5]},
                                     {"pol1", y1[7]},
                                     {"v2x", y2[3]},
                                     {"v2y", y2[4]},
                                     {"v2z", y2[5]},
                                     {"pol2", y2[7]},
                                     {"nx", normal[0]},
                                     {"ny", normal[1]},
                                     {"nz", normal[2]},
                                     {"solid1", static_cast<double>(leaving.ID)},
                                     {"solid2", static_cast<double>(entering.ID)}
                                    };

    vector<string> default_titles = {"jobnumber", "particle",
                                     "t", "x", "y", "z",
                                     "v1x", "v1y", "v1z", "pol1", "v2x", "v2y", "v2z", "pol2",
                                     "nx", "ny", "nz", "solid1", "solid2"};

    Log(p->GetName(), "hit", variables, default_titles);
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

    double B[3] = {0,0,0};
    state_type y(STATE_VARIABLES);
    trajectory_stepper.calc_state(x, y);
    field.BField(y[0],y[1],y[2],x,B);
    double Omega[3];
    p->SpinPrecessionAxis(x, trajectory_stepper, field, Omega[0], Omega[1], Omega[2]);

    state_type spin(spinstepper.current_state());
    if (x < spinstepper.current_time()){
        spinstepper.calc_state(x, spin);
    }

    map<string, double> variables = {{"jobnumber", static_cast<double>(jobnumber)},
                                     {"particle", static_cast<double>(p->GetParticleNumber())},
                                     {"t", x},
                                     {"x", y[0]},
                                     {"y", y[1]},
                                     {"z", y[2]},
                                     {"Sx", spin[0]},
                                     {"Sy", spin[1]},
                                     {"Sz", spin[2]},
                                     {"Wx", Omega[0]},
                                     {"Wy", Omega[1]},
                                     {"Wz", Omega[2]},
                                     {"Bx", B[0]},
                                     {"By", B[1]},
                                     {"Bz", B[2]}
                                    };

    vector<string> default_titles = {"jobnumber", "particle",
                                     "t", "x", "y", "z",
                                     "Sx", "Sy", "Sz", "Wx", "Wy", "Wz", "Bx", "By", "Bz"};

    Log(p->GetName(), "spin", variables, default_titles);
}

void TLogger::Log(const std::string &particlename, const std::string &suffix, const std::map<std::string, double> &variables, const std::vector<std::string> &default_titles){
    vector<string> titles;
    vector<double> vars;
    string filter;
    istringstream(config[particlename][suffix + "logfilter"]) >> filter;
    if (filter != "" and not EvalFormula(config, filter, variables)){
        return;
    }
    if (config[particlename][suffix + "logvars"] == ""){
        cout << suffix << "log for " << particlename << " is enabled but " << suffix << "logvars is empty. I will default to backward compatible output.\nSee example config on how to use the new logvars and logfilter options.\n";
        ostringstream os;
        copy(default_titles.begin(), default_titles.end(), ostream_iterator<string>(os, " "));
        config[particlename][suffix + "logvars"] = os.str();
    }
    istringstream varstr(config[particlename][suffix + "logvars"]);
    for (istream_iterator<string> var(varstr); var != istream_iterator<string>(); ++var){
        titles.push_back(*var);
        auto val = variables.find(*var);
        if (val != variables.end())
            vars.push_back(val->second);
        else
            vars.push_back(EvalFormula(config, *var, variables));
    }

    DoLog(particlename, suffix, titles, vars);
}


void TTextLogger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){
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

    TDirectory *rootdir = ROOTfile->mkdir("config", "config");
    for (auto &section: aconfig){
        TDirectory *dir = rootdir->mkdir(section.first.c_str(), section.first.c_str());
        dir->cd();
        for (auto &var: section.second){
            TObjString s(var.second.c_str());
            dir->WriteTObject(&s, var.first.c_str());
        }
    }
    ROOTfile->cd();
}

void TROOTLogger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){
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
