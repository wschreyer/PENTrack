#include "logger.h"

#include <atomic>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <thread>
#include <mutex>

#include <vector>
#include <memory>

using namespace std;

std::unique_ptr<TLogger> CreateLogger(TConfig& config){
  string logtype = "txt";
  auto ROOTlogvar = config["GLOBAL"].find("ROOTlog"); // check for old ROOTlog variable to stay backwards compatible
  if (ROOTlogvar != config["GLOBAL"].end()){
    bool ROOTlog;
    istringstream(ROOTlogvar->second) >> ROOTlog;
    if (ROOTlog) logtype = "ROOT";
  }
  istringstream(config["GLOBAL"]["logtype"]) >> logtype;

  if (logtype == "txt"){
    return std::unique_ptr<TLogger>(new TTextLogger(config));
  }
  else if (logtype == "ROOT"){
#ifdef USEROOT
    return std::unique_ptr<TLogger>(new TROOTLogger(config));
#else
    throw runtime_error("logtype ROOT is set but PENTrack was compiled without ROOT support!");
#endif
  }
  else if (logtype == "HDF5"){
#ifdef USEHDF5
    return std::unique_ptr<THDF5Logger>(new THDF5Logger(config));
#else
    throw runtime_error("logtype HDF5 is set but PENTrack was compiled without HDF5 support!");
#endif
  }
  else
    throw runtime_error("No valid logtype defined in config");
	    
}


void TLogger::Print(const std::unique_ptr<TParticle>& p, const value_type x, const state_type &y, const state_type &spin,
		    const TGeometry &geom, const TFieldManager &field, const std::string suffix){
  bool log = false;
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
  bool log = false;
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

void TLogger::PrintTrack(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, const value_type x, const state_type& y, const state_type &spin, const solid &sld, const TFieldManager &field, const bool force){
  bool log = false;
  double interval = 0.;
  istringstream(config[p->GetName()]["tracklog"]) >> log;
  istringstream(config[p->GetName()]["trackloginterval"]) >> interval;
  if (not log or interval <= 0)
    return;

  // // modified : Sly ToDo
  if (y[8] > 0 and int(y1[8]/interval) == int(y[8]/interval) and !force){ // if this is the first point or tracklength (y[8]) did cross an integer multiple of trackloginterval
    // std::cout << "exit print track" << std::endl;
    return;
  }

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
  bool log = false;
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
  bool log = false;
  double interval = 0.;
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

  // std::lock_guard<std::mutex> lock_log(log_mutex);
  DoLog(particlename, suffix, titles, vars);
}


void TTextLogger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){

  // std::cout << "\nThread " << std::this_thread::get_id() <<  " Do loggings " << std::endl << std::flush;
  // std::ostringstream oss;
  // oss << std::this_thread::get_id();
  // printf("Thread %s Do log \n", oss.str().c_str());

  // std::lock_guard<std::mutex> lock_log(mutexes[particlename + suffix]);
  std::lock_guard<std::mutex> lock_open(logstreams_mutex);

  try {
    std::ofstream &file = logstreams[particlename + suffix];

    if (!file.is_open()){
      std::ostringstream filename;
      filename << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << particlename << suffix << ".out";
      boost::filesystem::path outfile = outpath / filename.str();
      
      std::cout << "Thread " << std::this_thread::get_id() <<  " creating " << outfile << std::endl <<  std::flush;
      
      file.open(outfile.c_str());
      if(!file.is_open())
	{
	  throw std::runtime_error("Could not open " + outfile.native());
	}
      
      file << std::setprecision(std::numeric_limits<double>::digits10) << std::flush;
      // copy(titles.begin(), titles.end(), ostream_iterator<string>(file, " "));
      for(const auto& title : titles){
	file << title << " ";
      }
      file << '\n';
    }

    for(const auto& var : vars){
      file << var << " ";
    }
    // copy(vars.begin(), vars.end(), ostream_iterator<double>(file, " "));
    // file << '\n' << std::flush;
    file << '\n';
    // printf("Thread %s End log \n", oss.str().c_str());
  }
  catch (const std::length_error& le) {
    std::cerr << "Length error: " << le.what() << '\n';
  }
    
}



TTextLogger::~TTextLogger(){
  std::ostringstream oss;
  oss << std::this_thread::get_id();
  printf("Thread %s closing files \n", oss.str().c_str());
  for (auto &s: logstreams){if (s.second.is_open()) {s.second.close();}}
}


// // Original function : 
// void TTextLogger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars){

//   // Lock the mutex before accessing the text file
//   std::lock_guard<std::mutex> lock(text_mutex);
    
//   ofstream &file = logstreams[particlename + suffix];
//   if (!file.is_open()){
//     std::ostringstream filename;
//     filename << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << particlename << suffix << ".out";
//     boost::filesystem::path outfile = outpath / filename.str();
//     //		std::cout << "Creating " << outfile << '\n';
//     file.open(outfile.c_str());
//     if(!file.is_open())
//       {
// 	throw std::runtime_error("Could not open " + outfile.native());
//       }

//     file << std::setprecision(std::numeric_limits<double>::digits10);
//     copy(titles.begin(), titles.end(), ostream_iterator<string>(file, " "));
//     file << '\n';
//   }

//   copy(vars.begin(), vars.end(), ostream_iterator<double>(file, " "));
//   file << '\n';
// }


#ifdef USEROOT

#include "TObjString.h"

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




void TROOTLogger::DoLog(const std::string& particlename, const std::string& suffix, const std::vector<std::string>& titles, const std::vector<double>& vars) {
  const std::string name = particlename + suffix;
  
  // Check if tree already exists
  TNtupleD* tree = static_cast<TNtupleD*>(ROOTfile->Get(name.c_str()));
  
  if (!tree) {
    // Create new tree if it doesn't exist
    std::lock_guard<std::mutex> lock(root_mutex);
    tree = new TNtupleD(name.c_str(), name.c_str(), GetVarList(titles).c_str());
    ROOTfile->Add(tree);
    trees.push_back(std::unique_ptr<TNtupleD>(tree));
  }
  
  // Fill the tree with data
  std::lock_guard<std::mutex> lock(tree_mutex);
  tree->Fill(&vars[0]);
}

TROOTLogger::~TROOTLogger() {
  std::cout << "\nWriting root files & deleting them" <<std::flush;
  ROOTfile->Write();

}

#endif


#ifdef USEHDF5

THDF5Logger::THDF5Logger(TConfig& aconfig){
  config = aconfig;
  ostringstream filename;
  filename << setw(12) << std::setfill('0') << jobnumber << ".h5";
  boost::filesystem::path outfile = outpath / filename.str();
  HDF5file = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (HDF5file < 0) throw std::runtime_error("Could not create output file " + filename.str());

  std::string configGroupName = "config";
#if H5Gcreate_vers == 1
  auto ret = H5Gcreate1(HDF5file, configGroupName.c_str());
  if (ret == H5I_INVALID_HID) throw std::runtime_error("Could not create group in output file " + filename.str());
#elif H5Gcreate_vers == 2
  auto ret = H5Gcreate(HDF5file, configGroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (ret < 0) throw std::runtime_error("Could not create group in output file " + filename.str());
#else
  throw std::runtime_error("HDF5 library version not supported");
#endif
    
  const char *field_names[] = {"Variable", "Value"};
  struct configVariable{
    const char *variable;
    const char *value;
  };
  size_t offsets[2] = {HOFFSET(configVariable, variable), HOFFSET(configVariable, value)};
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, H5T_VARIABLE);
  hid_t field_types[2] = {string_type, string_type};
  for (auto &section: aconfig){
    auto N = section.second.size();
    configVariable data[N];
    int i = 0;
    for (auto &var: section.second){
      data[i++] = configVariable{var.first.c_str(), var.second.c_str()};
    }
    auto ret = H5TBmake_table(section.first.c_str(), HDF5file, (configGroupName + "/" + section.first).c_str(), 2, N, sizeof(configVariable), field_names, offsets, field_types, 10, nullptr, 0, data);
    if (ret < 0) throw std::runtime_error("Could not create table in output file " + filename.str());
  }
  H5Tclose(string_type);
  isClosed = false;

}



// Working, but slow : 
void THDF5Logger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) {
  // Check if interrupt signal is received
  if (quit.load()) {
    // Perform any necessary cleanup operations
    if (!isClosed) {
      H5Fclose(HDF5file);
      isClosed = true;
    }
    // Exit the function or program gracefully
    return;
  }
  
  string name = particlename + suffix;
  size_t Nfields = titles.size();
  size_t offsets[Nfields];
  for (size_t i = 0; i < Nfields; ++i) offsets[i] = i * sizeof(double);

  // std::cout << "\nsaving HDF5 table for thread " << std::this_thread::get_id() << name << std::endl;
  // Lock the mutex before accessing the HDF5 library
  std::lock_guard<std::mutex> lock(hdf5_mutex);
  // size_t mutexIndex = std::hash<std::string>{}(particlename) % hdf5_mutex_pool.size();
  // std::lock_guard<std::mutex> lock(hdf5_mutex_pool[mutexIndex]);

  
  if (H5Lexists(HDF5file, name.c_str(), H5P_DEFAULT) <= 0){
    const char *field_names[Nfields];
    hid_t field_types[Nfields];
    for (size_t i = 0; i < Nfields; ++i){
      field_names[i] = titles[i].c_str();
      field_types[i] = H5T_NATIVE_DOUBLE;
    }

    auto ret = H5TBmake_table(name.c_str(), HDF5file, name.c_str(), Nfields, 1, Nfields*sizeof(double), field_names, offsets, field_types, 10, nullptr, 1, vars.data());
    if (ret < 0) throw std::runtime_error("Could not create table " + name);
  }
  else{
    size_t sizes[Nfields];
    for (size_t i = 0; i < Nfields; ++i) sizes[i] = sizeof(double);
    auto ret = H5TBappend_records(HDF5file, name.c_str(), 1, Nfields*sizeof(double), offsets, sizes, vars.data());
    if (ret < 0) throw std::runtime_error("Could not write data to table " + name);
  }
}



// In construciton : localBuffer for threads
// void THDF5Logger::DoLog(const std::string &particlename, const std::string &suffix, const std::vector<std::string> &titles, const std::vector<double> &vars) {
//   string name = particlename + suffix;
//   size_t Nfields = titles.size();
//   size_t offsets[Nfields];
//   for (size_t i = 0; i < Nfields; ++i) offsets[i] = i * sizeof(double);

//   // Calculate the hash of the particlename to determine the lock index
//   std::hash<std::string> hash_fn;
//   size_t lockIndex = hash_fn(name) % hdf5_mutex_pool.size();

//   // Lock the corresponding mutex
//   std::lock_guard<std::mutex> lock(hdf5_mutex_pool[lockIndex]);

//   hid_t HDF5file = HDF5file_pool[lockIndex];

  
//   // Lock the mutex before accessing the HDF5 library
//   std::lock_guard<std::mutex> lock(hdf5_mutex);
//   // size_t mutexIndex = std::hash<std::string>{}(particlename) % hdf5_mutex_pool.size();
//   // std::lock_guard<std::mutex> lock(hdf5_mutex_pool[mutexIndex]);

  
//   if (H5Lexists(HDF5file, name.c_str(), H5P_DEFAULT) <= 0){
//     const char *field_names[Nfields];
//     hid_t field_types[Nfields];
//     for (size_t i = 0; i < Nfields; ++i){
//       field_names[i] = titles[i].c_str();
//       field_types[i] = H5T_NATIVE_DOUBLE;
//     }

//     auto ret = H5TBmake_table(name.c_str(), HDF5file, name.c_str(), Nfields, 1, Nfields*sizeof(double), field_names, offsets, field_types, 10, nullptr, 1, vars.data());
//     if (ret < 0) throw std::runtime_error("Could not create table " + name);
//   }
//   else{
//     size_t sizes[Nfields];
//     for (size_t i = 0; i < Nfields; ++i) sizes[i] = sizeof(double);
//     auto ret = H5TBappend_records(HDF5file, name.c_str(), 1, Nfields*sizeof(double), offsets, sizes, vars.data());
//     if (ret < 0) throw std::runtime_error("Could not write data to table " + name);
//   }
// }


THDF5Logger::~THDF5Logger() {
  if (!isClosed) {
    H5Fclose(HDF5file);
    isClosed = true;
  }
}
#endif
