/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#include "source.h"

#include <boost/format.hpp>

#include "neutron.h"
#include "proton.h"
#include "electron.h"
#include "mercury.h"
#include "xenon.h"
#include "geometry.h"
#include "globals.h"

using namespace std;

TParticleSource::TParticleSource(std::map<std::string, std::string> &sourceconf): fActiveTime(0), polarization(0), ParticleCounter(0){
  istringstream(sourceconf["particle"]) >> fParticleName;
  istringstream(sourceconf["ActiveTime"]) >> fActiveTime;
  istringstream(sourceconf["polarization"]) >> polarization;

  double rmin, rmax;
  istringstream(sourceconf["Emin"]) >> rmin;
  istringstream(sourceconf["Emax"]) >> rmax;
  spectrum = parse_distribution(sourceconf["spectrum"], rmin, rmax);

  istringstream(sourceconf["phi_v_min"]) >> rmin;
  istringstream(sourceconf["phi_v_max"]) >> rmax;
  phi_v = parse_distribution(sourceconf["phi_v"], rmin*conv, rmax*conv);

  istringstream(sourceconf["theta_v_min"]) >> rmin;
  istringstream(sourceconf["theta_v_max"]) >> rmax;
  theta_v = parse_distribution(sourceconf["theta_v"], rmin*conv, rmax*conv);

  if (sourceconf.count("pulseWidth") && sourceconf.count("pulseGap")) // Check for optional source pulsing params
    {
      istringstream(sourceconf["pulseWidth"]) >> pulseWidth;
      istringstream(sourceconf["pulseGap"]) >> pulseGap;
    } else {
    pulseWidth = 0;
    pulseGap = 0;
  }
	
  // Generate probability time start distribution for particle source
  std::vector<double> interval, weight;

  if (fActiveTime == 0){ // If all particles start at t = 0
    interval = {0, std::nextafter(0, std::numeric_limits<double>::max())};
    weight = {1};

  } else if (pulseGap > 0 && fActiveTime > 0) { // Pulsed particle source
    cout << "Pulsed particle source enabled\npulseWidth: " << pulseWidth << " pulseGap: " << pulseGap << "\n";
    for (double i = 0; i < fActiveTime; i += pulseGap) // Add intervals with source on and source off
      {
	interval.push_back(i);
	if (pulseWidth == 0) { // If user wants delta function pulses
	  interval.push_back( std::nextafter(i, std::numeric_limits<double>::max()) );
	  weight.push_back(1);
	} else if ( (i + pulseWidth) > fActiveTime){ // If last pulse extends over active time specification, shorten the pulse width
	  interval.push_back(fActiveTime); 
	  weight.push_back( (fActiveTime - i)/pulseWidth ); // adjust probability weighting accordingly
	} else {
	  interval.push_back(i + pulseWidth);
	  weight.push_back(1);
	}
	weight.push_back(0);
      }
    weight.pop_back();
  } else {	// Continuous source for duration of active time
    interval = {0, fActiveTime};
    weight = {1};
  }

  timedist = std::piecewise_constant_distribution<double> (interval.begin(), interval.end(), weight.begin());
	
  double mean = 2.97; // Utkarsh
  double sigma =  0.8113; // Utkarsh
  normtimedist = std::normal_distribution<double> (mean, sigma); // Utkarsh
}


TParticle* TParticleSource::CreateParticle(const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection, TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
  TParticle *p;
  if (fParticleName == NAME_NEUTRON)
    p = new TNeutron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, spinprojection, mc, geometry, field);
  else if (fParticleName == NAME_PROTON)
    p = new TProton(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, spinprojection, mc, geometry, field);
  else if (fParticleName == NAME_ELECTRON)
    p = new TElectron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, spinprojection, mc, geometry, field);
  else if (fParticleName == NAME_MERCURY)
    p = new TMercury(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, spinprojection, mc, geometry, field);
  else if (fParticleName == NAME_XENON) 
    p = new TXenon(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, spinprojection, mc, geometry, field);
  else{
    throw std::runtime_error("Could not create particle " + fParticleName);
  }
  return p;
}


TParticle* TSurfaceSource::CreateParticle(TMCGenerator &mc, TGeometry &geometry, const TFieldManager &field){
  CPoint p;
  CVector nv;
  unsigned ID;
  do{
    geometry.mesh.RandomPointOnSurface(p, nv, ID, mc, GetSourceVolumeBoundingBox());
  } while(!InSourceVolume(p[0], p[1], p[2]));
  p = p + nv*REFLECT_TOLERANCE; // move point slightly away from surface

  double Ekin = spectrum(mc);
  std::uniform_real_distribution<double> unidist(0., 2.*pi);
  double phi = unidist(mc);
  std::sincos_distribution<double> sincosdist(0., 0.5*pi);
  double theta = sincosdist(mc);
  if (Enormal > 0){
    double vnormal = sqrt(Ekin*cos(theta)*cos(theta) + Enormal); // add E_normal to component normal to surface
    double vtangential = sqrt(Ekin)*sin(theta);
    theta = atan2(vtangential, vnormal); // update angle
    Ekin = vnormal*vnormal + vtangential*vtangential; // update energy
  }

  double v[3] = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};
  double n[3] = {nv[0], nv[1], nv[2]};
  RotateVector(v, n);
  phi = atan2(v[1],v[0]);
  theta = acos(v[2]);

  std::polarization_distribution<double> pdist(polarization);
  int pol = pdist(mc); // randomly choose polarisation used for trajectory tracking

  return TParticleSource::CreateParticle(timedist(mc), p[0], p[1], p[2], Ekin, phi, theta, pol, polarization, mc, geometry, field);
}


void TVolumeSource::FindPotentialMinimum(TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
  cout << "Sampling phase space ";
  const int N = 100000;
  progress_display progress(N);
  std::polarization_distribution<double> pdist(polarization);
  for (int i = 0; i < N; i++){
    ++progress;
    double t = timedist(mc);
    double x, y, z;
    RandomPointInSourceVolume(x, y, z, mc); // dice point in source volume
    int pol = pdist(mc); // randomly choose polarisation used for trajectory tracking
    TParticle *p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, pol, polarization, mc, geometry, field); // create dummy particle with Ekin = 0
    double V = p->GetInitialTotalEnergy(geometry, field); // potential at particle position equals its total energy
    if (V < MinPot)
      MinPot = V; // remember minimal potential
    delete p;
    ParticleCounter--;
  }
  cout << " minimal potential = " << MinPot << "eV\n";
}

TParticle* TVolumeSource::CreateParticle(TMCGenerator &mc, TGeometry &geometry, const TFieldManager &field){
  std::polarization_distribution<double> pdist(polarization);
  if (fPhaseSpaceWeighting){ // if particle density should be weighted by available phase space
    if (MinPot == numeric_limits<double>::infinity()){ // if minimum potential energy has not yet been determined
      FindPotentialMinimum(mc, geometry, field); // find minimum potential energy
      if (MinPot > spectrum.max()) // abort program if spectrum completely out of potential range
	throw std::runtime_error( (boost::format("Error: your chosen spectrum is below the minimal potential energy in the source volume (%1% eV < %2% eV). Exiting!\n") % spectrum.max() % MinPot).str() );
    }

    if (MinPot > spectrum.min()){ // give warning if chosen spectrum contains energy ranges that are not possible
      cout << "Warning: your chosen spectrum contains energies below the minimal potential energy in the source volume (" << spectrum.min() << "eV < " << MinPot << "eV). The energy spectrum will be cut off!\n";
    }
    double H;
    do{
      H = spectrum(mc); // dice total(!) energy until one above minimum potential is found
    }while (H < MinPot);

    //		cout << "Trying to find starting point for particle with total energy " << H << "eV ...";
    for (int i = 0; true; i++){
      double t = timedist(mc); // dice start time
      double x, y, z;
      RandomPointInSourceVolume(x, y, z, mc); // dice point in source volume
      int pol = pdist(mc);
      TParticle *proposed_p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, pol, polarization, mc, geometry, field); // create new particle with Ekin = 0
      double V = proposed_p->GetInitialTotalEnergy(geometry, field); // potential at particle position equals its total energy
      delete proposed_p;
      ParticleCounter--; // delete particle and decrement particle counter

      if (H < V)
	continue; // if total energy < potential energy then particle is not possible at this point
      std::uniform_real_distribution<double> Hdist(0, sqrt(H - MinPot));
      if (sqrt(H - V) > Hdist(mc)){ // accept particle with probability sqrt(H-V)/sqrt(H-Vmin) (phase space weighting according to Golub)
	//				cout << " found after " << i+1 << " tries\n";
	return TParticleSource::CreateParticle(t, x, y, z, H - V, phi_v(mc), theta_v(mc), pol, polarization, mc, geometry, field); // if accepted, return new particle with correct Ekin
      }
    }
    assert(false); // this will never be reached
  }
  else{ // create particles uniformly distributed in volume
    double t;
    
    // do{ // Utkarsh
    //   t = std::exp(normtimedist(mc));
    // }while(!(0<t && t<fActiveTime));
 		
    t = timedist(mc); // original

    double x, y, z;
    RandomPointInSourceVolume(x, y, z, mc);
    int pol = pdist(mc);
    return TParticleSource::CreateParticle(t, x, y, z, spectrum(mc), phi_v(mc), theta_v(mc), pol, polarization, mc, geometry, field);
  }
}

TParticleSource* CreateParticleSource(TConfig &config, const TGeometry &geometry){
  std::map<std::string, std::string> &sc = config["SOURCE"];
  std::string sourcemode;
  std::istringstream(sc["sourcemode"]) >> sourcemode;

  TParticleSource *source = nullptr;
  if (sourcemode == "boxvolume"){
    source = new TCuboidVolumeSource(sc);
  }
  else if (sourcemode == "cylvolume"){
    source = new TCylindricalVolumeSource(sc);
  }
  else if (sourcemode == "cylvolumeext"){ // Utkarsh
    source = new TCylindricalVolumeSourceExt(sc);
  }
  else if (sourcemode == "STLvolume"){
    source = new TSTLVolumeSource(sc);
  }
  else if (sourcemode == "cylsurface"){
    source = new TCylindricalSurfaceSource(sc);
  }
  else if (sourcemode == "STLsurface"){
    source = new TSTLSurfaceSource(sc);
  }
  else
    throw std::runtime_error((boost::format("Could not load source %1%!") % sourcemode).str());

  return source;
}
