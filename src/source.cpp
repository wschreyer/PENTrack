/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#include "source.h"
#include "neutron.h"
#include "proton.h"
#include "electron.h"
#include "mercury.h"
#include "xenon.h"
#include "geometry.h"
#include "globals.h"

using namespace std;

TParticleSource::TParticleSource(std::map<std::string, std::string> &sourceconf, TMCGenerator &mc): fActiveTime(0), polarization(0), ParticleCounter(0){
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
}


TParticle* TParticleSource::CreateParticle(double t, double x, double y, double z, double E, double phi, double theta, double polarisation,
		TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	TParticle *p;
	if (fParticleName == NAME_NEUTRON)
		p = new TNeutron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, mc, geometry, field);
	else if (fParticleName == NAME_PROTON)
		p = new TProton(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, mc, geometry, field);
	else if (fParticleName == NAME_ELECTRON)
		p = new TElectron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, mc, geometry, field);
	else if (fParticleName == NAME_MERCURY)
		p = new TMercury(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, mc, geometry, field);
	else if (fParticleName == NAME_XENON) 
		p = new TXenon(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, mc, geometry, field);
	else{
		cout << "Could not create particle " << fParticleName << '\n';
		exit(-1);
	}
	return p;
}


TParticle* TSurfaceSource::CreateParticle(TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	CPoint p;
	CVector nv;
	do{
		std::uniform_real_distribution<double> areadist(0, area_sum.back());
		double area = areadist(mc);
		int index = std::distance(area_sum.begin(), std::lower_bound(area_sum.begin(), area_sum.end(), area));
		CTriangle tri = (geometry.mesh.GetTrianglesBegin() + index)->first;
		std::uniform_real_distribution<double> unidist(0, 1);
		double a = unidist(mc); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
		double b = unidist(mc);
		if (a+b > 1){
			a = 1 - a;
			b = 1 - b;
		}
		nv = tri.supporting_plane().orthogonal_vector();
		nv = nv/sqrt(nv.squared_length());
		p = tri[0] + a*(tri[1] - tri[0]) + b*(tri[2] - tri[0]);
	}while (!InSourceVolume(p[0], p[1], p[2]));
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

	std::uniform_real_distribution<double> timedist(0., fActiveTime);
	return TParticleSource::CreateParticle(timedist(mc), p[0], p[1], p[2], Ekin, phi, theta, polarization, mc, geometry, field);
}


void TVolumeSource::FindPotentialMinimum(TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	cout << "Sampling phase space ";
	const int N = 100000;
	int percent = 0;
	for (int i = 0; i < N; i++){
		PrintPercent((double)i/N, percent);
		std::uniform_real_distribution<double> timedist(0, fActiveTime);
		double t = timedist(mc);
		double x, y, z;
		RandomPointInSourceVolume(x, y, z, mc); // dice point in source volume
		TParticle *p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarization, mc, geometry, field); // create dummy particle with Ekin = 0
		double V = p->GetInitialTotalEnergy(geometry, field); // potential at particle position equals its total energy
		if (V < MinPot)
			MinPot = V; // remember minimal potential
		delete p;
		ParticleCounter--;
	}
	cout << " minimal potential = " << MinPot << "eV\n";
}

TParticle* TVolumeSource::CreateParticle(TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	std::uniform_real_distribution<double> timedist(0, fActiveTime);
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

		cout << "Trying to find starting point for particle with total energy " << H << "eV ...";
		for (int i = 0; true; i++){
			double t = timedist(mc); // dice start time
			double x, y, z;
			RandomPointInSourceVolume(x, y, z, mc); // dice point in source volume
			TParticle *proposed_p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarization, mc, geometry, field); // create new particle with Ekin = 0
			double V = proposed_p->GetInitialTotalEnergy(geometry, field); // potential at particle position equals its total energy
			delete proposed_p;
			ParticleCounter--; // delete particle and decrement particle counter

			if (H < V)
				continue; // if total energy < potential energy then particle is not possible at this point
			std::uniform_real_distribution<double> Hdist(0, sqrt(H - MinPot));
			if (sqrt(H - V) > Hdist(mc)){ // accept particle with probability sqrt(H-V)/sqrt(H-Vmin) (phase space weighting according to Golub)
				cout << " found after " << i+1 << " tries\n";
				return TParticleSource::CreateParticle(t, x, y, z, H - V, phi_v(mc), theta_v(mc), polarization, mc, geometry, field); // if accepted, return new particle with correct Ekin
			}
		}
		assert(false); // this will never be reached
	}
	else{ // create particles uniformly distributed in volume
		double t = timedist(mc);
		double x, y, z;
		RandomPointInSourceVolume(x, y, z, mc);
		cout << '\n';
		return TParticleSource::CreateParticle(t, x, y, z, spectrum(mc), phi_v(mc), theta_v(mc), polarization, mc, geometry, field);
	}
}

TParticleSource* CreateParticleSource(TConfig &config, TMCGenerator &mc, const TGeometry &geometry){
	std::map<std::string, std::string> &sc = config["SOURCE"];
	std::string sourcemode;
	std::istringstream(sc["sourcemode"]) >> sourcemode;

	TParticleSource *source = nullptr;
	if (sourcemode == "boxvolume"){
		source = new TCuboidVolumeSource(sc, mc);
	}
	else if (sourcemode == "cylvolume"){
		source = new TCylindricalVolumeSource(sc, mc);
	}
	else if (sourcemode == "STLvolume"){
		source = new TSTLVolumeSource(sc, mc);
	}
	else if (sourcemode == "cylsurface"){
		source = new TCylindricalSurfaceSource(sc, mc, geometry);
	}
	else if (sourcemode == "STLsurface"){
		source = new TSTLSurfaceSource(sc, mc, geometry);
	}
	else
		throw std::runtime_error((boost::format("Could not load source %1%!") % sourcemode).str());
	cout << '\n';

	return source;
}
