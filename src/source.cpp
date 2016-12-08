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

TParticleSource::TParticleSource(std::map<std::string, std::string> &sourceconf, const TMCGenerator &mc): fActiveTime(0), polarization(0), ParticleCounter(0){
	istringstream(sourceconf["particle"]) >> fParticleName;
	istringstream(sourceconf["ActiveTime"]) >> fActiveTime;
	istringstream(sourceconf["polarization"]) >> polarization;

	double rmin, rmax;
	istringstream(sourceconf["Emin"]) >> rmin;
	istringstream(sourceconf["Emax"]) >> rmax;
	spectrum = mc.ParseDist(sourceconf["spectrum"], rmin, rmax);

	istringstream(sourceconf["phi_v_min"]) >> rmin;
	istringstream(sourceconf["phi_v_max"]) >> rmax;
	phi_v = mc.ParseDist(sourceconf["phi_v"], rmin*conv, rmax*conv);

	istringstream(sourceconf["theta_v_min"]) >> rmin;
	istringstream(sourceconf["theta_v_max"]) >> rmax;
	theta_v = mc.ParseDist(sourceconf["theta_v"], rmin*conv, rmax*conv);
}


TParticle* TParticleSource::CreateParticle(double t, double x, double y, double z, double E, double phi, double theta, double polarisation,
		const TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
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


TParticle* TSurfaceSource::CreateParticle(const TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	CPoint p;
	CVector nv;
	do{
		double area = mc.UniformDist(0., area_sum.back());
		int index = std::distance(area_sum.begin(), std::lower_bound(area_sum.begin(), area_sum.end(), area));
		CTriangle tri = (geometry.mesh.GetTrianglesBegin() + index)->first;
		double a = mc.UniformDist(0, 1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
		double b = mc.UniformDist(0, 1);
		if (a+b > 1){
			a = 1 - a;
			b = 1 - b;
		}
		nv = tri.supporting_plane().orthogonal_vector();
		nv = nv/sqrt(nv.squared_length());
		p = tri[0] + a*(tri[1] - tri[0]) + b*(tri[2] - tri[0]);
	}while (!InSourceVolume(p[0], p[1], p[2]));
	p = p + nv*REFLECT_TOLERANCE; // move point slightly away from surface

	double Ekin = mc.SampleDist(spectrum);
	double phi_v = mc.UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
	double theta_v = mc.SinCosDist(0, 0.5*pi); // Lambert's law!
	if (Enormal > 0){
		double vnormal = sqrt(Ekin*cos(theta_v)*cos(theta_v) + Enormal); // add E_normal to component normal to surface
		double vtangential = sqrt(Ekin)*sin(theta_v);
		theta_v = atan2(vtangential, vnormal); // update angle
		Ekin = vnormal*vnormal + vtangential*vtangential; // update energy
	}

	double v[3] = {cos(phi_v)*sin(theta_v), sin(phi_v)*sin(theta_v), cos(theta_v)};
	double n[3] = {nv[0], nv[1], nv[2]};
	RotateVector(v, n);
	phi_v = atan2(v[1],v[0]);
	theta_v = acos(v[2]);

	return TParticleSource::CreateParticle(mc.UniformDist(0., fActiveTime), p[0], p[1], p[2], Ekin, phi_v, theta_v, polarization, mc, geometry, field);
}


void TVolumeSource::FindPotentialMinimum(const TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
	cout << "Sampling phase space ";
	const int N = 100000;
	int percent = 0;
	for (int i = 0; i < N; i++){
		PrintPercent((double)i/N, percent);
		double t = mc.UniformDist(0, fActiveTime); // dice start time
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

TParticle* TVolumeSource::CreateParticle(const TMCGenerator &mc, const TGeometry &geometry, const TFieldManager &field){
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
			H = mc.SampleDist(spectrum); // dice total(!) energy until one above minimum potential is found
		}while (H < MinPot);

		cout << "Trying to find starting point for particle with total energy " << H << "eV ...";
		for (int i = 0; true; i++){
			double t = mc.UniformDist(0, fActiveTime); // dice start time
			double x, y, z;
			RandomPointInSourceVolume(x, y, z, mc); // dice point in source volume
			TParticle *proposed_p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarization, mc, geometry, field); // create new particle with Ekin = 0
			double V = proposed_p->GetInitialTotalEnergy(geometry, field); // potential at particle position equals its total energy
			delete proposed_p;
			ParticleCounter--; // delete particle and decrement particle counter

			if (H < V)
				continue; // if total energy < potential energy then particle is not possible at this point
			if (sqrt(H - V) > mc.UniformDist(0, sqrt(H - MinPot))){ // accept particle with probability sqrt(H-V)/sqrt(H-Vmin) (phase space weighting according to Golub)
				cout << " found after " << i+1 << " tries\n";
				return TParticleSource::CreateParticle(t, x, y, z, H - V, mc.SampleDist(phi_v), mc.SampleDist(theta_v), polarization, mc, geometry, field); // if accepted, return new particle with correct Ekin
			}
		}
		assert(false); // this will never be reached
	}
	else{ // create particles uniformly distributed in volume
		double t = mc.UniformDist(0, fActiveTime);
		double x, y, z;
		RandomPointInSourceVolume(x, y, z, mc);
		cout << '\n';
		return TParticleSource::CreateParticle(t, x, y, z, mc.SampleDist(spectrum), mc.SampleDist(phi_v), mc.SampleDist(theta_v), polarization, mc, geometry, field);
	}
}

TParticleSource* CreateParticleSource(TConfig &config, const TMCGenerator &mc, const TGeometry &geometry){
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
