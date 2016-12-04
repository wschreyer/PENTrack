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

TParticle* TParticleSource::CreateParticle(double t, double x, double y, double z, double E, double phi, double theta, double polarisation){
	TParticle *p;
	if (fParticleName == NAME_NEUTRON)
		p = new TNeutron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_PROTON)
		p = new TProton(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_ELECTRON)
		p = new TElectron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_MERCURY)
		p = new TMercury(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_XENON) 
		p = new TXenon(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else{
		cout << "Could not create particle " << fParticleName << '\n';
		exit(-1);
	}
	return p;
}


TParticle* TSurfaceSource::CreateParticle(){
	CPoint p;
	CVector nv;
	do{
		double area = fmc->UniformDist(0., area_sum.back());
		int index = std::distance(area_sum.begin(), std::lower_bound(area_sum.begin(), area_sum.end(), area));
		CTriangle tri = (fgeom->mesh.GetTrianglesBegin() + index)->first;
		double a = fmc->UniformDist(0, 1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
		double b = fmc->UniformDist(0, 1);
		if (a+b > 1){
			a = 1 - a;
			b = 1 - b;
		}
		nv = tri.supporting_plane().orthogonal_vector();
		nv = nv/sqrt(nv.squared_length());
		p = tri[0] + a*(tri[1] - tri[0]) + b*(tri[2] - tri[0]);
	}while (!InSourceVolume(p[0], p[1], p[2]));
	p = p + nv*REFLECT_TOLERANCE; // move point slightly away from surface

	double Ekin = fmc->Spectrum(fParticleName);
	double phi_v = fmc->UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
	double theta_v = fmc->SinCosDist(0, 0.5*pi); // Lambert's law!
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
	double polarisation = fmc->pconfigs[fParticleName].polarization;

	return TParticleSource::CreateParticle(fmc->UniformDist(0., fActiveTime), p[0], p[1], p[2], Ekin, phi_v, theta_v, polarisation);
}


void TVolumeSource::FindPotentialMinimum(){
	cout << "Sampling phase space ";
	const int N = 100000;
	int percent = 0;
	for (int i = 0; i < N; i++){
		PrintPercent((double)i/N, percent);
		double t = fmc->UniformDist(0, fActiveTime); // dice start time
		double polarisation = fmc->pconfigs[fParticleName].polarization; // dice polarisation
		double x, y, z;
		RandomPointInSourceVolume(x, y, z); // dice point in source volume
		TParticle *p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarisation); // create dummy particle with Ekin = 0
		double V = p->GetInitialTotalEnergy(); // potential at particle position equals its total energy
		if (V < MinPot)
			MinPot = V; // remember minimal potential
		delete p;
		ParticleCounter--;
	}
	cout << " minimal potential = " << MinPot << "eV\n";
}

TParticle* TVolumeSource::CreateParticle(){
	if (fPhaseSpaceWeighting){ // if particle density should be weighted by available phase space
		TParticleConfig *mcconf = &fmc->pconfigs[fParticleName];
		if (MinPot == numeric_limits<double>::infinity()){ // if minimum potential energy has not yet been determined
			FindPotentialMinimum(); // find minimum potential energy
			if (MinPot > mcconf->spectrum.max()) // abort program if spectrum completely out of potential range
				throw std::runtime_error( (boost::format("Error: your chosen spectrum is below the minimal potential energy in the source volume (%1% eV < %2% eV). Exiting!\n") % mcconf->spectrum.max() % MinPot).str() );
		}

		if (MinPot > mcconf->spectrum.min()){ // give warning if chosen spectrum contains energy ranges that are not possible
			cout << "Warning: your chosen spectrum contains energies below the minimal potential energy in the source volume (" << mcconf->spectrum.min() << "eV < " << MinPot << "eV). The energy spectrum will be cut off!\n";
		}
		double H;
		do{
			H = fmc->Spectrum(fParticleName); // dice total(!) energy until one above minimum potential is found
		}while (H < MinPot);

		cout << "Trying to find starting point for particle with total energy " << H << "eV ...";
		for (int i = 0; true; i++){
			double t = fmc->UniformDist(0, fActiveTime); // dice start time
			double polarisation = mcconf->polarization; // dice polarisation
			double x, y, z;
			RandomPointInSourceVolume(x, y, z); // dice point in source volume
			TParticle *proposed_p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarisation); // create new particle with Ekin = 0
			double V = proposed_p->GetInitialTotalEnergy(); // potential at particle position equals its total energy
			delete proposed_p;
			ParticleCounter--; // delete particle and decrement particle counter

			if (H < V)
				continue; // if total energy < potential energy then particle is not possible at this point
			if (sqrt(H - V) > fmc->UniformDist(0, sqrt(H - MinPot))){ // accept particle with probability sqrt(H-V)/sqrt(H-Vmin) (phase space weighting according to Golub)
				cout << " found after " << i+1 << " tries\n";
				double phi_v, theta_v;
				fmc->AngularDist(fParticleName, phi_v, theta_v); // dice velocity direction
				return TParticleSource::CreateParticle(t, x, y, z, H - V, phi_v, theta_v, polarisation); // if accepted, return new particle with correct Ekin
			}
		}
		assert(false); // this will never be reached
	}
	else{ // create particles uniformly distributed in volume
		double t = fmc->UniformDist(0, fActiveTime);
		double E = fmc->Spectrum(fParticleName);
		double phi_v, theta_v;
		fmc->AngularDist(fParticleName, phi_v, theta_v);
		double polarisation = fmc->pconfigs[fParticleName].polarization;
		double x, y, z;
		RandomPointInSourceVolume(x, y, z);
		cout << '\n';
		return TParticleSource::CreateParticle(t, x, y, z, E, phi_v, theta_v, polarisation);
	}
}



TSource::TSource(TConfig &geometryconf, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: source(NULL){
	sourcemode = geometryconf["SOURCE"].begin()->first; // only first source in geometry.in is read in
	istringstream sourceconf(geometryconf["SOURCE"].begin()->second);
	string ParticleName;
	sourceconf >> ParticleName;

	double ActiveTime;
	bool PhaseSpaceWeighting;
	if (sourcemode == "boxvolume"){
		double x_min, x_max, y_min, y_max, z_min, z_max;
		sourceconf >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
		if (sourceconf)
			source = new TCuboidVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, x_min, x_max, y_min, y_max, z_min, z_max, mc, geometry, field);
	}
	else if (sourcemode == "cylvolume"){
		double r_min, r_max, phi_min, phi_max, z_min, z_max;
		sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
		if (sourceconf)
			source = new TCylindricalVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max, mc, geometry, field);
	}
	else if (sourcemode == "STLvolume"){
		string sourcefile;
		sourceconf >> sourcefile >> ActiveTime >> PhaseSpaceWeighting;
		if (sourceconf)
			source = new TSTLVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, sourcefile, mc, geometry, field);
	}
	else if (sourcemode == "cylsurface"){
		double r_min, r_max, phi_min, phi_max, z_min, z_max, E_normal;
		sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> E_normal;
		if (sourceconf)
			source = new TCylindricalSurfaceSource(ParticleName, ActiveTime, E_normal, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max, mc, geometry, field);
	}
	else if (sourcemode == "STLsurface"){
		string sourcefile;
		double E_normal;
		sourceconf >> sourcefile >> ActiveTime >> E_normal;
		if (sourceconf)
			source = new TSTLSurfaceSource(ParticleName, ActiveTime, sourcefile, E_normal, mc, geometry, field);
	}

	if (!source){
		cout << "\nCould not load source """ << sourcemode << """! Did you enter invalid parameters?\n";
		exit(-1);
	}
	cout << '\n';
}


TSource::~TSource(){
	if (source)
		delete source;
}


TParticle* TSource::CreateParticle(){
	return source->CreateParticle();
}
