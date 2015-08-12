/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#include "source.h"
#include "neutron.h"
#include "proton.h"
#include "electron.h"
#include "geometry.h"
#include "mc.h"
#include "trianglemesh.h"
#include "globals.h"

TParticleSource::TParticleSource(const string ParticleName, double ActiveTime, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: fActiveTime(ActiveTime), fParticleName(ParticleName), fmc(&mc), fgeom(&geometry), ffield(field), ParticleCounter(0){

}

TParticleSource::~TParticleSource(){

}


TParticle* TParticleSource::CreateParticle(double t, double x, double y, double z, double E, double phi, double theta, int polarisation){
	TParticle *p;
	if (fParticleName == NAME_NEUTRON)
		p = new TNeutron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_PROTON)
		p = new TProton(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else if (fParticleName == NAME_ELECTRON)
		p = new TElectron(++ParticleCounter, t, x, y, z, E, phi, theta, polarisation, *fmc, *fgeom, ffield);
	else{
		cout << "Could not create particle " << fParticleName << '\n';
		exit(-1);
	}
	return p;
}


TSurfaceSource::TSurfaceSource(const string ParticleName, double ActiveTime, double E_normal, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
		: TParticleSource(ParticleName, ActiveTime, mc, geometry, field), sourcearea(0), Enormal(E_normal){

}


TParticle* TSurfaceSource::CreateParticle(){
	double t = fmc->UniformDist(0, fActiveTime);
	double RandA = fmc->UniformDist(0,sourcearea);
	double SumA = 0;
	vector<TTriangle>::iterator i;
	for (i = sourcetris.begin(); i != sourcetris.end(); i++){
		SumA += i->area();
		if (RandA <= SumA) break;
	}
	double a = fmc->UniformDist(0,1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
	double b = fmc->UniformDist(0,1);
	if (a+b > 1){
		a = 1 - a;
		b = 1 - b;
	}
	CVector nv = i->normal();
	CPoint p = i->tri[0] + a*(i->tri[1] - i->tri[0]) + b*(i->tri[2] - i->tri[0]) + nv*REFLECT_TOLERANCE;

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
	int polarisation = fmc->DicePolarisation(fParticleName);

	return TParticleSource::CreateParticle(t, p[0], p[1], p[2], Ekin, phi_v, theta_v, polarisation);
}


TVolumeSource::TVolumeSource(std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
			:TParticleSource(ParticleName, ActiveTime, mc, geometry, field), MinPot(numeric_limits<double>::infinity()), fPhaseSpaceWeighting(PhaseSpaceWeighting){
}

void TVolumeSource::FindPotentialMinimum(){
	cout << "Sampling phase space ";
	const int N = 100000;
	int percent = 0;
	for (int i = 0; i < N; i++){
		PrintPercent((double)i/N, percent);
		double t = fmc->UniformDist(0, fActiveTime); // dice start time
		int polarisation = fmc->DicePolarisation(fParticleName); // dice polarisation
		double x, y, z;
		RandomPointInSourceVolume(x, y, z); // dice point in source volume
		TParticle *p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarisation); // create dummy particle with Ekin = 0
		double V = p->Hstart(); // potential at particle position equals its total energy
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
			if (MinPot > mcconf->Emax){ // abort program if spectrum completely out of potential range
				cout << "Error: your chosen spectrum is below the minimal potential energy in the source volume (" << mcconf->Emax << "eV < " << MinPot << "eV). Exiting!\n";
				exit(-1);
			}
		}

		if (MinPot > mcconf->Emin){ // give warning if chosen spectrum contains energy ranges that are not possible
			cout << "Warning: your chosen spectrum contains energies below the minimal potential energy in the source volume (" << mcconf->Emin << "eV < " << MinPot << "eV). The energy spectrum will be cut off!\n";
		}
		double H;
		do{
			H = fmc->Spectrum(fParticleName); // dice total(!) energy until one above minimum potential is found
		}while (H < MinPot);

		cout << "Trying to find starting point for particle with total energy " << H << "eV ...";
		for (int i = 0; true; i++){
			double t = fmc->UniformDist(0, fActiveTime); // dice start time
			int polarisation = fmc->DicePolarisation(fParticleName); // dice polarisation
			double x, y, z;
			RandomPointInSourceVolume(x, y, z); // dice point in source volume
			TParticle *proposed_p = TParticleSource::CreateParticle(t, x, y, z, 0, 0, 0, polarisation); // create new particle with Ekin = 0
			double V = proposed_p->Hstart(); // potential at particle position equals its total energy
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
		return NULL; // this will never be reached
	}
	else{ // create particles uniformly distributed in volume
		double t = fmc->UniformDist(0, fActiveTime);
		double E = fmc->Spectrum(fParticleName);
		double phi_v, theta_v;
		fmc->AngularDist(fParticleName, phi_v, theta_v);
		int polarisation = fmc->DicePolarisation(fParticleName);
		double x, y, z;
		RandomPointInSourceVolume(x, y, z);
		cout << '\n';
		return TParticleSource::CreateParticle(t, x, y, z, E, phi_v, theta_v, polarisation);
	}
}



TCuboidVolumeSource::TCuboidVolumeSource(const string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field), xmin(x_min), xmax(x_max), ymin(y_min), ymax(y_max), zmin(z_min), zmax(z_max){

}


void TCuboidVolumeSource::RandomPointInSourceVolume(double &x, double &y, double &z){
	x = fmc->UniformDist(xmin, xmax);
	y = fmc->UniformDist(ymin, ymax);
	z = fmc->UniformDist(zmin, zmax);
}



TCylindricalVolumeSource::TCylindricalVolumeSource(const string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){

}


void TCylindricalVolumeSource::RandomPointInSourceVolume(double &x, double &y, double &z){
	double r = fmc->LinearDist(rmin, rmax); // weighting because of the volume element and a r^2 probability outwards
	double phi_r = fmc->UniformDist(phimin,phimax);
	x = r*cos(phi_r);
	y = r*sin(phi_r);
	z = fmc->UniformDist(zmin,zmax);
}


bool TCylindricalSurfaceSource::InSourceVolume(CPoint p){
	double r = sqrt(p[0]*p[0] + p[1]*p[1]);
	double phi = atan2(p[1],p[0]);
	return (r >= rmin && r <= rmax &&
			phi >= phimin && phi <= phimax &&
			p[2] >= zmin && p[2] <= zmax); // check if point is in custom paramter range
}


TCylindricalSurfaceSource::TCylindricalSurfaceSource(const string ParticleName, double ActiveTime, double E_normal, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: TSurfaceSource(ParticleName, ActiveTime, E_normal, mc, geometry, field), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){
	for (CIterator i = geometry.mesh.triangles.begin(); i != geometry.mesh.triangles.end(); i++){
		if (InSourceVolume(i->tri[0]) && InSourceVolume(i->tri[1]) && InSourceVolume(i->tri[2])){
			sourcetris.push_back(*i);
			sourcearea += i->area();
		}
	}
	printf("Source Area: %g m^2\n",sourcearea);
}


TSTLVolumeSource::TSTLVolumeSource(const string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, string sourcefile, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field){
	kdtree.ReadFile(sourcefile.c_str(),0);
	kdtree.Init();
}


void TSTLVolumeSource::RandomPointInSourceVolume(double &x, double &y, double &z){
	double p[3];
	for(;;){
		p[0] = fmc->UniformDist(kdtree.tree.bbox().xmin(),kdtree.tree.bbox().xmax()); // random point
		p[1] = fmc->UniformDist(kdtree.tree.bbox().ymin(),kdtree.tree.bbox().ymax()); // random point
		p[2] = fmc->UniformDist(kdtree.tree.bbox().zmin(),kdtree.tree.bbox().zmax()); // random point
		if (kdtree.InSolid(p)){
			x = p[0];
			y = p[1];
			z = p[2];
			break;
		}
	}
}


TSTLSurfaceSource::TSTLSurfaceSource(const string ParticleName, double ActiveTime, string sourcefile, double E_normal, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field)
	: TSurfaceSource(ParticleName, ActiveTime, E_normal, mc, geometry, field){
	TTriangleMesh mesh;
	mesh.ReadFile(sourcefile.c_str(),0);
	mesh.Init();
	sourcearea = 0; // add triangles, whose vertices are all in the source volume, to sourcetris list
	for (CIterator i = fgeom->mesh.triangles.begin(); i != fgeom->mesh.triangles.end(); i++){
		if (mesh.InSolid(i->tri[0]) && mesh.InSolid(i->tri[1]) && mesh.InSolid(i->tri[2])){
			sourcetris.push_back(*i);
			sourcearea += i->area();
		}
	}
	printf("Source Area: %g m^2\n",sourcearea);
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
