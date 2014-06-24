/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "geometry.h"
#include "trianglemesh.h"
#include "globals.h"

/**
 * Virtual base class for all particle sources
 */
class TParticleSource{
protected:
	double fActiveTime; ///< Duration for which the source will be active
	const string fParticleName;
public:
	/**
	 * Constructor, should be called by every derived class
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 */
	TParticleSource(const string ParticleName, double ActiveTime): fActiveTime(ActiveTime), fParticleName(ParticleName){ };

	/**
	 * Destructor
	 */
	virtual ~TParticleSource(){ };

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. This value is usually not set in RandomSourcePoint, but it might be modified
	 * @param x Generated initial x coordinate
	 * @param y Generated initial y coordinate
	 * @param z Generated initial z coordinate
	 * @param phi_v Generated initial azimuth angle of velocity vector
	 * @param theta_v Generated initial polar angle of velocity vector
	 */
	virtual	void RandomSourcePoint(TMCGenerator &mc, long double &t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi_v, long double &theta_v) = 0;
};


/**
 * Virtual base class for surface source.
 *
 * It keeps a list of STL triangles on which initial coordinates are generated
 */
class TSurfaceSource: public TParticleSource{
protected:
	double sourcearea; ///< Net area of source surface
	double Enormal; ///< Boost given to particles starting from this surface
	vector<TTriangle> sourcetris; ///< List of triangles making up the source surface
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param E_normal Energy boost that should be given to particles starting from this surface source
	 */
	TSurfaceSource(const string ParticleName, double ActiveTime, double E_normal): TParticleSource(ParticleName, ActiveTime), Enormal(E_normal), sourcearea(0){ }

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. This value is modified when E_normal was given
	 * @param x Generated initial x coordinate
	 * @param y Generated initial y coordinate
	 * @param z Generated initial z coordinate
	 * @param phi_v Generated initial azimuth angle of velocity vector, distributed according to Lambert's law with an optional boost when E_normal was given
	 * @param theta_v Generated initial polar angle of velocity vector, distributed according to Lambert's law with an optional boost when E_normal was given
	 */
	void RandomSourcePoint(TMCGenerator &mc, long double &t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi_v, long double &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		long double p[3] = {x, y, z};
		long double n[3];
		double RandA = mc.UniformDist(0,sourcearea);
		double SumA = 0;
		vector<TTriangle>::iterator i;
		for (i = sourcetris.begin(); i != sourcetris.end(); i++){
			SumA += i->area();
			if (RandA <= SumA) break;
		}
		double a = mc.UniformDist(0,1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
		double b = mc.UniformDist(0,1);
		if (a+b > 1){
			a = 1 - a;
			b = 1 - b;
		}
		CVector nv = i->normal();
		CPoint v1 = i->tri[0] + a*(i->tri[1] - i->tri[0]) + b*(i->tri[2] - i->tri[0]) + nv*REFLECT_TOLERANCE;
		for (int j = 0; j < 3; j++){
			n[j] = nv[j]/sqrt(nv.squared_length());
			p[j] = v1[j];
		}

		x = p[0];
		y = p[1];
		z = p[2];
		phi_v = mc.UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
		theta_v = mc.SinCosDist(0, 0.5*pi); // Lambert's law!
		if (Enormal > 0){
			double vnormal = sqrt(Ekin*cos(theta_v)*cos(theta_v) + Enormal); // add E_normal to component normal to surface
			double vtangential = sqrt(Ekin)*sin(theta_v);
			theta_v = atan2(vtangential, vnormal); // update angle
			Ekin = vnormal*vnormal + vtangential*vtangential; // update energy
			long double v[3] = {cos(phi_v)*sin(theta_v), sin(phi_v)*sin(theta_v), cos(theta_v)};
			RotateVector(v,n);

			phi_v = atan2(v[1],v[0]);
			theta_v = acos(v[2]);
		}
	}
};

/**
 * Volume source generating points in a cylindrical coordinate range
 */
class TCylindricalVolumeSource: public TParticleSource{
private:
	double rmin, rmax, phimin, phimax, zmin, zmax;
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param r_min Minimal radial coordinate range
	 * @param r_max Maximal radial coordinate range
	 * @param phi_min Minimal azimuthal coordinate range
	 * @param phi_max Maximal azimuthal coordinate range
	 * @param z_min Minimal axial coordinate range
	 * @param z_max Maximal axial coordinate range
	 */
	TCylindricalVolumeSource(const string ParticleName, double ActiveTime, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max)
		: TParticleSource(ParticleName, ActiveTime), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){

	}

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. This value is modified when E_normal was given
	 * @param x Generated initial x coordinate
	 * @param y Generated initial y coordinate
	 * @param z Generated initial z coordinate
	 * @param phi_v Generated initial azimuth angle of velocity vector, isotropically distributed
	 * @param theta_v Generated initial polar angle of velocity vector, isotropically distributed
	 */
	void RandomSourcePoint(TMCGenerator &mc, long double &t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi_v, long double &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		double r = mc.LinearDist(rmin, rmax); // weighting because of the volume element and a r^2 probability outwards
		double phi_r = mc.UniformDist(phimin,phimax);
		x = r*cos(phi_r);
		y = r*sin(phi_r);
		z = mc.UniformDist(zmin,zmax);
		mc.AngularDist(fParticleName, phi_v, theta_v);
	}
};

/**
 * Surface source using a surface inside a cylindrial coordinate range
 */
class TCylindricalSurfaceSource: public TSurfaceSource{
private:
	double rmin, rmax, phimin, phimax, zmin, zmax;

	/**
	 * Check if a point is inside the cylindrical coordinate range
	 */
	bool InSourceVolume(CPoint p){
		double r = sqrt(p[0]*p[0] + p[1]*p[1]);
		double phi = atan2(p[1],p[0]);
		return (r >= rmin && r <= rmax &&
				phi >= phimin && phi <= phimax &&
				p[2] >= zmin && p[2] <= zmax); // check if point is in custom paramter range
	}

public:
	/**
	 * Constructor.
	 *
	 * Selects all triangles from TGeometry whose vertices are inside the cylindrical coordinate range
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param geometry TGeometry class from which the surface triangles shall be taken
	 * @param E_normal Optional energy boost for particles starting from this surface
	 * @param r_min Minimal radial coordinate range
	 * @param r_max Maximal radial coordinate range
	 * @param phi_min Minimal azimuthal coordinate range
	 * @param phi_max Maximal azimuthal coordinate range
	 * @param z_min Minimal axial coordinate range
	 * @param z_max Maximal axial coordinate range
	 */
	TCylindricalSurfaceSource(const string ParticleName, double ActiveTime, TGeometry &geometry, double E_normal, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max)
		: TSurfaceSource(ParticleName, ActiveTime, E_normal), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){
		for (CIterator i = geometry.mesh.triangles.begin(); i != geometry.mesh.triangles.end(); i++){
			if (InSourceVolume(i->tri[0]) && InSourceVolume(i->tri[1]) && InSourceVolume(i->tri[2])){
				sourcetris.push_back(*i);
				sourcearea += i->area();
			}
		}
		printf("Source Area: %g m^2\n",sourcearea);
	}
};

/**
 * Volume source.
 *
 * Starting points are generated within an STL solid
 */
class TSTLVolumeSource: public TParticleSource{
private:
	TTriangleMesh kdtree; ///< internal AABB tree storing the STL solid
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param sourcefile File from which the STL solid shall be read
	 */
	TSTLVolumeSource(const string ParticleName, double ActiveTime, string sourcefile): TParticleSource(ParticleName, ActiveTime){
		kdtree.ReadFile(sourcefile.c_str(),0);
		kdtree.Init();
	}

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. Not modified.
	 * @param x Generated initial x coordinate
	 * @param y Generated initial y coordinate
	 * @param z Generated initial z coordinate
	 * @param phi_v Generated initial azimuth angle of velocity vector, isotropically distributed
	 * @param theta_v Generated initial polar angle of velocity vector, isotropically distributed
	 */
	void RandomSourcePoint(TMCGenerator &mc, long double &t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi_v, long double &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		double p[3];
		for(;;){
			p[0] = mc.UniformDist(kdtree.tree.bbox().xmin(),kdtree.tree.bbox().xmax()); // random point
			p[1] = mc.UniformDist(kdtree.tree.bbox().ymin(),kdtree.tree.bbox().ymax()); // random point
			p[2] = mc.UniformDist(kdtree.tree.bbox().zmin(),kdtree.tree.bbox().zmax()); // random point
			if (kdtree.InSolid(p)){
				mc.AngularDist(fParticleName, phi_v, theta_v);
				x = p[0];
				y = p[1];
				z = p[2];
				break;
			}
		}
	}
};


/**
 * Surface source.
 *
 * Starting points are created on a surface of the experiment geometry whose triangles are all inside the given STL solid
 */
class TSTLSurfaceSource: public TSurfaceSource{
public:
	/**
	 * Constructor.
	 *
	 * Search for all triangles in geometry's mesh which are inside the STL solid given in sourcefile.
	 *
	 * @param ActiveTime Time for which the source is active.
	 * @param geometry Experiment geometry from which the surface is taken.
	 * @param sourcefile STL solid in which the surface should lie.
	 * @param E_normal Give particles starting at this source a velocity boost normal to the surface.
	 */
	TSTLSurfaceSource(const string ParticleName, double ActiveTime, TGeometry &geometry, string sourcefile, double E_normal): TSurfaceSource(ParticleName, ActiveTime, E_normal){
		TTriangleMesh mesh;
		mesh.ReadFile(sourcefile.c_str(),0);
		mesh.Init();
		sourcearea = 0; // add triangles, whose vertices are all in the source volume, to sourcetris list
		for (CIterator i = geometry.mesh.triangles.begin(); i != geometry.mesh.triangles.end(); i++){
			if (mesh.InSolid(i->tri[0]) && mesh.InSolid(i->tri[1]) && mesh.InSolid(i->tri[2])){
				sourcetris.push_back(*i);
				sourcearea += i->area();
			}
		}
		printf("Source Area: %g m^2\n",sourcearea);
	}
};

/**
 * Class which can produce random particle starting points from different sources.
 *
 * There are four source modes which can be specified in the [SOURCE] section of the configuration file:
 * "STLvolume" - random points inside a STL solid are created;
 * "STLsurface" - random points on triangles COMPLETELY surrounded by a STL solid are created;
 * "cylvolume" - random points in a cylindrical coordinate range are created;
 * "cylsurface"	- random points on triangles COMPLETELY inside a cylindrical coordinate range are created
 */
struct TSource{
public:
	string sourcemode; ///< volume/surface/customvol/customsurf
	TParticleSource *source; ///< TParticleSource contructed according to user chosen sourcemode

	/**
	 * Constructor, loads [SOURCE] section of configuration file, calculates TSource::Hmin_lfs and TSource::Hmin_hfs
	 *
	 * @param geometryconf TConfig struct containing [SOURCE] option map
	 * @param geom TGeometry class against which start points are checked and which contains source surfaces
	 * @param field TFieldManager class to calculate TSource::Hmin_lfs and TSource::Hmin_hfs
	 */
	TSource(TConfig &geometryconf, TGeometry &geom, TFieldManager &field): source(NULL){
		sourcemode = geometryconf["SOURCE"].begin()->first; // only first source in geometry.in is read in
		istringstream sourceconf(geometryconf["SOURCE"].begin()->second);
		string ParticleName;
		sourceconf >> ParticleName;

		double ActiveTime;
		if (sourcemode == "cylvolume"){
			double r_min, r_max, phi_min, phi_max, z_min, z_max;
			sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime;
			if (sourceconf)
				source = new TCylindricalVolumeSource(ParticleName, ActiveTime, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max);
		}
		else if (sourcemode == "STLvolume"){
			string sourcefile;
			sourceconf >> sourcefile >> ActiveTime;
			if (sourceconf)
				source = new TSTLVolumeSource(ParticleName, ActiveTime, sourcefile);
		}
		else if (sourcemode == "cylsurface"){
			double r_min, r_max, phi_min, phi_max, z_min, z_max, E_normal;
			sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> E_normal;
			if (sourceconf)
				source = new TCylindricalSurfaceSource(ParticleName, ActiveTime, geom, E_normal, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max);
		}
		else if (sourcemode == "STLsurface"){
			string sourcefile;
			double E_normal;
			sourceconf >> sourcefile >> ActiveTime >> E_normal;
			if (sourceconf)
				source = new TSTLSurfaceSource(ParticleName, ActiveTime, geom, sourcefile, E_normal);
		}

		if (!source){
			cout << "\nCould not load source """ << sourcemode << """! Did you enter invalid parameters?\n";
			exit(-1);
		}
		cout << '\n';
	};


	/**
	 * Destructor. Delete TSource::source.
	 */
	~TSource(){
		if (source)
			delete source;
	}


	/**
	 * Create random point in source volume.
	 *
	 * "volume": Random points with isotropic velocity distribution inside the bounding box of TSource::kdtree are produced until TSource::InSourceVolume is true for a point.
	 * "customvol": Random points with isotropic velocity deistribution in the coordinate range ([TSource::r_min..TSource::r_max]:[TSource::phi_min..TSource::phi_max]:[TSource::z_min..TSource::z_max]) are produced.
	 * "surface/customsurf": A random point on a random triangle from the TSource::sourcetris list is chosen (weighted by triangle area).
	 * 						Starting angles are cosine-distributed to the triangle normal and then angles and Ekin are modified according to TSource::E_normal.
	 * This is repeated for all source modes until the point does not lie inside a solid in TSource::geometry.
	 *
	 * @param mc TMCGenerator to produce random number distributions
	 * @param t Start time of particle
	 * @param Ekin Kinetic start energy of particle, might be modified by E_normal
	 * @param x Returns starting x coordinate of particle
	 * @param y Returns starting y coordinate of particle
	 * @param z Returns starting z coordinate of particle
	 * @param phi_v Returns azimuth of velocity vector
	 * @param theta_v Returns polar angle of velocity vector
	 */
	void RandomPointInSourceVolume(TMCGenerator &mc, long double &t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi_v, long double &theta_v){
		source->RandomSourcePoint(mc, t, Ekin, x, y, z, phi_v, theta_v);
	};
};



#endif /* SOURCE_H_ */
