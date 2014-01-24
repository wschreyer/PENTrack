/*
 * source.h
 *
 *  Created on: 10.12.2013
 *      Author: wolfgang
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "geometry.h"
#include "kdtree.h"
#include "globals.h"

/**
 * Virtual base class for all particle sources
 */
class TParticleSource{
protected:
	Doub fActiveTime; ///< Duration for which the source will be active
public:
	/**
	 * Constructor, should be called by every derived class
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 */
	TParticleSource(Doub ActiveTime): fActiveTime(ActiveTime){ };

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
	 * @param x Generated initial y coordinate
	 * @param x Generated initial z coordinate
	 * @param x Generated initial azimuth angle of velocity vector
	 * @param x Generated initial polar angle of velocity vector
	 */
	virtual	void RandomSourcePoint(TMCGenerator &mc, Doub &t, Doub &Ekin, Doub &x, Doub &y, Doub &z, Doub &phi_v, Doub &theta_v) = 0;
};


/**
 * Virtual base class for surface source.
 *
 * It keeps a list of STL triangles on which initial coordinates are generated
 */
class TSurfaceSource: public TParticleSource{
protected:
	Doub sourcearea; ///< Net area of source surface
	Doub Enormal; ///< Boost given to particles starting from this surface
#ifndef USE_CGAL
	vector<Triangle*> sourcetris; ///< list of triangles (part of TSource::geometry) inside the source volume (only relevant if source mode is "surface" or "customsurf")
#else
	vector<TTriangle*> sourcetris;
#endif
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param E_normal Energy boost that should be given to particles starting from this surface source
	 */
	TSurfaceSource(Doub ActiveTime, Doub E_normal): TParticleSource(ActiveTime), Enormal(E_normal), sourcearea(0){ }

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. This value is modified when E_normal was given
	 * @param x Generated initial x coordinate
	 * @param x Generated initial y coordinate
	 * @param x Generated initial z coordinate
	 * @param x Generated initial azimuth angle of velocity vector, distributed according to Lambert's law with an optional boost when E_normal was given
	 * @param x Generated initial polar angle of velocity vector, distributed according to Lambert's law with an optional boost when E_normal was given
	 */
	void RandomSourcePoint(TMCGenerator &mc, Doub &t, Doub &Ekin, Doub &x, Doub &y, Doub &z, Doub &phi_v, Doub &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		Doub p[3] = {x, y, z};
		Doub n[3];
		Doub RandA = mc.UniformDist(0,sourcearea);
		Doub CurrA = 0, SumA = 0;
#ifndef USE_CGAL
		vector<Triangle*>::iterator i;
		for (i = sourcetris.begin(); i != sourcetris.end(); i++){
			n[0] = (*i)->normal[0];
			n[1] = (*i)->normal[1];
			n[2] = (*i)->normal[2];
			CurrA = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			SumA += CurrA;
			if (RandA <= SumA) break; // select random triangle, weighted by triangle area
		}
#else
		vector<TTriangle*>::iterator i;
		for (i = sourcetris.begin(); i != sourcetris.end(); i++){
			SumA += sqrt((*i)->squared_area());
			if (RandA <= SumA) break;
		}
#endif
		Doub a = mc.UniformDist(0,1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
		Doub b = mc.UniformDist(0,1);
		if (a+b > 1){
			a = 1 - a;
			b = 1 - b;
		}
#ifndef USE_CGAL
		for (int j = 0; j < 3; j++){
			n[j] /= CurrA; // normalize normal vector
			p1[j] = (*i)->vertex[0][j] + a*((*i)->vertex[1][j] - (*i)->vertex[0][j]) + b*((*i)->vertex[2][j] - (*i)->vertex[0][j]);
			p1[j] += REFLECT_TOLERANCE*n[j]; // add some tolerance to avoid starting inside solid
		}
#else
		K::Vector_3 nv = (*i)->supporting_plane().orthogonal_vector();
		CPoint v1 = (**i)[0] + a*((**i)[1] - (**i)[0]) + b*((**i)[2] - (**i)[0]) + nv*REFLECT_TOLERANCE;
		for (int j = 0; j < 3; j++){
			n[j] = nv[j]/sqrt(nv.squared_length());
			p[j] = v1[j];
		}
#endif

		x = p[0];
		y = p[1];
		z = p[2];
		phi_v = mc.UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
		theta_v = mc.SinCosDist(0, 0.5*pi); // Lambert's law!
		if (Enormal > 0){
			Doub vnormal = sqrt(Ekin*cos(theta_v)*cos(theta_v) + Enormal); // add E_normal to component normal to surface
			Doub vtangential = sqrt(Ekin)*sin(theta_v);
			theta_v = atan2(vtangential, vnormal); // update angle
			Ekin = vnormal*vnormal + vtangential*vtangential; // update energy
			Doub v[3] = {cos(phi_v)*sin(theta_v), sin(phi_v)*sin(theta_v), cos(theta_v)};
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
	Doub rmin, rmax, phimin, phimax, zmin, zmax;
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
	TCylindricalVolumeSource(Doub ActiveTime, Doub r_min, Doub r_max, Doub phi_min, Doub phi_max, Doub z_min, Doub z_max)
		: TParticleSource(ActiveTime), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){

	}

	/**
	 * Generate a random time, point and velocity direction
	 *
	 * @param mc TMCGenerator class, used to generate different random distributions
	 * @param t Generated start time, usually within [0..ActiveTime]
	 * @param Ekin Initial kinetic energy. This value is modified when E_normal was given
	 * @param x Generated initial x coordinate
	 * @param x Generated initial y coordinate
	 * @param x Generated initial z coordinate
	 * @param x Generated initial azimuth angle of velocity vector, isotropically distributed
	 * @param x Generated initial polar angle of velocity vector, isotropically distributed
	 */
	void RandomSourcePoint(TMCGenerator &mc, Doub &t, Doub &Ekin, Doub &x, Doub &y, Doub &z, Doub &phi_v, Doub &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		Doub r = mc.LinearDist(rmin, rmax); // weighting because of the volume element and a r^2 probability outwards
		Doub phi_r = mc.UniformDist(phimin,phimax);
		x = r*cos(phi_r);
		y = r*sin(phi_r);
		z = mc.UniformDist(zmin,zmax);
		mc.IsotropicDist(phi_v, theta_v);
	}
};

/**
 * Surface source using a surface inside a cylindrial coordinate range
 */
class TCylindricalSurfaceSource: public TSurfaceSource{
private:
	Doub rmin, rmax, phimin, phimax, zmin, zmax;

	/**
	 * Check if a point is inside the cylindrical coordinate range
	 */
	bool InSourceVolume(const float p[3]){
		Doub r = sqrt(p[0]*p[0] + p[1]*p[1]);
		Doub phi = atan2(p[1],p[0]);
		return (r >= rmin && r <= rmax &&
				phi >= phimin && phi <= phimax &&
				p[2] >= zmin && p[2] <= zmax); // check if point is in custom paramter range
	}

#ifdef USE_CGAL
	bool InSourceVolume(CPoint p){
		float pp[3] = {p[0],p[1],p[2]};
		return InSourceVolume(pp);
	}
#endif

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
	TCylindricalSurfaceSource(Doub ActiveTime, TGeometry &geometry, Doub E_normal, Doub r_min, Doub r_max, Doub phi_min, Doub phi_max, Doub z_min, Doub z_max)
		: TSurfaceSource(ActiveTime, E_normal), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){
#ifndef USE_CGAL
		for (vector<Triangle>::iterator i = geometry.kdtree->alltris.begin(); i != geometry.kdtree->alltris.end(); i++){
			if (InSourceVolume(i->vertex[0]) && InSourceVolume(i->vertex[1]) && InSourceVolume(i->vertex[2])){
				sourcetris.push_back(&*i);
				long double *n = i->normal;
				sourcearea += sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			}
		}
#else
		for (CIterator i = geometry.kdtree->triangles.begin(); i != geometry.kdtree->triangles.end(); i++){
			if (InSourceVolume((*i)[0]) && InSourceVolume((*i)[1]) && InSourceVolume((*i)[2])){
				sourcetris.push_back(&*i);
				sourcearea += sqrt(i->squared_area());
			}
		}
#endif
		printf("Source Area: %LG m^2\n",sourcearea);
	}
};

/**
 * Volume source.
 *
 * Starting points are generated within an STL solid
 */
class TSTLVolumeSource: public TParticleSource{
private:
	KDTree kdtree; ///< internal KDTree storing the STL solid
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param sourcefile File from which the STL solid shall be read
	 */
	TSTLVolumeSource(Doub ActiveTime, string sourcefile): TParticleSource(ActiveTime){
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
	 * @param x Generated initial y coordinate
	 * @param x Generated initial z coordinate
	 * @param x Generated initial azimuth angle of velocity vector, isotropically distributed
	 * @param x Generated initial polar angle of velocity vector, isotropically distributed
	 */
	void RandomSourcePoint(TMCGenerator &mc, Doub &t, Doub &Ekin, Doub &x, Doub &y, Doub &z, Doub &phi_v, Doub &theta_v){
		t = mc.UniformDist(0, fActiveTime);
		Doub p[3];
		for(;;){
#ifndef USE_CGAL
			p[0] = mc.UniformDist(kdtree->lo[0],kdtree->hi[0]); // random point
			p[1] = mc.UniformDist(kdtree->lo[1],kdtree->hi[1]); // random point
			p[2] = mc.UniformDist(kdtree->lo[2],kdtree->hi[2]); // random point
#else
			p[0] = mc.UniformDist(kdtree.tree.bbox().xmin(),kdtree.tree.bbox().xmax()); // random point
			p[1] = mc.UniformDist(kdtree.tree.bbox().ymin(),kdtree.tree.bbox().ymax()); // random point
			p[2] = mc.UniformDist(kdtree.tree.bbox().zmin(),kdtree.tree.bbox().zmax()); // random point
#endif
			if (kdtree.InSolid(p)){
				mc.IsotropicDist(phi_v, theta_v);
				x = p[0];
				y = p[1];
				z = p[2];
				break;
			}
		}
	}
	virtual void CalcHmin(TFieldManager *field = NULL){

	}
};

class TSTLSurfaceSource: public TSurfaceSource{
public:
	TSTLSurfaceSource(Doub ActiveTime, TGeometry &geometry, string sourcefile, Doub E_normal): TSurfaceSource(ActiveTime, E_normal){
		KDTree kdtree;
		kdtree.ReadFile(sourcefile.c_str(),0);
		kdtree.Init();
		sourcearea = 0; // add triangles, whose vertices are all in the source volume, to sourcetris list
#ifndef USE_CGAL
		for (vector<Triangle>::iterator i = geometry.kdtree->alltris.begin(); i != geometry.kdtree->alltris.end(); i++){
			if (kdtree.InSolid(i->vertex[0]) && kdtree.InSolid(i->vertex[1]) && kdtree.InSolid(i->vertex[2])){
				sourcetris.push_back(&*i);
				long double *n = i->normal;
				sourcearea += sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			}
		}
#else
		for (CIterator i = geometry.kdtree->triangles.begin(); i != geometry.kdtree->triangles.end(); i++){
			if (kdtree.InSolid((*i)[0]) && kdtree.InSolid((*i)[1]) && kdtree.InSolid((*i)[2])){
				sourcetris.push_back(&*i);
				sourcearea += sqrt(i->squared_area());
			}
		}
#endif
		printf("Source Area: %LG m^2\n",sourcearea);
	}
	virtual void CalcHmin(TFieldManager *field = NULL){

	}
};

/**
 * Class which can produce random particle starting points from different sources.
 *
 * There are four source modes which can be specified in the [SOURCE] section of the configuration file:
 * "volume" - random points inside a STL solid are created;
 * "surface" - random points on triangles COMPLETELY surrounded by a STL solid are created;
 * "customvol" - random points in a cylindrical coordinate range are created;
 * "customsurf"	- random points on triangles COMPLETELY inside a cylindrical coordinate range are created
 */
struct TSource{
public:
	string sourcemode; ///< volume/surface/customvol/customsurf
	TParticleSource *source;

	/**
	 * Constructor, loads [SOURCE] section of configuration file, calculates TSource::Hmin_lfs and TSource::Hmin_hfs
	 *
	 * @param geometryin Configuration file containing [SOURCE] section
	 * @param geom TGeometry class against which start points are checked and which contains source surfaces
	 * @param field TFieldManager class to calculate TSource::Hmin_lfs and TSource::Hmin_hfs
	 */
	TSource(const char *geometryin, TGeometry &geom, TFieldManager &field){
		TConfig geometryconf;
		ReadInFile(geometryin, geometryconf);
		sourcemode = geometryconf["SOURCE"].begin()->first; // only first source in geometry.in is read in
		istringstream sourceconf(geometryconf["SOURCE"].begin()->second);

		Doub ActiveTime;
		if (sourcemode == "customvol"){
			Doub r_min, r_max, phi_min, phi_max, z_min, z_max;
			sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime;
			source = new TCylindricalVolumeSource(ActiveTime, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max);
		}
		else if (sourcemode == "volume"){
			string sourcefile;
			sourceconf >> sourcefile >> ActiveTime;
			source = new TSTLVolumeSource(ActiveTime, sourcefile);
		}
		else if (sourcemode == "customsurf"){
			Doub r_min, r_max, phi_min, phi_max, z_min, z_max, E_normal;
			sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> E_normal;
			source = new TCylindricalSurfaceSource(ActiveTime, geom, E_normal, r_min, r_max, phi_min*conv, phi_max*conv, z_min, z_max);
		}
		else if (sourcemode == "surface"){
			string sourcefile;
			Doub E_normal;
			sourceconf >> sourcefile >> ActiveTime >> E_normal;
			source = new TSTLSurfaceSource(ActiveTime, geom, sourcefile, E_normal);
		}
		else{
			cout << "Uknown source " << sourcemode << "! Stopping!\n";
			exit(-1);
		}
	};


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
	void RandomPointInSourceVolume(TMCGenerator &mc, Doub &t, Doub &Ekin, Doub &x, Doub &y, Doub &z, Doub &phi_v, Doub &theta_v){
		source->RandomSourcePoint(mc, t, Ekin, x, y, z, phi_v, theta_v);
	};
};



#endif /* SOURCE_H_ */
