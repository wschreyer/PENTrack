/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include <string>
#include <limits>

#include "particle.h"
#include "mc.h"

/**
 * Virtual base class for all particle sources
 */
class TParticleSource{
protected:
	double fActiveTime; ///< Duration for which the source will be active
	const std::string fParticleName; ///< Name of particle that the source should create
	TMCGenerator *fmc; ///< TMCGenerator class passed in constructor
	TGeometry *fgeom; ///< TGeometry class passed in constructor
	TFieldManager *ffield; ///< TFieldManager class passed in constructor
public:
	int ParticleCounter; ///< Count number of particles created by source
	/**
	 * Constructor, should be called by every derived class
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TParticleSource(const std::string ParticleName, double ActiveTime, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		fActiveTime(ActiveTime), fParticleName(ParticleName), fmc(&mc), fgeom(&geometry), ffield(field), ParticleCounter(0){

	}

	/**
	 * Destructor
	 */
	virtual ~TParticleSource(){

	}

	/**
	 * Basic creation of a new particle. Maps particle name to the corresponding class
	 *
	 * @param t Starting time
	 * @param x x coordinate of creation point
	 * @param y y coordinate of creation point
	 * @param z z coordinate of creation point
	 * @param E Initial kinetic energy
	 * @param phi Azimuthal angle of initial velocity vector
	 * @param theta Polar angle of initial velocity vector
	 * @param polarisation Initial polarisation of particle (-1, 0, 1)
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle(double t, double x, double y, double z, double E, double phi, double theta, int polarisation);


	/**
	 * Virtual routine that has to be implemented by every derived source class
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	virtual TParticle* CreateParticle() = 0;
};


/**
 * Virtual base class for surface source.
 *
 * It keeps a list of STL triangles on which initial coordinates are generated
 */
class TSurfaceSource: public TParticleSource{
protected:
	double Enormal; ///< Boost given to particles starting from this surface
	std::vector<double> area_sum; ///< list of areas of corresponding triangles in TGeometry class, each entry accumulates areas of all preceding triangles that are contained in the source volume

	/**
	 * Check if point is inside the source volume.
	 *
	 * Abstract function, has to be implemented by every derived class.
	 */
	virtual bool InSourceVolume(const double x, const double y, const double z) const = 0;
public:
	/**
	 * Constructor.
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param E_normal Energy boost that should be given to particles starting from this surface source
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TSurfaceSource(const std::string ParticleName, double ActiveTime, double E_normal, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TParticleSource(ParticleName, ActiveTime, mc, geometry, field), Enormal(E_normal){

	}

	/**
	 * Create new particle on surface
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle();
};



/**
 * Virtual base class for volume source.
 *
 * Can generate a particle in volume source, optionally weighted by available phase space.
 * Derived sources need to implement the RandomPointInSourceVolume routine.
 */
class TVolumeSource: public TParticleSource{
private:
	/**
	 * find potential minimum in source volume
	 */
	void FindPotentialMinimum();
protected:
	double MinPot; ///< minimal potential energy in source volume
	bool fPhaseSpaceWeighting; ///< Tells source to weight particle density according to available phase space.

	/**
	 * Produce random point in the source volume
	 *
	 * Has to be implemented by every derived class
	 *
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	virtual void RandomPointInSourceVolume(double &x, double &y, double &z) const = 0;
public:
	/**
	 * Constructor
	 *
	 * Create generic volume source
	 *
	 * @param ParticleName Name of particle that the source should create
	 * @param ActiveTime Duration for which the source shall be active
	 * @param PhaseSpaceWeighting If this is set true, the source will weight the particle density by available phase space
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TVolumeSource(std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TParticleSource(ParticleName, ActiveTime, mc, geometry, field), MinPot(std::numeric_limits<double>::infinity()), fPhaseSpaceWeighting(PhaseSpaceWeighting){

	}

	/**
	 * Create particle in source volume
	 *
	 * Particle density distribution can be weighted by available phase space
	 *
	 */
	TParticle* CreateParticle();
};

/**
 * Volume source generating points in a cuboid coordinate range
 */
class TCuboidVolumeSource: public TVolumeSource{
private:
	double xmin, xmax, ymin, ymax, zmin, zmax;

	/**
	 * Produce random point in the source volume
	 *
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	virtual void RandomPointInSourceVolume(double &x, double &y, double &z) const final{
		x = fmc->UniformDist(xmin, xmax);
		y = fmc->UniformDist(ymin, ymax);
		z = fmc->UniformDist(zmin, zmax);
	}
public:
	/**
	 * Constructor.
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param PhaseSpaceWeighting If this is set true, the source will weight the particle density by available phase space
	 * @param x_min Minimal radial coordinate range
	 * @param x_max Maximal radial coordinate range
	 * @param y_min Minimal azimuthal coordinate range
	 * @param y_max Maximal azimuthal coordinate range
	 * @param z_min Minimal axial coordinate range
	 * @param z_max Maximal axial coordinate range
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TCuboidVolumeSource(const std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field), xmin(x_min), xmax(x_max), ymin(y_min), ymax(y_max), zmin(z_min), zmax(z_max){

	}
};


/**
 * Volume source generating points in a cylindrical coordinate range
 */
class TCylindricalVolumeSource: public TVolumeSource{
private:
	double rmin, rmax, phimin, phimax, zmin, zmax;

	/**
	 * Produce random point in the source volume
	 *
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	void RandomPointInSourceVolume(double &x, double &y, double &z) const final{
		double r = fmc->LinearDist(rmin, rmax); // weighting because of the volume element and a r^2 probability outwards
		double phi_r = fmc->UniformDist(phimin,phimax);
		x = r*cos(phi_r);
		y = r*sin(phi_r);
		z = fmc->UniformDist(zmin,zmax);
	}
public:
	/**
	 * Constructor.
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param PhaseSpaceWeighting If this is set true, the source will weight the particle density by available phase space
	 * @param r_min Minimal radial coordinate range
	 * @param r_max Maximal radial coordinate range
	 * @param phi_min Minimal azimuthal coordinate range
	 * @param phi_max Maximal azimuthal coordinate range
	 * @param z_min Minimal axial coordinate range
	 * @param z_max Maximal axial coordinate range
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TCylindricalVolumeSource(const std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){

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
	bool InSourceVolume(const double x, const double y, const double z) const final{
		double r = sqrt(x*x + y*y);
		double phi = atan2(y, x);
		return r > rmin && r < rmax && phi > phimin && phi < phimax && z > zmin && z < zmax;
	}

	/**
	 * Check if a triangle is inside the cylidrical coordinate range
	 */
	bool InSourceVolume(const CTriangle &tri) const{
		for (int i = 0; i < 3; i++){
			if (!InSourceVolume(tri[i][0], tri[i][1], tri[i][2]))
				return false;
		}
		return true;
	}

public:
	/**
	 * Constructor.
	 *
	 * Selects all triangles from TGeometry whose vertices are inside the cylindrical coordinate range
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param E_normal Optional energy boost for particles starting from this surface
	 * @param r_min Minimal radial coordinate range
	 * @param r_max Maximal radial coordinate range
	 * @param phi_min Minimal azimuthal coordinate range
	 * @param phi_max Maximal azimuthal coordinate range
	 * @param z_min Minimal axial coordinate range
	 * @param z_max Maximal axial coordinate range
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TCylindricalSurfaceSource(const std::string ParticleName, double ActiveTime, double E_normal, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
			TSurfaceSource(ParticleName, ActiveTime, E_normal, mc, geometry, field), rmin(r_min), rmax(r_max), phimin(phi_min), phimax(phi_max), zmin(z_min), zmax(z_max){
		std::transform(	geometry.mesh.GetTrianglesBegin(),
						geometry.mesh.GetTrianglesEnd(),
						std::back_inserter(area_sum),
						[this](std::pair<CTriangle, double> tri){
							if (InSourceVolume(tri.first))
								return sqrt(tri.first.squared_area());
							return 0.;
						}
		);
		std::partial_sum(area_sum.begin(), area_sum.end(), area_sum.begin());
		std::cout << "Source area " << area_sum.back() << "m^2\n";
	}
};

/**
 * Volume source.
 *
 * Starting points are generated within an STL solid
 */
class TSTLVolumeSource: public TVolumeSource{
private:
	TTriangleMesh sourcevol; ///< internal AABB tree storing the STL solid

	/**
	 * Produce random point in the source volume
	 *
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	void RandomPointInSourceVolume(double &x, double &y, double &z) const final{
		CBox bbox = sourcevol.GetBoundingBox();
		do{
			x = fmc->UniformDist(bbox.xmin(), bbox.xmax()); // random point
			y = fmc->UniformDist(bbox.ymin(), bbox.ymax()); // random point
			z = fmc->UniformDist(bbox.zmin(), bbox.zmax()); // random point
		}while (!sourcevol.InSolid(x, y, z));
	}
public:
	/**
	 * Constructor.
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Duration for which the source shall be active
	 * @param PhaseSpaceWeighting If this is set true, the source will weight the particle density by available phase space
	 * @param sourcefile File from which the STL solid shall be read
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TSTLVolumeSource(const std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, std::string sourcefile, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TVolumeSource(ParticleName, ActiveTime, PhaseSpaceWeighting, mc, geometry, field){
		sourcevol.ReadFile(sourcefile, 0);
	}
};


/**
 * Surface source.
 *
 * Starting points are created on a surface of the experiment geometry whose triangles are all inside the given STL solid
 */
class TSTLSurfaceSource: public TSurfaceSource{
private:
	TTriangleMesh sourcevol;
	bool InSourceVolume(const double x, const double y, const double z) const{
		return sourcevol.InSolid(x, y, z);
	}
	bool InSourceVolume(const CTriangle &tri) const{
		for (int i = 0; i < 3; i++){
			if (!sourcevol.InSolid(tri[i]))
				return false;
		}
		return true;
	}
public:
	/**
	 * Constructor.
	 *
	 * Search for all triangles in geometry's mesh which are inside the STL solid given in sourcefile.
	 *
	 * @param ParticleName Name of particle type that the source should produce
	 * @param ActiveTime Time for which the source is active.
	 * @param sourcefile STL solid in which the surface should lie.
	 * @param E_normal Give particles starting at this source a velocity boost normal to the surface.
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TSTLSurfaceSource(const std::string ParticleName, double ActiveTime, std::string sourcefile, double E_normal, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field):
		TSurfaceSource(ParticleName, ActiveTime, E_normal, mc, geometry, field){
		sourcevol.ReadFile(sourcefile, 0);

		std::transform(	geometry.mesh.GetTrianglesBegin(),
						geometry.mesh.GetTrianglesEnd(),
						std::back_inserter(area_sum),
						[this](std::pair<CTriangle, double> tri){
							if (InSourceVolume(tri.first))
								return sqrt(tri.first.squared_area());
							return 0.;
						}
		);
		std::partial_sum(area_sum.begin(), area_sum.end(), area_sum.begin());
		std::cout << "Source area " << area_sum.back() << "m^2\n";
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
	std::string sourcemode; ///< volume/surface/customvol/customsurf
	TParticleSource *source; ///< TParticleSource contructed according to user chosen sourcemode


	/**
	 * Constructor, loads [SOURCE] section of configuration file, calculates TSource::Hmin_lfs and TSource::Hmin_hfs
	 *
	 * @param geometryconf TConfig struct containing [SOURCE] option map
	 * @param mc TMCGenerator that will be used to produce random numbers
	 * @param geometry TGeometry in which particles will be created
	 * @param field TFieldManager fields in which particles will be created
	 */
	TSource(TConfig &geometryconf, TMCGenerator &mc, TGeometry &geometry, TFieldManager *field);


	/**
	 * Destructor. Delete TSource::source.
	 */
	~TSource();


	/**
	 * Create new particle emerging from srouce
	 *
	 * "STLvolume": Random points with isotropic velocity distribution inside the bounding box of TSource::kdtree are produced until TSource::InSourceVolume is true for a point.
	 * "cylvolume": Random points with isotropic velocity deistribution in the coordinate range ([TSource::r_min..TSource::r_max]:[TSource::phi_min..TSource::phi_max]:[TSource::z_min..TSource::z_max]) are produced.
	 * "STLsurface/cylsurface": A random point on a random triangle from the TSource::sourcetris list is chosen (weighted by triangle area).
	 * 						Starting angles are cosine-distributed to the triangle normal and then angles and Ekin are modified according to TSource::E_normal.
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle();

};



#endif /* SOURCE_H_ */
