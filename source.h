/**
 * \file
 * Contains base class TParticleSource and several dervied particle source classes.
 * Class TSource creates one of these according to user input.
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include <string>

#include "particle.h"
#include "mc.h"

using namespace std;

/**
 * Virtual base class for all particle sources
 */
class TParticleSource{
protected:
	double fActiveTime; ///< Duration for which the source will be active
	const string fParticleName; ///< Name of particle that the source should create
public:
	int ParticleCounter; ///< Count number of particles created by source
	/**
	 * Constructor, should be called by every derived class
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 */
	TParticleSource(const string ParticleName, double ActiveTime);

	/**
	 * Destructor
	 */
	virtual ~TParticleSource();

	/**
	 * Basic creation of a new particle. Maps particle name to the corresponding class
	 *
	 * @param mc random number generator
	 * @param t Starting time
	 * @param x x coordinate of creation point
	 * @param y y coordinate of creation point
	 * @param z z coordinate of creation point
	 * @param E Initial kinetic energy
	 * @param phi Azimuthal angle of initial velocity vector
	 * @param theta Polar angle of initial velocity vector
	 * @param polarisation Initial polarisation of particle (-1, 0, 1)
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle(TMCGenerator &mc, double t, double x, double y, double z, double E, double phi, double theta, int polarisation, TGeometry &geometry, TFieldManager *field);


	/**
	 * Virtual routine that has to be implemented by every derived source class
	 *
	 * @param mc Random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional field (can be NULL)
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	virtual TParticle* CreateParticle(TMCGenerator &mc, TGeometry &geometry, TFieldManager *field) = 0;
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
	TSurfaceSource(const string ParticleName, double ActiveTime, double E_normal);

	/**
	 * Create new particle on surface
	 *
	 * @param mc random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle(TMCGenerator &mc, TGeometry &geometry, TFieldManager *field);
};


class TVolumeSource: public TParticleSource{
protected:
	bool fPhaseSpaceWeighting;

	/**
	 * Produce random point in the source volume
	 *
	 * Has to be implemented by every derived class
	 *
	 * @param mc Random number generator
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	virtual void RandomPointInSourceVolume(TMCGenerator &mc, double &x, double &y, double &z) = 0;
public:
	/**
	 * Constructor
	 *
	 * Create generic volume source
	 *
	 * @param ParticleName Name of particle that the source should create
	 * @param ActiveTime Duration for which the source shall be active
	 * @param PhaseSpaceWeighting If this is set true, the source will weight the particle density by available phase space
	 */
	TVolumeSource(std::string ParticleName, double ActiveTime, bool PhaseSpaceWeighting);

	/**
	 * Create particle in source volume
	 *
	 * Particle density distribution can be weighted by available phase space
	 *
	 * @param mc Random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 */
	TParticle* CreateParticle(TMCGenerator &mc, TGeometry &geometry, TFieldManager *field);
};

/**
 * Volume source generating points in a cylindrical coordinate range
 */
class TCylindricalVolumeSource: public TVolumeSource{
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
	TCylindricalVolumeSource(const string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max);


	/**
	 * Produce random point in the source volume
	 *
	 * @param mc Random number generator
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	virtual void RandomPointInSourceVolume(TMCGenerator &mc, double &x, double &y, double &z);
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
	bool InSourceVolume(CPoint p);

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
	TCylindricalSurfaceSource(const string ParticleName, double ActiveTime, TGeometry &geometry, double E_normal, double r_min, double r_max, double phi_min, double phi_max, double z_min, double z_max);
};

/**
 * Volume source.
 *
 * Starting points are generated within an STL solid
 */
class TSTLVolumeSource: public TVolumeSource{
private:
	TTriangleMesh kdtree; ///< internal AABB tree storing the STL solid
public:
	/**
	 * Constructor.
	 *
	 * @param ActiveTime Duration for which the source shall be active
	 * @param sourcefile File from which the STL solid shall be read
	 */
	TSTLVolumeSource(const string ParticleName, double ActiveTime, bool PhaseSpaceWeighting, string sourcefile);


	/**
	 * Produce random point in the source volume
	 *
	 * @param mc Random number generator
	 * @param x Returns x coordinate
	 * @param y Returns y coordinate
	 * @param z Returns z coordinate
	 */
	virtual void RandomPointInSourceVolume(TMCGenerator &mc, double &x, double &y, double &z);
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
	TSTLSurfaceSource(const string ParticleName, double ActiveTime, TGeometry &geometry, string sourcefile, double E_normal);
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
	TSource(TConfig &geometryconf, TGeometry &geom, TFieldManager &field);


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
	 * @param mc random number generator
	 * @param geometry Experiment geometry
	 * @param field Optional fields (can be NULL)
	 *
	 * @return Returns newly created particle, memory has to be freed by user
	 */
	TParticle* CreateParticle(TMCGenerator &mc, TGeometry &geometry, TFieldManager *field);

};



#endif /* SOURCE_H_ */
