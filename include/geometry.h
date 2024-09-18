/**
 * \file
 * Contains class to include experiment geometry.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "config.h"
#include "pathing.h"
#include "mesh.h"
#include "globals.h"


static const double REFLECT_TOLERANCE = 1e-8;  ///< max distance of reflection point to actual surface collision point

/// Struct to store material properties (read from geometry.in, right now only for neutrons)
struct material{
	std::string name;
	double FermiReal; ///< Real part of Fermi potential
	double FermiImag; ///< Imaginary part of Fermi potential
	double DiffProb; ///< Diffuse reflection probability
	double SpinflipProb; ///< Probability for spin flip on reflection
	double RMSRoughness; ///< RMS roughness of surface, for MicroRoughness model reflections
	double CorrelLength; ///< Correlation length of surface roughness, for MicroRoughness model reflections
	double InternalBField; ///< Internal magnetic field, for magnetized materials
	double ModifiedLambertProb; ///< Probability of diffuse reflection according to Modified Lambert model
	double LossPerBounce; ///< Probability of loss when hitting the wall (independent of neutron energy)
	double microfacetDistributionWidth; ///< Width of Beckmann microfacet distribution
	double microfacetDistributionWidthExponent; ///< Exponent for wavelength-dependent microfacet distribution = microfacetDistributionWidth * (wavelength[nm]/cos(incident angle))^microfacetDistributionWidthExponent)
};

///Read material properties, except name, from input stream
std::istream& operator>>(std::istream &str, material &mat);


///Write material properties into output stream
std::ostream& operator<<(std::ostream &str, const material &mat);



struct solid{
	unsigned ID;
	material mat;
};

/**
 * Structure returned by TGeometry::GetCollisions.
 */
struct TCollision{
	double s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
	std::array<double, 3> normal; ///< normal (length = 1) of intersected surface
	std::array<double, 3> surfaceVelocity; ///< velocity vector of surface
	unsigned ID; ///< ID of solid the intersected surface belongs to

	/**
	 * Overloaded operator, needed for sorting
	 * 
	 * Ascending distance along segment, descending ID if distance equal
	 */
	inline bool operator < (const TCollision c) const {
		if (s == c.s)
			return ID > c.ID;
		else
			return s < c.s;
	};
};


/**
 * Rigid body described by a surface mesh, with material properties and an optional path along which it is moving
 */
template<class Material, class Path, class Mesh>
struct TBody{
	unsigned ID; ///< Unique ID number
	boost::filesystem::path filename; ///< File the surface mesh is read from
	Material material; ///< Material properties
	std::optional<Path> path; ///< Optional path along which the mesh is moving
	Mesh mesh; ///< Class describing the surface mesh
};


/**
 * Construct a TBody class from a parameter string, a list of materials, and a list of path descriptions.
 * 
 * @param ID Unique ID number
 * @param parameters Parameter string containing filename to read surface mesh from, material name, and path name
 * @param materials List of materials
 * @param paths List of path descriptions
 */
template<class Material, class Path>
TBody<Material, Path, TMesh<double> > constructBody(const unsigned ID, const std::string &parameters, const std::map<std::string, Material> &materials, const std::map<std::string, Path> &paths){
	boost::filesystem::path filename;
	std::string materialName;
	std::istringstream ss(parameters);
	ss >> filename >> materialName;
	if (not ss){
		throw std::runtime_error("Invalid parameters "" + parameters + "" for solid " + ID);
	}
	std::string pathName;
	ss >> pathName;
	std::optional<Path> path;
	if (ss){
		path = paths.at(pathName);
	}
	return {ID, filename, materials.at(materialName), std::move(path), constructMesh(boost::filesystem::absolute(filename, configpath.parent_path()).native())};
}

/**
 * Class describing experiment geometry.
 *
 * Contains a list of bodies, provides intersection tests and random point generators
 */
struct TGeometry{
	std::vector<TBody<material, TPathInterpolator<std::vector<double> >, TMesh<double> > > solids; ///< List of bodies
	material defaultMaterial; ///< "vacuum", these material properties are provided if a point is not inside any other solid
		
	/**
	 * Constructor, reads geometry configuration parameters to generate list of TBody classes.
	 *
	 * @param geometryin TConfig struct containing geometry configuration paramters in [GEOMETRY], [MATERIALS], and [TRANSFORMATIONS] sections
	 */
	TGeometry(TConfig &geometryin);


	/**
	 * Check if segment is intersecting with geometry bounding box.
	 *
	 * @param tstart Time of segment start
	 * @param start Position vector of segment start
	 * @param tend Time of segment end
	 * @param end Position vector of segment end
	 *
	 * @return Returns true if segment is intersecting bounding box
	 */
	bool CheckSegment(const double tstart, const std::array<double, 3> &start, const double tend, const std::array<double, 3> &end) const;
	

	/**
	 * Checks if line segment p1->p2 collides with a surface.
	 *
	 * Calls KDTree::Collision to check for collisions and removes all collisions
	 * which should be ignored (given by ignore times in geometry configuration file).
	 *
	 * @param tstart Start time of line segment
	 * @param start Start point of line segment
	 * @param tend End time of line segment
	 * @param end End point of line segment
	 * @param colls List of collisions
	 *
	 * @return Returns true if line segment collides with a surface
	 */
	bool GetCollisions(const double tstart, const std::array<double, 3> &start, const double tend, const std::array<double, 3> &end, std::vector<TCollision> &colls) const;
	
		
	/**
	 * Get solids in which the point p lies
	 *
	 * @param t Time
	 * @param p Point to test
	 *
	 * @return List of solids in which the point is inside
	 */
	std::vector<solid> GetSolids(const double t, const std::array<double, 3> &p) const;


	/**
	 * Get solid with highest priority in which the point p lies
	 *
	 * @param t Time
	 * @param p Point to test
	 *
	 * @return Returns solid with highest priority
	 */
	solid GetSolid(const double t, const std::array<double, 3> &p) const{
		auto solids = GetSolids(t, p);
		return solids.back();
	}


	/**
	 * Get solid with given ID
	 * 
	 * @param ID ID
	 * 
	 * @return Returns solid with given ID
	 */
	solid GetSolid(const unsigned ID) const{
		for (auto &body: solids){
			if (body.ID == ID){
				return {ID, body.material};
			}
		}
		throw std::runtime_error((boost::format("Could not find solid with ID %s") % ID).str());
	}


	/**
	 * Generate random points homogeneously distributed on all surfaces of the geometry, optionally within a bounding box
	 * 
	 * @param rand Random number generator class
	 * @param boundingBox Optional bounding box the point should be contained in
	 * 
	 * @returns Tuple of two vectors (point on a surface and surface normal at that point) and ID number of the surface
	 */
	template <class RandomGenerator>
	std::tuple<std::array<double, 3>, std::array<double, 3>, unsigned int> RandomPointOnSurface(RandomGenerator &rand, const std::pair<std::array<double, 3>, std::array<double, 3>> &boundingBox){
		std::vector<double> weights;
		for (auto &body: solids){
			if (overlappingBoundingBoxes(body.mesh, boundingBox)){
				if (not body.path){
					weights.push_back(body.mesh.area);
				}
				else{
					//std::cout << "Generating points on moving surfaces is not implemented!\n";
					weights.push_back(0.);
				}
			}
			else{
				weights.push_back(0.);
			}
		}
		std::discrete_distribution<size_t> d(weights.begin(), weights.end());
		while (true){
			auto i = d(rand);
			auto [p, n] = randomPointOnSurface(solids[i].mesh, rand);
			if (p[0] > boundingBox.first[0] and p[1] > boundingBox.first[1] and p[2] > boundingBox.first[2] and p[0] < boundingBox.second[0] and p[1] < boundingBox.second[1] and p[2] < boundingBox.second[2]){
				return std::make_tuple(p, n, solids[i].ID);
			}
		}
	}


	/**
	 * Generate random points homogeneously distributed on all surfaces of the geometry
	 * 
	 * @param rand Random number generator class
	 * 
	 * @returns Tuple of two vectors (point on a surface and surface normal at that point) and ID number of the surface
	 */
	template <class RandomGenerator>
	std::tuple<std::array<double, 3>, std::array<double, 3>, unsigned int> RandomPointOnSurface(double time, RandomGenerator &rand){
		std::vector<double> weights;
		for (auto &body: solids){
			weights.push_back(body.mesh.area);
		}
		std::discrete_distribution<size_t> d(weights.begin(), weights.end());
		auto i = d(rand);
		auto [p, n] = solids[i].path ? randomPointOnMovingSurface(solids[i].mesh, time, rand, *solids[i].path) : randomPointOnSurface(solids[i].mesh, rand);
		return std::make_tuple(p, n, solids[i].ID);
	}

};



#endif /*GEOMETRY_H_*/
