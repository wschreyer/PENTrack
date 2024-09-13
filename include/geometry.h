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
 * Structure returned by TTriangleMesh::Collision.
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


template<class Material, class Path, class Mesh>
struct TBody{
	unsigned ID;
	boost::filesystem::path filename;
	Material material;
	std::optional<Path> path;
	Mesh mesh;
};

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
 * Class to include experiment geometry.
 *
 * Loads solids and materials from geometry.in, maintains solids list, checks for collisions.
 */
struct TGeometry{
	std::vector<TBody<material, TPathInterpolator<std::vector<double> >, TMesh<double> > > solids; ///< solids list
	material defaultMaterial; ///< "vacuum", this solid's properties are used when the particle is not inside any other solid
		
	/**
	 * Constructor, reads geometry configuration file, loads triangle meshes.
	 *
	 * @param geometryin TConfig struct containing MATERIALS and GEOMETRY config section
	 */
	TGeometry(TConfig &geometryin);


	/**
	 * Check if segment is intersecting with geometry bounding box.
	 *
	 * @param y1 Position vector of segment start
	 * @param y2 Position vector of segment end
	 *
	 * @return Returns true if segment is intersecting bounding box
	 */
	bool CheckSegment(const double t1, const std::array<double, 3> &y1, const double t2, const std::array<double, 3> &y2) const;
	

	/**
	 * Checks if line segment p1->p2 collides with a surface.
	 *
	 * Calls KDTree::Collision to check for collisions and removes all collisions
	 * which should be ignored (given by ignore times in geometry configuration file).
	 *
	 * @param x1 Start time of line segment
	 * @param p1 Start point of line segment
	 * @param x2 End time of line segment
	 * @param p2 End point of line segment
	 * @param colls List of collisions, paired with bool indicator it it should be ignored
	 *
	 * @return Returns true if line segment collides with a surface
	 */
	bool GetCollisions(const double x1, const std::array<double, 3> &p1, const double x2, const std::array<double, 3> &p2, std::multimap<TCollision, bool> &colls) const;
	
		
	/**
	 * Get solids in which the point p lies
	 *
	 * @param t Time
	 * @param p Point to test
	 *
	 * @return Map of solids in which the point is inside paired with information if it was ignored or not
	 */
	std::vector<std::pair<solid, bool> > GetSolids(const double t, const std::array<double, 3> &p) const;


	/**
	 * Get solid with highest priority in which the point p lies
	 *
	 * @param t Time
	 * @param p Point to test
	 *
	 * @return Returns solid with highest priority, that was not ignored at time t
	 */
	solid GetSolid(const double t, const std::array<double, 3> &p) const{
		auto solids = GetSolids(t, p);
		return solids.back().first;
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

	template <class RandomGenerator>
	std::tuple<std::array<double, 3>, std::array<double, 3>, unsigned int> RandomPointOnSurface(RandomGenerator &rand, const std::optional<std::pair<std::array<double, 3>, std::array<double, 3>>> &boundingBox = std::nullopt){
		std::vector<double> weights;
		for (auto &body: solids){
			if (not body.path and (not boundingBox or overlappingBoundingBoxes(body.mesh, *boundingBox))){
				weights.push_back(body.mesh.area);
			}
			else{
				weights.push_back(0.);
			}
		}
		std::discrete_distribution<size_t> d(weights.begin(), weights.end());
		auto i = d(rand);
		auto [p, n] = randomPointOnSurface(solids[i].mesh, rand);
		return std::make_tuple(p, n, solids[i].ID);
	}
};



#endif /*GEOMETRY_H_*/
