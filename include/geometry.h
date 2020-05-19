/**
 * \file
 * Contains class to include experiment geometry.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <vector>
#include <map>

#include "trianglemesh.h"
#include "config.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

/// Struct to store material properties (read from geometry.in, right now only for neutrons)
struct material{
	std::string name; ///< Material name
	double FermiReal; ///< Real part of Fermi potential
	double FermiImag; ///< Imaginary part of Fermi potential
	double DiffProb; ///< Diffuse reflection probability
	double SpinflipProb; ///< Probability for spin flip on reflection
	double RMSRoughness; ///< RMS roughness of surface, for MicroRoughness model reflections
	double CorrelLength; ///< Correlation length of surface roughness, for MicroRoughness model reflections
    double InternalBField; ///< Internal magnetic field, for magnetized materials
    double ModifiedLambertProb; ///< Probability of diffuse reflection according to Modified Lambert model
    double LossPerBounce; ///< Probability of loss when hitting the wall (independent of neutron energy)
};

///Read material properties, except name, from input stream
std::istream& operator>>(std::istream &str, material &mat);
///Write material properties into output stream
std::ostream& operator<<(std::ostream &str, const material &mat);


/// Struct to store solid information (read from geometry.in)
struct solid{
	boost::filesystem::path filename; ///< name of file containing STL mesh
	std::string name; ///< name of solid
	material mat; ///< material of solid
	unsigned ID; ///< ID of solid
	std::vector<std::pair<double, double> > ignoretimes; ///< pairs of times, between which the solid should be ignored

	/**
	 * Comparison operator used to sort solids by priority (descending)
	 *
	 * @param s solid struct compared to this solid
	 *
	 * @return Returns true if this solid's ID is larger than the other one's
	 */
	bool operator< (const solid s) const { return ID > s.ID; };

	/**
	 * Check if solid is ignored at a certain time
	 * 
	 * @param t Time
	 * 
	 * @return Returns true if solid is ignored at time t
	 */
	bool is_ignored(const double t) const{
	    return std::any_of(ignoretimes.begin(), ignoretimes.end(),
	            [&t](const std::pair<double, double> &its){ return t >= its.first && t < its.second; }
	            ); // check if collision time lies between any pair of ignore times};
	}
};

///Read solid properties, except ID, from input stream
std::istream& operator>>(std::istream &str, solid &model);
///Write solid properties into output stream
std::ostream& operator<<(std::ostream &str, const solid &sld);


/**
 * Class to include experiment geometry.
 *
 * Loads solids and materials from geometry.in, maintains solids list, checks for collisions.
 */
struct TGeometry{
	private:
		std::vector<solid> solids; ///< solids list
	public:
		TTriangleMesh mesh; ///< kd-tree structure containing triangle meshes from STL-files
		solid defaultsolid; ///< "vacuum", this solid's properties are used when the particle is not inside any other solid
		
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
		bool CheckSegment(const double y1[3], const double y2[3]) const{
			return mesh.InBoundingBox(CSegment(CPoint(y1[0], y1[1], y1[2]), CPoint(y2[0], y2[1], y2[2])));
		};
		

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
		bool GetCollisions(const double x1, const double p1[3], const double x2, const double p2[3], std::multimap<TCollision, bool> &colls) const;
		
			
		/**
		 * Get solids in which the point p lies
		 *
		 * @param t Time
		 * @param p Point to test
		 *
		 * @return Map of solids in which the point is inside paired with information if it was ignored or not
		 */
		std::vector<std::pair<solid, bool> > GetSolids(const double t, const double p[3]) const;


		/**
		 * Get solid with highest priority in which the point p lies
		 *
		 * @param t Time
		 * @param p Point to test
		 *
		 * @return Returns solid with highest priority, that was not ignored at time t
		 */
		solid GetSolid(const double t, const double p[3]) const;


		/**
		 * Get solid with given ID
		 * 
		 * @param ID ID
		 * 
		 * @return Returns solid with given ID
		 */
		 solid GetSolid(const unsigned ID) const;
};

#endif /*GEOMETRY_H_*/
