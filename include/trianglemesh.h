/**
 * \file
 * This algorithm uses the CGAL AABB_tree structure to search
 * for collisions with a surface consisting of a list
 * of triangles.
 * Initially, the triangles are read from a set of STL-files
 * (http://www.ennex.com/~fabbers/StL.asp)	via
 * ReadFile(filename,surfacetype) and stored in the AABB_tree
 * via Init().
 * You can define a surfacetype for each file which is
 * returned on collision tests to identify different surfaces
 * during runtime.
 * During runtime segments point1->point2 can be checked for
 * intersection with the surface via
 * Collision(point1,point2,list of TCollision). Collision returns
 * true if an intersection occurred and gives the parametric
 * coordinate s of the intersection point (I=p1+s*(p2-p1)),
 * the normal n and the surfacetype of the intersected surface.
 *
 */

#ifndef TRIANGLEMESH_H_
#define TRIANGLEMESH_H_

#include <vector>
#include <memory>
#include <random>

#include <algorithm>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

static const double REFLECT_TOLERANCE = 1e-8;  ///< max distance of reflection point to actual surface collision point

typedef CGAL::Simple_cartesian<double> CKernel; ///< Geometric Kernel used for CGAL types
typedef CKernel::Segment_3 CSegment; ///< CGAL segment type
typedef CKernel::Point_3 CPoint; ///< CGAL point type
typedef CKernel::Vector_3 CVector; ///< CGAL vector type
typedef CKernel::Iso_cuboid_3 CCuboid; ///< CGAL cuboid type

typedef CGAL::Surface_mesh<CPoint> CMesh; ///< CGAL triangle mesh type
typedef CGAL::AABB_face_graph_triangle_primitive<CMesh> CPrimitive; ///< CGAL triangle type contained in AABB tree
typedef CGAL::AABB_traits<CKernel, CPrimitive> CTraits; ///< CGAL triangle traits type
typedef CGAL::AABB_tree<CTraits> CTree; ///< CGAL AABB tree type containing CPrimitives
typedef boost::optional< CTree::Intersection_and_primitive_id<CSegment>::Type > CIntersection; ///< CGAL segment-triangle intersection type


/**
 * Structure returned by TTriangleMesh::Collision.
 */
struct TCollision{
	double s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
	double normal[3]; ///< normal (length = 1) of intersected surface
	unsigned ID; ///< ID of solid the intersected surface belongs to
	double distnormal; ///< distance between start- and endpoint of colliding segment, projected onto normal direction

	/**
	 * Create TCollision object
	 *
	 * @param segment Segment that collided with mesh
	 * @param n Normal vector of hit surface
	 * @param point Collision point
	 * @param aID ID of hit surface
	 */
	TCollision(const CSegment &segment, const CVector &n, const CPoint &point, const unsigned aID){
      s = /*std::min(1., std::max(0.,*/ (point - segment.start())*segment.to_vector()/segment.squared_length()/*))*/;
      ID = aID;
      normal[0] = n[0];
      normal[1] = n[1];
      normal[2] = n[2];
      distnormal = segment.to_vector()*n;
    };

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
 * Class to hold your STL geometry and do intersection tests.
 */
class TTriangleMesh{
private:
	/**
	 * Class containing triangle mesh and AABB tree for each loaded StL file
	 */
    struct CTriangleMesh{
        std::unique_ptr<CMesh> mesh; ///< Triangle mesh
        std::unique_ptr<CTree> tree; ///< Axis-aligned bounding-box tree for fast intersection search
        int ID; ///< unique ID for each StL file
        std::discrete_distribution<size_t> triangle_sampler; ///< Probability distribution to randomly sample triangles from mesh weighted by their areas.
    };
	std::vector<CTriangleMesh> meshes; ///< List of triangle meshes from all loaded StL files
	std::discrete_distribution<size_t> mesh_sampler; ///< Probability distribution to randomly sample meshes weighted by their areas

public:
	/**
	 * Read STL-file.
	 *
	 * @param filename Filename of STL file
	 * @param ID ID of solid assigned to this STL file
	 *
	 * @return Returns name of mesh in file
	 */
	std::string ReadFile(const std::string &filename, const int ID);

	/**
	 * Test line segment p1->p2 for collision with all triangles in previously read files.
	 *
	 * @param p1 Line start point
	 * @param p2 Line end point
	 *
	 * @return Returns vector containing collisions
	 */
	std::vector<TCollision> Collision(const std::vector<double> &p1, const std::vector<double> &p2) const;

	/**
	 * Test if point is inside the mesh
	 *
	 * @param p Point
	 *
	 * @return Returns true if point is inside the mesh
	 */
	template<typename T> bool InSolid(const T p[3]) const{
		return InSolid(p[0], p[1], p[2]);
	}

	/**
	 * Get overall bounding box containing all meshes
	 * 
	 * @return Overall bounding box.
	 */
	CCuboid GetBoundingBox() const{
	    std::vector<CCuboid> b;
	    std::transform(meshes.begin(), meshes.end(), std::back_inserter(b), [](const CTriangleMesh &m){ return m.tree->bbox(); });
	    return CGAL::bbox_3(b.begin(), b.end());
	}

	/**
	 * Test if point is inside the mesh
	 *
	 * @param x X coordinate of point
	 * @param y Y coordinate of point
	 * @param z Z coordinate of point
	 *
	 * @return Returns true if point inside the mesh
	 */
	bool InSolid(const double x, const double y, const double z) const;
	/**
	  * Test if point is inside the mesh
	  *
	  * @param p Point
	  *
	  * @return Returns true if point is inside the mesh
	  */
	template<class Point> bool InSolid(Point p) const{
		return InSolid(p[0], p[1], p[2]);
	}

	/**
	 * Return list of solids the point is inside of
	 * @param p Point
	 * @return List of solid IDs
	 */
    template<class Point> std::vector<unsigned> GetSolids(Point p) const{
        std::vector<unsigned> solids;
        for (auto &m: meshes){
            if (m.tree->number_of_intersected_primitives(CKernel::Ray_3(CPoint(p[0],p[1],p[2]), CVector(0., 0., 1.))) % 2 != 0)
                solids.push_back(m.ID);
        }
        return solids;
    }

	/**
	 * Check if point is contained in bounding box
	 * 
	 * @param p Point
	 * 
	 * @return Returns true if point is contained in bounding box
	 */
	template<class Object> bool InBoundingBox(Object p) const{
        return std::any_of(meshes.begin(), meshes.end(), [&p](const CTriangleMesh &mesh){ return CGAL::do_intersect(p, mesh.tree->bbox()); });
	}

	/**
	 * Return random point on surface
	 * 
	 * @param p Returned point
	 * @param n Returned normal vector of surface at point
	 * @param ID Returned ID of surface at point
	 * @param rand Random number generator
	 * @param bbox Bounding box that point should be contained in
	 */
	template<class Point, class Vector, class RandomGenerator, class BoundingBox> void RandomPointOnSurface(Point &p, Vector &n, unsigned &ID, RandomGenerator &rand, BoundingBox bbox){
        size_t meshidx;
        do{
            meshidx = mesh_sampler(rand);
        }while (not CGAL::do_intersect(meshes[meshidx].tree->bbox(), bbox));
        ID = meshes[meshidx].ID;
        CMesh::Face_index faceidx(meshes[meshidx].triangle_sampler(rand));
        std::vector<CPoint> vertices;
        for (auto v: meshes[meshidx].mesh->vertices_around_face(meshes[meshidx].mesh->halfedge(faceidx))) {
            vertices.push_back(meshes[meshidx].mesh->point(v));
        }
        std::uniform_real_distribution<double> unidist(0, 1);
        double a = unidist(rand); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
        double b = unidist(rand);
        if (a+b > 1){
            a = 1 - a;
            b = 1 - b;
        }
        CPoint pp = vertices[0] + a*(vertices[1] - vertices[0]) + b*(vertices[2] - vertices[0]);
        CVector nv = CGAL::Polygon_mesh_processing::compute_face_normal(faceidx, *meshes[meshidx].mesh);
        p = {pp.x(), pp.y(), pp.z()};
        n = {nv.x(), nv.y(), nv.z()};
	}

	/**
	 * Return random point in volume bounded by mesh
	 * 
	 * @param rand Random number generator
	 * 
	 * @return Point
	 */
	template<class RandomGenerator> std::array<double, 3> RandomPointInVolume(RandomGenerator &rand) const{
        std::array<double, 3> p;
        do{
            p = RandomPointInBoundingBox(rand);
        }while (!InSolid(p));
        return p;
    }

	/**
	 * Return random point in bounding box
	 * 
	 * @param rand Random number generator
	 * 
	 * @return Point
	 */
    template<class RandomGenerator> std::array<double, 3> RandomPointInBoundingBox(RandomGenerator &rand) const{
        std::vector<double> bvols;
        std::transform(meshes.begin(), meshes.end(), std::back_inserter(bvols), [](const CTriangleMesh &mesh){ return CCuboid(mesh.tree->bbox()).volume(); });
        std::discrete_distribution<unsigned> dist(bvols.begin(), bvols.end());
        CCuboid bbox = meshes[dist(rand)].tree->bbox();
        std::uniform_real_distribution<double> unidist(0, 1);
        return {bbox.xmin() + unidist(rand)*(bbox.xmax() - bbox.xmin()),
                bbox.ymin() + unidist(rand)*(bbox.ymax() - bbox.ymin()),
                bbox.zmin() + unidist(rand)*(bbox.zmax() - bbox.zmin())};
    }
};

#endif // TRIANGLEMESH_H_
