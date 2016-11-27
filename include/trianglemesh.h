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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> CKernel; ///< Geometric Kernel used for CGAL types
typedef CKernel::Segment_3 CSegment; ///< CGAL segment type
typedef CKernel::Point_3 CPoint; ///< CGAL point type
typedef CKernel::Triangle_3 CTriangle; ///< CGAL triangle type
typedef CKernel::Vector_3 CVector; ///< CGAL vector type
typedef CGAL::Bbox_3 CBox; ///< CGAL Box type
typedef std::vector< std::pair<CTriangle, int> > CTriangleList; ///< list used to store all triangles
typedef CTriangleList::const_iterator CIterator; ///< Iterator of triangle list. This is stored in the AABB tree.

/**
 * This primitive provides the conversion facilities between CTriangle and the types needed by CGAL AABB_tree
 */
struct CPrimitive {
public:
	typedef CIterator Id; ///< Type returned by CPrimitive::id().
	typedef CPoint Point; ///< Type returned by CPrimitive::reference_point().
	typedef CTriangle Datum; ///< Type returned by CPrimitive::datum().
private:
	Id m_pt; ///< this is what the AABB tree stores internally
public:
	/**
	 * Needed default constructor
	 */
	CPrimitive(): m_pt() {}

	/**
	 * Constructor
	 *
	 * this constructor is the one that receives the iterators from the
	 * iterator range given as input to the AABB_tree
	 */
	CPrimitive(CIterator it): m_pt(it) {}

	/**
	 * Return internal iterator.
	 */
	const Id& id() const { return m_pt; }

	/**
	 * Return the CGAL Primitive
	 */
	Datum datum() const{ return m_pt->first; }
	/**
	 * Return a reference point.
	 *
	 * returns a reference point which must be on the primitive
	 */
	Point reference_point() const{ return m_pt->first.vertex(0); }
};



typedef CGAL::AABB_traits<CKernel, CPrimitive> CTraits; ///< CGAL triangle traits type
typedef CGAL::AABB_tree<CTraits> CTree; ///< CGAL AABB tree type containing CPrimitives
#if CGAL_VERSION_NR<1040301000
	typedef boost::optional< CTree::Object_and_primitive_id > CIntersection; ///< CGAL 4.2 or older segment-triangle intersection type
#else
	typedef boost::optional< CTree::Intersection_and_primitive_id<CSegment>::Type > CIntersection; ///< CGAL 4.3 segment-triangle intersection type
#endif



/**
 * Structure returned by TTriangleMesh::Collision.
 */
struct TCollision{
	double s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
	double normal[3]; ///< normal (length = 1) of intersected surface
	int sldindex; ///< index of solid to which the intersected surface belongs
	double distnormal; ///< distance between start- and endpoint of colliding segment, projected onto normal direction

	/**
	 * Create TCollision object
	 *
	 * @param segment Segment that collided with mesh
	 * @param tri Triangle the segment collided with
	 * @param point Collision point
	 */
	TCollision(const CSegment &segment, const CIterator &tri, const CPoint &point){
		s = sqrt((point - segment.start()).squared_length()/segment.squared_length());
		sldindex = tri->second;
		CVector n = tri->first.supporting_plane().orthogonal_vector();
		n = n/sqrt(n.squared_length());
		normal[0] = n[0];
		normal[1] = n[1];
		normal[2] = n[2];
		distnormal = segment.to_vector()*n;
	}

	/**
	 * Overloaded operator, needed for sorting
	 */
	inline bool operator < (const TCollision c) const {
		if (s == c.s)
			return sldindex > c.sldindex;
		else
			return s < c.s;
	};
};



/**
 * Class to hold your STL geometry and do intersection tests.
 */
class TTriangleMesh{
private:
	CTriangleList triangles; ///< list of triangles
	CTree tree; ///< AABB tree

public:
	/**
	 * Read STL-file.
	 *
	 * @param filename Filename of STL file
	 * @param sldindex Index of solid properties assigned to this STL file
	 *
	 * @return Returns name of mesh in file
	 */
	std::string ReadFile(const std::string &filename, const int sldindex);

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
	 * Test if point is inside the mesh
	 *
	 * @param x X coordinate of point
	 * @param y Y coordinate of point
	 * @param z Z coordinate of point
	 *
	 * @return Returns true if point inside the mesh
	 */
	bool InSolid(const double x, const double y, const double z) const{
		std::vector<TCollision> colls = Collision(std::vector<double>{x, y, z}, std::vector<double>{x, y, tree.bbox().zmin() - 1});
		return (colls.size() & 1); // return true if collided with surface and number of collisions is odd
	}

	/**
	  * Test if point is inside the mesh
	  *
	  * @param p Point
	  *
	  * @return Returns true if point is inside the mesh
	  */
	bool InSolid(const CPoint &p) const{
		return InSolid(p[0], p[1], p[2]);
	}

	/**
	 * Return iterator to triangle list
	 */
	CIterator GetTrianglesBegin() const{
		return triangles.begin();
	}

	/**
	 * Return iterator to triangle list
	 */
	CIterator GetTrianglesEnd() const{
		return triangles.end();
	}

	/**
	 * Return bounding box of mesh
	 */
	CBox GetBoundingBox() const{
		return tree.bbox();
	}
};

#endif // TRIANGLEMESH_H_
