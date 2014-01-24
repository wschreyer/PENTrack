/**
 * \file
 * This algorithm uses a kd-tree structure to search
 * for collisions with a surface consisting of a list
 * of triangles.
 * Initially, the triangles are read from a set of STL-files
 * (http://www.ennex.com/~fabbers/StL.asp)	via
 * ReadFile(filename,surfacetype) and stored in the kd-tree
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

#ifndef KDTREE_H_
#define KDTREE_H_

#include <vector>
#include <set>
#include <list>

#define USE_CGAL

#ifdef USE_CGAL

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::Segment_3 CSegment;
typedef K::Point_3 CPoint;
typedef K::Triangle_3 CTriangle;
class TTriangle: public K::Triangle_3{
public:
	int ID;
	TTriangle(K::Point_3 a, K::Point_3 b, K::Point_3 c, int aID): K::Triangle_3(a, b, c), ID(aID) { };
};
typedef std::list<TTriangle>::iterator CIterator;
typedef CGAL::AABB_triangle_primitive<K,CIterator> CPrimitive;
typedef CGAL::AABB_traits<K, CPrimitive> CTraits;
typedef CGAL::AABB_tree<CTraits> CTree;

#else

/**
 * Triangle class.
 *
 * Stores vertices, normal vector and an ID.
 */
struct Triangle{
    const float *vertex[3]; ///< the three vertices of the triangle (in counterclockwise order according to STL standard)
    long double normal[3]; ///< normal vector of triangle plane (length = area)
    unsigned ID;          ///< user defined surface number

    /**
     * Constructor.
     *
     * Calculates ::normal
     *
     * @param vertices Array of three vertices with three components
     */
    Triangle(const float *vertices[3]);

    /**
     * Does the line segment p1 -> p2 intersect the triangle?
     *
     * http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#Segment-Triangle
     *
     * @param p1 Start point of line segment
     * @param p2 End point of line segment
     * @param s Length fraction of the intersection point (P = p1 + s*(p2 - p1))
     *
     * @return Returns true if they intersect
     */
    bool intersect(const long double p1[3], const long double p2[3], long double &s);
};
#endif

/**
 * Structure that is returned by KDTree::Collision.
 */
struct TCollision{
	long double s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
	long double normal[3]; ///< normal (length = 1) of intersected surface
	unsigned ID; ///< ID of intersected surface
	inline bool operator < (const TCollision c) const { return s < c.s; }; ///< overloaded operator, needed for sorting
};

/**
 * General kd-tree class.
 *
 * Encapsulates all the nodes.
 */
class KDTree{
#ifndef USE_CGAL
	private:
    	/**
    	 * Triangle vertex class.
    	 *
    	 * Each vertex is only stored once in a set to save memory.
    	 */
        struct TVertex{
            float vertex[3]; ///< coordinates of the vertex

            /**
             * Comparison operator.
             *
             * Needed to keep only one unique vertex in the set.
             *
             * @param v Vertex to which shall be compared
             */
            inline bool operator == (const TVertex v) const { 
            	return  std::abs(vertex[0] - v.vertex[0]) <= 0 &&
            			std::abs(vertex[1] - v.vertex[1]) <= 0 &&
            			std::abs(vertex[2] - v.vertex[2]) <= 0;
            };

            /**
             * Sort operator.
             *
             * Needed to keep only one unique vertex in the set.
             *
             * @param v Vertex to which shall be compared
             */
        	inline bool operator < (const TVertex v) const {
        		if (vertex[0] == v.vertex[0]){
        			if (vertex[1] == v.vertex[1])
        				return vertex[2] < v.vertex[2];
        			else
        				return vertex[1] < v.vertex[1];
        		}
        		else return vertex[0] < v.vertex[0];
        	};
        };

        /**
         * Kd-tree node class.
         */
        class KDNode{
            private:
                KDNode *hichild; ///< First leaf of this node
                KDNode *lochild; ///< Second leaf of this node
                KDNode *parent; ///< Root of this node
                float lo[3]; ///< The three largest components of the node vertices
                float hi[3]; ///< The three lowest components of the node vertices
                int splitdir; ///< Direction in which this node is split (0=x, 1=y, 2=z)
                int depth; ///< Depth of the node in the tree
                unsigned tricount; ///< Count of triangles stored in this node
                std::vector<Triangle*> tris; ///< List of triangles stored in this node

                /**
                 * Test if triangle intersects with this node.
                 *
                 * @param tri Triangle that shall be tested for intersection.
                 *
                 * @return Returns true if triangle intersects this node.
                 */
                bool TriangleInBox(Triangle *tri);

                /**
                 * Test all triangles in this tree and his leaves for intersection with line segment p1->p2.
                 *
                 * @param p1 Line start point
                 * @param p2 Line end point
                 * @param colls Found intersections are added to this list
                 *
                 * @return Returns true if at least one intersection was found
                 */
                bool TestCollision(const long double p1[3], const long double p2[3], std::set<TCollision> &colls);
            public:
                /**
                 * Constructor.
                 *
                 * @param boxlo Three lowest components of node vertices
                 * @param boxhi Three highes components of node vertices
                 * @param depth Depth of this node in the tree
                 * @param aparent Root of this node
                 */
                KDNode(const float boxlo[3], const float boxhi[3], const int depth, KDNode *aparent); // constructor
                ~KDNode(); ///< Destructor

                /**
                 * Add triangle to this node.
                 *
                 * @param tri Triangle to add
                 */
                void AddTriangle(Triangle *tri);

                /**
                 * Split this node into two leaves.
                 *
                 * Splits the node into two new leaves if no stop criterion is met.
                 */
                void Split();   // split node in two leaves

                /**
                 * Find the smallest node which contains the line segment p1->p2 and call TestCollision there.
                 *
                 * @param p1 Line start point
                 * @param p2 Line end point
                 * @param lastnode Returns smallest node which contains line segment so next search can be started there
                 * @param colls Found collisions are added to this list
                 *
                 * @return Returns true if at least one collision was found
                 */
                bool Collision(const long double p1[3], const long double p2[3], KDNode* &lastnode, std::set<TCollision> &colls);

                /**
                 * Check if point is inside this node
                 *
                 * @param p Point
                 *
                 * @return Returns true if point is contained in node
                 */
                template <typename coord> bool PointInBox(const coord p[3]) {
                	return ((p[0] <= hi[0]) && (p[0] >= lo[0]) && (p[1] <= hi[1]) && (p[1] >= lo[1]) && (p[2] <= hi[2]) && (p[2] >= lo[2]));
                }

                /**
                 *  Test if line segment p1->p2 intersects this node.
                 *
                 *  @param p1 Line start point
                 *  @param p2 Line end point
                 *
                 *  @return Returns true if line intersects this node
                 */
                template <typename coord> bool SegmentInBox(const coord p1[3], const coord p2[3]);

        };

        std::set<TVertex> allvertices; ///< list of all triangle vertices
        KDNode *root; ///< root node
        KDNode *lastnode; ///< remember last collision-tested node to speed up search when segments are adjacent

    public:
        KDTree(); ///< constructor
        ~KDTree(); ///< destructor
        float lo[3],hi[3]; ///< bounding box of root node
        std::vector<Triangle> alltris; ///< list of all triangles in the tree
#else
    public:
        std::list<TTriangle> triangles;
        CTree tree;
#endif
    public:
        /**
         * Read STL-file.
         *
         * @param filename Filename of STL file
         * @param ID Assign this ID to all triangles in this file
         * @param name Returns name of file
         */
        void ReadFile(const char *filename, const unsigned ID, char name[80] = NULL);
        void Init(); ///< create KD-tree

        /**
         * Test line segment p1->p2 for collision with all triangles in previously read files.
         *
         * @param p1 Line start point
         * @param p2 Line end point
         * @param colls Found collisions are added to this list
         *
         * @return Returns true if at least one collision was found
         */
        bool Collision(const long double p1[3], const long double p2[3], std::set<TCollision> &colls);

#ifndef USE_CGAL
       /**
         * Test if point is inside root node.
         *
         * @param p Point
         */
        bool SegmentInBox(const long double p1[3], const long double p2[3]){ return (root && root->SegmentInBox(p1,p2)); };
#endif
};

#endif // KDTREE_H_
