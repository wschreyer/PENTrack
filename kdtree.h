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

#include <cmath>
#include <vector>
#include <set>
#include <list>

using namespace std;

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

/**
 * Structure that is returned by KDTree::Collision.
 */
struct TCollision{
	long double s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
	long double normal[3]; ///< normal (length = 1) of intersected surface
	unsigned ID; ///< ID of intersected surface
	Triangle *tri; ///< pointer to intersected triangle
	inline bool operator < (const TCollision c) const { return s < c.s; }; ///< overloaded operator, needed for sorting
	inline bool operator == (const TCollision c) const { return c.tri == tri; }; ///< overloaded operator, to get rid of two equal collisions
};

/**
 * General kd-tree class.
 *
 * Encapsulates all the nodes.
 */
class KDTree{
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
            	return  abs(vertex[0] - v.vertex[0]) <= 0 &&
            			abs(vertex[1] - v.vertex[1]) <= 0 &&
            			abs(vertex[2] - v.vertex[2]) <= 0;
            };

            /**
             * Sort operator.
             *
             * Needed to keep only one unique vertex in the set.
             *
             * @param v Vertex to which shall be compared
             */
        	inline bool operator < (const TVertex v) const { return !(*this == v); };
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
                short splitdir; ///< Direction in which this node is split (0=x, 1=y, 2=z)
                short depth; ///< Depth of the node in the tree
                unsigned tricount; ///< Count of triangles stored in this node
                vector<Triangle*> tris; ///< List of triangles stored in this node

                /**
                 *  Test if line segment p1->p2 intersects this node.
                 *
                 *  @param p1 Line start point
                 *  @param p2 Line end point
                 *
                 *  @return Returns true if line intersects this node
                 */
                template <typename coord> bool SegmentInBox(const coord p1[3], const coord p2[3]);

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
                bool TestCollision(const long double p1[3], const long double p2[3], list<TCollision> &colls);
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
                bool Collision(const long double p1[3], const long double p2[3], KDNode* &lastnode, list<TCollision> &colls);

                /**
                 * Check if point is inside this node
                 *
                 * @param p Point
                 *
                 * @return Returns true if point is contained in node
                 */
                template <typename coord> bool PointInBox(const coord p[3]) {
                	return ((p[0] <= hi[0]) && (p[0] >= lo[0]) && (p[1] <= hi[1]) && (p[1] >= lo[1]) && (p[2] <= hi[2]) && (p[2] >= lo[2]));
                };
        };

        set<TVertex> allvertices; ///< list of all triangle vertices
        KDNode *root; ///< root node
        KDNode *lastnode; ///< remember last collision-tested node to speed up search when segments are adjacent
    public:
        KDTree(); ///< constructor
        ~KDTree(); ///< destructor
        float lo[3],hi[3]; ///< bounding box of root node
        vector<Triangle> alltris; ///< list of all triangles in the tree

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
        bool Collision(const long double p1[3], const long double p2[3], list<TCollision> &colls);

        /**
         * Test if point is inside root node.
         *
         * @param p Point
         */
        bool PointInBox(const long double p[3]){ return (root && root->PointInBox(p)); };
};
