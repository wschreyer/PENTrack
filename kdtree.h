/**************************************************************
	This algorithm uses a kd-tree structure to search
	for collisions with a surface consisting of a list
	of triangles.
	Initially, the triangles are read from a set of STL-files
	(http://www.ennex.com/~fabbers/StL.asp)	via
	ReadFile(filename,surfacetype) and stored in the kd-tree
	via Init().
	You can define a surfacetype for each file which is
	returned on collision tests to identify different surfaces
	during runtime.
	During runtime segments point1->point2 can be checked for
	intersection with the surface via
	Collision(point1,point2,list of TCollision). Collision returns
	true if an intersection occurred and gives the parametric
	coordinate s of the intersection point (I=p1+s*(p2-p1)),
	the normal n and the surfacetype of the intersected surface.
*****************************************************************/


#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <cmath>

using namespace std;

// Triangle class
struct Triangle{
    const float *vertex[3]; // the three vertices of the triangle (in counterclockwise order according to STL standard)
    unsigned ID;          // user defined surface number
    void CalcNormal(long double normal[3]);
    bool intersect(const long double p1[3], const long double p2[3], long double &s);    // does the segment defined by p1,p2 intersect the triangle?
};

// structure that is returned by KDTree::Collision
struct TCollision{
	long double s, normal[3], Area; // parametric coordinate of intersection point, normal and area of intersected surface
	unsigned ID;	// ID of intersected surface
	Triangle *tri;	// pointer to intersected triangle
	bool operator < (TCollision c){ return s < c.s; };	// overloaded operator, needed for sorting
	bool operator == (TCollision c){ return tri == c.tri; };	// overloaded operator, needed for unique()
};

// root of kd-Tree
class KDTree{
    private:
    	// triangle vertex class, each vertex is only stored once in a set to save memory
        struct TVertex{
            float vertex[3];
            inline bool operator == (const TVertex v) const { 
            	return abs(vertex[0] - v.vertex[0]) < 1e-10 && 
            			abs(vertex[1] - v.vertex[1]) < 1e-10 && 
            			abs(vertex[2] - v.vertex[2]) < 1e-10; 
            };
        	inline bool operator < (const TVertex v) const { return !(*this == v); };
        };

        // kd-Tree node class
        class KDNode{
            private:
                KDNode *hichild, *lochild, *parent; // leaves and root of this node
                float lo[3], hi[3]; // node box defined by the three highest and three lowest coordinates
                short splitdir, depth; // direction in which this node will be split to form new leaves (0 = x, 1 = y, 2 = z), depth of this node
                unsigned tricount;  // count of stored triangles
                list<Triangle*> tris;   // list of triangles in this box
                template <typename coord> bool SegmentInBox(const coord p1[3], const coord p2[3]);    // test if segment p1->p2 cuts through box
                bool TriangleInBox(Triangle *tri); // test if triangle cuts through box
                bool TestCollision(const long double p1[3], const long double p2[3], list<TCollision> &colls); // test all triangles for intersection in this tree and his leaves
            public:
                KDNode(const float boxlo[3], const float boxhi[3], const int depth, KDNode *aparent); // constructor
                ~KDNode();  // destructor
                void AddTriangle(Triangle *tri);    // add triangle to node
                void Split();   // split node in two leaves
                bool Collision(const long double p1[3], const long double p2[3], KDNode* &lastnode, list<TCollision> &colls);  // find the smallest box which contains the segment p1->p2 and call TestCollision there
                template <typename coord> bool PointInBox(const coord p[3]) {
                	return ((p[0] <= hi[0]) && (p[0] >= lo[0]) && (p[1] <= hi[1]) && (p[1] >= lo[1]) && (p[2] <= hi[2]) && (p[2] >= lo[2]));
                }; // test if point is inside box
        };

        set<TVertex> allvertices;	// list of triangle vertices
        KDNode *root;   // root node
        KDNode *lastnode;    // remember last collision-tested node to speed up search when segments are adjacent
        int (*FLog)(const char*, ...);	// function to print log messages (e.g. printf), given by constructor parameter
    public:
        KDTree(int (*ALog)(const char*, ...));   // constructor, parameter: see FLog above
        ~KDTree();  // destructor
        float lo[3],hi[3];    // bounding box of root node
        list<Triangle*> alltris; // list of all triangles in the tree
        void ReadFile(const char *filename, const unsigned ID, char name[80] = NULL);    // read STL-file
        void Init();    // create KD-tree
        bool Collision(const long double p1[3], const long double p2[3], list<TCollision> &colls);  // test segment p1->p2 for collision with triangle
        bool PointInBox(const long double p[3]){ return (root && root->PointInBox(p)); };  // test if point is inside root node
};
