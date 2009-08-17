/**************************************************************
	This algorithm uses a kd-tree structure to search 
	for collisions with a surface consisting of a list
	of triangles.
	Initially, the triangles are read from a set of STL-files
	(http://www.ennex.com/~fabbers/StL.asp)	via
	ReadFile(filename,surfacetype) and stored in the kd-tree
	via Init(ControlPoint).
	You can define a surfacetype for each file which is
	returned on collision tests to identify different surfaces
	during runtime.
	By passing a control point to Init, the algorithm tests
	if the point is contained inside a closed(!) surface and
	can then check via PointInVolume(point) if another point
	is inside the same volume.
	During runtime segments point1->point2 can be checked for
	intersection with the surface via
	Collision(point1,point2,s,n,surfacetype). Collision returns
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

struct TCollision{
	long double s, normalx,normaly,normalz;
	unsigned ID;	
	void *tri;
	bool operator < (TCollision c){ return s < c.s; };
};

// root of kd-Tree
class KDTree{
    private:

        // Triangle class
        class Triangle{
            private:
                float v[3],w[3];    // vectors from first vertex to the other two vertices
                long double vw,ww,vv;     // dotproducts of v and w (precalculated for intersection algorithm)
            public:
                float vertex[3][3]; // the three vertices of the triangle (in counterclockwise order according to STL standard)
                float lo[3],hi[3];  // bounding box defined by the three lowest and three highest coordinates
                long double normal[3];    // normal of triangle-plane
                unsigned ID;          // user defined surface number
//                short normalIO;           // normal points into volume (1) / out of volume (-1) / not specified (0)
                Triangle* neighbours[3];    // neighbouring triangles
                Triangle(ifstream &f, const unsigned aID);  // constructor, read vertices and normal from STL-file
                bool intersect(const long double p1[3], const long double p2[3], long double &s);    // does the segment defined by p1,p2 intersect the triangle?
//                void SetNormal(const short anormalIO);
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
                KDNode(const float boxlo[3], const float boxhi[3], const int depth, const short splitdir, KDNode *aparent); // constructor
                ~KDNode();  // destructor
                void AddTriangle(Triangle *tri);    // add triangle to list or to leaves
                void Split();   // split node in two leaves
                void FindNeighbour(Triangle* tri, const short vertexnumber); // find neighbours of a triangle
                bool Collision(const long double p1[3], const long double p2[3], KDNode* &lastnode, list<TCollision> &colls);  // find the smallest box which contains the segment p1->p2 and call TestCollision there
                template <typename coord> bool PointInBox(const coord p[3]) {
                	return ((p[0] <= hi[0]) && (p[0] >= lo[0]) &&
							(p[1] <= hi[1]) && (p[1] >= lo[1]) &&
							(p[2] <= hi[2]) && (p[2] >= lo[2]));
                }; // test if point is inside box
                unsigned facecount();   // count triangles in this tree and his leaves
        };

        list<Triangle*> alltris; // list of all triangles in the tree
        KDNode *root;   // root node
        KDNode *lastnode;    // remember last collision-tested node to speed up search when segments are adjacent
    public:
        KDTree();   // constructor
        ~KDTree();  // destructor
        float lo[3],hi[3];    // bounding box of root node
        void ReadFile(const char *filename, const unsigned ID, char name[80] = NULL);    // read STL-file
        void Init(/*const long double PointInVolume[3] = NULL*/);    // create tree
        bool Collision(const long double p1[3], const long double p2[3], list<TCollision> &colls);  // test segment p1->p2 for collision with triangle
        bool PointInBox(const long double p[3]){ return (root && root->PointInBox(p)); };  // test if point is inside root node
//        bool PointInVolume(const long double p[3]);
};
