#include <iostream>
#include <fstream>
#include <list>
#include <cmath>

using namespace std;

// root of kd-Tree
class KDTree{
    private:

        // Triangle class
        class Triangle{
            private:
                long double v[3],w[3];    // vectors from first vertex to the other two vertices
                long double vw,ww,vv;     // dotproducts of v and w (precalculated for intersection algorithm)
            public:
                long double vertex[3][3]; // the three vertices of the triangle (in counterclockwise order according to STL standard)
                long double lo[3],hi[3];  // bounding box defined by the three lowest and three highest coordinates
                long double normal[3];    // normal of triangle-plane
                int surfacetype;          // user defined surface number
                short normalIO;           // normal points into volume (1) / out of volume (-1) / not specified (0)
                Triangle* neighbours[3];    // neighbouring triangles
                Triangle(ifstream &f, const int asurfacetype);  // constructor, read vertices and normal from STL-file
                bool intersect(const long double p1[3], const long double p2[3], long double &s);    // does the segment defined by p1,p2 intersect the triangle?
                void SetNormal(const short anormalIO);
        };

        // kd-Tree node class
        class KDNode{
            private:
                KDNode *hichild, *lochild, *parent; // leaves and root of this node
                long double lo[3], hi[3]; // node box defined by the three highest and three lowest coordinates
                short splitdir, depth; // direction in which this node will be split to form new leaves (0 = x, 1 = y, 2 = z), depth of this node
                unsigned tricount;  // count of stored triangles, count of stored "big" triangles (volume > 1/8 of box volume)
                list<Triangle*> tris;   // list of triangles in this box
                bool SegmentInBox(const long double p1[3], const long double p2[3]);    // test if segment p1->p2 cuts through box
                bool TriangleInBox(Triangle *tri); // test if triangle cuts through box
                bool TestCollision(const long double p1[3], const long double p2[3], long double &s, Triangle* &tri); // test all triangles for intersection in this tree and his leaves
            public:
                KDNode(const long double boxlo[3], const long double boxhi[3], const int depth, const short splitdir, KDNode *aparent); // constructor
                ~KDNode();  // destructor
                void AddTriangle(Triangle *tri);    // add triangle to list or to leaves
                void Split();   // split node in two leaves
                void FindNeighbour(Triangle* tri, const short vertexnumber); // find neighbours of a triangle
                bool Collision(const long double p1[3], const long double p2[3], long double &s, KDNode* &lastnode, Triangle* &tri);  // find the smallest box which contains the segment p1->p2 and call TestCollision there
                inline bool PointInBox(const long double p[3]); // test if point is inside box
                unsigned facecount();   // count triangles in this tree and his leaves
        };

        list<Triangle*> alltris; // list of all triangles in the tree
        KDNode *root;   // root node
        KDNode *lastnode;    // remember last collision-tested node to speed up search when segments are adjacent
    public:
        KDTree();   // constructor
        ~KDTree();  // destructor
        long double lo[3],hi[3];    // bounding box of root node
        void ReadFile(const char *filename, const int surfacetype);    // read STL-file
        void Init(const long double PointInVolume[3] = NULL);    // create tree
        bool Collision(const long double p1[3], const long double p2[3], long double &s, long double normal[3], int &surfacetype);  // test segment p1->p2 for collision with triangle
        bool PointInBox(const long double p[3]);  // test if point is inside root node
        bool PointInVolume(const long double p[3]);
};
