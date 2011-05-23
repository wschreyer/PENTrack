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


#include <cstring>
#include <fstream>
#include <algorithm>

#include "kdtree.h"

#define MIN_SIZE 0.02		// when the size of a node is smaller than that, it will be splitted no more
#define MAX_FACECOUNT 10        // when the number of faces in one node exceeds this value the node is split up in two new nodes
#define TOLERANCE 0         // kd-nodes will overlap this far

// misc functions

//Dot-Product of two vectors
template <typename coord1, typename coord2> inline long double DotProduct(const coord1 x[3], const coord2 y[3]){
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline void CrossProduct(float v1[3], float v2[3], long double v[3]){
    for (short i = 0; i < 3; i++)
        v[i] = v1[(i+1)%3]*v2[(i+2)%3] - v1[(i+2)%3]*v2[(i+1)%3];  // recalculate normal
}

//-------------------------------------------------------------------------------------------------------------
// Triangle class definition

Triangle::Triangle(const float *vertices[3]){
	vertex[0] = vertices[0];
	vertex[1] = vertices[1];
	vertex[2] = vertices[2];
    float v[3] = {vertex[1][0] - vertex[0][0], vertex[1][1] - vertex[0][1], vertex[1][2] - vertex[0][2]}; // triangle edge vectors
    float w[3] = {vertex[2][0] - vertex[0][0], vertex[2][1] - vertex[0][1], vertex[2][2] - vertex[0][2]};
    CrossProduct(v,w,normal);
    normal[0] *= 0.5;
    normal[1] *= 0.5;
    normal[2] *= 0.5; // length = triangle area
}


// does the segment p1->p2 intersect the triangle?
// http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#Segment-Triangle
bool Triangle::intersect(const long double p1[3], const long double p2[3], long double &s){
    long double u[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};   // segment direction vector
    long double un = DotProduct(u,normal);
    if (un == 0)   // direction vector parallel to triangle plane?
        return false;
    long double x[3] = {p1[0] - vertex[0][0], p1[1] - vertex[0][1], p1[2] - vertex[0][2]};   // vector from first vertex to first segment point
    s = -DotProduct(x,normal) / un;  // parametric coordinate of intersection point (i = p1 + s*u)
    if ((s <= 0) || (s > 1)) // intersection point lies outside of the segment
        return false;
    x[0] += s*u[0];
    x[1] += s*u[1];
    x[2] += s*u[2]; // vector from first vertex to intersection point
    
    // calculate barycentric coordinates a,b of intersection point
    float v[3] = {vertex[1][0] - vertex[0][0], vertex[1][1] - vertex[0][1], vertex[1][2] - vertex[0][2]}; // triangle edge vectors
    float w[3] = {vertex[2][0] - vertex[0][0], vertex[2][1] - vertex[0][1], vertex[2][2] - vertex[0][2]};
    long double vw = DotProduct(v,w), vv = DotProduct(v,v), ww = DotProduct(w,w);
    long double parametric_factor = 1/(vw*vw - vv*ww);
    long double xw = DotProduct(x,w)*parametric_factor, xv = DotProduct(x,v)*parametric_factor;
    long double a = vw*xw - ww*xv;
    if ((a < 0) || (a > 1))
        return false;
    long double b = vw*xv - vv*xw;
    if ((b >= 0) && (a + b <= 1)){  // intersection point lies inside the triangle?
        return true;
    }
    return false;
}


//----------------------------------------------------------------------------------------------------------------------------------
// kd-tree node class definition

unsigned facecount;
unsigned nodecount;
unsigned emptynodecount;

// constructor
KDTree::KDNode::KDNode(const float boxlo[3], const float boxhi[3], const int adepth, KDNode *aparent){
    for (int i = 0; i < 3; i++){
        lo[i] = boxlo[i];   // copy box coordinates into instance variables
        hi[i] = boxhi[i];
    }
    depth = adepth;         // set instance variables according to parameters
    if (hi[0] - lo[0] > hi[1] - lo[1]){
        if (hi[0] - lo[0] > hi[2] - lo[2]) splitdir = 0;
        else splitdir = 2;
    }
    else{
        if (hi[1] - lo[1] > hi[2] - lo[2]) splitdir = 1;
        else splitdir = 2;
    }
    parent = aparent;
    hichild = lochild = NULL;   // initialize the remaining variables
    tricount = 0;
}

// destructor
KDTree::KDNode::~KDNode(){
    if (hichild) delete hichild;    // free leaves
    if (lochild) delete lochild;
}

// test if a segment goes through the box
template <typename coord> bool KDTree::KDNode::SegmentInBox(const coord p1[3], const coord p2[3]){
    if (PointInBox(p1) || PointInBox(p2))   // one of the segment end points in box?
       return true;
    long double s, a, b, w;
    short j,k;
    for (short i = 0; i < 3; i++){  // segment cuts one of the six box faces?
        j = (i+1)%3;
        k = (i+2)%3;
		w = p2[i] - p1[i];
        s = (lo[i] - p1[i])/w;  // ith coordinate of intersection point
        if (s >= 0 && s <= 1){
            a = p1[j] + s*(p2[j] - p1[j]);
            b = p1[k] + s*(p2[k] - p1[k]);  // other two coordinates of intersection point on box face?
            if (((a >= lo[j]) && (a <= hi[j])) &&
				((b >= lo[k]) && (b <= hi[k]))) return true;
        }
        s = (hi[i] - p1[i])/w;
        if (s >= 0 && s <= 1){
            a = p1[j] + s*(p2[j] - p1[j]);
            b = p1[k] + s*(p2[k] - p1[k]);
            if (((a >= lo[j]) && (a <= hi[j])) &&
				((b >= lo[k]) && (b <= hi[k]))) return true;
        }
    }
    return false;
}

// test if a triangle is contained in the box
bool KDTree::KDNode::TriangleInBox(Triangle *tri){
    float trilo[3], trihi[3];
    for (short i = 0; i < 3; i++){
        trilo[i] = min(min(tri->vertex[0][i], tri->vertex[1][i]), tri->vertex[2][i]);
        trihi[i] = max(max(tri->vertex[0][i], tri->vertex[1][i]), tri->vertex[2][i]);
    }
    if ((trilo[0] <= hi[0]) && (trilo[1] <= hi[1]) && (trilo[2] <= hi[2]) && // do the bounding boxes intersect?
        (trihi[0] >= lo[0]) && (trihi[1] >= lo[1]) && (trihi[2] >= lo[2])){
		if (SegmentInBox(tri->vertex[0],tri->vertex[1]) || SegmentInBox(tri->vertex[1], tri->vertex[2]) || SegmentInBox(tri->vertex[2], tri->vertex[0]))
			return true;    // one of the triangle sides cuts through the box?

		// one of the 12 box edges cuts through the triangle?
		long double s = INFINITY;
		long double p1[3] = {lo[0],lo[1],lo[2]};	// lololo -> hilolo
		long double p2[3] = {hi[0],lo[1],lo[2]};
		if (tri->intersect(p1,p2,s))
			return true;
		p1[0] = hi[0];	// hilolo -> hihilo
		p2[1] = hi[1];
		if (tri->intersect(p1,p2,s))
			return true;
		p1[1] = hi[1];	// hihilo -> lohilo
		p2[0] = lo[0];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[0] = lo[0];	// lohilo -> lololo
		p2[1] = lo[1];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[1] = lo[1];	// lololo -> lolohi
		p2[2] = hi[2];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[2] = hi[2];	// lolohi -> hilohi
		p2[0] = hi[0];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[0] = hi[0];	// hilohi -> hihihi
		p2[1] = hi[1];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[1] = hi[1];	// hihihi -> lohihi
		p2[0] = lo[0];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[0] = lo[0];	// lohihi -> lolohi
		p2[1] = lo[1];
		if (tri->intersect(p1,p2,s))
			return true;

		p2[1] = hi[1];	// lohihi -> lohilo
		p2[2] = lo[2];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[0] = hi[0];	// hihihi -> hihilo
		p2[0] = hi[0];
		if (tri->intersect(p1,p2,s))
			return true;

		p1[1] = lo[1];	// hilohi -> hilolo
		p2[1] = lo[1];
		if (tri->intersect(p1,p2,s))
			return true;

    }
    return false;
}

// add triangle to node
void KDTree::KDNode::AddTriangle(Triangle *tri){
    tris.push_back(tri);    // add triangle to list
    tricount++;    // increase triangle counter
}

// split node into two new leaves
void KDTree::KDNode::Split(){
    if ((hi[splitdir]-lo[splitdir] > MIN_SIZE) && (tricount > MAX_FACECOUNT)){ // only split if node not too deep and contains enough triangles
        float newlo[3], newhi[3];
        int newdepth = depth + 1;
        for (short i = 0; i < 3; i++){
            newlo[i] = lo[i];
            newhi[i] = hi[i];
        }

        newlo[splitdir] += (hi[splitdir] - lo[splitdir])/2 - TOLERANCE; // split this node in half in splitdirection (with small overlap)
        hichild = new KDNode(newlo,newhi,newdepth,this);    // create new leaves
        newhi[splitdir] = newlo[splitdir] + 2*TOLERANCE;
        newlo[splitdir] = lo[splitdir];
        lochild = new KDNode(newlo,newhi,newdepth,this);

        vector<Triangle*> lotris, hitris;
        vector<Triangle*>::iterator i;
        for (i = tris.begin(); i != tris.end(); i++){
            bool inhi = hichild->TriangleInBox(*i);
            bool inlo = lochild->TriangleInBox(*i);
            if (inhi)
                hitris.push_back(*i);
            if (inlo)
                lotris.push_back(*i);
            else if (!inhi && !inlo)
                printf("Gap in tree!\n");
        }
        if (hitris != lotris){
            tris.clear();
            for (i = hitris.begin(); i != hitris.end(); i++)
                hichild->AddTriangle(*i);
            for (i = lotris.begin(); i != lotris.end(); i++)
                lochild->AddTriangle(*i);
            hitris.clear();
            lotris.clear();
            hichild->Split();   // split leaves
            lochild->Split();
        }
        else{
            hitris.clear();
            lotris.clear();
            delete hichild;
            hichild = NULL;
            delete lochild;
            lochild = NULL;
        }
    }
    else{
        facecount += tricount;
        nodecount++;
        if (tricount == 0) emptynodecount++;
    }
}

// test all triangles in this node and his leaves for intersection with segment p1->p2
bool KDTree::KDNode::TestCollision(const long double p1[3], const long double p2[3], list<TCollision> &colls){
    if (tricount == 0) return false;
    bool result = false;
    long double s_loc;
    for (vector<Triangle*>::iterator i = tris.begin(); i != tris.end(); i++){    // iterate through triangles stored in node
        if ((*i)->intersect(p1,p2,s_loc)){
        	long double n = sqrt(DotProduct((*i)->normal,(*i)->normal));
        	TCollision c;
        	c.s = s_loc;
        	c.normal[0] = (*i)->normal[0]/n;	// return normalized normal vector
        	c.normal[1] = (*i)->normal[1]/n;
        	c.normal[2] = (*i)->normal[2]/n;
        	c.ID = (*i)->ID;
        	c.tri = *i;
        	colls.push_back(c);
        	result = true;
        }
    }
    bool inhi = false, inlo = false;
    if (hichild && (inhi = hichild->SegmentInBox(p1,p2)))
        result |= hichild->TestCollision(p1,p2,colls);   // test in leaves
    if (lochild && (inlo = lochild->SegmentInBox(p1,p2)))
        result |= lochild->TestCollision(p1,p2,colls);
    else if ((hichild || lochild) && !inhi && !inlo)
        printf("Gap in tree!\n");
    return result;
}

// find smallest box which contains segment and call TestCollision there
bool KDTree::KDNode::Collision(const long double p1[3], const long double p2[3], KDNode* &lastnode, list<TCollision> &colls){
    colls.clear();
    if (hichild && hichild->PointInBox(p1) && hichild->PointInBox(p2))    // if both segment points are contained in the first leave, test collision there
        return hichild->Collision(p1,p2,lastnode,colls);
    else if (lochild && lochild->PointInBox(p1) && lochild->PointInBox(p2))   // if both segment points are contained in the second leave, test collision there
        return lochild->Collision(p1,p2,lastnode,colls);
    else if (PointInBox(p1) && PointInBox(p2)){ // else if both segment point are contained in this box, test collision here
        lastnode = this;    // remember this node to speed up search when segments are adjacent
        return TestCollision(p1,p2,colls);
    }
    else if (parent) // else if parent exists, test collision there
        return parent->Collision(p1,p2,lastnode,colls);
    else if (SegmentInBox(p1,p2)) // else if a part of the segment is inside the (root) box
        return TestCollision(p1,p2,colls);
    return false;
}

//-----------------------------------------------------------------------------------------------------------
// kd-tree class definition

// constructor
KDTree::KDTree(){
    root = lastnode = NULL;
    lo[0] = lo[1] = lo[2] = INFINITY;
    hi[0] = hi[1] = hi[2] = -INFINITY;
}

// destructor
KDTree::~KDTree(){
    if (root) delete root;  // delete tree
}

// read triangles from STL-file
void KDTree::ReadFile(const char *filename, const unsigned ID, char name[80]){
	if (root) return; // stop if Init was called already
    ifstream f(filename, fstream::binary);
    if (f.is_open()){
        char header[80];
        unsigned int filefacecount, i;
        f.read((char*)header,80);   // read header
        for (i = 79; i >= 0; i--){
            if (header[i] == ' ') header[i] = 0;    // trim trailing whitespaces
            else break;
        }
        if (name) memcpy(name,header,80);
        f.read((char*)&filefacecount,4);
        printf("Reading '%.80s' from '%s' containing %u triangles ... ",header,filename,filefacecount);    // print header
		fflush(stdout);

        for (i = 0; i < filefacecount && !f.eof(); i++){
            f.seekg(3*4,fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
            TVertex v;
            const float *vertices[3];
            for (short j = 0; j < 3; j++){
                f.read((char*)&v.vertex[0],4);
                f.read((char*)&v.vertex[1],4);
                f.read((char*)&v.vertex[2],4);
                vertices[j] = allvertices.insert(v).first->vertex; // insert vertex into vertex list and save pointer to inserted vertex into triangle
            }
            f.seekg(2,fstream::cur);    // 2 attribute bytes, not used in the STL standard (http://www.ennex.com/~fabbers/StL.asp)

            Triangle tri(vertices);
            tri.ID = ID;
            for (short j = 0; j < 3; j++){
                lo[j] = min(lo[j], min(min(tri.vertex[0][j], tri.vertex[1][j]), tri.vertex[2][j]));  // save lowest and highest vertices to get the size for the root node
                hi[j] = max(hi[j], max(max(tri.vertex[0][j], tri.vertex[1][j]), tri.vertex[2][j]));
            }
            alltris.push_back(tri); // add triangle to list
        }
        f.close();
        printf("Read %u triangles\n",i);
    }
    else{
        printf("Could not open '%s'!\n",filename);
        exit(-1);
    }
}

// build search tree
void KDTree::Init(){
    printf("Edges are (%f %f %f),(%f %f %f)\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);  // print the size of the root node

    root = new KDNode(lo,hi,0,NULL);  // create root node
    facecount = 0;
    nodecount = 1;
    emptynodecount = 0;

	printf("Building KD-tree ... ");
	fflush(stdout);
    for (vector<Triangle>::iterator it = alltris.begin(); it != alltris.end(); it++)
        root->AddTriangle(&(*it)); // add triangles to root node
    root->Split();  // split root node
    printf("Wrote %u triangles in %u nodes (%u empty)\n\n",facecount,nodecount,emptynodecount);   // print number of triangles contained in tree
    fflush(stdout);
}

// test segment p1->p2 for collision with a triangle and return a list of all found collisions
bool KDTree::Collision(const long double p1[3], const long double p2[3], list<TCollision> &colls){
    if (!root) return false; // stop if Init was not called yet
    if (!lastnode) lastnode = root; // start search in last tested node, when last node is not known start in root node
    if (lastnode->Collision(p1,p2,lastnode,colls)){
    	colls.sort();
    	colls.unique();
    	return true;
    }
    else return false;

}

