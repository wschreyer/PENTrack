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


#include "kdtree.h"

#define MAX_DEPTH 19            // max. depth of tree
#define MAX_FACECOUNT 10        // when the number of faces in one node exceeds this value the node is split up in two new nodes
#define TOLERANCE 1e-10         // points that are closer than this distance will be considered equal

// misc functions

// Read vector from STL-file
// STL format description: http://www.ennex.com/~fabbers/StL.asp
void ReadVector(ifstream &f, float v[3]){
    for (short i = 0; i < 3; i++){
        f.read((char*)&v[i],4);
    }
}

// test if two points are equal (closer than TOLERANCE)
inline bool PointsEqual(float p1[3], float p2[3]){
    //return p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2];
	float d1 = p1[0] - p2[0], d2 = p1[1] - p2[1], d3 = p1[2] - p2[2];
	if (abs(d1) < TOLERANCE && abs(d2) < TOLERANCE && abs(d3) < TOLERANCE)
    //if (d1*d1 + d2*d2 + d3*d3 < TOLERANCE*TOLERANCE)
        return true;
    return false;
}

//Dot-Product of two vectors
template <typename coord1, typename coord2> inline long double DotProduct(const coord1 x[3], const coord2 y[3]){
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//-------------------------------------------------------------------------------------------------------------
// Triangle class definition

// constructor
KDTree::Triangle::Triangle(ifstream &f, const unsigned aID){
    ID = aID;
//    normalIO = 0;
    f.seekg(3*4,fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
    short i;
    for (i = 0; i < 3; i++){
//        neighbours[i] = NULL;
        ReadVector(f, vertex[i]);   // read vertices
    }
    f.seekg(2,fstream::cur);    // 2 attribute bytes, not used in the STL standard (http://www.ennex.com/~fabbers/StL.asp)
    for (i = 0; i < 3; i++){
        hi[i] = max(max(vertex[0][i],vertex[1][i]),vertex[2][i]);   // calculate bounding box
        lo[i] = min(min(vertex[0][i],vertex[1][i]),vertex[2][i]);
        v[i] = vertex[1][i] - vertex[0][i]; // calculate edge vectors
        w[i] = vertex[2][i] - vertex[0][i];
    }
    for (i = 0; i < 3; i++)
        normal[i] = v[(i+1)%3]*w[(i+2)%3] - v[(i+2)%3]*w[(i+1)%3];  // recalculate normal

    vw = DotProduct(v,w); vv = DotProduct(v,v); ww = DotProduct(w,w);   // precalculate dotproducts
    long double parametric_factor = 1/(vw*vw - vv*ww);    // needed for parametric coordinates
    vw *= parametric_factor;
    vv *= parametric_factor;
    ww *= parametric_factor;
}


// does the segment p1->p2 intersect the triangle?
// http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#Segment-Triangle
bool KDTree::Triangle::intersect(const long double p1[3], const long double p2[3], long double &s){
    long double u[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};   // segment direction vector
    long double un = DotProduct(u,normal);
    if (un == 0)   // direction vector parallel to triangle plane?
        return false;
    long double x[3] = {p1[0] - vertex[0][0], p1[1] - vertex[0][1], p1[2] - vertex[0][2]};   // vector from first vertex to first segment point
    s = -DotProduct(x,normal) / un;  // parametric coordinate of intersection point (i = p1 + new_s*u)
    if ((s <= 0) || (s > 1)) // intersection point lies outside of the segment
        return false;
    x[0] += s*u[0];
    x[1] += s*u[1];
    x[2] += s*u[2]; // vector from first vertex to intersection point
    long double xw = DotProduct(x,w), xv = DotProduct(x,v);   // calculate parametric coordinates a,b of intersection point
    long double a = vw*xw - ww*xv;
    if ((a < 0) || (a > 1))
        return false;
    long double b = vw*xv - vv*xw;
    if ((b >= 0) && (a + b <= 1)){  // intersection point lies inside the triangle?
        return true;
    }
    return false;
}

/*void KDTree::Triangle::SetNormal(const short anormalIO){
    if (normalIO == 0)
        normalIO = anormalIO;
    else if (normalIO == anormalIO)
        return;
    else{
        printf("Normal direction switched!\n");
        return;
    }
    neighbours[0]->SetNormal(anormalIO);
    neighbours[1]->SetNormal(anormalIO);
    neighbours[2]->SetNormal(anormalIO);
}*/

//----------------------------------------------------------------------------------------------------------------------------------
// kd-tree node class definition

// constructor
KDTree::KDNode::KDNode(const float boxlo[3], const float boxhi[3], const int adepth, const short asplitdir, KDNode *aparent){
    for (int i = 0; i < 3; i++){
        lo[i] = boxlo[i];   // copy box coordinates into instance variables
        hi[i] = boxhi[i];
    }
    depth = adepth;         // set instance variables according to parameters
    splitdir = asplitdir;
    parent = aparent;
    hichild = lochild = NULL;   // initialize the remaining variables
    tricount = 0;
}

// destructor
KDTree::KDNode::~KDNode(){
    if (!tris.empty())   // empty triangle list
        tris.clear();
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
    if ((tri->lo[0] <= hi[0]) && (tri->lo[1] <= hi[1]) && (tri->lo[2] <= hi[2]) && // do the bounding boxes intersect?
        (tri->hi[0] >= lo[0]) && (tri->hi[1] >= lo[1]) && (tri->hi[2] >= lo[2])){
		if (SegmentInBox(tri->vertex[0],tri->vertex[1]) || SegmentInBox(tri->vertex[1], tri->vertex[2]) || SegmentInBox(tri->vertex[2], tri->vertex[0]))
			return true;    // one of the triangle sides cuts through the box?

		// one of the six box edges cuts through the triangle?
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
    if ((depth < MAX_DEPTH) && (tricount > MAX_FACECOUNT)){ // only split if node not too deep and contains enough triangles
        float newlo[3], newhi[3];
        int newdepth = depth + 1, newsplitdir = (splitdir+1)%3;
        for (short i = 0; i < 3; i++){
            newlo[i] = lo[i];
            newhi[i] = hi[i];
        }

        newlo[splitdir] += (hi[splitdir] - lo[splitdir])/2 - TOLERANCE; // split this node in half in splitdirection (with small overlap)
        hichild = new KDNode(newlo,newhi,newdepth,newsplitdir,this);    // create new leaves
        newhi[splitdir] = newlo[splitdir] + 2*TOLERANCE;
        newlo[splitdir] = lo[splitdir];
        lochild = new KDNode(newlo,newhi,newdepth,newsplitdir,this);

        Triangle *tri;
        while (!tris.empty()){  // empty list and add all triangles to leaves
            tri = tris.front();
            tris.pop_front();
            bool inhi = hichild->TriangleInBox(tri);
            bool inlo = lochild->TriangleInBox(tri);
            if (inhi)
                hichild->AddTriangle(tri);
            if (inlo)
                lochild->AddTriangle(tri);
            else if (!inhi && !inlo)
                printf("Gap in tree!\n");
        }
        hichild->Split();   // split leaves
        lochild->Split();
    }
}

// find neighbours of triangle
/*void KDTree::KDNode::FindNeighbour(Triangle* tri, const short vertexnumber){
    // go down to smallest node containing the triangle vertex
    if (hichild && tri->vertex[vertexnumber][splitdir] >= hichild->lo[splitdir])
        hichild->FindNeighbour(tri,vertexnumber);
    if (lochild && tri->vertex[vertexnumber][splitdir] <= lochild->hi[splitdir])
        lochild->FindNeighbour(tri,vertexnumber);
    else if (!hichild && !lochild){ // when reached smallest node, find neighbours
        for (list<Triangle*>::iterator i = tris.begin(); i != tris.end(); i++){ // search triangles sharing this vertex
            if (tri == (*i))
                continue;
            for (short k = 0; k < 3; k++){
                if (!(*i)->neighbours[(k+2)%3] && PointsEqual(tri->vertex[vertexnumber],(*i)->vertex[k])){
                    if (PointsEqual(tri->vertex[(vertexnumber+1)%3],(*i)->vertex[(k+2)%3])){
                        tri->neighbours[vertexnumber] = *i;
                        (*i)->neighbours[(k+2)%3] = tri;
                        return;
                    }
//                    else if (PointsEqual(tri->vertex[(vertexnumber+1)%3],(*i)->vertex[(k+1)%3])){
//                        tri->AddNeighbour(*i,vertexnumber);
//                        (*i)->AddNeighbour(tri,k);
//                        return;
                    }
                }
            }
        }
    }
}*/

// test all triangles in this node and his leaves for intersection with segment p1->p2
bool KDTree::KDNode::TestCollision(const long double p1[3], const long double p2[3], list<TCollision> &colls){
    if (tricount == 0) return false;
    bool result = false;
    long double s_loc;
    for (list<Triangle*>::iterator i = tris.begin(); i != tris.end(); i++){    // iterate through triangles stored in node
        if ((*i)->intersect(p1,p2,s_loc)){
        	long double n = sqrt(DotProduct((*i)->normal,(*i)->normal));
        	TCollision c;
        	c.s = s_loc;
        	c.normalx = (*i)->normal[0]/n;	// return normalized normal vector
        	c.normaly = (*i)->normal[1]/n;
        	c.normalz = (*i)->normal[2]/n;    
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
    else if (SegmentInBox(p1,p2)){
        return TestCollision(p1,p2,colls);
    }
    else
        printf("Segment outside of box!");
    return false;
}

// count triangles in this node and his leaves
unsigned KDTree::KDNode::facecount(){
    if (!hichild && !lochild) return tricount;
    unsigned c = 0;
    if (hichild) c += hichild->facecount();
    if (lochild) c += lochild->facecount();
    return c;
}

//-----------------------------------------------------------------------------------------------------------
// kd-tree class definition

// constructor
KDTree::KDTree(void (*ALog)(const char*, ...)){
    root = lastnode = NULL;
    lo[0] = lo[1] = lo[2] = INFINITY;
    hi[0] = hi[1] = hi[2] = -INFINITY;
    FLog = ALog;
}

// destructor
KDTree::~KDTree(){
    alltris.clear(); // empty triangle list
    if (root) delete root;  // delete tree
}

// read triangles from STL-file
void KDTree::ReadFile(const char *filename, const unsigned ID, char name[80]){
	if (root) return; // stop if Init was called already
    ifstream f(filename, fstream::binary);
    if (f.is_open()){
        unsigned i,j;
        char header[80];
        unsigned long int filefacecount;
        f.read((char*)header,80);   // read header
        for (i = 79; i >= 0; i--){
            if (header[i] == ' ') header[i] = 0;    // trim trailing whitespaces
            else break;
        }
        if (name) memcpy(name,header,80);
        f.read((char*)&filefacecount,4);
        FLog("Reading '%.80s' from '%s' containing %lu triangles ... ",header,filename,filefacecount);    // print header
        flush(cout);

        Triangle *tri;
        for (i = 0; i < filefacecount && !f.eof(); i++){
            tri = new Triangle(f,ID);  // read triangles from file
            for (j = 0; j < 3; j++){
                lo[j] = min(lo[j],tri->lo[j]);  // save lowest and highest vertices to get the size for the root node
                hi[j] = max(hi[j],tri->hi[j]);
            }
            alltris.push_back(tri); // add triangles to list
        }
        f.close();
        FLog("Read %d triangles\n",i);
    }
    else{
        FLog("Could not open '%s'!\n",filename);
        exit(-1);
    }
}

// build search tree
void KDTree::Init(/*const long double PointInVolume[3]*/){
    FLog("Edges are (%f %f %f),(%f %f %f)\n",lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);  // print the size of the root node

    root = new KDNode(lo,hi,0,2,NULL);  // create root node

	FLog("Building tree ... ");
	flush(cout);
    for (list<Triangle*>::iterator it = alltris.begin(); it != alltris.end(); it++)
        root->AddTriangle(*it); // add triangles to root node
    root->Split();  // split root node
    FLog("Wrote %u triangles\n",root->facecount());   // print number of triangles contained in tree

/*    if (PointInVolume){
		cout << "Finding neighbours...\n";
	    for (list<Triangle*>::iterator it = alltris.begin(); it != alltris.end(); it++){    // find neighbours for every triangle
	        for (unsigned i = 0; i < 3; i++){   // iterate vertices
	            if (!(*it)->neighbours[i])
	                root->FindNeighbour(*it,i);
	            if (!(*it)->neighbours[i]) // if a triangle has less than three neighbours the model is broken
	                printf("%p missing neighbour %u (%f %f %f),(%f %f %f)!\n",*it,i,(*it)->vertex[i][0],(*it)->vertex[i][1],(*it)->vertex[i][2],
	                                                                                (*it)->vertex[(i+1)%3][0],(*it)->vertex[(i+1)%3][1],(*it)->vertex[(i+1)%3][2]);
	        }
	    }

		cout << "Defining interior...\n";
        long double p2[3] = {PointInVolume[0],PointInVolume[1],lo[2]};
        list<long double> s;
        list<Triangle*> tri;
        if (root->Collision(PointInVolume,p2,lastnode,s,tri)){
            if (tri.front()->normal[2] < 0)
                tri.front()->SetNormal(-1);
            else
                tri.front()->SetNormal(1);
        }
        else
            cout << "Controlpoint not inside closed surface!\n";
    }
    cout << "Geometry loaded!\n";
*/
}

// test segment p1->p2 for collision with a triangle and return a list of all found collisions
bool KDTree::Collision(const long double p1[3], const long double p2[3], list<TCollision> &colls){
    if (!root) return false; // stop if Init was not called yet
    if (!lastnode) lastnode = root; // start search in last tested node, when last node is not known start in root node
    if (lastnode->Collision(p1,p2,lastnode,colls)){
    	colls.sort();
    	return true;
    }
    else return false;
    
}

// test if a point is located inside a closed(!) volume (volume defined by Init)
/*bool KDTree::PointInVolume(const long double p[3]){
	if (PointInBox(p)){
	    long double p2[3] = {p[0],p[1],lo[2]};
	    list<TCollision> colls;*/
//	    return (root && root->Collision(p,p2,lastnode,colls) && ((colls.front().normalz < 0)/* == (tri.front()->normalIO > 0)*/));
/*	}
	return false;
}*/

