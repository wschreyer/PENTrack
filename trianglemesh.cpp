#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <fstream>

#include "trianglemesh.h"


// read triangles from STL-file
void TTriangleMesh::ReadFile(const char *filename, const unsigned ID, char name[80]){
	std::ifstream f(filename, std::fstream::binary);
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

        for (i = 0; i < filefacecount && !f.eof(); i++){
            f.seekg(3*4,std::fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
            float v[3][3];
            for (short j = 0; j < 3; j++){
                 f.read((char*)v[j],12);
            }
            f.seekg(2,std::fstream::cur);    // 2 attribute bytes, not used in the STL standard (http://www.ennex.com/~fabbers/StL.asp)
            triangles.push_back(TTriangle(CPoint(v[0][0], v[0][1], v[0][2]),
            							  CPoint(v[1][0], v[1][1], v[1][2]),
            							  CPoint(v[2][0], v[2][1], v[2][2]), ID));
        }
        f.close();
        printf("Read %u triangles\n",i);
    }
    else{
        printf("Could not open '%s'!\n",filename);
        std::exit(-1);
    }
}

// build search tree
void TTriangleMesh::Init(){
	tree.rebuild(triangles.begin(), triangles.end());
	printf("Edges are (%f %f %f),(%f %f %f)\n",tree.bbox().min(0),tree.bbox().min(1),tree.bbox().min(2),
												tree.bbox().max(0),tree.bbox().max(1),tree.bbox().max(2));  // print the size of the root node
}

// test segment p1->p2 for collision with triangles and return a list of all found collisions
bool TTriangleMesh::Collision(const long double p1[3], const long double p2[3], std::set<TCollision> &colls){
	CPoint point1(p1[0], p1[1], p1[2]);
	CPoint point2(p2[0], p2[1], p2[2]);
	CSegment segment(point1, point2);

	std::list<CIntersection> out;
	tree.all_intersections(segment, std::back_inserter(out));
	for (std::list<CIntersection>::iterator i = out.begin(); i != out.end(); i++){
		if (*i){
#if CGAL_VERSION_NR<1040301000
			const CPoint *collp = CGAL::object_cast<CPoint>(&(*i)->first);
#else
			CPoint *collp = boost::get<CPoint>(&((*i)->first));
#endif
			if (collp){
				TCollision coll;
				coll.s = sqrt(pow(collp->x() - p1[0], 2) + pow(collp->y() - p1[1], 2) + pow(collp->z() - p1[2], 2))
						/sqrt(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2) + pow(p2[2] - p1[2], 2));
				CIterator tri = (*i)->second;
				CVector n = tri->tri.supporting_plane().orthogonal_vector();
				coll.ID = tri->ID;
				coll.normal[0] = n[0]/sqrt(n.squared_length());
				coll.normal[1] = n[1]/sqrt(n.squared_length());
				coll.normal[2] = n[2]/sqrt(n.squared_length());
				colls.insert(coll);
			}
		}
	}
	return !(colls.empty());
}

