#include "trianglemesh.h"

#include <fstream>
#include "boost/format.hpp"


// read triangles from STL-file
std::string TTriangleMesh::ReadFile(const std::string &filename, const int sldindex){
	std::ifstream f(filename, std::fstream::binary);
	if (!f.is_open())
		throw std::runtime_error( (boost::format("Could not open %1%") % filename).str() );

	char header[80];
	f.read(header, 80); // read 80-byte header
	std::string sldname(std::move(header), 80);
	sldname.erase(sldname.find_last_not_of(" ") + 1); // strip trailing whitespace from header

	unsigned int filefacecount;
	unsigned int i = triangles.size();
	f.read((char*)&filefacecount,4);
	if (filefacecount == 0)
		throw std::runtime_error( (boost::format("%1% contains no triangles") % filename).str() );
	std::cout << "Reading '" << sldname << "' from '" << filename << "' containing " << filefacecount << " triangles ... ";    // print header

	float v[3][3];
	while (f){
		f.seekg(3*4,std::fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
		for (short j = 0; j < 3; j++){
			 f.read((char*)v[j],12);
		}
		f.seekg(2,std::fstream::cur);    // 2 attribute bytes, not used in the STL standard (http://www.ennex.com/~fabbers/StL.asp)
		if (f){
			CTriangle tri(CPoint(v[0][0], v[0][1], v[0][2]), CPoint(v[1][0], v[1][1], v[1][2]), CPoint(v[2][0], v[2][1], v[2][2]));
			if (tri.is_degenerate())
				throw std::runtime_error( (boost::format("%1% contains a degenerate triangle!") % filename).str() );
			triangles.push_back(std::make_pair(tri, sldindex));
		}
	}
	f.close();

	i = triangles.size() - i;
	if (i != filefacecount)
		throw std::runtime_error( (boost::format("%1% should contain %2% triangles but read %3%") % filename % filefacecount % i).str() );
	std::cout << "Read " << i << " triangles\n";

	tree.rebuild(triangles.begin(), triangles.end()); // rebuild AABB tree
	return sldname;
}


// test segment p1->p2 for collision with triangles and return a list of all found collisions
std::vector<TCollision> TTriangleMesh::Collision(const std::vector<double> &p1, const std::vector<double> &p2) const{
	CSegment segment(CPoint(p1[0], p1[1], p1[2]), CPoint(p2[0], p2[1], p2[2]));
	std::vector<CIntersection> out;
	std::vector<TCollision> colls;
	tree.all_intersections(segment, std::back_inserter(out)); // search intersections of segment with mesh
	for (std::vector<CIntersection>::iterator i = out.begin(); i != out.end(); i++){
		#if CGAL_VERSION_NR<1040301000
			const CPoint *collp = CGAL::object_cast<CPoint>(&(*i)->first);
		#else
			const CPoint *collp = boost::get<CPoint>(&((*i)->first));
		#endif
		if (collp) // if intersection is a point
			colls.push_back(TCollision(segment, (*i)->second, *collp)); // add collision to list
	}

	std::sort(	colls.begin(),
				colls.end(),
				[](const TCollision &c1, const TCollision &c2){
					if (c1.s == c2.s){
						assert(c1.sldindex != c2.sldindex);
						return c1.sldindex > c2.sldindex;
					}
					else
						return c1.s < c2.s;
				}
	);
	return colls;
}

