#include "trianglemesh.h"

#include <fstream>
#include <random>
#include <boost/format.hpp>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

// read triangles from STL-file
std::string TTriangleMesh::ReadFile(const std::string &filename, const int ID){
	std::ifstream f(filename, std::fstream::binary);
	if (!f.is_open())
		throw std::runtime_error( (boost::format("Could not open %1%") % filename).str() );

	char header[80];
	f.read(header, 80); // read 80-byte header
	std::string sldname(std::move(header), 80);
	sldname.erase(sldname.find_last_not_of(" ") + 1); // strip trailing whitespace from header

	unsigned int filefacecount;
	f.read((char*)&filefacecount,4);
	if (filefacecount == 0)
		throw std::runtime_error( (boost::format("%1% contains no triangles") % filename).str() );
	std::cout << "Reading '" << sldname << "' from '" << filename << "' containing " << filefacecount << " triangles ... ";    // print header

	std::vector<CPoint> vertices;
	std::vector<std::vector<size_t> > faces;
	while (f){
		f.seekg(3*4,std::fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
        std::vector<size_t> vidx;
		for (short j = 0; j < 3; j++){
		    float v[3];
		    f.read((char*)v,12);
			if (f) {
                CPoint p(v[0], v[1], v[2]);
                auto vertex = std::find_if(vertices.begin(), vertices.end(), [&p](const CPoint &p2) { return CSegment(p, p2).squared_length() < REFLECT_TOLERANCE * REFLECT_TOLERANCE; });
                vidx.push_back(std::distance(vertices.begin(), vertex));
                if (vertex == vertices.end())
                    vertices.push_back(p);
            }
		}
		if (f) {
            faces.push_back(vidx);
            f.seekg(2, std::fstream::cur);    // 2 attribute bytes, not used in the STL standard (http://www.ennex.com/~fabbers/StL.asp)
        }
	}
	f.close();

	if (faces.size() != filefacecount)
		throw std::runtime_error( (boost::format("%1% should contain %2% triangles but read %3%") % filename % filefacecount % faces.size()).str() );

    CGAL::Polygon_mesh_processing::repair_polygon_soup(vertices, faces);
    if (not CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(faces))
        throw(std::runtime_error("Triangles do not form a mesh"));
    std::unique_ptr<CMesh> mesh(new CMesh());
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(vertices, faces, *mesh);
    double A = CGAL::Polygon_mesh_processing::area(*mesh)*1e4;
    double V = CGAL::Polygon_mesh_processing::volume(*mesh)*1e6;
    std::cout << "Built mesh with " << faces.size() << " triangles (" << A << "cm2, " << V << "cm3)\n";
    if (not CGAL::is_closed(*mesh))
        throw std::runtime_error("Mesh is not closed");
        //std::cout << "Mesh is not closed!\n";
    if (not CGAL::Polygon_mesh_processing::does_bound_a_volume(*mesh))
        throw std::runtime_error("Mesh does not bound a volume");
        //std::cout << "Mesh does not bound a volume!\n";
    if (CGAL::Polygon_mesh_processing::does_self_intersect(*mesh))
        throw std::runtime_error("Mesh is self-intersecting");
        //std::cout << "Mesh is self-intersecting\n";

    std::vector<double> areas;
    std::transform(mesh->faces_begin(), mesh->faces_end(), std::back_inserter(areas), [&mesh](const CMesh::Face_index &fi){ return CGAL::Polygon_mesh_processing::face_area(fi, *mesh); });
    std::discrete_distribution<size_t> triangle_sampler(areas.begin(), areas.end());

    std::unique_ptr<CTree> tree(new CTree(mesh->faces_begin(), mesh->faces_end(), *mesh));
    tree->accelerate_distance_queries();

    std::vector<double> total_areas;
    std::transform(meshes.begin(), meshes.end(), std::back_inserter(total_areas), [](const CTriangleMesh &m){ return CGAL::Polygon_mesh_processing::area(*m.mesh); });
    mesh_sampler = std::discrete_distribution<size_t>(total_areas.begin(), total_areas.end());

    meshes.push_back({std::move(mesh), std::move(tree), ID, triangle_sampler});

	return sldname;
}


// test segment p1->p2 for collision with triangles and return a list of all found collisions
std::vector<TCollision> TTriangleMesh::Collision(const std::vector<double> &p1, const std::vector<double> &p2) const{
	CSegment segment(CPoint(p1[0], p1[1], p1[2]), CPoint(p2[0], p2[1], p2[2]));
	std::vector<TCollision> colls;
	for (auto &it: meshes) {
        std::vector<CIntersection> out;
        it.tree->all_intersections(segment, std::back_inserter(out)); // search intersections of segment with mesh
        for (auto &i: out){
            const CPoint *collp = boost::get<CPoint>(&(i->first));
            if (collp) { // if intersection is a point
                CVector n = CGAL::Polygon_mesh_processing::compute_face_normal(i->second, *it.mesh);
                colls.push_back(TCollision(segment, n, *collp, it.ID)); // add collision to list
            }
            else
                throw std::runtime_error("Segment-triangle intersection happened to not be a point");
        }
    }

	std::sort(	colls.begin(),
				colls.end(),
				[](const TCollision &c1, const TCollision &c2){
					if (c1.s == c2.s){
                        //std::cout << "Coincident collision between solids " << c1.ID << " and " << c2.ID << std::endl;
						return c1.ID > c2.ID;
					}
					else
						return c1.s < c2.s;
				}
	);
	return colls;

}


bool TTriangleMesh::InSolid(const double x, const double y, const double z) const{
    return std::any_of(meshes.begin(), meshes.end(), [x,y,z](const CTriangleMesh &mesh){
        return mesh.tree->number_of_intersected_primitives(CKernel::Ray_3(CPoint(x,y,z), CVector(0.,0.,1.))) % 2 != 0;
    });
}
