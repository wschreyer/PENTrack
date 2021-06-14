#include "trianglemesh.h"

#include <fstream>
#include <random>
#include <boost/format.hpp>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/repair.h>


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
	std::cout << "Reading '" << filename << "' containing " << filefacecount << " triangles ... ";    // print header

	std::vector<CPoint> vertices;
	std::vector<std::vector<size_t> > faces;
	while (f){
		f.seekg(3*4,std::fstream::cur);  // skip normal in STL-file (will be calculated from vertices)
        std::vector<size_t> vidx;
		for (short j = 0; j < 3; j++){
		    float v[3];
		    f.read((char*)v,12);
			if (f) {
                CPoint p(std::abs(v[0]) < REFLECT_TOLERANCE ? 0. : v[0], std::abs(v[1]) < REFLECT_TOLERANCE ? 0. : v[1], std::abs(v[2]) < REFLECT_TOLERANCE ? 0. : v[2]);
                auto vertex = vertices.end(); //std::find_if(vertices.begin(), vertices.end(), [&p](const CPoint &p2) { return CGAL::squared_distance(p, p2) < std::pow(REFLECT_TOLERANCE, 2); });
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

    namespace PMP = CGAL::Polygon_mesh_processing;
    typedef boost::graph_traits<CMesh>::face_descriptor fd;
    auto cerr_precision = std::cerr.precision(3);
    auto cout_precision = std::cout.precision(3);
    PMP::repair_polygon_soup(vertices, faces/*, CGAL::parameters::require_same_orientation(true)*/);
    PMP::orient_polygon_soup(vertices, faces);
    std::unique_ptr<CMesh> mesh(new CMesh());

    if (not PMP::is_polygon_soup_a_polygon_mesh(faces))
        //throw(std::runtime_error("Triangles do not form a mesh"));
        std::cerr << "Triangles in " << filename << " do not form a mesh\n";
    PMP::polygon_soup_to_polygon_mesh(vertices, faces, *mesh);
//    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(*mesh);
    double A = PMP::area(*mesh)*1e4;
    double V = PMP::volume(*mesh)*1e6;

    auto fccmap = mesh->add_property_map<fd, boost::graph_traits<CMesh>::faces_size_type>("f:CC").first;
    auto num = PMP::connected_components(*mesh, fccmap);
    std::cout << "built mesh with " << mesh->number_of_faces() << " triangles and " << num << " components (" << A << "cm2, " << V << "cm3)\n";
    int affected_components = 0;
    double border_length = 0.;
    double self_intersecting_area = 0.;
    for (size_t i = 0; i < num; ++i) {
        CGAL::Face_filtered_graph<CMesh> ffg(*mesh, i, fccmap);
        bool not_closed = not CGAL::is_closed(ffg);
        bool not_bounding = not PMP::does_bound_a_volume(ffg);
        bool self_intersecting = PMP::does_self_intersect(ffg);
        if (not_closed){
            std::vector<boost::graph_traits<CMesh>::halfedge_descriptor> border_edges;
            PMP::border_halfedges(ffg, std::back_inserter(border_edges));
            for (auto edge: border_edges)
                border_length += PMP::edge_length(edge, *mesh);
        }
        if (self_intersecting){
            std::vector<std::pair<fd, fd> > self_intersecting_face_pairs;
            PMP::self_intersections(ffg, std::back_inserter(self_intersecting_face_pairs));
            std::set<fd> self_intersecting_faces;
            for (auto face_pair: self_intersecting_face_pairs) {
                self_intersecting_faces.insert(face_pair.first);
                self_intersecting_faces.insert(face_pair.second);
            }
            for (auto face: self_intersecting_faces)
                self_intersecting_area += PMP::face_area(face, *mesh);
        }
        if (not_closed or not_bounding or self_intersecting) {
            ++affected_components;
        }
    }
    if (affected_components > 0) {
        std::cerr << "\nWarning: " << affected_components << " of " << num << " components in "
                  << filename << " have holes with total circumference "
                  << border_length * 1e2 << "cm and " << self_intersecting_area * 1e4
                  << "cm2 of their area is self-intersecting!\n\n";
    }
    std::cerr.precision(cout_precision);
    std::cerr.precision(cerr_precision);

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
