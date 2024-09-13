#include "vectormath.h"
#include <boost/format.hpp>

#include "mesh.h"

TMesh<double> constructMesh(const std::string &filename){
    namespace PMP = CGAL::Polygon_mesh_processing;
    std::cout << "Reading '" << filename << " ... ";
    TMesh<double> mesh;
    mesh.name = filename;
    mesh.mesh = std::unique_ptr<TMesh<double>::CMesh>(new TMesh<double>::CMesh());
    if (not PMP::IO::read_polygon_mesh(filename, *mesh.mesh, PMP::parameters::verbose(true))){
        throw std::runtime_error("Could not read polygon mesh from " + filename);
    }
    PMP::merge_reversible_connected_components(*mesh.mesh);
    mesh.fccmap = mesh.mesh->add_property_map<TMesh<double>::CFace, std::size_t>("f:CC").first;
    mesh.normalmap = mesh.mesh->add_property_map<TMesh<double>::CFace, TMesh<double>::CVector>("f:normals").first;
    PMP::compute_face_normals(*mesh.mesh, mesh.normalmap);
    mesh.area = PMP::area(*mesh.mesh);
    mesh.volume = PMP::volume(*mesh.mesh);

    auto cerr_precision = std::cerr.precision(3);
    auto cout_precision = std::cout.precision(3);

    auto num = PMP::connected_components(*mesh.mesh, mesh.fccmap);
    std::cout << "built mesh with " << mesh.mesh->number_of_faces() << " triangles and " << num << " components (" << mesh.area*1e4 << "cm2, " << mesh.volume*1e6 << "cm3).\n";
    if (not CGAL::is_closed(*mesh.mesh))
        std::cerr << "Mesh does not appear to be closed!\n";
    if (not PMP::does_bound_a_volume(*mesh.mesh))
        std::cerr << "Mesh does not appear to bound a volume!\n";
    if (PMP::does_self_intersect(*mesh.mesh))
        std::cerr << "Mesh appears to self intersect!\n";
    if (not PMP::is_outward_oriented(*mesh.mesh))
        std::cerr << "Mesh appears to be incorrectly oriented!\n";
    std::cout.precision(cout_precision);
    std::cerr.precision(cerr_precision);

    std::vector<double> areas;
    for (auto face: mesh.mesh->faces()){
        areas.push_back(PMP::face_area(face, *mesh.mesh));
    }
    mesh.triangle_sampler = std::discrete_distribution<size_t>(areas.begin(), areas.end());

    mesh.tree = std::unique_ptr<TMesh<double>::CTree>(new TMesh<double>::CTree(mesh.mesh->faces_begin(), mesh.mesh->faces_end(), *mesh.mesh));
    mesh.tree->build(*mesh.mesh);

//    mesh.tree.accelerate_distance_queries();

    return mesh;
}
