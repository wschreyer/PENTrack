#ifndef MESH_H_
#define MESH_H_

#include <random>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/iterator/function_output_iterator.hpp>

/**
 * Structure returned by TMovingMesh::findCollisions.
 */
template<class Vector, typename Real = typename Vector::value_type>
struct TCollision{
    Vector collisionPoint;
    Real s; ///< parametric coordinate of intersection point (P = p1 + s*(p2 - p1))
    Vector normal; ///< normal (length = 1) of intersected surface
    Vector velocity; ///< velocity of surface
};


template<typename Real>
struct TMesh{
    typedef typename CGAL::Simple_cartesian<Real> CKernel; ///< Geometric Kernel used for CGAL types
    typedef typename CKernel::Point_3 CPoint; ///< CGAL point type
    typedef typename CKernel::Vector_3 CVector; ///< CGAL vector type
    typedef typename CKernel::Segment_3 CSegment; ///< CGAL segment type
    typedef typename CKernel::Ray_3 CRay;
    typedef typename CGAL::Surface_mesh<CPoint> CMesh; ///< CGAL triangle mesh type
    typedef typename CGAL::AABB_face_graph_triangle_primitive<CMesh> CPrimitive; ///< CGAL triangle type contained in AABB tree
    typedef typename CGAL::AABB_traits<CKernel, CPrimitive> CTraits; ///< CGAL triangle traits type
    typedef typename CGAL::AABB_tree<CTraits> CTree; ///< CGAL AABB tree type containing CPrimitives
    typedef typename boost::optional<typename CTree::Intersection_and_primitive_id<CSegment>::Type > CIntersection;
    typedef typename boost::graph_traits<CMesh>::face_descriptor CFace;

    CMesh mesh; ///< Triangle mesh
    CTree tree; ///< Axis-aligned bounding-box tree for fast intersection search
    typename CMesh::Property_map<CFace, std::size_t> fccmap;
    std::discrete_distribution<size_t> triangle_sampler; ///< Probability distribution to randomly sample triangles from mesh weighted by their areas.
    std::string name;
};

template<class Mesh>
void printMeshInfo(const Mesh &mesh){
    namespace PMP = CGAL::Polygon_mesh_processing;
    auto cerr_precision = std::cerr.precision(3);
    auto cout_precision = std::cout.precision(3);
    auto A = PMP::area(mesh.mesh)*1e4;
    auto V = PMP::volume(mesh.mesh)*1e6;

    auto num = PMP::connected_components(mesh.mesh, mesh.fccmap);
    std::cout << "built mesh with " << mesh.mesh.number_of_faces() << " triangles and " << num << " components (" << A << "cm2, " << V << "cm3)\n";
    int affected_components = 0;
    double border_length = 0.;
    double self_intersecting_area = 0.;
    for (size_t i = 0; i < num; ++i) {
        CGAL::Face_filtered_graph<typename Mesh::CMesh> ffg(mesh.mesh, i, mesh.fccmap);
        bool not_closed = not CGAL::is_closed(ffg);
        bool not_bounding = not PMP::does_bound_a_volume(ffg);
        bool self_intersecting = PMP::does_self_intersect(ffg);
        if (not_closed){
            PMP::border_halfedges(ffg, boost::make_function_output_iterator([&border_length, &mesh](auto edge){ border_length += PMP::edge_length(edge, mesh.mesh); }));
        }
        if (self_intersecting){
            PMP::self_intersections(ffg, boost::make_function_output_iterator([&self_intersecting_area, &mesh](auto faces){ self_intersecting_area += PMP::face_area(faces.first, mesh.mesh) + PMP::face_area(faces.second, mesh.mesh); }));
        }
        if (not_closed or not_bounding or self_intersecting) {
            ++affected_components;
        }
    }
    if (affected_components > 0) {
        std::cerr << "\nWarning: " << affected_components << " of " << num << " components in "
                << mesh.name << " have holes with total circumference "
                << border_length * 1e2 << "cm and " << self_intersecting_area * 1e4
                << "cm2 of their area is self-intersecting!\n\n";
    }
    std::cout.precision(cout_precision);
    std::cerr.precision(cerr_precision);
}


TMesh<double> createMesh(const std::string &filename){
    namespace PMP = CGAL::Polygon_mesh_processing;
    std::cout << "Reading '" << filename << " ... ";
    TMesh<double> mesh;
    mesh.name = filename;
    if (not PMP::IO::read_polygon_mesh(filename, mesh.mesh)){
        throw std::runtime_error("Could not read polygon mesh from " + filename);
    }
    mesh.fccmap = mesh.mesh.add_property_map<TMesh<double>::CFace, std::size_t>("f:CC").first;
    printMeshInfo(mesh);

    std::vector<double> areas;
    for (auto face: mesh.mesh.faces()){
        areas.push_back(PMP::face_area(face, mesh.mesh));
    }
    mesh.triangle_sampler = std::discrete_distribution<size_t>(areas.begin(), areas.end());

    mesh.tree = TMesh<double>::CTree(mesh.mesh.faces_begin(), mesh.mesh.faces_end(), mesh.mesh);
    mesh.tree.accelerate_distance_queries();

    return mesh;
}


template<class Mesh, typename Real, class Vector>
std::vector<std::pair<Vector, Vector> > findCollisionsWithMesh(const Mesh &mesh, Real timeStart, const Vector &segmentStart, Real timeEnd, const Vector &segmentEnd){
    typename Mesh::CSegment segment(typename Mesh::CPoint(boost::qvm::X(segmentStart), boost::qvm::Y(segmentStart), boost::qvm::Z(segmentStart)),
                                    typename Mesh::CPoint(boost::qvm::X(segmentEnd), boost::qvm::Y(segmentEnd), boost::qvm::Z(segmentEnd)));
    std::vector<std::pair<Vector, Vector> > colls;
    std::vector<typename Mesh::CIntersection> out;
    mesh.tree.all_intersections(segment, std::back_inserter(out)); // search intersections of segment with mesh
    for (auto i: out){
        auto collisionPoint = boost::get<typename Mesh::CPoint>(&(i->first));
        if (collisionPoint) { // if intersection is a point
            auto normal = CGAL::Polygon_mesh_processing::compute_face_normal(i->second, mesh.mesh);
            colls.push_back(std::make_pair(Vector{collisionPoint->x(), collisionPoint->y(), collisionPoint->z()}, Vector{normal.x(), normal.y(), normal.z()})); // add collision to list
        }
        else
            throw std::runtime_error("Segment-triangle intersection happened to not be a point");
    }
    return colls;
}


template<class Mesh, typename Real = typename Mesh::Real, class Vector, class Path>
std::vector<TCollision<Vector> > findCollisionsWithMovingMesh(const Mesh &mesh, const Real timeStart, const Vector &segmentStart, const Real timeEnd, const Vector &segmentEnd, const Path &path){
    auto segmentStartLocal = passiveTransform(segmentStart, interpolateTransformation(path, timeStart));
    auto segmentEndLocal = passiveTransform(segmentEnd, interpolateTransformation(path, timeEnd));
    auto collisions = findCollisionsWithMesh(mesh, timeStart, segmentStartLocal, timeEnd, segmentEndLocal);
    std::vector<TCollision<Vector> > movingCollisions;
    for (auto collision: collisions){
        auto s = boost::qvm::dot(collision.first - segmentStartLocal, segmentEndLocal - segmentStartLocal)/boost::qvm::mag(segmentEndLocal - segmentStartLocal);
        auto collisionTime = timeStart + s*(timeEnd - timeStart);
        auto collisionPoint = activeTransform(collision.first, interpolateTransformation(path, collisionTime));
        movingCollisions.push_back({collisionPoint, s, activeRotate(collision.second, interpolateRotation(path, collisionTime)), interpolateRelativeVelocity(path, collisionTime, collisionPoint)});
    }
    return movingCollisions;
}

template<class Mesh, typename Real, class Vector, class Path>
bool inMovingSolid(const Mesh &mesh, const Real time, const Vector &position, const Path &path){
    return inSolid(mesh, passiveTransform(position, interpolateTransformation(path, time)));
}

template<class Mesh, class Vector>
bool inSolid(const Mesh &mesh, const Vector &position){
    return mesh.tree.number_of_intersected_primitives(typename Mesh::CRay(typename Mesh::CPoint(position[0], position[1], position[2]), typename Mesh::CVector(0.,0.,1.))) % 2 != 0;
}

template<class Mesh, class RandomGenerator, typename Real = double>
std::pair<std::array<Real, 3>, std::array<Real, 3> > randomPointOnSurface(Mesh &mesh, RandomGenerator &rand){
    auto faceidx = *(mesh.mesh.faces_begin() + mesh.triangle_sampler(rand));
    Real u(std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rand));
    Real v(std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rand)*(1 - u));
    auto pp = CGAL::Polygon_mesh_processing::construct_point(std::make_pair(faceidx, std::array<Real, 3>{u, v, 1 - u - v}), mesh.mesh);
    auto nv = CGAL::Polygon_mesh_processing::compute_face_normal(faceidx, mesh.mesh);
    return {{pp.x(), pp.y(), pp.z()}, {nv.x(), nv.y(), nv.z()}};
}

template<class Mesh, typename Real = typename Mesh::Real, class Path, class RandomGenerator>
std::pair<std::array<Real, 3>, std::array<Real, 3> > randomPointOnMovingSurface(Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto [position, normal] = randomPointOnSurface(mesh, rand);
    position = activeTransform(position, interpolateTransformation(path, time));
    normal = activeRotate(normal, interpolateRotation(path, time));
    return std::make_pair(position, normal);
}

template<class Mesh, class RandomGenerator> 
std::array<double, 3> randomPointInVolume(const Mesh &mesh, RandomGenerator &rand){
    std::array<double, 3> p;
    do{
        p = randomPointInBoundingBox(mesh, rand);
    }while(not inSolid(mesh, p));
    return p;
}

template<class Mesh, class RandomGenerator>
std::array<double, 3> randomPointInBoundingBox(const Mesh &mesh, RandomGenerator &rand){
    auto bbox = mesh.tree.bbox();
    std::uniform_real_distribution<double> unidist(0, 1);
    return {bbox.xmin() + unidist(rand)*(bbox.xmax() - bbox.xmin()),
            bbox.ymin() + unidist(rand)*(bbox.ymax() - bbox.ymin()),
            bbox.zmin() + unidist(rand)*(bbox.zmax() - bbox.zmin())};
}

template<class Mesh, typename Real = typename Mesh::Real, class RandomGenerator, class Path> 
std::array<Real, 3> randomPointInMovingVolume(const Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto bbox = mesh.tree.bbox();
    std::array<Real, 3> p;
    do{
        p = randomPointInBoundingBox(mesh, rand);
    }while(not inSolid(mesh, p));
    return activeTransform(p, interpolateTransformation(path, time));
}

template<class Mesh, typename Real = typename Mesh::Real, class RandomGenerator, class Path>
std::array<Real, 3> randomPointInMovingBoundingBox(const Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto p = randomPointInLocalBoundingBox(mesh, rand);
    return activeTransform(p, interpolateTransformation(path, time));
}



#endif // MESH_H_