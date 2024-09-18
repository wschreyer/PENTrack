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

/**
 * Structure returned by TMovingMesh::findCollisions.
 */
template<class Vector, typename Real = typename Vector::value_type>
struct collision{
    Real collisionTime;
    Vector collisionPoint;
    Vector surfaceNormal; ///< normal (length = 1) of intersected surface
    Vector surfaceVelocity; ///< velocity of surface
};


/**
 * Class containing a triangle mesh and an AABB tree for fast intersection tests.
 */
template<typename Real>
struct TMesh{
    typedef Real RealType;
    typedef typename CGAL::Simple_cartesian<Real> CKernel; ///< Geometric Kernel used for CGAL types
    typedef typename CKernel::Point_3 CPoint; ///< CGAL point type
    typedef typename CKernel::Vector_3 CVector; ///< CGAL vector type
    typedef typename CKernel::Segment_3 CSegment; ///< CGAL segment type
    typedef typename CKernel::Ray_3 CRay; ///< CGAL ray type
    typedef typename CGAL::Surface_mesh<CPoint> CMesh; ///< CGAL triangle mesh type
    typedef typename CGAL::AABB_face_graph_triangle_primitive<CMesh> CPrimitive; ///< CGAL triangle type contained in AABB tree
    typedef typename CGAL::AABB_traits<CKernel, CPrimitive> CTraits; ///< CGAL triangle traits type
    typedef typename CGAL::AABB_tree<CTraits> CTree; ///< CGAL AABB tree type containing CPrimitives
    typedef typename boost::optional<typename CTree::Intersection_and_primitive_id<CSegment>::Type > CIntersection; ///< CGAL type returned by intersection test between mesh and a segment
    typedef typename boost::graph_traits<CMesh>::face_descriptor CFace; ///< CGAL mesh face type

    std::unique_ptr<CMesh> mesh; ///< Triangle mesh
    std::unique_ptr<CTree> tree; ///< Axis-aligned bounding-box tree for fast intersection search
    typename CMesh::Property_map<CFace, std::size_t> fccmap; ///< Property map assigning each mesh face to a connected component
    typename CMesh::Property_map<CFace, CVector> normalmap; ///< Property map assigning each mesh face a normal
    std::discrete_distribution<size_t> triangle_sampler; ///< Probability distribution to randomly sample triangles from mesh weighted by their areas.
    std::string name; ///< Name of the mesh
    Real area; ///< Area of the mesh
    Real volume; ///< Volume of the mesh
};


/**
 * Construct a TMesh by reading from a file.
 * 
 * Supports STL and several other file types, see
 * https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__IO__grp.html
 * 
 * @param filename File name to read from
 */
TMesh<double> constructMesh(const std::string &filename);


/**
 * Find collisions between a path segment and a mesh
 * 
 * @param mesh Class describing the mesh
 * @param segmentStart Starting point of segment
 * @param segmentEnd End point of segment
 * 
 * @returns List of a pair of vectors (collision point and normal of the collision surface)
 */
template<class Mesh, class Vector>
std::vector<std::pair<Vector, Vector> > findCollisionsWithMesh(const Mesh &mesh, const Vector &segmentStart,const Vector &segmentEnd){
    typename Mesh::CSegment segment(typename Mesh::CPoint(boost::qvm::X(segmentStart), boost::qvm::Y(segmentStart), boost::qvm::Z(segmentStart)),
                                    typename Mesh::CPoint(boost::qvm::X(segmentEnd), boost::qvm::Y(segmentEnd), boost::qvm::Z(segmentEnd)));
    std::vector<std::pair<Vector, Vector> > colls;
    std::vector<typename Mesh::CIntersection> out;
    mesh.tree->all_intersections(segment, std::back_inserter(out)); // search intersections of segment with mesh
    for (auto i: out){
        auto collisionPoint = boost::get<typename Mesh::CPoint>(&(i->first));
        if (collisionPoint) { // if intersection is a point
            auto normal = mesh.normalmap[i->second];//CGAL::Polygon_mesh_processing::compute_face_normal(i->second, *mesh.mesh);
            colls.push_back(std::make_pair(Vector{collisionPoint->x(), collisionPoint->y(), collisionPoint->z()}, Vector{normal.x(), normal.y(), normal.z()})); // add collision to list
        }
        else
            throw std::runtime_error("Segment-triangle intersection happened to not be a point");
    }
    return colls;
}


/**
 * Find collisions between a path segment and a mesh moving along a path
 * 
 * @param mesh Class describing the mesh
 * @param timeStart Starting time of segment
 * @param segmentStart Starting point of segment
 * @param timeEnd End time of segment
 * @param segmentEnd End point of segment
 * @param path Class describing the path along which the mesh is moving
 * 
 * @returns List of collision classes
 */
template<class Mesh, typename Real = typename Mesh::Real, class Vector, class Path>
std::vector<collision<Vector, Real> > findCollisionsWithMovingMesh(const Mesh &mesh, const Real timeStart, const Vector &segmentStart, const Real timeEnd, const Vector &segmentEnd, const Path &path){
    auto segmentStartLocal = passiveTransform(segmentStart, interpolateTransformation(path, timeStart));
    auto segmentEndLocal = passiveTransform(segmentEnd, interpolateTransformation(path, timeEnd));
    auto collisions = findCollisionsWithMesh(mesh, segmentStartLocal, segmentEndLocal);
    std::vector<collision<Vector, Real> > movingCollisions;
    for (auto coll: collisions){
        auto s = boost::qvm::dot(coll.first - segmentStartLocal, segmentEndLocal - segmentStartLocal)/boost::qvm::mag(segmentEndLocal - segmentStartLocal);
        auto collisionTime = timeStart + s*(timeEnd - timeStart);
        auto collisionPoint = activeTransform(coll.first, interpolateTransformation(path, collisionTime));
        movingCollisions.push_back({collisionTime, collisionPoint, activeRotate(coll.second, interpolateRotation(path, collisionTime)), interpolateRelativeVelocity(path, collisionTime, collisionPoint)});
    }
    return movingCollisions;
}


/**
 * Test if point is inside a mesh moving along a path
 * 
 * @param mesh Class describing the mesh
 * @param time Point in time
 * @param position Point in space
 * @param path Class descirbing the path along which the mesh is moving
 * 
 * @returns True if point is contained within mesh
*/
template<class Mesh, typename Real, class Vector, class Path>
bool inMovingSolid(const Mesh &mesh, const Real time, const Vector &position, const Path &path){
    return inSolid(mesh, passiveTransform(position, interpolateTransformation(path, time)));
}


/**
 * Test if point is inside a mesh
 * 
 * @param mesh Class describing the mesh
 * @param position Point in space
 * 
 * @returns True if point is contained within mesh
 */
template<class Mesh, class Vector>
bool inSolid(const Mesh &mesh, const Vector &position){
    return mesh.tree->number_of_intersected_primitives(typename Mesh::CRay(typename Mesh::CPoint(position[0], position[1], position[2]), typename Mesh::CVector(0.,0.,1.))) % 2 != 0;
}


/**
 * Generate random points homogeneously distributed on a surface mesh
 * 
 * @param mesh Class describing the mesh
 * @param rand Random number generator class
 * 
 * @returns A pair of vectors (random position on mesh surface and surface normal at this point)
 */
template<class Mesh, class RandomGenerator, typename Real = typename Mesh::RealType>
std::pair<std::array<Real, 3>, std::array<Real, 3> > randomPointOnSurface(Mesh &mesh, RandomGenerator &rand){
    auto faceidx = *(mesh.mesh->faces_begin() + mesh.triangle_sampler(rand));
    Real u(std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rand));
    Real v(std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rand));
    if (u + v >= 1){
        u = 1 - u;
        v = 1 - v;
    }
    auto pp = CGAL::Polygon_mesh_processing::construct_point(std::make_pair(faceidx, std::array<Real, 3>{u, v, 1 - u - v}), *mesh.mesh);
    auto nv = mesh.normalmap[faceidx]; //CGAL::Polygon_mesh_processing::compute_face_normal(faceidx, *mesh.mesh);
    return {{pp.x(), pp.y(), pp.z()}, {nv.x(), nv.y(), nv.z()}};
}


/**
 * Generate random points homogeneously distributed on a surface mesh moving along a path
 * 
 * @param mesh Class describing the mesh
 * @param time Time when the point should be generated
 * @param rand Random number generator class
 * @param path Class describing the path along which the mesh is moving
 * 
 * @returns A pair of vectors (random position on mesh surface and surface normal at this point)
 */
template<class Mesh, typename Real = typename Mesh::RealType, class Path, class RandomGenerator>
std::pair<std::array<Real, 3>, std::array<Real, 3> > randomPointOnMovingSurface(Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto [position, normal] = randomPointOnSurface(mesh, rand);
    position = activeTransform(position, interpolateTransformation(path, time));
    normal = activeRotate(normal, interpolateRotation(path, time));
    return std::make_pair(position, normal);
}


/**
 * Generate random points homogeneously distributed inside a mesh
 * 
 * @param mesh Class describing the mesh
 * @param rand Random number generator class
 * 
 * @returns Random point contained within the mesh
 */
template<class Mesh, class RandomGenerator> 
std::array<double, 3> randomPointInVolume(const Mesh &mesh, RandomGenerator &rand){
    std::array<double, 3> p;
    do{
        p = randomPointInBoundingBox(mesh, rand);
    }while(not inSolid(mesh, p));
    return p;
}


/**
 * Generate random points homogeneously distributed inside the bounding box of a mesh
 * 
 * @param mesh Class describing the mesh
 * @param rand Random number generator class
 * 
 * @returns Random point contained within the ounding box of the mesh
 */
template<class Mesh, class RandomGenerator>
std::array<double, 3> randomPointInBoundingBox(const Mesh &mesh, RandomGenerator &rand){
    auto bbox = mesh.tree->bbox();
    std::uniform_real_distribution<double> unidist(0, 1);
    return {bbox.xmin() + unidist(rand)*(bbox.xmax() - bbox.xmin()),
            bbox.ymin() + unidist(rand)*(bbox.ymax() - bbox.ymin()),
            bbox.zmin() + unidist(rand)*(bbox.zmax() - bbox.zmin())};
}


/**
 * Generate random points homogeneously distributed inside a mesh moving along a path
 * 
 * @param mesh Class describing the mesh
 * @param time Time when the point should be generated
 * @param rand Random number generator class
 * @param path Class describing the path along which the mesh is moving
 * 
 * @returns Random point contained within the mesh
 */
template<class Mesh, typename Real = typename Mesh::RealType, class RandomGenerator, class Path> 
std::array<Real, 3> randomPointInMovingVolume(const Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto bbox = mesh.tree->bbox();
    std::array<Real, 3> p;
    do{
        p = randomPointInBoundingBox(mesh, rand);
    }while(not inSolid(mesh, p));
    return activeTransform(p, interpolateTransformation(path, time));
}


/**
 * Generate random points homogeneously distributed inside the bounding box of a mesh
 * 
 * @param mesh Class describing the mesh
 * @param time Time when the point should be generated
 * @param rand Random number generator class
 * @param path Class describing the path along which the mesh is moving
 * 
 * @returns Random point contained within the ounding box of the mesh
 */
template<class Mesh, typename Real = typename Mesh::RealType, class RandomGenerator, class Path>
std::array<Real, 3> randomPointInMovingBoundingBox(const Mesh &mesh, const Real time, RandomGenerator &rand, const Path &path){
    auto p = randomPointInLocalBoundingBox(mesh, rand);
    return activeTransform(p, interpolateTransformation(path, time));
}


/**
 * Check if segment is (partially) contained in the bounding box of a mesh
 * 
 * @param mesh Class describing the mesh
 * @param segmentStart Starting point of segment
 * @param segmentEnd End point of segment
 * 
 * @returns True if segment is (partially) contained in the bounding box of the mesh
 */
template<class Mesh, class Vector>
bool inBoundingBox(const Mesh &mesh, const Vector &segmentStart, const Vector &segmentEnd){
    return CGAL::do_intersect(mesh.tree->bbox(), typename Mesh::CSegment(typename Mesh::CPoint(boost::qvm::X(segmentStart), boost::qvm::Y(segmentStart), boost::qvm::Z(segmentStart)),
                                                                         typename Mesh::CPoint(boost::qvm::X(segmentEnd), boost::qvm::Y(segmentEnd), boost::qvm::Z(segmentEnd))));
}


/**
 * Check if segment is (partially) contained in the bounding box of a mesh moving along a path
 * 
 * @param mesh Class describing the mesh
 * @param timeStart Starting time of segment
 * @param segmentStart Starting point of segment
 * @param timeEnd End time of segment
 * @param segmentEnd End point of segment
 * @param path Class describing the path along which the mesh is moving
 * 
 * @returns True if segment is (partially) contained in the bounding box of the mesh
 */
template<class Mesh, typename Real = typename Mesh::RealType, class Vector, class Path>
bool inMovingBoundingBox(const Mesh &mesh, const Real timeStart, const Vector &segmentStart, const Real timeEnd, const Vector &segmentEnd, const Path &path){
    return inBoundingBox(mesh, passiveTransform(segmentStart, interpolateTransformation(path, timeStart)), passiveTransform(segmentEnd, interpolateTransformation(path, timeEnd)));
}


/**
 * Check if a bounding box overlaps the bounding box of a mesh
 * 
 * @param mesh Class describing the mesh
 * @param boundingBox Pair of vectors defining opposite corners (minimum and maximum coordinates) of a bounding box 
 * 
 * @returns True if the bounding boxes are overlapping
 */
template<class Mesh, class Vector>
bool overlappingBoundingBoxes(const Mesh &mesh, const std::pair<Vector, Vector> &boundingBox){
    CGAL::Bbox_3 bbox(boundingBox.first[0], boundingBox.first[1], boundingBox.first[2], boundingBox.second[0], boundingBox.second[1], boundingBox.second[2]);
    return CGAL::do_overlap(mesh.tree->bbox(), bbox);
}


#endif // MESH_H_