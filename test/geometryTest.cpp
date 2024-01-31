#include <boost/test/unit_test.hpp>

#include <random>

#include "pathing.h"
#include "mesh.h"
#include "vectormath.h"

using namespace std;

template<class Vector1, class Vector2 = Vector1> void compareVector(Vector1 v1, Vector2 v2){
    BOOST_TEST_CONTEXT("v1 = " << v1[0] << " " << v1[1] << " " << v1[2] << ", v2 = " << v2[0] << " " << v2[1] << " " << v2[2]){
        BOOST_CHECK_SMALL(boost::qvm::mag(v2 - v1), 1e-13);
    }
}

BOOST_AUTO_TEST_CASE(PathTest){
    vector<double> times{0, 1, 2};
    vector<array<double, 3> > translations{{0,0,0}, {1,0,0}, {2,0,0}};
    vector<array<double, 3> > velocities{{1,0,0}, {0,0,0}, {1,0,0}};
    vector<array<double, 3> > rotations{{0,0,0}, {0,0,1}, {0,0,2}};
    vector<array<double, 3> > angularvelocities{{0,0,1}, {0,0,0}, {0,0,1}};
    array<double, 3> position{1, 0, 0};
    auto path = createPathInterpolator("0 (0 0 0) (1 0 0) (0 0 0) (0 0 57.295779513082320876798154814105) " 
                                       "1 (1 0 0) (0 0 0) (0 0 57.295779513082320876798154814105) (0 0 0) "
                                       "2 (2 0 0) (1 0 0) (0 0 114.59155902616464175359630962821) (0 0 57.295779513082320876798154814105)");

    // we start at origin with velocity in x direction and angular velocity around z
    compareVector(interpolateTranslation(path, 0.),                 {0,0,0});
    compareVector(interpolateRotation(path, 0.),                    {0,0,0});
    compareVector(interpolateVelocity(path, 0.),                    {1,0,0});
    compareVector(interpolateAngularVelocity(path, 0.),             {0,0,1});
    compareVector(interpolateRelativeVelocity(path, 0.,  position), {1,1,0});
    compareVector(activeRotate(    position, interpolateRotation(path, 0.)),       {1,0,0});
    compareVector(passiveRotate(   position, interpolateRotation(path, 0.)),       {1,0,0});
    compareVector(activeTransform( position, interpolateTransformation(path, 0.)), {1,0,0});
    compareVector(passiveTransform(position, interpolateTransformation(path, 0.)), {1,0,0});

    // velocities initially increase to reach new position at t = 1
    compareVector(interpolateTranslation(path, 0.5),                {0.625,0,0});
    compareVector(interpolateRotation(path, 0.5),                   {0,0,0.625});
    compareVector(interpolateVelocity(path, 0.5),                   {1.25,0,0});
    compareVector(interpolateAngularVelocity(path, 0.5),            {0,0,1.25});
    compareVector(interpolateRelativeVelocity(path, 0.5, position), {1.25,0.375*1.25,0});
    compareVector(activeRotate(    position,  interpolateRotation(path, 0.5)),       {std::cos(0.625),std::sin(0.625),0});
    compareVector(passiveRotate(   position,  interpolateRotation(path, 0.5)),       {std::cos(-0.625),std::sin(-0.625),0});
    compareVector(activeTransform( position,  interpolateTransformation(path, 0.5)), {std::cos(0.625) + 0.625,std::sin(0.625), 0}); // rotate, then translate
    compareVector(passiveTransform(position,  interpolateTransformation(path, 0.5)), {0.375*std::cos(-0.625),0.375*std::sin(-0.625),0}); // translate (1,0,0) -> (0.375,0,0), then rotate

    // at t = 1 we reach position (1,0,0) and rotation (0,0,1) with 0 velocities
    compareVector(interpolateTranslation(path, 1.),                 {1,0,0});
    compareVector(interpolateRotation(path, 1.),                    {0,0,1});
    compareVector(interpolateVelocity(path, 1.),                    {0,0,0});
    compareVector(interpolateAngularVelocity(path, 1.),             {0,0,0});
    compareVector(interpolateRelativeVelocity(path, 1.,  position), {0,0,0});
    compareVector(activeRotate(    position,  interpolateRotation(path, 1.)),       {std::cos(1),std::sin(1),0});
    compareVector(passiveRotate(   position,  interpolateRotation(path, 1.)),       {std::cos(-1),std::sin(-1),0});
    compareVector(activeTransform( position,  interpolateTransformation(path, 1.)), {std::cos(1) + 1,std::sin(1),0}); // rotate, then translate
    compareVector(passiveTransform(position,  interpolateTransformation(path, 1.)), {0,0,0}); // translate (1,0,0) -> (0,0,0), then rotate

    // at t = 2 we reach position (2,0,0) and rotation (0,0,2) with same velocity as in the beginning
    compareVector(interpolateTranslation(path, 2.),                 {2,0,0});
    compareVector(interpolateRotation(path, 2.),                    {0,0,2});
    compareVector(interpolateVelocity(path, 2.),                    {1,0,0});
    compareVector(interpolateAngularVelocity(path, 2.),             {0,0,1});
    compareVector(interpolateRelativeVelocity(path, 2.,  position), {1,-1,0});
    compareVector(activeRotate(    position,  interpolateRotation(path, 2.)),       {std::cos(2),std::sin(2),0});
    compareVector(passiveRotate(   position,  interpolateRotation(path, 2.)),       {std::cos(-2),std::sin(-2),0});
    compareVector(activeTransform( position,  interpolateTransformation(path, 2.)), {std::cos(2) + 2,std::sin(2),0}); // rotate, then translate
    compareVector(passiveTransform(position,  interpolateTransformation(path, 2.)), {-1*std::cos(-2),-1*std::sin(-2),0}); // translate (1,0,0) -> (-1,0,0), then rotate

    // path should be extrapolated with constant velocities beyond the time range (0..2)
    compareVector(interpolateTranslation(path, 3.),                 {3,0,0});
    compareVector(interpolateRotation(path, 3.),                    {0,0,3});
    compareVector(interpolateVelocity(path, 3.),                    {1,0,0});
    compareVector(interpolateAngularVelocity(path, 3.),             {0,0,1});
    compareVector(interpolateRelativeVelocity(path, 3.,  position), {1,-2,0});
    compareVector(activeRotate(    position,  interpolateRotation(path, 3.)),       {std::cos(3),std::sin(3),0});
    compareVector(passiveRotate(   position,  interpolateRotation(path, 3.)),       {std::cos(-3),std::sin(-3),0});
    compareVector(activeTransform( position,  interpolateTransformation(path, 3.)), {std::cos(3) + 3,std::sin(3),0}); // rotate, then translate
    compareVector(passiveTransform(position,  interpolateTransformation(path, 3.)), {-2*std::cos(-3),-2*std::sin(-3),0}); // translate (1,0,0) -> (-2,0,0), then rotate

    compareVector(interpolateTranslation(path, -1.),                 {-1,0,0});
    compareVector(interpolateRotation(path, -1.),                    {0,0,-1});
    compareVector(interpolateVelocity(path, -1.),                    {1,0,0});
    compareVector(interpolateAngularVelocity(path, -1.),             {0,0,1});
    compareVector(interpolateRelativeVelocity(path, -1.,  position), {1,2,0});
    compareVector(activeRotate(    position,  interpolateRotation(path, -1.)),       {std::cos(-1),std::sin(-1),0});
    compareVector(passiveRotate(   position,  interpolateRotation(path, -1.)),       {std::cos(1),std::sin(1),0});
    compareVector(activeTransform( position,  interpolateTransformation(path, -1.)), {std::cos(-1) - 1,std::sin(-1),0}); // rotate, then translate
    compareVector(passiveTransform(position,  interpolateTransformation(path, -1.)), {2*std::cos(1),2*std::sin(1),0}); // translate (1,0,0) -> (2,0,0), then rotate

    compareVector(activeTransform(passiveTransform({1, 0, 0}, interpolateTransformation(path, 2.)), interpolateTransformation(path, 2.)), {1,0,0}); // doing both active and passive transformation should give back the inital vector
}


BOOST_AUTO_TEST_CASE(movingMeshTest){
    random_device rd;
    mt19937 rng(rd());
    std::uniform_real_distribution<double> unidist(0, 1);
    auto path = createPathInterpolator("0 (0 0 0) (0 0 0) (0 0 0) (0 0 0) 1 (1 0 0) (0 0 0) (0 0 90) (0 0 0)");
    auto mesh = createMesh("test/UnitCubeCenteredOn(0,0,0).STL");
    for (int i = 0; i < 10; ++i){
        auto t = unidist(rng);
        auto pIn = randomPointInMovingVolume(mesh, t, rng, path);
        BOOST_CHECK(inMovingSolid(mesh, t, pIn, path));
        auto pInLocal = passiveTransform(pIn, interpolateTransformation(path, t));
        BOOST_CHECK(abs(pInLocal[0]) < 0.5 and abs(pInLocal[1]) < 0.5 and abs(pInLocal[2]) < 0.5);
        auto [pOn, n] = randomPointOnMovingSurface(mesh, t, rng, path);
        auto pOnLocal = passiveTransform(pOn, interpolateTransformation(path, t));
        BOOST_CHECK(abs(abs(pOnLocal[0]) - 0.5) < 1e-13 or abs(abs(pOnLocal[1]) - 0.5) < 1e-13 or abs(abs(pOnLocal[2]) - 0.5) < 1e-13);
        BOOST_CHECK(abs(pOnLocal[0]) <= 0.5 + 1e-13 and abs(pOnLocal[1]) <= 0.5 + 1e-13 and abs(pOnLocal[2]) <= 0.5 + 1e-13);
        auto nLocal = passiveRotate(n, interpolateRotation(path, t));
        BOOST_CHECK(abs(abs(nLocal[0]) - 1) < 1e-13 or abs(abs(nLocal[1]) - 1) < 1e-13 or abs(abs(nLocal[2]) - 1) < 1e-13);
        BOOST_CHECK_SMALL(boost::qvm::mag(n) - 1, 1e-13);

        for (int dim = 0; dim < 3; ++dim){
            std::array<double, 3> point{0,0,0};
            point[dim] += 1;
            point = activeTransform(point, interpolateTransformation(path, t));
            BOOST_CHECK(not inMovingSolid(mesh, t, point, path));
            auto collisions = findCollisionsWithMovingMesh(mesh, t, pIn, t, point, path);
            BOOST_CHECK_EQUAL(std::size(collisions), 1);
            if (not std::empty(collisions)){
                auto collpLocal = passiveTransform(collisions[0].collisionPoint, interpolateTransformation(path, t));
                BOOST_CHECK_SMALL(collpLocal[dim] - 0.5, 1e-13);
                BOOST_CHECK_SMALL(boost::qvm::mag(collisions[0].normal) - 1, 1e-13);
                compareVector(passiveRotate(collisions[0].normal, interpolateRotation(path, t)), passiveTransform(point, interpolateTransformation(path, t)));
                compareVector(collisions[0].velocity, interpolateRelativeVelocity(path, t, collisions[0].collisionPoint));
            }
        }

    }

}

