#ifndef PATHING_H_
#define PATHING_H_

#include <boost/math/interpolators/cubic_hermite.hpp>

#include "vectormath.h"

/**
 * Initialize interpolators for a set of times, vectors, and the vectors' time derivatives
 * 
 * @param times List of times
 * @param vectors List of vectors
 * @param vectorDerivatives List of vector time derivatives
 * 
 * @return List of twice-differentiable interpolators for each of the three translation or orientation vector components
*/
template<class RandomAccessContainer, class VectorContainer>
std::array<boost::math::interpolators::cubic_hermite<RandomAccessContainer>, 3> constructVectorInterpolator(RandomAccessContainer times, VectorContainer vectors, VectorContainer vectorDerivatives){
    std::array<RandomAccessContainer, 3> x, v;
    if (std::empty(times)){
        times.push_back(0);
        for (int i = 0; i < 3; ++i){
            x[i].push_back(0);
            v[i].push_back(0);
        }
    }
    for (int i = 0; i < 3; ++i){ // transform list of 3D vectors into 3 lists of vector components
        for (auto vector: vectors) x[i].push_back(vector[i]);
        for (auto vectorDerivative: vectorDerivatives) v[i].push_back(vectorDerivative[i]);
    }
    if (std::size(times) == 1){
        times.push_back(times.back() + 1);
        for (int i = 0; i < 3; ++i){
            x[i].push_back(x[i].back() + v[i].back());
            v[i].push_back(v[i].back());
        }
    }
    return {boost::math::interpolators::cubic_hermite<RandomAccessContainer>(RandomAccessContainer(times), std::move(x[0]), std::move(v[0])),
            boost::math::interpolators::cubic_hermite<RandomAccessContainer>(RandomAccessContainer(times), std::move(x[1]), std::move(v[1])),
            boost::math::interpolators::cubic_hermite<RandomAccessContainer>(RandomAccessContainer(times), std::move(x[2]), std::move(v[2]))};
}


/**
 * Class to describe the path of a rigid body in 3D.
 * It smoothly interpolates between 3D translations, velocities, orientations, and angular velocities at specified times.
 * Rotations are applied in the local coordinate system (rotations are applied before translations).
*/
template<class RandomAccessContainer>
struct TPathInterpolator{
    using Interpolator = typename boost::math::interpolators::cubic_hermite<RandomAccessContainer>;
    std::array<Interpolator, 3> translationInterpolator; ///< twice-differentiable interpolators for 3D translation
    std::array<Interpolator, 3> rotationInterpolator; ///< twice-differentiable interpolators for 3D orientation

    /**
     * Creates translation and orientation interpolator between given times.
     * Rotations are applied in the local coordinate system (rotations applied first, then translation).
     * 
     * @param times List of times
     * @param translations List of translation vectors
     * @param velocities List of velocity vectors
     * @param rotations List of rotation (axis-angle) vectors
     * @param angularvelocities List of angular velocity vectors
    */
    template<class VectorContainer>
    TPathInterpolator(RandomAccessContainer times, VectorContainer translations, VectorContainer velocities, VectorContainer rotations, VectorContainer angularvelocities):
        translationInterpolator(constructVectorInterpolator(times, translations, velocities)),
        rotationInterpolator(constructVectorInterpolator(times, rotations, angularvelocities)){

    };



};

inline
TPathInterpolator<std::vector<double> > constructPathInterpolator(const std::string &pathParameters){
    std::istringstream str(pathParameters);
    std::vector<double> times;
    std::array<std::vector<std::array<double, 3> >, 4> vectors;
    str >> std::ws;
    while (not str.eof()){
        double time;
        str >> time;
        if (str){
            times.push_back(time);
        }
        for (int i = 0; i < 4; ++i){
            double x, y, z;
            char c1, c2;
            str >> c1 >> x >> y >> z >> c2;
            if (str and c1 == '(' and c2 == ')'){
                vectors[i].push_back({x, y, z});
                if (i > 1){
                    vectors[i].back() *= M_PI/180; // convert rotation and angular velocity from degrees to radians
                }
            }
            else{
                throw std::runtime_error("Invalid trajectory parameters '" + pathParameters + "'");
            }
        }
        str >> std::ws;
        
    }
    return TPathInterpolator(times, vectors[0], vectors[1], vectors[2], vectors[3]);
}




/**
 * Provide translation at specified time
 * 
 * @param time Time to interpolate
 * 
 * @return Translation vector
*/
template<class VectorInterpolator, typename Real>
std::array<Real, 3> interpolateVector(const VectorInterpolator &interpolator, Real time){
    auto t = std::clamp(time, interpolator[0].domain().first, interpolator[0].domain().second);
    std::array<Real, 3> x;
    std::transform(std::begin(interpolator), std::end(interpolator), std::begin(x), [t](auto interp){ return interp(t); });
    if (t != time){
        for (int i = 0; i < 3; ++i){
            x[i] += interpolator[i].prime(t)*(time - t);
        }
    }
    return x;
}


template<class VectorInterpolator, typename Real>
std::array<Real, 3> interpolateVectorDerivative(const VectorInterpolator &interpolator, Real time){
    auto t = std::clamp(time, interpolator[0].domain().first, interpolator[0].domain().second);
    std::array<Real, 3> v;
    std::transform(std::begin(interpolator), std::end(interpolator), std::begin(v), [t](auto interp){ return interp.prime(t); });
    return v;
}


/**
 * Provide translation at specified time
 * 
 * @param time Time to interpolate
 * 
 * @return Translation vector
*/
template<class Path, typename Real>
std::array<Real, 3> interpolateTranslation(const Path &path, Real time){
    return interpolateVector(path.translationInterpolator, time);
}


/**
 * Provide rotation in local coordinate system at specified time
 * 
 * @param time Time to interpolate
 * 
 * @return Rotation vector (axis-angle)
*/
template<class Path, typename Real>
std::array<Real, 3> interpolateRotation(const Path &path, const Real time){
    return interpolateVector(path.rotationInterpolator, time);
}

/**
 * Provide velocity at specified time
 * 
 * @param time Time to interpolate
 * 
 * @return Velocity vector
*/
template<class Path, typename Real>
std::array<Real, 3> interpolateVelocity(const Path &path, const Real time){
    return interpolateVectorDerivative(path.translationInterpolator, time);
}

/**
 * Provide angular velocity in local coordinate system at specified time
 * 
 * @param time Time to interpolate
 * 
 * @return Angular velocity vector
*/
template<class Path, typename Real>
std::array<Real, 3> interpolateAngularVelocity(const Path &path, Real time){
    return interpolateVectorDerivative(path.rotationInterpolator, time);
}

/**
 * Provide velocity relative to fixed position at specified time
 * 
 * @param time Time to interpolate
 * @param position Position at which to calculate rotational velocity
 * 
 * @return Velocity vector
*/
template<class Path, typename Real, class Vector>
Vector interpolateRelativeVelocity(const Path &path, Real time, const Vector &position){
    Vector x = interpolateTranslation(path, time);
    Vector v = interpolateVelocity(path, time);
    Vector a = interpolateAngularVelocity(path, time);
    return v + boost::qvm::cross(a, position - x);
}


template<class Path, typename Real>
std::pair<std::array<Real, 3>, std::array<Real, 3> > interpolateTransformation(const Path &path, Real time){
    return {interpolateTranslation(path, time), interpolateRotation(path, time)};
}


#endif // PATHING_H_
