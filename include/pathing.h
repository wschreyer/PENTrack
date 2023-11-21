#ifndef PATHING_H_
#define PATHING_H_

#include <boost/math/interpolators/cubic_hermite.hpp>

#include "vectormath.h"

/**
 * Class to describe the path of a rigid body in 3D.
 * It smoothly interpolates between 3D translations, velocities, orientations, and angular velocities at specified times.
 * Rotations are applied in the local coordinate system (rotations are applied before translations).
*/
template<class RandomAccessContainer, class VectorContainer>
class TPath{
private:
    using Real = typename RandomAccessContainer::value_type;
    using Vector = typename VectorContainer::value_type;
    using interpolator = typename boost::math::interpolators::cubic_hermite<RandomAccessContainer>;
    std::array<interpolator, 3> translationInterpolator; ///< twice-differentiable interpolators for 3D translation
    std::array<interpolator, 3> rotationInterpolator; ///< twice-differentiable interpolators for 3D orientation

    /**
     * Initialize interpolators for translation or orientation
     * 
     * @param times List of times
     * @param translations List of translation or orientation (axis-angle) vectors for each time
     * @param velocities List of velocity or angular velocity vectors for each time
     * 
     * @return List of twice-differentiable interpolators for each of the three translation or orientation vector components
    */
    std::array<interpolator, 3> initInterpolator(RandomAccessContainer times, VectorContainer translations, VectorContainer velocities){
        if (std::empty(times)){
            return {interpolator({0}, {0}, {0}),
                    interpolator({0}, {0}, {0}),
                    interpolator({0}, {0}, {0})};
        }
        std::array<RandomAccessContainer, 3> x, v;
        for (int i = 0; i < 3; ++i){ // transform list of 3D vectors into 3 lists of vector components
            std::transform(translations.begin(), translations.end(), std::back_inserter(x[i]), [i](Vector ts){ return ts[i]; });
            std::transform(velocities.begin(), velocities.end(), std::back_inserter(v[i]), [i](Vector vs){ return vs[i]; });
        }
        return {interpolator(RandomAccessContainer(times), std::move(x[0]), std::move(v[0])), // we need to pass a copy of times because the interpolator constructor expects moveable parameters
                interpolator(RandomAccessContainer(times), std::move(x[1]), std::move(v[1])),
                interpolator(RandomAccessContainer(times), std::move(x[2]), std::move(v[2]))};
    }
public:
    /**
     * Constructor
     * 
     * Creates translation and orientation interpolator between given times.
     * Rotations are applied in the local coordinate system (rotations applied first, then translation).
     * 
     * @param times List of times
     * @param translations List of translation vectors
     * @param velocities List of velocity vectors
     * @param rotations List of rotation (axis-angle) vectors
     * @param angularvelocities List of angular velocity vectors
    */
    TPath(RandomAccessContainer times, VectorContainer translations, VectorContainer velocities, VectorContainer rotations, VectorContainer angularvelocities):
        translationInterpolator(initInterpolator(times, translations, velocities)), rotationInterpolator(initInterpolator(times, rotations, angularvelocities)){

    }

    /**
     * Provide translation at specified time
     * 
     * @param time Time to interpolate
     * 
     * @return Translation vector
    */
    Vector translation(Real time){
        auto t = std::clamp(time, translationInterpolator[0].domain().first, translationInterpolator[0].domain().second);
        Vector x{translationInterpolator[0](t), translationInterpolator[1](t), translationInterpolator[2](t)};
        if (t != time){
            Vector v{translationInterpolator[0].prime(t), translationInterpolator[1].prime(t), translationInterpolator[2].prime(t)};
            x = x + v*(time - t);
        }
        return x;
    }


    /**
     * Provide rotation in local coordinate system at specified time
     * 
     * @param time Time to interpolate
     * 
     * @return Rotation vector (axis-angle)
    */
    Vector rotation(Real time){
        auto t = std::clamp(time, rotationInterpolator[0].domain().first, rotationInterpolator[0].domain().second);
        Vector r{rotationInterpolator[0](t), rotationInterpolator[1](t), rotationInterpolator[2](t)};
        if (t != time){
            Vector a{rotationInterpolator[0].prime(t), rotationInterpolator[1].prime(t), rotationInterpolator[2].prime(t)};
            r = r + a*(time - t);
        }
        return r;
    }

    /**
     * Provide velocity at specified time
     * 
     * @param time Time to interpolate
     * 
     * @return Velocity vector
    */
    Vector velocity(Real time){
        auto t = std::clamp(time, translationInterpolator[0].domain().first, translationInterpolator[0].domain().second);        
        return Vector{translationInterpolator[0].prime(t), translationInterpolator[1].prime(t), translationInterpolator[2].prime(t)};
    }

    /**
     * Provide angular velocity in local coordinate system at specified time
     * 
     * @param time Time to interpolate
     * 
     * @return Angular velocity vector
    */
    Vector angularVelocity(Real time){
        auto t = std::clamp(time, rotationInterpolator[0].domain().first, rotationInterpolator[0].domain().second);        
        return Vector{rotationInterpolator[0].prime(t), rotationInterpolator[1].prime(t), rotationInterpolator[2].prime(t)};
    }

    /**
     * Provide velocity relative to fixed position at specified time
     * 
     * @param time Time to interpolate
     * @param position Position at which to calculate rotational velocity
     * 
     * @return Velocity vector
    */
    Vector relativeVelocity(Real time, Vector position){
        Vector x = translation(time);
        Vector v = velocity(time);
        Vector a = angularVelocity(time);
        return v + boost::qvm::cross(a, position - x);
    }

    /**
     * Transform position vector into rest frame of body moving along path (passive transformation)
     * 
     * @param time Time
     * @param position Position vector
    */
    Vector passive_transform(Real time, Vector position){
        auto x = translation(time);
        auto r = rotation(time);
        auto angle = boost::qvm::mag(r);
        if (angle == 0){
            return position - x;
        }
        else{
            return boost::qvm::rot_quat(r/angle, -angle)*(position - x);
        }
    }
};

#endif // PATHING_H_
