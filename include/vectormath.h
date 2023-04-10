#ifndef VECTORMATH_H_
#define VECTORMATH_H_

#include <vector>
#include <array>

#include <boost/qvm/all.hpp>

namespace boost { namespace qvm {

  /**
   * Implement boost::qvm vector traits for std::array<3> so we can use them with boost::qvm vector math methods
  */
  template <typename T>
  struct vec_traits<std::array<T, 3> >: vec_traits_defaults<std::array<T, 3>, T, 3>
  {

    template <int I>
    static inline T & write_element( std::array<T, 3> & v ) {
      return v[I];
    }

    static inline T & write_element_idx( int i, std::array<T, 3> & v ) {
      return v[i];
    } //optional

  };


  /**
   * Implement boost::qvm vector traits for std::array<4> so we can use them with boost::qvm vector math methods
  */
  template <typename T>
  struct vec_traits<std::array<T, 4> >: vec_traits_defaults<std::array<T, 4>, T, 4>
  {

    template <int I>
    static inline T & write_element( std::array<T, 4> & v ) {
      return v[I];
    }

    static inline T & write_element_idx( int i, std::array<T, 4> & v ) {
      return v[i];
    } //optional

  };

  /**
   * Implement boost::qvm vector traits for std::vector (assuming 3 elements) so we can use them with boost::qvm vector math methods
  */
  template <typename T>
  struct vec_traits<std::vector<T> >: vec_traits_defaults<std::vector<T>, T, 3>
  {

    template <int I>
    static inline T & write_element( std::vector<T> & v ) {
      return v[I];
    }

    static inline T & write_element_idx( int i, std::vector<T> & v ) {
      return v[i];
    } //optional

  };

} }


using namespace boost::qvm;

/**
 * Construct a vector orthogonal to given vector.
 * Works with any class that has boost::qvm vector traits implemented.
 * 
 * @param v A vector
 * 
 * @returns A vector orthogonal to v
*/
template<class Vector>
Vector constructoOrthogonalVector(const Vector &v){
    return {std::copysign(v[2], v[0]), std::copysign(v[2], v[1]), -std::copysign(v[0], v[2]) - std::copysign(v[1], v[2])};
}


/**
 * Rotate a vector.
 *
 * Rotate vector into new coordinate system defined by new z axis "z" and a vector "xz" in the new x-z plane (active transformation).
 * Works with any class that has boost::qvm vector traits implemented.
 * 
 * @param v Vector to be rotated
 * @param z Z axis of new coordinate system
 * @param xz Vector in x-z plane of new coordinate system. If this vector is parallel to z an arbitrary vector orthogonal z is used instead.
 * 
 * @returns Vector rotated into new coordinate system (active transformation)
 */
template<class Vector1, class Vector2 = Vector1, class Vector3 = Vector1>
void RotateVector(Vector1 &v, const Vector2 z, const Vector3 xz){
    if (z[0] == 0 and z[1] == 0 and xz[1] == 0){
        if (z[2] == 0)
            throw std::runtime_error("RotateVector called with null vector as parameter!");
        return; // if vector z is parallel to z axis and xz is in x-z plane we don't have to do anything
    }
    typedef typename vec_traits<Vector2>::scalar_type T;
    mat<T, 3, 3> R;
    col<1>(R) = boost::qvm::cross(z, xz); // vector parallel to new y axis
    if (col<1>(R) == zero_vec<T, 3>()){ // if z and xz are parallel their cross product is zero, and we can just use any vector orthogonal to z as a new y axis
        col<1>(R) = constructoOrthogonalVector(z); // construct a vector orthogonal to z
    }
    normalize(col<1>(R));
    col<0>(R) = normalized(cross(col<1>(R), z));
    col<2>(R) = normalized(z);
    vref(v) = R*v;
}

/**
 * Rotate a vector.
 *
 * Rotate vector into new coordinate system defined by new z axis "z" (active transformation).
 * Works with any class that has boost::qvm vector traits implemented.
 * 
 * @param v Vector to be rotated
 * @param z Z axis of new coordinate system
 * 
 * @returns Vector rotated into new coordinate system (active transformation)
 */
template<class Vector1, class Vector2 = Vector1>
void RotateVector(Vector1 &v, const Vector2 z){
    RotateVector(v, z, z);
}


/**
 * Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta.
 *
 * Adapted from ROOT
 * Works with any class that has boost::qvm vector traits implemented
 *
 * @param beta Vector v/c of moving reference frame
 * @param p Four-vector
 */
template<class Vector3, class Vector4>
void BOOST(const Vector3 &beta, Vector4 &p){
   //Boost this Lorentz vector (copy&paste from ROOT)
   auto b2 = mag_sqr(beta);
   auto gamma = 1.0 / std::sqrt(1.0 - b2);
   auto bp = dot(beta, YZW(p));
   auto gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

   YZW(p) = YZW(p) + gamma2*bp*beta + gamma*beta*X(p);
   X(p) = gamma*(X(p) + bp);
}


#endif /*VECTORMATH_H_*/