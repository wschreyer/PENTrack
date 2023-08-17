/**
 * \file
 * Analytical B fields
 *
 */

#ifndef ANALYTICFIELDS_H_
#define ANALYTICFIELDS_H_

// #define PY_SSIZE_T_CLEAN
// #include <Python.h>
// #include <boost/python.hpp>

#include "field.h"


/**
 * B field described by:
 *
 * B_x = a1 * exp(- a2* x + a3) + c1
 *
 * B_y = y * a1 * a2 / 2 * exp(- a2* x + a3) + c2
 *
 * B_z = z * a1 * a2 / 2 * exp(- a2* x + a3) + c2
 *
 */
class TExponentialFieldX: public TField{
private:
  double a1, a2, a3, c1, c2;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla]
   * @param _a2
   * @param _a3
   * @param _c1 [Tesla]
   * @param _c2 [Tesla]
   */
  TExponentialFieldX(const double _a1, const double _a2, const double _a3, const double _c1, const double _c2);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};

/**
 * B field described by:
 *
 * B_z = a1*x + a2
 *
 */
class TLinearFieldZ: public TField{
private:
  double a1, a2;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla/meter]
   * @param _a2	[Tesla]
   */
  TLinearFieldZ(const double _a1, const double _a2);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};

/**
 * B field with quadratic gradient described by:
 *
 * B_z = a_1/2 z^2 + a_2 z + z0
 *
 * dBdz =  a_1 z + a_2
 *
 * where Bx, By, and derivatives satisfy Maxwell's equation
 */
class TB0GradZ: public TField{
private:
  double a1, a2, z0;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla/meter^2]
   * @param _a2	[Tesla/meter]
   * @param _z0 [Tesla]
   */
  TB0GradZ(const double _a1, const double _a2, const double _z0);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};

/**
 * B field with cubic gradient dependent on x described by:
 *
 * B_z = (a_1 x^2 + a_2 x + a3) z + z0
 *
 * dBdz = a_1 x^2 + a_2 x + a3
 *
 * where Bx, By, and derivatives satisfy Maxwell's equation
 */
class TB0GradX2: public TField{
private:
  double a1, a2, a3, z0;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla/meter^3]
   * @param _a2	[Tesla/meter^2]
   * @param _a3	[Tesla/meter]
   * @param _z0 [Tesla]
   */
  TB0GradX2(const double _a1, const double _a2, const double _a3, const double _z0);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};

/**
 * B field with quadratic gradient described by:
 *
 * B_z = a_1 xyz + a_2 z + z0
 *
 * dBdz =  a_1 xy + a_2
 *
 * where Bx, By, and derivatives satisfy Maxwell's equation
 */
class TB0GradXY: public TField{
private:
  double a1, a2, z0;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla/meter^3]
   * @param _a2	[Tesla/meter]
   * @param _z0 [Tesla]
   */
  TB0GradXY(const double _a1, const double _a2, const double _z0);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};

/**
 * B0 field with a component dependent on XY:
 *
 * B_z = a_1 xy + z0
 *
 * B_y = a_1 xz
 *
 * B_x = a_1 yz
 *
 * where Bx, By, and derivatives satisfy Maxwell's equation
 */
class TB0_XY: public TField{
private:
  double a1, z0;
public:
  /**
   *
   * Constructor requires:
   *
   * @param _a1 [Tesla/meter^2]
   * @param _z0 [Tesla]
   */
  TB0_XY(const double _a1, const double _z0);

  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const override;

  /**
   * Adds no electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override  {};
};


/**
 * Class calculating magnetic field from user-defined formulas
 * 
 * Field gradients are calculated numerically with a five-point stencil method
 */
class TCustomBField: public TField{
private:
  std::unique_ptr<double> tvar, xvar, yvar, zvar; ///< Variables used to evaluate formulas, they need to be pointers to make sure references stored in the exprtk expression do not get invalidated when copying
  std::array<exprtk::expression<double>, 3> Bexpr; ///< Formula interpreters, one for each field component
public:
  /**
   * Constructor
   * 
   * @param _Bx String containing formula for x component of field
   * @param _By String containing formula for y component of field
   * @param _Bz String containing formula for z component of field
   */
  TCustomBField(const std::string &_Bx, const std::string &_By, const std::string &_Bz);


  /**
   * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z and time t
   * 
   * Spatial derivatives are calculated numerically using a five-point stencil method
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t the time
   * @param B Returns magnetic-field components
   * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
   **/
  void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = nullptr) const override;

  /**
   * Returns empty electric field.
   * Required because we are inheriting from abstract class TField which has virtual void EField function.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Electric potential
   * @param Ei Electric field components
   **/
  void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override {};
};


// /**
//  * Class calculating magnetic field from magpylib
//  * 
//  * Field gradients are calculated numerically with a five-point stencil method
//  */
// class TMagpy: public TField{
// private:
//   // PyObject *pMagnetObject;
//   // PyObject *pBFieldFunc;
//   boost::python::object bpBFieldFunc;
//   boost::python::object bpMagnetObject;

// public:
//   /**
//    * Constructor
//    * 
//    * @ft string of the file name containing a function that returns a magnet magpy object.
//    */
//   TMagpy(const std::string ft);


//   /**
//    * Calculates B field B[3] and the derivatives dBidxj[3][3] for a given point x,y,z and time t
//    * 
//    * Spatial derivatives are calculated numerically using a five-point stencil method
//    *
//    * @param x Cartesian x coordinate
//    * @param y Cartesian y coordinate
//    * @param z Cartesian z coordinate
//    * @param t the time
//    * @param B Returns magnetic-field components
//    * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
//    **/
//   void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = nullptr) const override;

//   /**
//    * Returns empty electric field.
//    * Required because we are inheriting from abstract class TField which has virtual void EField function.
//    *
//    * @param x Cartesian x coordinate
//    * @param y Cartesian y coordinate
//    * @param z Cartesian z coordinate
//    * @param t Time
//    * @param V Electric potential
//    * @param Ei Electric field components
//    **/
//   void EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const override {};
// };


#endif /*ANALYTICFIELDS_H_*/
