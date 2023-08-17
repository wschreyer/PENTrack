/**
 * \file
 * Python B fields
 *
 */

#ifndef PENTRACK_PYTHONFIELDS_H_
#define PENTRACK_PYTHONFIELDS_H_


// // #define PY_SSIZE_T_CLEAN
// // #include <Python.h>
#include <boost/python.hpp>

#include "field.h"


/**
 * Class calculating magnetic field from magpylib
 * 
 * Field gradients are calculated numerically with a five-point stencil method
 */
class TPythonField: public TField{
private:
  // PyObject *pMagnetObject;
  // PyObject *pBFieldFunc;
  boost::python::object bpBFieldFunc;
  boost::python::object buildSourceFunc;
  boost::python::object bpSourceObject;
  bool temporal;
public:
  /**
   * Constructor
   * 
   * @ft string of the file name containing a function that returns a magnet magpy object.
   */
  TPythonField(const std::string ft, const bool tp);


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

#endif /*PENTRACK_PYTHONFIELDS_H_*/
