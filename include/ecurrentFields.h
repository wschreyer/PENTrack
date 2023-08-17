/**
 * \file
 * Electric current B fields
 *
 */


// #include <cmath>
#include <vector>

#include "field.h"


#ifndef PENTRACK_ELECTRICCURRENTFIELDS_H_
#define PENTRACK_ELECTRICCURRENTFIELDS_H_


class TECurrentField: public TField{

public:
  struct WireSegment {
    double x1, y1, z1;  // Starting point of the wire segment
    double x2, y2, z2;  // Ending point of the wire segment
    // double current;     // Current flowing through the wire segment
  };

  /**
   * Constructor
   * 
   * @ft string of the file name containing a function that returns a magnet magpy object.
   */
  
  TECurrentField(const std::string ft, const std::string &_I);
  
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
  
private:
  std::vector<WireSegment> wireSegments;
  std::unique_ptr<double> tvar; ///< Variables used to evaluate formulas, they need to be pointers to make sure references stored in the exprtk expression do not get invalidated when copying
  exprtk::expression<double> Iexpr; ///< Formula interpreters for the electrical current
};

#endif /*PENTRACK_ELECTRICCURRENTFIELDS_H_*/
