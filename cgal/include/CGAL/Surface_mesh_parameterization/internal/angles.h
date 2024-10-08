// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Surface_mesh_parameterization/include/CGAL/Surface_mesh_parameterization/internal/angles.h $
// $Id: angles.h f2b1124 %aI Mael Rouxel-Labbé
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ANGLES_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ANGLES_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_type_config.h>

#include <cmath>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

//                                                                           ->  ->
// Returns the cotangent of the corner (P,Q,R) (i.e. the cotan of the angle (QP, QR) ).
template<typename K>
typename K::FT cotangent(const typename K::Point_3& P,
                         const typename K::Point_3& Q,
                         const typename K::Point_3& R)
{
  typedef typename K::FT                   NT;
  typedef typename K::Vector_3             Vector_3;

  Vector_3 u = P - Q;
  Vector_3 v = R - Q;
  NT dot = (u * v);
  Vector_3 cross_vector = CGAL::cross_product(u, v);
  NT cross_norm = CGAL::sqrt(cross_vector * cross_vector);
  if(cross_norm != NT(0))
    return (dot / cross_norm);
  else
    return 0; // undefined
}

//                                                                       ->  ->
// Returns the tangent of the corner (P,Q,R) (i.e. the tangent of angle (QP, QR) ).
template<typename K>
typename K::FT tangent(const typename K::Point_3& P,
                       const typename K::Point_3& Q,
                       const typename K::Point_3& R)
{
  typedef typename K::FT                 NT;
  typedef typename K::Vector_3           Vector_3;

  Vector_3 u = P - Q;
  Vector_3 v = R - Q;
  NT dot = (u * v);
  Vector_3 cross_vector = CGAL::cross_product(u, v);
  NT cross_norm = CGAL::sqrt(cross_vector * cross_vector);
  if(dot != NT(0))
    return (cross_norm / dot);
  else
    return 0; // undefined
}

// Fixes the sine to be within [-1;1].
template<typename K>
typename K::FT fix_sine(typename K::FT sine)
{
  if(sine >= 1)
    return 1;
  else if(sine <= -1)
    return -1;
  else
    return sine;
}

//                                                                       ->  ->
// Returns the angle (in radians) of the corner (P,Q,R) (i.e. the angle (QP, QR) ).
template<typename K>
typename K::FT compute_angle_rad(const typename K::Point_3& P,
                                 const typename K::Point_3& Q,
                                 const typename K::Point_3& R)
{
  typedef typename K::FT                           NT;
  typedef typename K::Vector_3                     Vector_3;

  Vector_3 u = P - Q;
  Vector_3 v = R - Q;

  NT product = CGAL::sqrt(u * u) * CGAL::sqrt(v * v);
  if(product == NT(0))
    return 0;

  // cosine
  NT dot = (u * v);
  NT cosine = dot / product;

  // sine
  Vector_3 w = CGAL::cross_product(u, v);
  NT abs_sine = CGAL::sqrt(w * w) / product;

  if(cosine >= NT(0))
    return std::asin(fix_sine<K>(abs_sine));
  else
    return CGAL_PI - std::asin(fix_sine<K>(abs_sine));
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ANGLES_H
