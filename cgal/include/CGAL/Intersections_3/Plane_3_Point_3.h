// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Intersections_3/include/CGAL/Intersections_3/Plane_3_Point_3.h $
// $Id: Plane_3_Point_3.h 057f4ea %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_3_PLANE_3_POINT_3_PLANE_3_H
#define CGAL_INTERSECTIONS_3_PLANE_3_POINT_3_PLANE_3_H

#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Point_3 &pt,
             const typename K::Plane_3 &plane,
             const K&)
{
  return plane.has_on(pt);
}

template <class K>
inline bool
do_intersect(const typename K::Plane_3 &plane,
             const typename K::Point_3 &pt,
             const K&)
{
  return plane.has_on(pt);
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_3, typename K::Plane_3>::result_type
intersection(const typename K::Point_3 &pt,
             const typename K::Plane_3 &plane,
             const K& k)
{
  if (do_intersect(pt,plane,k))
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Plane_3>(pt);
  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Plane_3>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Plane_3, typename K::Point_3>::result_type
intersection(const typename K::Plane_3 &plane,
             const typename K::Point_3 &pt,
             const K& k)
{
  return internal::intersection(pt, plane,k);
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Point_3, Plane_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Point_3, Plane_3, 3)

} //namespace CGAL
#endif // CGAL_INTERSECTIONS_3_PLANE_3_POINT_3_PLANE_3_H
