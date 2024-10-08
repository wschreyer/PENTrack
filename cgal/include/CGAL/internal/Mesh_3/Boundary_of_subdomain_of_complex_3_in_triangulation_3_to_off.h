// Copyright (c) 2008-2014  GeometryFactory (France).
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Mesh_3/include/CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h $
// $Id: Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h e5cbe7e %aI Mael Rouxel-Labbé
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H
#define CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>

namespace CGAL {

namespace internal {

template<typename PointContainer, typename FaceContainer>
std::ostream&
output_polygon_soup_to_off(const PointContainer& points,
                           const FaceContainer& faces,
                           std::ostream& out)
{
  typedef typename FaceContainer::value_type             Face;

  const std::size_t np = points.size();
  const std::size_t nf = faces.size();

  out << "OFF\n" << np << " " << nf << " 0\n";

  for(std::size_t i=0; i<np; ++i)
    out << points[i] << '\n';

  for(std::size_t i=0; i<nf; ++i)
  {
    const Face& f = faces[i];
    const std::size_t fs = f.size();
    out << fs;
    for(std::size_t j=0; j<fs; ++j)
      out << " " << f[j];
    out << '\n';
  }

  return out;
}

template <typename C3T3>
std::ostream&
output_boundary_of_c3t3_to_off(const C3T3& c3t3,
                               typename C3T3::Subdomain_index sd_index,
                               std::ostream& out,
                               bool normals_point_outside_of_the_subdomain = true)
{
  typedef typename C3T3::Triangulation::Geom_traits::Point_3             Point;
  typedef std::vector<std::size_t>                                       Face;

  std::vector<Point> points;
  std::vector<Face> faces;

  CGAL::Mesh_3::internal::facets_in_complex_3_to_triangle_soup(c3t3, sd_index, points, faces, normals_point_outside_of_the_subdomain);

  return output_polygon_soup_to_off(points, faces, out);
}

template <typename C3T3>
std::ostream&
output_facets_in_complex_to_off(const C3T3& c3t3,
                                std::ostream& out)
{
  typedef typename C3T3::Triangulation::Geom_traits::Point_3             Point;
  typedef std::vector<std::size_t>                                       Face;

  std::vector<Point> points;
  std::vector<Face> faces;

  CGAL::Mesh_3::internal::facets_in_complex_3_to_triangle_soup(c3t3, points, faces);

  return output_polygon_soup_to_off(points, faces, out);
}

} } // end of namespace CGAL::internal

#endif // CGAL_INTERNAL_MESH_3_BOUNDARY_OF_SUDDOMAIN_OF_COMPLEX_3_IN_TRIANGULATION_3_TO_OFF_H
