// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Stream_support/include/CGAL/IO/File_writer_OFF_impl.h $
// $Id: File_writer_OFF_impl.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/IO/File_writer_OFF.h>

namespace CGAL {

CGAL_INLINE_FUNCTION
void
File_writer_OFF::
write_header( std::ostream& o,
              std::size_t   vertices,
              std::size_t   halfedges,
              std::size_t   facets,
              bool          normals) {
    m_out = &o;
    m_header.set_vertices(  vertices);
    // Don't. This halfdges aren't trusted:
    // m_header.set_halfedges( halfedges);
    (void)halfedges;
    m_header.set_facets(  facets);
    m_header.set_normals( normals);
    // Print header.
    out() << m_header;
}
} //namespace CGAL
// EOF //
