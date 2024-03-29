// Copyright (c) 1999,2001  
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/STL_Extension/include/CGAL/Fourtuple.h $
// $Id: Fourtuple.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_FOURTUPLE_H
#define CGAL_FOURTUPLE_H

#include <CGAL/config.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Fourtuple.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS "Please use `CGAL::cpp11::array` instead."
#include <CGAL/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {

template < class T >
struct CGAL_DEPRECATED Fourtuple
{
  typedef T value_type;

  T  e0, e1, e2, e3;

  Fourtuple()
  {}

  Fourtuple(const T & a0, const T & a1, const T & a2, const T & a3)
    : e0(a0), e1(a1), e2(a2), e3(a3)
  {}
};

} //namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_FOURTUPLE_H
