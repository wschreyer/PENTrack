// Copyright (c) 2001  
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Kernel_d/include/CGAL/Kernel_d/PointHd_impl.h $
// $Id: PointHd_impl.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_POINTHD_C
#define CGAL_POINTHD_C
namespace CGAL {
#define PointHd PointHd2

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Vector res = t(vector_rep());
  return PointHd<RT,LA>(dimension(),res.begin(),res.end()); }

template <class RT, class LA>
VectorHd<RT,LA> PointHd<RT,LA>::operator-(const Origin&) const 
{ return VectorHd<RT,LA>(Base(*this)); }

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator+(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr()->homogeneous_add(ptr(), v.ptr()); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA> PointHd<RT,LA>::operator-(const VectorHd<RT,LA> &v) const
{ PointHd<RT,LA> res(dimension()); 
  res.ptr()->homogeneous_sub(ptr(), v.ptr()); 
  return res; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator+= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr()->homogeneous_add(old.ptr(), v.ptr()); 
  return *this; 
}

template <class RT, class LA>
PointHd<RT,LA>& PointHd<RT,LA>::operator-= (const VectorHd<RT,LA>& v)
{ int d = dimension(); 
  PointHd<RT,LA> old(*this); 
  *this = PointHd<RT,LA>(d); 
  ptr()->homogeneous_sub(old.ptr(), v.ptr()); 
  return *this; 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, PointHd<RT,LA>& p)
{ p.copy_on_write(); p.ptr()->read(I); 
  CGAL_assertion_msg((p.homogeneous(p.dimension()) > 0), 
  "operator>>: denominator of point must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const PointHd<RT,LA>& p)
{ p.ptr()->print(O,"PointHd"); return O; } 

#undef PointHd 
} //namespace CGAL
#endif // CGAL_POINTHD_C
