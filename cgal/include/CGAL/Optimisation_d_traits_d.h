// Copyright (c) 1997-2001  
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Optimisation_basic/include/CGAL/Optimisation_d_traits_d.h $
// $Id: Optimisation_d_traits_d.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_D_TRAITS_D_H
#define CGAL_OPTIMISATION_D_TRAITS_D_H

// includes
#  include <CGAL/Optimisation/Access_dimension_d.h>
#  include <CGAL/Optimisation/Access_coordinates_begin_d.h>
#  include <CGAL/Optimisation/Construct_point_d.h>

namespace CGAL {

// Class declaration
// =================
template < class K_, class ET_ = typename K_::RT,
                     class NT_ = typename K_::RT >
class Optimisation_d_traits_d;

// Class interface
// ===============
template < class K_, class ET_, class NT_>
class Optimisation_d_traits_d {
  public:
    // self
    typedef  K_                         K;
    typedef  ET_                        ET;
    typedef  NT_                        NT;
    typedef  Optimisation_d_traits_d<K,ET,NT>
                                        Self;

    // types
    typedef  typename K::Point_d        Point_d;

    typedef  typename K::Rep_tag        Rep_tag;

    typedef  typename K::RT             RT;
    typedef  typename K::FT             FT;

    typedef  CGAL::Access_dimension_d<K>      Access_dimension_d;
    typedef  CGAL::Access_coordinates_begin_d<K>
                                        Access_coordinates_begin_d;

    typedef  CGAL::_Construct_point_d<K>       Construct_point_d;

    // creation
    Optimisation_d_traits_d( ) { }
    Optimisation_d_traits_d( const Optimisation_d_traits_d<K_,ET_,NT_>&) {}

    // operations
    Access_dimension_d
    access_dimension_d_object( ) const
        { return Access_dimension_d(); }

    Access_coordinates_begin_d
    access_coordinates_begin_d_object( ) const
        { return Access_coordinates_begin_d(); }

    Construct_point_d
    construct_point_d_object( ) const
        { return Construct_point_d(); }
};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_D_TRAITS_D_H

// ===== EOF ==================================================================
