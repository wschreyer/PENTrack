// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/HalfedgeDS/include/CGAL/HalfedgeDS_halfedge_max_base_with_id.h $
// $Id: HalfedgeDS_halfedge_max_base_with_id.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_halfedge_base.h>

namespace CGAL {

template < class Refs, class ID>
class HalfedgeDS_halfedge_max_base_with_id : public HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true >
{
public:
    typedef HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true> Base ;
    
    typedef typename Base::Base_base Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type mID ;
    
public:

    HalfedgeDS_halfedge_max_base_with_id( size_type i = size_type(-1) ) : mID(i) {}
    
    size_type&       id()       { return mID; }
    size_type const& id() const { return mID; }
};

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H //
// EOF //
