// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Surface_sweep_2/include/CGAL/Surface_sweep_2/Default_visitor.h $
// $Id: Default_visitor.h 0b5353c %aI Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_DEFAULT_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_DEFAULT_VISITOR_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the Default_visitor class-template.
 */

#include <CGAL/Surface_sweep_2/Default_visitor_base.h>
#include <CGAL/Surface_sweep_2/Default_event.h>
#include <CGAL/Surface_sweep_2/Default_subcurve.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Default_visitor
 *
 * An empty surface-sweep visitor that does little. It is used as a base-class
 * for other concrete visitors that produce some output. It is also used to
 * set default values for the allocator, event, and subcurve types.
 */
template <typename Visitor_,
          typename GeometryTraits_2,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          typename Event_ = Default_event<GeometryTraits_2, Allocator_>,
          typename Subcurve_ = Default_subcurve<GeometryTraits_2, Event_,
                                                Allocator_> >
class Default_visitor : public Default_visitor_base<GeometryTraits_2, Event_,
                                                    Subcurve_, Allocator_,
                                                    Visitor_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Allocator_                                    Allocator;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
