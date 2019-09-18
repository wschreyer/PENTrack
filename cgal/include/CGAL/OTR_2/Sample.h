// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14/Optimal_transportation_reconstruction_2/include/CGAL/OTR_2/Sample.h $
// $Id: Sample.h ee57fc2 %aI Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando de Goes, Pierre Alliez, Ivo Vigan, Clément Jamin

#ifndef CGAL_OTR2_SAMPLE_H
#define CGAL_OTR2_SAMPLE_H

#include <CGAL/license/Optimal_transportation_reconstruction_2.h>



/// \cond SKIP_IN_MANUAL

namespace CGAL {
namespace OTR_2 {

template <class Traits_>
class Sample
{
public:
  typedef typename Traits_::FT FT;
  typedef typename Traits_::Point_2 Point;

private:
  Point m_point;
  FT m_mass;

  FT m_dist2_to_edge;
  FT m_coordinate;

  FT m_backup_dist2;
  FT m_backup_coord;

public:
  Sample(const Point& point,
      const FT mass = FT(1))
  : m_point(point),
    m_mass(mass),
    m_dist2_to_edge(0),
    m_coordinate(0),
    m_backup_dist2(0),
    m_backup_coord(0)
  {
  }

  Sample(const Sample& sample)
  : m_point(sample.point()),
    m_mass(sample.mass()),
    m_dist2_to_edge(0),
    m_coordinate(0),
    m_backup_dist2(0),
    m_backup_coord(0)
  {
  }

  ~Sample() { }

  const Point& point() const { return m_point; }
  Point& point() { return m_point; }

  const FT& mass() const { return m_mass; }
  FT& mass() { return m_mass; }

  const FT& distance2() const { return m_dist2_to_edge; }
  FT& distance2() { return m_dist2_to_edge; }

  const FT& coordinate() const { return m_coordinate; }
  FT& coordinate() { return m_coordinate; }

  void backup()
  {
    m_backup_dist2 = m_dist2_to_edge;
    m_backup_coord = m_coordinate;
  }

  void restore()
  {
    m_dist2_to_edge = m_backup_dist2;
    m_coordinate = m_backup_coord;
  }
};

template <class Sample_>
class Sample_with_priority
{
public:
  typedef typename Sample_::FT FT;

private:
  Sample_* m_sample;
  FT m_priority;

public:
  Sample_with_priority(Sample_* sample, const FT priority = FT(0))
  {
    m_sample   = sample;
    m_priority = priority;
  }

  Sample_with_priority(const Sample_with_priority& psample)
  {
    m_sample   = psample.sample();
    m_priority = psample.priority();
  }

  ~Sample_with_priority() { }

  Sample_with_priority& operator = (const Sample_with_priority& psample)
  {
    m_sample   = psample.sample();
    m_priority = psample.priority();
    return *this;
  }

  Sample_* sample() const { return m_sample; }

  const FT priority() const { return m_priority; }
};


template <class T>
struct greater_priority
{
  bool operator() (const T& a, const T& b) const
  {
    return ( a.priority() > b.priority() );
  }
};

} } //end namespaces

/// \endcond

#endif // CGAL_OTR2_SAMPLE_H