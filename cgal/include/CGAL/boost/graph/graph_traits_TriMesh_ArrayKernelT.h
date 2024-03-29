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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/BGL/include/CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h $
// $Id: graph_traits_TriMesh_ArrayKernelT.h 150bde7 %aI Maxime Gimeno
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Andreas Fabri, Philipp Moeller

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_TRIMESH_ARRAYKERNELT_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_TRIMESH_ARRAYKERNELT_H

// http://openmesh.org/Documentation/OpenMesh-Doc-Latest/classOpenMesh_1_1Concepts_1_1KernelT.html
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <CGAL/boost/graph/properties_TriMesh_ArrayKernelT.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#define OPEN_MESH_CLASS OpenMesh::TriMesh_ArrayKernelT<K>
#include <CGAL/boost/graph/graph_traits_OpenMesh.h>

#endif // CGAL_BOOST_GRAPH_TRAITS_TRIMESH_ARRAYKERNELT_H
