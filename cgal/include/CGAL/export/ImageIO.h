// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
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
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/CGAL_ImageIO/include/CGAL/export/ImageIO.h $
// $Id: ImageIO.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IMAGEIO_EXPORT_H
#define CGAL_IMAGEIO_EXPORT_H

#include <CGAL/config.h>
#include <CGAL/export/helpers.h>

#if defined(CGAL_BUILD_SHARED_LIBS) && ! defined(CGAL_HEADER_ONLY)

#  if defined(CGAL_ImageIO_EXPORTS) // defined by CMake or in cpp files of the dll

#    define CGAL_IMAGEIO_EXPORT CGAL_DLL_EXPORT
#    define CGAL_IMAGEIO_EXPIMP_TEMPLATE

#  else // not CGAL_ImageIO_EXPORTS

#    define CGAL_IMAGEIO_EXPORT CGAL_DLL_IMPORT
#    define CGAL_IMAGEIO_EXPIMP_TEMPLATE extern

#  endif // not CGAL_IMAGEIO_EXPORTS

#else // not CGAL_BUILD_SHARED_LIBS

#  define CGAL_IMAGEIO_EXPORT
#  define CGAL_IMAGEIO_EXPIMP_TEMPLATE

#endif // not CGAL_BUILD_SHARED_LIBS

#endif //  CGAL_IMAGEIO_EXPORT_H


