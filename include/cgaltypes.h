#ifndef INCLUDE_CGALTYPES_H_
#define INCLUDE_CGALTYPES_H_

#include <CGAL/Simple_cartesian.h>


typedef CGAL::Simple_cartesian<double> CKernel; ///< Geometric Kernel used for CGAL types
typedef CKernel::Segment_3 CSegment; ///< CGAL segment type
typedef CKernel::Point_3 CPoint; ///< CGAL point type
typedef CKernel::Vector_3 CVector; ///< CGAL vector type
typedef CKernel::Triangle_3 CTriangle; ///< CGAL triangle type
typedef CGAL::Bbox_3 CBox; ///< CGAL Box type


#endif /* INCLUDE_CGALTYPES_H_ */
