
/**
 * \file
 * Implementation of a magnetic field determined by 
 * coefficients, provided as inputs by the user, of an expansion in
 * terms of harmonic polynomials.
*/

#include "harmonicfields.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <cmath>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/constants/constants.hpp>

HarmonicExpandedBField::HarmonicExpandedBField(const double _xoff, const double _yoff, const double _zoff, 
		const double _axis_x, const double _axis_y, const double _axis_z, const double _angle, 
		const double G0, const double G1, const double G2, const double G3, const double G4, const double G5, const double G6, const double G7, const double G8, 
		const double G9, const double G10, const double G11, const double G12, const double G13, const double G14, const double G15, const double G16, const double G17, 
		const double G18, const double G19, const double G20, const double G21, const double G22, const double G23){

	// the offset values
	xoff = _xoff;
	yoff = _yoff;
	zoff = _zoff;

	// parameters defining the axis-angle quaternion rotation
	axis_x = _axis_x;
	axis_y = _axis_y;
	axis_z = _axis_z;
	angle  = _angle;

	// G parameters
	G[0]  = G0;
	G[1]  = G1;
	G[2]  = G2;
	G[3]  = G3;
	G[4]  = G4;
	G[5]  = G5;
	G[6]  = G6;
	G[7]  = G7;
	G[8]  = G8;
	G[9]  = G9;
	G[10] = G10;
	G[11] = G11;
	G[12] = G12;
	G[13] = G13;
	G[14] = G14;
	G[15] = G15;
	G[16] = G16;
	G[17] = G17;
	G[18] = G18;
	G[19] = G19;
	G[20] = G20;
	G[21] = G21;
	G[22] = G22;
	G[23] = G23;

	}

void HarmonicExpandedBField::BField(const double _x, const double _y, const double _z, const double t, double B[3], double dBidxj[3][3]) const{

	// Updating the x, y, z values with the given offset values
	double x = _x + xoff;
	double y = _y + yoff;
	double z = _z + zoff;

	/* 
	Information about Calculating the Harmonic Polynomial Expansion

	The calculations below are hard-coded using the cartesian form of the 
	harmonic polynomial expansion. Alternatively, this can be carried out 
	in an iterative fashion using the BOOST library for "Legendre and 
	associated polynomials":
	
	boost.org/doc/libs/1_65_0/libs/math/doc/html/math_toolkit/sf_poly/legendre.html 

	This was not chosen for the implementation here because it was believed
	that this hard-coded representation displayed the physical form of the 
	functions more explicitly. Perhaps more importantly, this form is much
	easier to debug. 

	If in future, simulations wish to explore harmonic expansions above third 
	order, returning to an implementation using BOOST's legendre library may
	be preferrable. 
	*/

	// Calculating B_x
	B[0] = 	  G[2] 	\
			+ G[3]  * (y) \
			+ G[5]  * (-x / 2) \
			+ G[6]  * (z) \
			+ G[7]  * (x) \
			+ G[8]  * (2 * x * y) \
			+ G[9]  * (2 * y * z) \
			+ G[10] * (-x * y / 2) \
			+ G[11] * (-x * z) \
			+ G[12] * (-(3 * std::pow(x,2) + std::pow(y,2) - 4 * std::pow(z,2)) / 4) \
			+ G[13] * (2 * x * z) \
			+ G[14] * ((std::pow(x,2) - std::pow(y,2))) \
			+ G[15] * ((3 * std::pow(x,2) * y - std::pow(y,3))) \
			+ G[16] * (6 * x * y * z) \
			+ G[17] * (-(3 * std::pow(x,2) * y + std::pow(y,3) - 6 * y * std::pow(z,2)) / 2) \
			+ G[18] * (-(3 / 2) * (x * y * z)) \
			+ G[19] * ((3 / 8) * (std::pow(x,3) + x * std::pow(y,2) - 4 * x * std::pow(z,2))) \
			+ G[20] * ((-1 / 4) * (9 * std::pow(x,2) * z + 3 * std::pow(y,2) * z - 4 * std::pow(z,3))) \
			+ G[21] * (-std::pow(x,3) + 3 * x * std::pow(z,2)) \
			+ G[22] * (3 * (std::pow(x,2) * z - std::pow(y,2) * z)) \
			+ G[23] * (std::pow(x,3) - 3 * x * std::pow(y,2));
	
	// Calculating B_y
	B[1] =	  G[0] 	\
			+ G[3]  * (x) \
			+ G[4]  * (z) \
			+ G[5]  * (-y / 2) \
			+ G[7]  * (-y) \
			+ G[8]  * (std::pow(x,2) - std::pow(y,2)) \
			+ G[9]  * (2 * x * z) \
			+ G[10] * ((std::pow(x,2) + 3 * std::pow(y,2) - 4 * std::pow(z,2)) / 4) \
			+ G[11] * (-y * z) \
			+ G[12] * (-x * y / 2) \
			+ G[13] * (-2 * y * z) \
			+ G[14] * (-2 * x * y) \
			+ G[15] * (std::pow(x,3) - 3 * x * std::pow(y,2)) \
			+ G[16] * (3 * (std::pow(x,2) * z - std::pow(y,2) * z)) \
			+ G[17] * (-(std::pow(x,3) + 3 * x * std::pow(y,2) - 6 * x * std::pow(z,2)) / 2) \
			+ G[18] * ( (-1 / 4) * (3 * std::pow(x,2) * z + 9 * std::pow(y,2) * z - 4 * std::pow(z,3))) \
			+ G[19] * ((3 / 8) * (std::pow(x,2) * y + std::pow(y,3) - 4 * y * std::pow(z,2))) \
			+ G[20] * ((-3 / 2) * x * y * z) \
			+ G[21] * (-3 * y * std::pow(z,2) + std::pow(y,3)) \
			+ G[22] * (-6 * x * y * z) \
			+ G[23] * (-3 * std::pow(x,2) * y + std::pow(y,3));

	// Calculating B_z
	B[2] =	  G[1] 	\
			+ G[4]  * (y) \
			+ G[5]  * (z) \
			+ G[6]  * (x) \
			+ G[9]  * (2 * x * y) \
			+ G[10] * (2 * y * z) \
			+ G[11] * (std::pow(z,2) - (1 / 2) * (std::pow(x,2) + std::pow(y,2))) \
			+ G[12] * (2 * x * z) \
			+ G[13] * (std::pow(x,2) - std::pow(y,2)) \
			+ G[16] * (3 * std::pow(x,2) * y - std::pow(y,3))\
			+ G[17] * (6 * x * y * z)  \
			+ G[18] * (3 * y * std::pow(z,2) - (3/4) * (std::pow(x,2) * y + std::pow(y,3))) \
			+ G[19] * (std::pow(z,3) - (3 / 2) * z * (std::pow(x,2) + std::pow(y,2))) \
			+ G[20] * (3 * x * std::pow(z,2) - (3/4) * (std::pow(x,3) + x * std::pow(y,2))) \
			+ G[21] * (3 * (std::pow(x,2) * z - std::pow(y,2) * z)) \
			+ G[22] * (std::pow(x,3) - 3 * x * std::pow(y,2));


	if (dBidxj != nullptr){

		// Calculating dBxdx
		dBidxj[0][0] = -G[5]/2 + G[7] + 2*G[8]*y - G[10]*y/2 - G[11]*z - \
					3*G[12]*x/2 + 2*G[13]*z + 2*G[14]*x + 6*G[15]*x*y + \
					6*G[16]*y*z - 3*G[17]*x*y - 1.5*G[18]*y*z + \
					G[19]*(1.125*std::pow(x,2) + 0.375*std::pow(y,2) - 1.5*std::pow(z,2)) - \
					4.5*G[20]*x*z + G[21]*(-3*std::pow(x,2) + 3*std::pow(z,2)) + \
					6*G[22]*x*z + G[23]*(3*std::pow(x,2) - 3*std::pow(y,2));

		// Calculating dBxdy
		dBidxj[0][1] = G[3] + 2*G[8]*x + 2*G[9]*z - G[10]*x/2 - G[12]*y/2 - \
					2*G[14]*y + G[15]*(3*std::pow(x,2) - 3*std::pow(y,2)) + 6*G[16]*x*z+\
					G[17]*(-3*std::pow(x,2)/2 - 3*std::pow(y,2)/2 + 3*std::pow(z,2)) - \
					1.5*G[18]*x*z + 0.75*G[19]*x*y - 1.5*G[20]*y*z - \
					6*G[22]*y*z - 6*G[23]*x*y;

		// Calculating dBxdz
		dBidxj[0][2] = G[6] + 2*G[9]*y - G[11]*x + 2*G[12]*z + 2*G[13]*x + \
					6*G[16]*x*y + 6*G[17]*y*z - 1.5*G[18]*x*y - 3.0*G[19]*x*z+\
					G[20]*(-2.25*std::pow(x,2) - 0.75*std::pow(y,2) + 3.0*std::pow(z,2)) + \
					6*G[21]*x*z + G[22]*(3*std::pow(x,2) - 3*std::pow(y,2));

		// Calculating dBydx
		dBidxj[1][0] = G[3] + 2*G[8]*x + 2*G[9]*z + G[10]*x/2 - G[12]*y/2 - \
					2*G[14]*y + G[15]*(3*std::pow(x,2) - 3*std::pow(y,2)) + 6*G[16]*x*z+\
					G[17]*(-3*std::pow(x,2)/2 - 3*std::pow(y,2)/2 + 3*std::pow(z,2)) - \
					1.5*G[18]*x*z + 0.75*G[19]*x*y - 1.5*G[20]*y*z - \
					6*G[22]*y*z - 6*G[23]*x*y;

		// Calculating dBydy
		dBidxj[1][1] = -G[5]/2 - G[7] - 2*G[8]*y + 3*G[10]*y/2 - G[11]*z - \
					G[12]*x/2 - 2*G[13]*z - 2*G[14]*x - 6*G[15]*x*y - \
					6*G[16]*y*z - 3*G[17]*x*y - 4.5*G[18]*y*z + \
					G[19]*(0.375*std::pow(x,2) + 1.125*std::pow(y,2) - 1.5*std::pow(z,2)) - \
					1.5*G[20]*x*z + G[21]*(3*std::pow(y,2) - 3*std::pow(z,2)) - \
					6*G[22]*x*z + G[23]*(-3*std::pow(x,2) + 3*std::pow(y,2));

		// Calculating dBydz
		dBidxj[1][2] = G[4] + 2*G[9]*x - 2*G[10]*z - G[11]*y - 2*G[13]*y + \
					G[16]*(3*std::pow(x,2) - 3*std::pow(y,2)) + 6*G[17]*x*z + \
					G[18]*(-0.75*std::pow(x,2) - 2.25*std::pow(y,2) + 3.0*std::pow(z,2)) - \
					3.0*G[19]*y*z - 1.5*G[20]*x*y - 6*G[21]*y*z - 6*G[22]*x*y;

		// Calculating dBzdx
		dBidxj[2][0] = G[6] + 2*G[9]*y - 1.0*G[11]*x + 2*G[12]*z + 2*G[13]*x + \
					6*G[16]*x*y + 6*G[17]*y*z - 1.5*G[18]*x*y - 3.0*G[19]*x*z+\
					G[20]*(-2.25*std::pow(x,2) - 0.75*std::pow(y,2) + 3*std::pow(z,2)) + \
					6*G[21]*x*z + G[22]*(3*std::pow(x,2) - 3*std::pow(y,2));

		// Calculating dBzdy
		dBidxj[2][1] = G[4] + 2*G[9]*x + 2*G[10]*z - 1.0*G[11]*y - 2*G[13]*y + \
					G[16]*(3*std::pow(x,2) - 3*std::pow(y,2)) + 6*G[17]*x*z + \
					G[18]*(-0.75*std::pow(x,2) - 2.25*std::pow(y,2) + 3*std::pow(z,2)) - \
					3.0*G[19]*y*z - 1.5*G[20]*x*y - 6*G[21]*y*z - 6*G[22]*x*y;

		// Calculating dBzdz
		dBidxj[2][2] = G[5] + 2*G[10]*y + 2*G[11]*z + 2*G[12]*x + 6*G[17]*x*y + \
					6*G[18]*y*z + G[19]*(-1.5*std::pow(x,2) - 1.5*std::pow(y,2) + \
					3*std::pow(z,2)) + 6*G[20]*x*z + G[21]*(3*std::pow(x,2) -3*std::pow(y,2));
	}

	/* 
	Information about Rotations with Quaternions

	from the double-specific constructor documentation of BOOST
	boost.org/doc/libs/1_70_0/libs/math/doc/html/math_toolkit/quat_mem_fun.html

	explicit quaternion(double const & requested_a = 0.0, 
						double const & requested_b = 0.0, 
						double const & requested_c = 0.0, 
						double const & requested_d = 0.0);

	This constructs a quartention. From the useful wikipedia page on rotations
	and quartenions... wikipedia.org/wiki/Quaternions_and_spatial_rotation

	"In a programmatic implementation, this is achieved by constructing a 
	quaternion whose vector part is p and real part equals zero and then
	performing the quaternion multiplication. The vector part of the resulting 
	quaternion is the desired vector p′."

	here they refer to the original vector, which we want to rotate, as p. The 
	resulting, rotated, vector is called p'. Refer to the referenced wikipedia
	page for details on how to construct the desired quartenion q. On that 
	wikipedia page, an example rotation in 3D space is provided. 

	(!!!) due to the limitation of boost version within the UCN cluster, 
	we are not able to use the QVM library (earliest v1.62 of boost) which 
	contains more updated functionality specific to the rotations 
	implemented here. See http://boostorg.github.io/qvm/#quat_rotate
	and https://github.com/boostorg/qvm
	*/

	// The B vector is rotated, as well as the three dBxdxj, dBydxj dBzdxj 
	// directional gradient vectors, using the same quaternion conjugation.
	// This achieves a full rotation of the magnetic field and its spatial 
	// derivatives.

	// if the user provides values of 0 for every input to the axis of rotation
	// then no rotation is performed.
	if (axis_x != 0 || axis_y != 0 || axis_z != 0){

        // construct the quarternion, p, which represents the vector to be 
		// rotated
        // here the 0th component is equal to zero, and the remaining three 
        // components are those of the earlier computed magnetic field
        boost::math::quaternion<double> p(0, B[0], B[1], B[2]);

        // the components of the axis of rotation (member variables of this 
		// class) are used to construct a normalized vector.
        boost::numeric::ublas::vector<double> axis(3, 0);
        axis(0) = axis_x;
        axis(1) = axis_y;
        axis(2) = axis_z;

        // the rotational axis vector is normalized
        axis = (axis / boost::numeric::ublas::norm_2(axis));

        // calculate quarternion, q, whih represents the rotation
        boost::math::quaternion<double> q(std::cos(angle / 2), 
                                          std::sin(angle / 2) * axis(0), 
                                          std::sin(angle / 2) * axis(1), 
                                          std::sin(angle / 2) * axis(2));

        // get the conjugate of the quaternion, q'
        boost::math::quaternion<double> q_prime = boost::math::conj(q);

        // perform multiplication to rotate
        boost::math::quaternion<double> p_prime = q * p * q_prime;

        // extract the resulting vector components, post-rotation
        B[0] = p_prime.R_component_2();
        B[1] = p_prime.R_component_3();
        B[2] = p_prime.R_component_4();

	if (dBidxj != nullptr){
	        // We create three quaternions, one for the gradient of each of the 
	        // gradient vectors

      	  boost::math::quaternion<double> p_x(0, 
      	                                      dBidxj[0][0],
      	                                      dBidxj[0][1],
      	                                      dBidxj[0][2]);

      	  boost::math::quaternion<double> p_y(0, 
      	                                      dBidxj[1][0], 
      	                                      dBidxj[1][1], 
      	                                      dBidxj[1][2]);

      	  boost::math::quaternion<double> p_z(0, 
      	                                      dBidxj[2][0], 
      	                                      dBidxj[2][1], 
      	                                      dBidxj[2][2]);

      	  // we use the same quaternion, q, and conjugate q' for these rotations
      	  boost::math::quaternion<double> p_x_prime = q * p_x * q_prime;
      	  boost::math::quaternion<double> p_y_prime = q * p_y * q_prime;
      	  boost::math::quaternion<double> p_z_prime = q * p_z * q_prime;

      	  // extract the resulting vector components, post-rotation
      	  // assign the components to the dBidxj matrix
      	  dBidxj[0][0] = p_x_prime.R_component_2();
      	  dBidxj[0][1] = p_x_prime.R_component_3();
      	  dBidxj[0][2] = p_x_prime.R_component_4();

      	  dBidxj[1][0] = p_y_prime.R_component_2();
      	  dBidxj[1][1] = p_y_prime.R_component_3();
      	  dBidxj[1][2] = p_y_prime.R_component_4();

      	  dBidxj[2][0] = p_z_prime.R_component_2();
      	  dBidxj[2][1] = p_z_prime.R_component_3();
      	  dBidxj[2][2] = p_z_prime.R_component_4();
	}

    }

	// }
}

