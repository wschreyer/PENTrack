/**
 * \file
 * Tricubic interpolation of 3D field tables.
 * Also a linked list merge sort utilized to
 */

#ifndef COMSOL_3D_H_
#define COMSOL_3D_H_

#include "field.h"

#include <vector>

/**
 * Class for tricubic field interpolation, create one for every table file you want to use.
 *
 * Reads in a text file where the first 3 columns are x y z coordinates
 * and the next 3 are Bx By Bz. (e.g. a COMSOL arrow plot)
 *
 * Lines beginning with % are ignored
 *
 * Units assumed to be meters and Tesla. You can change this by altering
 * alengthconv = __, double aBconv = ___, in this header file
 *
 *
 */
class comsolField3: public TField{
	private:
		int xl; ///< size of the table in x direction (i.e. the number of x values not counting duplicates)
		int yl; ///< size of the table in y direction (i.e. the number of y values not counting duplicates)
		int zl; ///< size of the table in z direction (i.e. the number of z values not counting duplicates)
		double xdist; ///< distance between grid points in x direction
		double ydist; ///< distance between grid points in y direction
		double zdist; ///< distance between grid points in z direction
		double x_mi; ///< lower x coordinate of cuboid grid
		double y_mi; ///< lower y coordinate of cuboid grid
		double z_mi; ///< lower z coordinate of cuboid grid
		std::vector<std::vector<double> > Bxc; ///< interpolation coefficients for magnetic x component
		std::vector<std::vector<double> > Byc; ///< interpolation coefficients for magnetic y component
		std::vector<std::vector<double> > Bzc; ///< interpolation coefficients for magnetic z component
		double BoundaryWidth; ///< if this is larger 0, the field will be smoothly reduced to 0 in this boundary around the tabulated field cuboid
		double lengthconv; ///< Factor to convert length units in file to PENTrack units
		double Bconv; ///< Factor to convert magnetic field units in file to PENTrack units

		/**
		 *
		 * Read a text file containing x,y, z B_x, B_y ,B_z columns, delineated by space, comma, or tab
		 * Sets comsolField3::xl, comsolField3::yl, comsolField3::zl, comsolField3::xdist, comsolField3::ydist, comsolField3::zdist, comsolField3::x_mi, comsolField3::y_mi, comsolField3::z_mi according to the values in the table file which are used to determine the needed indeces on interpolation.
		 *
		 * @param tabfile Path to table file
		 * @param Bscale Magnetic field is always scaled by this factor
		 * @param Escale Electric field is always scaled by this factor
		 * @param BxTab Returns magnetic field x components at grid points
		 * @param ByTab Returns magnetic field y components at grid points
		 * @param BzTab Returns magnetic field z components at grid points
		 */
		void ReadTabFile(const std::string &tabfile, std::vector<double> &BxTab,
				std::vector<double> &ByTab, std::vector<double> &BzTab);


		/**
		 * Print some information for each table column
		 *
		 * @param BxTab B_x column
		 * @param ByTab B_y column
		 * @param BzTab B_z column
		 */
		void CheckTab(const std::vector<double> &BxTab, const std::vector<double> &ByTab,
				const std::vector<double> &BzTab);


		/**
		 * Calculate spatial derivatives of a table column using 1D cubic spline interpolations.
		 *
		 * @param Tab Table column of which derivatives shall be calculated
		 * @param Tab1 Returns derivatives with respect to x
		 * @param Tab2 Returns derivatives with respect to y
		 * @param Tab3 Returns derivatives with respect to z
		 * @param Tab12 Returns second derivatives with respect to x and y
		 * @param Tab13 Returns second derivatives with respect to x and z
		 * @param Tab23 Returns second derivatives with respect to y and z
		 * @param Tab123 Returns third derivatives with respect to x, y and z
		 */
		void CalcDerivs(const std::vector<double> &Tab, std::vector<double> &Tab1, std::vector<double> &Tab2, std::vector<double> &Tab3,
						std::vector<double> &Tab12, std::vector<double> &Tab13, std::vector<double> &Tab23, std::vector<double> &Tab123) const;


		/**
		 * Calculate tricubic interpolation coefficients for a table column
		 *
		 * Calls comsolField3::CalcDerivs and determines the interpolation coefficients with ::tricubic_get_coeff
		 *
		 * @param coeff Returns coefficients
		 * @param Tab Table column
		 */
		void PreInterpol(std::vector<std::vector<double> > &coeff, const std::vector<double> &Tab) const;


		/**
		 * Smoothly reduce the field at the edges of the tabulated region
		 *
		 * If coordinates are within BoundaryWidth of the edges of the tabulated field,
		 * the field and its derivatives are scaled by the SmthrStp and SmthrStpDer functions.
		 *
		 * @param x x coordinate
		 * @param y y coordinate
		 * @param z z coordinate
		 * @param F Field value at (x,y,z)
		 * @param dFdxi Field derivatives at (x,y,z)
		 */
		void FieldSmthr(const double x, const double y, const double z, double &F, double dFdxi[3]) const;

		/**
		 * Smooth function used to scale the field at the edges
		 *
		 * @param x function parameter (x = 0..1)
		 *
		 * @return Returns number between 0 and 1, smoothly rising with x
		 */
		double SmthrStp(const double x) const;

		/**
		 * Derivative of SmthStpDer
		 *
		 * @param x function parameter (x = 0..1)
		 *
		 * @return Returns derivative of SmthrStp at parameter x
		 */
		double SmthrStpDer(const double x) const;
	public:
		/**
		 * Constructor.
		 *
		 * Calls comsolField3::ReadTabFile, comsolField3::CheckTab and for each column comsolField3::PreInterpol
		 *
		 * @param tabfile Path of table file
		 * @param Bscale Time-dependent scaling formula for magnetic field
		 * @param Escale Time-dependent scaling formula for electric field
		 * @param aBoundaryWidth Sets comsolField3::BoundaryWidth
		 * @param alengthconv Factor to convert length units in file to PENTrack units (default: expect m, convert to m)
		 * @param aBconv Factor to convert magnetic field units in file to PENTrack units (default: expect Tesla, convert to Tesla)
		 */
		comsolField3(const std::string &tabfile, const std::string &Bscale,
				const double aBoundaryWidth, const double alengthconv = 1, const double aBconv = 1);


		/**
		 * Get magnetic field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from comsolField3::x_mi, comsolField3::xdist, comsolField3::y_mi, comsolField3::ydist, comsolField3::z_mi, comsolField3::zdist
		 * and evaluates the interpolation polynom tricubic.h#tricubic_eval.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param B Returns magnetic-field components
		 * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
		 */
		void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = NULL) const;


		/**
		 * Adds no electric field.
		 * Required because we are inheriting from abstract class TField which has virtual void EField function.
		 *
		 * @param x Cartesian x coordinate
		 * @param y Cartesian y coordinate
		 * @param z Cartesian z coordinate
		 * @param t Time
		 * @param V Electric potential
		 * @param Ei Electric field components
		 * @param dEidxj Spatial derivatives of electric field components
		 **/
		void EField(const double x, const double y, const double z, const double t,
				double &V, double Ei[3], double dEidxj[3][3] = NULL) const {};
};



#endif // COMSOL_3D_H_
