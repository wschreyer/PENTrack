/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */

#ifndef FIELD_3D_H_
#define FIELD_3D_H_

#include "field.h"

#include <vector>

#include "boost/multi_array.hpp"

/**
 * Class for tricubic field interpolation, create one for every table file you want to use.
 *
 * This class loads a tabulated magnetic and electric field on a rectilinear, three-dimensional grid and
 * calculates tricubic interpolation coefficients (4x4x4 = 64 for each grid point) to allow fast evaluation of the fields at arbitrary points.
 *
 */
class TabField3: public TField{
private:
        std::array<std::vector<double>, 3> xyz; ///< coordinates of points on interpolation grid
        typedef boost::multi_array<double, 3> array3D;
        typedef std::array<double, 64> tricubic_coeff; ///< interpolation coefficients for one grid cell
        typedef boost::multi_array<tricubic_coeff, 3> field_type; ///< interpolation coefficients for all grid cells
        std::array<field_type, 3> Bc; ///< interpolation coefficients for magnetix x,y, and z components
        field_type Vc; ///< interpolation coefficients for electric potential
private:
		/**
		 * Print some information for each table column
		 *
         * @param B Lists of Bx, By, and Bz magnetic field components on grid points
         * @param V List of electric potentials on grid points
		 */
        void CheckTab(const std::array<std::vector<double>, 3> &B, const std::vector<double> &V);


		/**
         * Calculate spatial derivatives of a table column along one dimension using 1D cubic spline interpolations.
		 *
         * @param Tab 3D array of field components on grid points
         * @param diff_dim Coordinate dimension to differentiate (0, 1, or 2 for x, y, or z)
         * @param DiffTab Returns 3D array of field components differentiated with respect to dimension diff_dim
         */
        void CalcDerivs(const array3D &Tab, const unsigned long diff_dim, array3D &DiffTab) const;


		/**
		 * Calculate tricubic interpolation coefficients for a table column
		 *
		 * Calls TabField3::CalcDerivs and determines the interpolation coefficients with ::tricubic_get_coeff
		 *
         * @param Tab 3D array of field components on grid
         * @param coeff Returns 3D array of tricubic interpolation coefficients for each grid cell
         */
        void PreInterpol(const array3D &Tab, field_type &coeff) const;


		/**
		 * Interpolate field component at a specific point.
		 *
		 * Uses a binary search to find the grid cell that contains the point (x,y,z) and and calculates the tricubic interpolation using the coefficients belonging to the grid cell.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param coeffs 3D array of tricubic interpolation coefficients
		 * @param F Returns interpolated field component
		 * @param dFdxi Returns spatial derivatives of field component (can be null)
		 */
		void Interpolate(const double x, const double y, const double z,
                            const field_type& coeffs, double &F, double dFdxi[3]) const;
	public:
		/**
		 * Constructor.
		 *
		 * Calls TabField3::ReadTabFile, TabField3::CheckTab and for each column TabField3::PreInterpol
		 *
         * @param xyzTab Lists of x, y, and z coordinates of grid points
         * @param BTab Lists of Bx, By, and Bz magnetic field components on grid points
         * @param VTab List of electric potentials on grid points
		 */
        TabField3(const std::array<std::vector<double>, 3> &xyzTab, const std::array<std::vector<double>, 3> &BTab, const std::vector<double> &VTab);


		/**
		 * Get magnetic field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField3::x_mi, TabField3::xdist, TabField3::y_mi, TabField3::ydist, TabField3::z_mi, TabField3::zdist
		 * and evaluates the interpolation polynom tricubic.h#tricubic_eval.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param B Returns magnetic-field components
		 * @param dBidxj Returns spatial derivatives of magnetic-field components (optional)
		 */
		void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const override;


		/**
		 * Get electric field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField3::x_mi, TabField3::xdist, TabField3::y_mi, TabField3::ydist, TabField3::z_mi, TabField3::zdist
		 * and evaluates the interpolation polynom tricubic.h#tricubic_eval.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param V Returns electric potential
		 * @param Ei Returns electric field (negative spatial derivatives of V)
		 */
		void EField(const double x, const double y, const double z, const double t,
				double &V, double Ei[3]) const override;
};

/**
 * Read 3D table file exported from OPERA
 * @param params String containing parameters defined in config.in. Should contain field type "3Dtable", file name, magnetic field scaling formula, electric field scaling formula, and boundary width
 * @return Pointer to created class, derived from TField
 */
TFieldContainer ReadOperaField3(const std::string &params, const std::map<std::string, std::string> &formulas);

/**
* Read generic file containing table of magnetic field mapped on list of points, e.g. exported from COMSOL
* @param params String containing parameters defined in config.in. Should contain field type "COMSOL", file name, magnetic field scaling formula, and boundary width
* @return Pointer to created class, derived from TField
*/
TFieldContainer ReadComsolField(const std::string &params, const std::map<std::string, std::string> &formulas);

#endif // FIELD_3D_H_
