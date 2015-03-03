/**
 * \file
 * Bicubic interpolation of axisymmetric field tables.
 */

#ifndef FIELD_2D_H_
#define FIELD_2D_H_

#include <vector>

#include "interpolation.h"

#include "field.h"

/**
 * Class for bicubic field interpolation, create one for every table file you want to use.
 *
 * This class loads a special file format from "Vectorfields Opera" containing a regular, rectangular table of magnetic and electric fields and
 * calculates bicubic interpolation coefficients (4x4 matrix for each grid point) to allow fast evaluation of the fields at arbitrary points.
 * Therefore it assumes that the fields are axisymmetric around the z axis.
 *
 */
class TabField: public TField{
	private:
		int m; ///< radial size of the table file
		int n; ///< axial size of the arrays
		alglib::real_1d_array rind, zind, BrTab, BphiTab, BzTab, ErTab, EphiTab, EzTab, VTab;
		alglib::spline2dinterpolant Brc, Bphic, Bzc, Erc, Ephic, Ezc, Vc;
		double NullFieldTime; ///< Time before magnetic field is ramped (passed by constructor)
		double RampUpTime; ///< field is ramped linearly from 0 to 100% in this time (passed by constructor)
		double FullFieldTime; ///< Time the field stays at 100% (passed by constructor)
		double RampDownTime; ///< field is ramped down linearly from 100% to 0 in this time (passed by constructor)


		/**
		 * Reads an Opera table file.
		 *
		 * File has to contain x and z coordinates, it may contain B_x, B_y ,B_z, E_x, E_y, E_z and V columns. If V is present, E_i are ignored.
		 * Sets TabField::m, TabField::n, TabField::rdist, TabField::zdist, TabField::r_mi, TabField::z_mi according to the values in the table file which are used to determine the needed indeces on interpolation.
		 *
		 * @param tabfile Path to table file
		 * @param Bscale Magnetic field is always scaled by this factor
		 * @param Escale Electric field is always scaled by this factor
		 * @param BrTab Returns radial magnetic field components at (r,z) = (x,z)
		 * @param BphiTab Returns azimuthal magnetic field components at (r,z) = (x,z)
		 * @param BzTab Returns axial magnetic field components at (r,z) = (x,z)
		 * @param ErTab Returns radial electric field components at (r,z) = (x,z)
		 * @param EphiTab Returns azimuthal electric field components at (r,z) = (x,z)
		 * @param EzTab Returns axial electric field components at (r,z) = (x,z)
		 * @param VTab Returns electric potential at (r,z) = (x,z)
		 */
		void ReadTabFile(const char *tabfile, double Bscale, double Escale);


		/**
		 * Print some information for each table column
		 *
		 * @param BrTab B_r column
		 * @param BphiTab B_phi column
		 * @param BzTab B_z column
		 * @param ErTab E_r column
		 * @param EphiTab E_phi column
		 * @param EzTab E_z column
		 * @param VTab	V column
		 */
		void CheckTab();


		/**
		 * Get magnetic field scale factor for a specific time.
		 *
		 * Determined by TabField::NullFieldTime, TabField::RampUpTime, TabField::FullFieldTime, TabField::RampDownTime
		 *
		 * @param t Time
		 *
		 * @return Returns magnetic field scale factor
		 */
		double BFieldScale(double t);
	public:
		/**
		 * Constructor.
		 *
		 * Calls TabField::ReadTabFile, TabField::CheckTab and for each column TabField::PreInterpol
		 *
		 * @param tabfile Path of table file
		 * @param Bscale Magnetic field is always scaled by this factor
		 * @param Escale Electric field is always scaled by this factor
		 * @param aNullFieldTime Sets TabField::NullFieldTime
		 * @param aRampUpTime Sets TabField::RampUpTime
		 * @param aFullFieldTime Sets TabField::FullFieldTime
		 * @param aRampDownTime Set TabField::RampDownTime
		 */
		TabField(const char *tabfile, double Bscale, double Escale,
				double aNullFieldTime, double aRampUpTime, double aFullFieldTime, double aRampDownTime);

		~TabField();

		/**
		 * Get magnetic field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField::r_mi, TabField::rdist, TabField::z_mi, TabField::zdist
		 * and evaluates the interpolation polynom ::bcuint.
		 * These radial, axial und azimuthal components have to be rotated into cartesian coordinate system.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param B Adds magnetic field components B[0..2][0], their derivatives B[0..2][1..3], the absolute value B[3][0] and its derivatives B[3][1..3]
		 */
		void BField(double x, double y, double z, double t, double B[4][4]);


		/**
		 * Get electric field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField::r_mi, TabField::rdist, TabField::z_mi, TabField::zdist
		 * and evaluates the interpolation polynom ::bcuint.
		 * These radial, axial und azimuthal components have to be rotated into cartesian coordinate system.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param V Returns electric potential
		 * @param Ei Return electric field (negative spatial derivatives of V)
		 *
		 * @return Returns true if electric field could be evaluated at this point
		 */
		void EField(double x, double y, double z, double t, double &V, double Ei[3]);
};


#endif // FIELD_2D_H_
