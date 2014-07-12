/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */

#ifndef FIELD_3D_H_
#define FIELD_3D_H_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "libtricubic/tricubic.h"
#include "field.h"

int INDEX_3D(int xi, int yi, int zi, int xsize, int ysize, int zsize){
	return zi*xsize*ysize + yi*xsize + xi;
}

/**
 * Class for tricubic field interpolation, create one for every table file you want to use.
 *
 * This class loads a special file format from "Vectorfields Opera" containing a regular, cuboid table of magnetic and electric fields and
 * calculates tricubic interpolation coefficients (4x4x4 = 64 for each grid point) to allow fast evaluation of the fields at arbitrary points.
 *
 */
class TabField3: public TField{
	private:
		int xl; ///< size of the table in x direction
		int yl; ///< size of the table in y direction
		int zl; ///< size of the table in z direction
		double xdist; ///< distance between grid points in x direction
		double ydist; ///< distance between grid points in y direction
		double zdist; ///< distance between grid points in z direction
		double x_mi; ///< lower x coordinate of cuboid grid
		double y_mi; ///< lower y coordinate of cuboid grid
		double z_mi; ///< lower z coordinate of cuboid grid
		vector<vector<double> > Bxc; ///< interpolation coefficients for magnetic x component
		vector<vector<double> > Byc; ///< interpolation coefficients for magnetic y component
		vector<vector<double> > Bzc; ///< interpolation coefficients for magnetic z component
		vector<vector<double> > Vc; ///< interpolation coefficients for electric potential
		double NullFieldTime; ///< Time before magnetic field is ramped (passed by constructor)
		double RampUpTime; ///< field is ramped linearly from 0 to 100% in this time (passed by constructor)
		double FullFieldTime; ///< Time the field stays at 100% (passed by constructor)
		double RampDownTime; ///< field is ramped down linearly from 100% to 0 in this time (passed by constructor)


		/**
		 * Reads an Opera table file.
		 *
		 * File has to contain x,y and z coordinates, it may contain B_x, B_y ,B_z and V columns.
		 * Sets TabField3::xl, TabField3::yl, TabField3::zl, TabField3::xdist, TabField3::ydist, TabField3::zdist, TabField3::x_mi, TabField3::y_mi, TabField3::z_mi according to the values in the table file which are used to determine the needed indeces on interpolation.
		 *
		 * @param tabfile Path to table file
		 * @param Bscale Magnetic field is always scaled by this factor
		 * @param Escale Electric field is always scaled by this factor
		 * @param BxTab Returns magnetic field x components at grid points
		 * @param ByTab Returns magnetic field y components at grid points
		 * @param BzTab Returns magnetic field z components at grid points
		 * @param VTab Returns electric potential at grid points
		 */
		void ReadTabFile(const char *tabfile, double Bscale, double Escale,
				vector<double> &BxTab, vector<double> &ByTab, vector<double> &BzTab, vector<double> &VTab){
			ifstream FIN(tabfile, ifstream::in);
			if (!FIN.is_open()){
				printf("\nCould not open %s!\n",tabfile);
				exit(-1);
			}
			printf("\nReading %s ",tabfile);
			string line;
			FIN >> xl >> yl >> zl;

			getline(FIN,line);
			getline(FIN,line);
			getline(FIN,line);
			getline(FIN,line);
			getline(FIN,line);

			if (line.find("BX") != string::npos){
				BxTab.resize(xl*yl*zl);
				getline(FIN,line);
			}
			if (line.find("BY") != string::npos){
				ByTab.resize(xl*yl*zl);
				getline(FIN,line);
			}
			if (line.find("BZ") != string::npos){
				BzTab.resize(xl*yl*zl);
				getline(FIN,line);
			}

			if (line.find("V") != string::npos){	// file contains potential?
				VTab.resize(xl*yl*zl);
				getline(FIN,line);
			}

			if (!FIN || line.substr(0,2) != " 0"){
				printf("%s not found or corrupt! Exiting...\n",tabfile);
				exit(-1);
			}

			vector<double> xind(xl), yind(yl), zind(zl);
			int xi = 0, yi = 0, zi = -1, perc = 0;
			double x, y, z, val;
			while (FIN.good()){
				FIN >> x;
				FIN >> y;
				FIN >> z;
				if (!FIN) break;
				x *= lengthconv;
				y *= lengthconv;
				z *= lengthconv;
				if (zi >= 0 && z < zind[zi]){
					if (yi >= 0 && y < yind[yi]){
						xi++;
						yi = 0;
					}
					else yi++;
					zi = 0;
				}
				else zi++;

				int i3 = INDEX_3D(xi, yi, zi, xl, yl, zl);
				// status if read is displayed
				PrintPercent((float)i3/(xl*yl*zl), perc);

				xind[xi] = x;
				yind[yi] = y;
				zind[zi] = z;
				if (BxTab.size() > 0){
					FIN >> val;
					BxTab[i3] = val*Bconv*Bscale;
				}
				if (ByTab.size() > 0){
					FIN >> val;
					ByTab[i3] = val*Bconv*Bscale;
				}
				if (BzTab.size() > 0){
					FIN >> val;
					BzTab[i3] = val*Bconv*Bscale;
				}
				if (VTab.size() > 0){
					FIN >> val;
					VTab[i3] = val*Escale;
				}
				FIN >> ws;

			}

			printf("\n");
			if (xi+1 != xl || yi + 1 != yl || zi+1 != zl){
				printf("The header says the size is %i by %i by %i, actually it is %i by %i by %i! Exiting...\n", xl, yl, zl, xi+1, yi+1, zi+1);
				exit(-1);
			}
			FIN.close();

			if (Bscale == 0){
				BxTab.clear();
				ByTab.clear();
				BzTab.clear();
			}
			if (Escale == 0){
				VTab.clear();
			}

			xdist = xind[1] - xind[0];
			ydist = yind[1] - yind[0];
			zdist = zind[1] - zind[0];
			x_mi = xind[0];
			y_mi = yind[0];
			z_mi = zind[0];

		};


		/**
		 * Print some information for each table column
		 *
		 * @param BxTab B_x column
		 * @param ByTab B_y column
		 * @param BzTab B_z column
		 * @param VTab V column
		 */
		void CheckTab(vector<double> &BxTab, vector<double> &ByTab, vector<double> &BzTab, vector<double> &VTab){
			//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
			printf("The arrays are %d by %d by %d.\n",xl,yl,zl);
			printf("The x values go from %G to %G\n",x_mi,x_mi + xdist*(xl-1));
			printf("The y values go from %G to %G\n",y_mi,y_mi + ydist*(yl-1));
			printf("The z values go from %G to %G.\n",z_mi,z_mi + zdist*(zl-1));
			printf("xdist = %G ydist = %G zdist = %G\n",xdist,ydist,zdist);

			double Babsmax = 0, Babsmin = 9e99, Babs;
			double Vmax = 0, Vmin = 9e99;
			for (int j=0; j < xl; j++)
			{
				for (int k=0; k < yl; k++)
				{
					for (int l=0; l < zl; l++)
					{
						Babs = 0;
						int i3 = INDEX_3D(j, k, l, xl, yl, zl);
						if (BxTab.size() > 0)
							Babs += BxTab[i3] * BxTab[i3];
						if (ByTab.size() > 0)
							Babs += ByTab[i3] * ByTab[i3];
						if (BzTab.size() > 0)
							Babs += BzTab[i3] * BzTab[i3];
						Babsmax = max(sqrt(Babs),Babsmax);
						Babsmin = min(sqrt(Babs),Babsmin);
						if (VTab.size() > 0)
						{
							Vmax = max(VTab[i3],Vmax);
							Vmin = min(VTab[i3],Vmin);
						}

					}
				}
			}

			printf("The input table file has values of |B| from %G T to %G T and values of V from %G V to %G V\n",Babsmin,Babsmax,Vmin,Vmax);
		};


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
		void CalcDerivs(vector<double> &Tab, vector<double> &Tab1, vector<double> &Tab2, vector<double> &Tab3,
						vector<double> &Tab12, vector<double> &Tab13, vector<double> &Tab23, vector<double> &Tab123)
		{
			vector<double> x(xl), y(xl);
			for (int zi = 0; zi < zl; zi++){
				for (int yi = 0; yi < yl; yi++){
					for (int xi = 0; xi < xl; xi++){
						x[xi] = x_mi + xi*xdist;
						y[xi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
				    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, xl);
				    gsl_spline_init(spline, &x[0], &y[0], xl);
					for (int xi = 0; xi < xl; xi++)
						Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[xi], NULL); // get derivatives dF/dx from spline interpolation
					gsl_spline_free(spline);
				}
			}

			x.resize(yl);
			y.resize(yl);
			for (int xi = 0; xi < xl; xi++){
				for (int zi = 0; zi < zl; zi++){
					for (int yi = 0; yi < yl; yi++){
						x[yi] = y_mi + yi*ydist;
						y[yi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, yl);
					gsl_spline_init(spline, &x[0], &y[0], yl);
					for (int yi = 0; yi < yl; yi++)
						Tab2[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[yi], NULL); // get derivatives dF/dy from spline interpolation
					gsl_spline_free(spline);
				}
			}

			x.resize(zl);
			y.resize(zl);
			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zl);
					gsl_spline_init(spline, &x[0], &y[0], zl);
					for (int zi = 0; zi < zl; zi++)
						Tab3[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[zi], NULL); // get derivatives dF/dz from spline interpolation
					gsl_spline_free(spline);
				}
			}

			x.resize(yl);
			y.resize(yl);
			for (int xi = 0; xi < xl; xi++){
				for (int zi = 0; zi < zl; zi++){
					for (int yi = 0; yi < yl; yi++){
						x[yi] = y_mi + yi*ydist;
						y[yi] = Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, yl);
					gsl_spline_init(spline, &x[0], &y[0], yl);
					for (int yi = 0; yi < yl; yi++)
						Tab12[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[yi], NULL); // get cross derivatives d2F/dxdy from spline coefficients
					gsl_spline_free(spline);
				}
			}

			x.resize(zl);
			y.resize(zl);
			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab1[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zl);
					gsl_spline_init(spline, &x[0], &y[0], zl);
					for (int zi = 0; zi < zl; zi++)
						Tab13[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[zi], NULL); // get cross derivatives d2F/dxdz from spline coefficients
					gsl_spline_free(spline);
				}
			}

			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab2[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zl);
					gsl_spline_init(spline, &x[0], &y[0], zl); // splineinterpolate dF/dy in z direction
					for (int zi = 0; zi < zl; zi++)
						Tab23[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[zi], NULL); // get cross derivatives d2F/dydz from spline coefficients
					gsl_spline_free(spline);
				}
			}

			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab12[INDEX_3D(xi, yi, zi, xl, yl, zl)];
					}
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zl);
					gsl_spline_init(spline, &x[0], &y[0], zl); // splineinterpolate d2F/dxdy in z direction
					for (int zi = 0; zi < zl; zi++)
						Tab123[INDEX_3D(xi, yi, zi, xl, yl, zl)] = gsl_spline_eval_deriv(spline, x[zi], NULL); // get cross derivatives d3F/dxdydz from spline coefficients
					gsl_spline_free(spline);
				}
			}

		};


		/**
		 * Calculate tricubic interpolation coefficients for a table column
		 *
		 * Calls TabField3::CalcDerivs and determines the interpolation coefficients with ::tricubic_get_coeff
		 *
		 * @param coeff Returns coefficients
		 * @param Tab Table column
		 */
		void PreInterpol(vector<vector<double> > &coeff, vector<double> &Tab){
			vector<double> Tab1(xl*yl*zl), Tab2(xl*yl*zl), Tab3(xl*yl*zl);
			vector<double> Tab12(xl*yl*zl), Tab13(xl*yl*zl), Tab23(xl*yl*zl), Tab123(xl*yl*zl);
			CalcDerivs(Tab,Tab1,Tab2,Tab3,Tab12,Tab13,Tab23,Tab123);

			// allocating space for the preinterpolation
			// The B*c are 5D arrays with xl x yl x zl x 4 x 4 fields
			coeff.resize((xl-1)*(yl-1)*(zl-1));

			int indx, indy, indz;
			double yyy[8], yyy1[8], yyy2[8], yyy3[8], yyy12[8], yyy13[8], yyy23[8], yyy123[8]; //rectangle with values at the 4 corners for interpolation

			for (indx=0; indx < xl-1; indx++){
				for (indy=0; indy < yl-1; indy++){
					for (indz=0; indz < zl-1; indz++){
						// fill cube with values and derivatives
						// order: see Lekien, Marsden: "Tricubic interpolation in three dimensions"
						int i3[8] ={INDEX_3D(indx  , indy  , indz  , xl, yl, zl),
									INDEX_3D(indx+1, indy  , indz  , xl, yl, zl),
									INDEX_3D(indx  , indy+1, indz  , xl, yl, zl),
									INDEX_3D(indx+1, indy+1, indz  , xl, yl, zl),
									INDEX_3D(indx  , indy  , indz+1, xl, yl, zl),
									INDEX_3D(indx+1, indy  , indz+1, xl, yl, zl),
									INDEX_3D(indx  , indy+1, indz+1, xl, yl, zl),
									INDEX_3D(indx+1, indy+1, indz+1, xl, yl, zl)};

						for (int i = 0; i < 8; i++){
							yyy[i]    = Tab[i3[i]];
							yyy1[i]   = Tab1[i3[i]]*xdist;
							yyy2[i]   = Tab2[i3[i]]*ydist;
							yyy3[i]   = Tab3[i3[i]]*zdist;
							yyy12[i]  = Tab12[i3[i]]*xdist*ydist;
							yyy13[i]  = Tab13[i3[i]]*xdist*zdist;
							yyy23[i]  = Tab23[i3[i]]*ydist*zdist;
							yyy123[i] = Tab123[i3[i]]*xdist*ydist*zdist;
						}

						// determine coefficients of interpolation
						int ind3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
						coeff[ind3].resize(64);
						tricubic_get_coeff(&coeff[ind3][0], yyy, yyy1, yyy2, yyy3, yyy12, yyy13, yyy23, yyy123);
					}
				}
			}
		};


		/**
		 * Get magnetic field scale factor for a specific time.
		 *
		 * Determined by TabField3::NullFieldTime, TabField3::RampUpTime, TabField3::FullFieldTime, TabField3::RampDownTime
		 *
		 * @param t Time
		 *
		 * @return Returns magnetic field scale factor
		 */
		double BFieldScale(double t){
			if (t < NullFieldTime || t >= NullFieldTime + RampUpTime + FullFieldTime + RampDownTime)
				return 0;
			else if (t >= NullFieldTime && t < NullFieldTime + RampUpTime){
				// ramping up field smoothly with cosine
				//result = (0.5 - 0.5*cos(pi*(t - CleaningTime - FillingTime)/RampUpTime)) * BFeldSkalGlobal;

				// linear ramp
				return (t - NullFieldTime)/RampUpTime;
			}
			else if (t >= NullFieldTime + RampUpTime && t < NullFieldTime + RampUpTime + FullFieldTime)
				return 1;
			else if (t >= NullFieldTime + RampUpTime + FullFieldTime && t < NullFieldTime + RampUpTime + FullFieldTime + RampDownTime){
				// ramping down field smoothly with cosine
				//result = (0.5 + 0.5*cos(pi*(t - (RampUpTime + CleaningTime + FillingTime + FullFieldTime)) / RampDownTime)) * BFeldSkalGlobal;

				// linear ramp
				return (1 - (t - RampUpTime - NullFieldTime - FullFieldTime)/RampDownTime);
			}
			else return 0;
		};
	public:
		/**
		 * Constructor.
		 *
		 * Calls TabField3::ReadTabFile, TabField3::CheckTab and for each column TabField3::PreInterpol
		 *
		 * @param tabfile Path of table file
		 * @param Bscale Magnetic field is always scaled by this factor
		 * @param Escale Electric field is always scaled by this factor
		 * @param aNullFieldTime Sets TabField3::NullFieldTime
		 * @param aRampUpTime Sets TabField3::RampUpTime
		 * @param aFullFieldTime Sets TabField3::FullFieldTime
		 * @param aRampDownTime Set TabField3::RampDownTime
		 */
		TabField3(const char *tabfile, double Bscale, double Escale,
				double aNullFieldTime, double aRampUpTime, double aFullFieldTime, double aRampDownTime){
			NullFieldTime = aNullFieldTime;
			RampUpTime = aRampUpTime;
			FullFieldTime = aFullFieldTime;
			RampDownTime = aRampDownTime;

			vector<double> BxTab, ByTab, BzTab;	// Bx/By/Bz values
			vector<double> VTab; // potential values
			ReadTabFile(tabfile,Bscale,Escale,BxTab,ByTab,BzTab,VTab); // open tabfile and read values into arrays

			CheckTab(BxTab,ByTab,BzTab,VTab); // print some info

			printf("Starting Preinterpolation ... ");
			float size = 0;
			if (BxTab.size() > 0){
				printf("Bx ... ");
				fflush(stdout);
				PreInterpol(Bxc,BxTab); // precalculate interpolation coefficients for B field
				size += float(Bxc.size()*64*sizeof(double)/1024/1024);
				BxTab.clear();
			}
			if (ByTab.size() > 0){
				printf("By ... ");
				fflush(stdout);
				PreInterpol(Byc,ByTab);
				size += float(Byc.size()*64*sizeof(double)/1024/1024);
				ByTab.clear();
			}
			if (BzTab.size() > 0){
				printf("Bz ... ");
				fflush(stdout);
				PreInterpol(Bzc,BzTab);
				size += float(Bzc.size()*64*sizeof(double)/1024/1024);
				BzTab.clear();
			}
			if (VTab.size() > 0){
				printf("V ... ");
				fflush(stdout);
				PreInterpol(Vc,VTab);
				size += float(Vc.size()*64*sizeof(double)/1024/1024);
				VTab.clear();
			}
			printf("Done (%.f MB)\n",size);
		}

		/**
		 * Get magnetic field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField3::x_mi, TabField3::xdist, TabField3::y_mi, TabField3::ydist, TabField3::z_mi, TabField3::zdist
		 * and evaluates the interpolation polynom ::tricubic_eval.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param B Returns magnetic field components B[0..2][0], their derivatives B[0..2][1..3], the absolute value B[3][0] and its derivatives B[3][1..3]
		 */
		void BField(long double x, long double y, long double z, long double t, long double B[4][4]){
			long double Bscale = BFieldScale(t);
			// get coordinate index
			int indx = (int)floor((x - x_mi)/xdist);
			int indy = (int)floor((y - y_mi)/ydist);
			int indz = (int)floor((z - z_mi)/zdist);
			if (Bscale != 0 && indx >= 0 && indx < xl - 1 && indy >= 0 && indy < yl - 1 && indz >= 0 && indz < zl - 1){
				// scale coordinates to unit cube
				x = (x - x_mi - indx*xdist)/xdist;
				y = (y - y_mi - indy*ydist)/ydist;
				z = (z - z_mi - indz*zdist)/zdist;
				int i3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
				// tricubic interpolation
				if (Bxc.size() > 0){
					B[0][0] += Bscale*tricubic_eval(&Bxc[i3][0], x, y, z);
					B[0][1] += Bscale*tricubic_eval(&Bxc[i3][0], x, y, z, 1, 0, 0)/xdist;
					B[0][2] += Bscale*tricubic_eval(&Bxc[i3][0], x, y, z, 0, 1, 0)/ydist;
					B[0][3] += Bscale*tricubic_eval(&Bxc[i3][0], x, y, z, 0, 0, 1)/zdist;
				}
				if (Byc.size() > 0){
					B[1][0] += Bscale*tricubic_eval(&Byc[i3][0], x, y, z);
					B[1][1] += Bscale*tricubic_eval(&Byc[i3][0], x, y, z, 1, 0, 0)/xdist;
					B[1][2] += Bscale*tricubic_eval(&Byc[i3][0], x, y, z, 0, 1, 0)/ydist;
					B[1][3] += Bscale*tricubic_eval(&Byc[i3][0], x, y, z, 0, 0, 1)/zdist;
				}
				if (Bzc.size() > 0){
					B[2][0] += Bscale*tricubic_eval(&Bzc[i3][0], x, y, z);
					B[2][1] += Bscale*tricubic_eval(&Bzc[i3][0], x, y, z, 1, 0, 0)/xdist;
					B[2][2] += Bscale*tricubic_eval(&Bzc[i3][0], x, y, z, 0, 1, 0)/ydist;
					B[2][3] += Bscale*tricubic_eval(&Bzc[i3][0], x, y, z, 0, 0, 1)/zdist;
				}
			}
		};

		/**
		 * Get electric field at a specific point.
		 *
		 * Searches the right interpolation coefficients by determining the indices from TabField3::x_mi, TabField3::xdist, TabField3::y_mi, TabField3::ydist, TabField3::z_mi, TabField3::zdist
		 * and evaluates the interpolation polynom ::tricubic_eval.
		 *
		 * @param x X coordinate where the field shall be evaluated
		 * @param y Y coordinate where the field shall be evaluated
		 * @param z Z coordinate where the field shall be evaluated
		 * @param t Time
		 * @param V Returns electric potential
		 * @param Ei Returns electric field (negative spatial derivatives of V)
		 */
		void EField(long double x, long double y, long double z, long double t, long double &V, long double Ei[3]){
			if (Vc.size() > 0 &&
				(x - x_mi)/xdist > 0 && (x - x_mi - xl*xdist)/xdist < 0 &&
				(y - y_mi)/ydist > 0 && (y - y_mi - yl*ydist)/ydist < 0 &&
				(z - z_mi)/zdist > 0 && (z - z_mi - zl*zdist)/zdist < 0){
				// get coordinate index
				int indx = (int)floor((x - x_mi)/xdist);
				int indy = (int)floor((y - y_mi)/ydist);
				int indz = (int)floor((z - z_mi)/zdist);
				// scale coordinates to unit cube
				x = (x - x_mi - indx*xdist)/xdist;
				y = (y - y_mi - indy*ydist)/ydist;
				z = (z - z_mi - indz*zdist)/zdist;
				int i3 = INDEX_3D(indx, indy, indz, xl-1, yl-1, zl-1);
				// tricubic interpolation
				V += tricubic_eval(&Vc[i3][0], x, y, z);
				Ei[0] += -tricubic_eval(&Vc[i3][0], x, y, z, 1, 0, 0)/xdist;
				Ei[1] += -tricubic_eval(&Vc[i3][0], x, y, z, 0, 1, 0)/ydist;
				Ei[2] += -tricubic_eval(&Vc[i3][0], x, y, z, 0, 0, 1)/zdist;
			}
		};
};

#endif // FIELD_3D_H_
