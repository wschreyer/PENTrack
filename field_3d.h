/**
 * \file
 * Tricubic interpolation of 3D field tables.
 */

#ifndef FIELD_3D_H_
#define FIELD_3D_H_

#include "libtricubic/tricubic.h"
#include "field.h"

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
		long double xdist; ///< distance between grid points in x direction
		long double ydist; ///< distance between grid points in y direction
		long double zdist; ///< distance between grid points in z direction
		long double x_mi; ///< lower x coordinate of cuboid grid
		long double y_mi; ///< lower y coordinate of cuboid grid
		long double z_mi; ///< lower z coordinate of cuboid grid
		NRMat3d<double[64]> *Bxc; ///< interpolation coefficients for magnetic x component
		NRMat3d<double[64]> *Byc; ///< interpolation coefficients for magnetic y component
		NRMat3d<double[64]> *Bzc; ///< interpolation coefficients for magnetic z component
		NRMat3d<double[64]> *Vc; ///< interpolation coefficients for electric potential
		long double NullFieldTime; ///< Time before magnetic field is ramped (passed by constructor)
		long double RampUpTime; ///< field is ramped linearly from 0 to 100% in this time (passed by constructor)
		long double FullFieldTime; ///< Time the field stays at 100% (passed by constructor)
		long double RampDownTime; ///< field is ramped down linearly from 100% to 0 in this time (passed by constructor)


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
		void ReadTabFile(const char *tabfile, long double Bscale, long double Escale,
				Mat3DDoub_O *&BxTab, Mat3DDoub_O *&ByTab, Mat3DDoub_O *&BzTab, Mat3DDoub_O *&VTab){
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
				BxTab = new Mat3DDoub(xl,yl,zl);
				getline(FIN,line);
			}
			if (line.find("BY") != string::npos){
				ByTab = new Mat3DDoub(xl,yl,zl);
				getline(FIN,line);
			}
			if (line.find("BZ") != string::npos){
				BzTab = new Mat3DDoub(xl,yl,zl);
				getline(FIN,line);
			}

			if (line.find("V") != string::npos){	// file contains potential?
				VTab = new Mat3DDoub(xl,yl,zl);
				getline(FIN,line);
			}

			if (!FIN || line.substr(0,2) != " 0"){
				printf("%s not found or corrupt! Exiting...\n",tabfile);
				exit(-1);
			}

			VecDoub xind(xl), yind(yl), zind(zl);
			int xi = 0, yi = 0, zi = -1, perc = 0;
			long double x, y, z, val;
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

				// status if read is displayed
				percent(xi*yi*zi ,0, xl*yl*zl, perc);

				xind[xi] = x;
				yind[yi] = y;
				zind[zi] = z;
				if (BxTab){
					FIN >> val;
					(*BxTab)[xi][yi][zi] = val*Bconv*Bscale;
				}
				if (ByTab){
					FIN >> val;
					(*ByTab)[xi][yi][zi] = val*Bconv*Bscale;
				}
				if (BzTab){
					FIN >> val;
					(*BzTab)[xi][yi][zi] = val*Bconv*Bscale;
				}
				if (VTab){
					FIN >> val;
					(*VTab)[xi][yi][zi] = val*Escale;
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
				delete BxTab;
				BxTab = NULL;
				delete ByTab;
				ByTab = NULL;
				delete BzTab;
				BzTab = NULL;
			}
			if (Escale == 0){
				delete VTab;
				VTab = NULL;
			}

			xdist = xind[1] - xind[0];
			ydist = yind[1] - yind[0];
			zdist = zind[1] - zind[0];
			x_mi = xind[0];
			y_mi = yind[0];
			z_mi = zind[0];

/*
			x_mi = -xind[xl-1];
			y_mi = -yind[yl-1];
			z_mi = zind[0];
			xl = 2*(xind[xl-1])/xdist + 1;
			yl = 2*(yind[yl-1])/ydist + 1;
			int zeroindx = -xind[0]/xdist;
			int zeroindy = -yind[0]/ydist;
			int newzeroindx = (xl - 1)/2;
			int newzeroindy = (yl - 1)/2;

			Mat3DDoub *temp;
			if (BxTab){
				temp = BxTab;
				BxTab = new Mat3DDoub(xl,yl,zl);
				for (xi = 0; xi < (xl - 1)/2; xi++){
					for (yi = 0; yi < (yl - 1)/2; yi++){
						for (zi = 0; zi < zl; zi++){
							(*BxTab)[newzeroindx+xi][newzeroindy+yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper right quadrant
							(*BxTab)[newzeroindx+xi][newzeroindy-yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*BxTab)[newzeroindx-xi][newzeroindy-yi][zi] = -(*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*BxTab)[newzeroindx-xi][newzeroindy+yi][zi] = -(*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper left quadrant
						}
					}
				}
				delete temp;
			}
			if (ByTab){
				temp = ByTab;
				ByTab = new Mat3DDoub(xl,yl,zl);
				for (xi = 0; xi < (xl - 1)/2; xi++){
					for (yi = 0; yi < (yl - 1)/2; yi++){
						for (zi = 0; zi < zl; zi++){
							(*ByTab)[newzeroindx+xi][newzeroindy+yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper right quadrant
							(*ByTab)[newzeroindx+xi][newzeroindy-yi][zi] = -(*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*ByTab)[newzeroindx-xi][newzeroindy-yi][zi] = -(*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*ByTab)[newzeroindx-xi][newzeroindy+yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper left quadrant
						}
					}
				}
				delete temp;
			}
			if (BzTab){
				temp = BzTab;
				BzTab = new Mat3DDoub(xl,yl,zl);
				for (xi = 0; xi < (xl - 1)/2; xi++){
					for (yi = 0; yi < (yl - 1)/2; yi++){
						for (zi = 0; zi < zl; zi++){
							(*BzTab)[newzeroindx+xi][newzeroindy+yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper right quadrant
							(*BzTab)[newzeroindx+xi][newzeroindy-yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*BzTab)[newzeroindx-xi][newzeroindy-yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // lower right quadrant
							(*BzTab)[newzeroindx-xi][newzeroindy+yi][zi] = (*temp)[zeroindx + xi][zeroindy + yi][zi]; // upper left quadrant
						}
					}
				}
				delete temp;
			}
*/
		};


		/**
		 * Print some information for each table column
		 *
		 * @param BxTab B_x column
		 * @param ByTab B_y column
		 * @param BzTab B_z column
		 * @param VTab V column
		 */
		void CheckTab(Mat3DDoub_I *BxTab, Mat3DDoub_I *ByTab, Mat3DDoub_I *BzTab, Mat3DDoub_I *VTab){
			//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
			printf("The arrays are %d by %d by %d.\n",xl,yl,zl);
			printf("The x values go from %LG to %LG\n",x_mi,x_mi + xdist*(xl-1));
			printf("The y values go from %LG to %LG\n",y_mi,y_mi + ydist*(yl-1));
			printf("The z values go from %LG to %LG.\n",z_mi,z_mi + zdist*(zl-1));
			printf("xdist = %LG ydist = %LG zdist = %LG\n",xdist,ydist,zdist);

			long double Babsmax = 0, Babsmin = 9e99, Babs;
			long double Vmax = 0, Vmin = 9e99;
			for (int j=0; j < xl; j++)
			{
				for (int k=0; k < yl; k++)
				{
					for (int l=0; l < zl; l++)
					{
						Babs = 0;
						if (BxTab) Babs += (*BxTab)[j][k][l] * (*BxTab)[j][k][l];
						if (ByTab) Babs += (*ByTab)[j][k][l] * (*ByTab)[j][k][l];
						if (BzTab) Babs += (*BzTab)[j][k][l] * (*BzTab)[j][k][l];
						Babsmax = max(sqrt(Babs),Babsmax);
						Babsmin = min(sqrt(Babs),Babsmin);
						if (VTab)
						{
							Vmax = max((*VTab)[j][k][l],Vmax);
							Vmin = min((*VTab)[j][k][l],Vmin);
						}

					}
				}
			}

			printf("The input table file has values of |B| from %LG T to %LG T and values of V from %LG V to %LG V\n",Babsmin,Babsmax,Vmin,Vmax);
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
		void CalcDerivs(Mat3DDoub_I &Tab, Mat3DDoub_O &Tab1, Mat3DDoub_O &Tab2, Mat3DDoub_O &Tab3,
						Mat3DDoub_O &Tab12, Mat3DDoub_O &Tab13, Mat3DDoub_O &Tab23, Mat3DDoub_O &Tab123)
		{
			VecDoub x(xl), y(xl);
			for (int zi = 0; zi < zl; zi++){
				for (int yi = 0; yi < yl; yi++){
					for (int xi = 0; xi < xl; xi++){
						x[xi] = x_mi + xi*xdist;
						y[xi] = Tab[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate field values F in x-direction
					for (int xi = 0; xi < xl-1; xi++)
						Tab1[xi][yi][zi] = (y[xi+1] - y[xi])/xdist - (spli.y2[xi+1] - spli.y2[xi])*xdist/6 - spli.y2[xi]*xdist/2; // get dF/dx from spline coefficients spli.y2
					Tab1[xl-1][yi][zi] = (y[xl-1] - y[xl-2])/xdist - (spli.y2[xl-1] - spli.y2[xl-2])*xdist/6 + spli.y2[xl-1]*xdist/2;
				}
			}

			x.resize(yl);
			y.resize(yl);
			for (int xi = 0; xi < xl; xi++){
				for (int zi = 0; zi < zl; zi++){
					for (int yi = 0; yi < yl; yi++){
						x[yi] = y_mi + yi*ydist;
						y[yi] = Tab[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate field values F in y-diretion
					for (int yi = 0; yi < yl-1; yi++)
						Tab2[xi][yi][zi] = (y[yi+1] - y[yi])/ydist - (spli.y2[yi+1] - spli.y2[yi])*ydist/6 - spli.y2[yi]*ydist/2; // get dF/dy from spline coefficients spli.y2
					Tab2[xi][yl-1][zi] = (y[yl-1] - y[yl-2])/ydist - (spli.y2[yl-1] - spli.y2[yl-2])*ydist/6 + spli.y2[yl-1]*ydist/2;
				}
			}

			x.resize(zl);
			y.resize(zl);
			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate field values F in z-diretion
					for (int zi = 0; zi < zl-1; zi++)
						Tab3[xi][yi][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get dF/dz from spline coefficients spli.y2
					Tab3[xi][yi][zl-1] = (y[zl-1] - y[zl-2])/zdist - (spli.y2[zl-1] - spli.y2[zl-2])*zdist/6 + spli.y2[zl-1]*zdist/2;
				}
			}

			x.resize(yl);
			y.resize(yl);
			for (int xi = 0; xi < xl; xi++){
				for (int zi = 0; zi < zl; zi++){
					for (int yi = 0; yi < yl; yi++){
						x[yi] = y_mi + yi*ydist;
						y[yi] = Tab1[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate dF/dx in y direction
					for (int yi = 0; yi < yl-1; yi++)
						Tab12[xi][yi][zi] = (y[yi+1] - y[yi])/ydist - (spli.y2[yi+1] - spli.y2[yi])*ydist/6 - spli.y2[yi]*ydist/2; // get cross derivatives d2F/dxdy from spline coefficients
					Tab12[xi][yl-1][zi] = (y[yl-1] - y[yl-2])/ydist - (spli.y2[yl-1] - spli.y2[yl-2])*ydist/6 + spli.y2[yl-1]*ydist/2;
				}
			}

			x.resize(zl);
			y.resize(zl);
			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab1[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate dF/dx in z direction
					for (int zi = 0; zi < zl-1; zi++)
						Tab13[xi][yi][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get cross derivatives d2F/dxdz from spline coefficients
					Tab13[xi][yi][zl-1] = (y[zl-1] - y[zl-2])/zdist - (spli.y2[zl-1] - spli.y2[zl-2])*zdist/6 + spli.y2[zl-1]*zdist/2;
				}
			}

			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab2[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate dF/dy in z direction
					for (int zi = 0; zi < zl-1; zi++)
						Tab23[xi][yi][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get cross derivatives d2F/dydz from spline coefficients
					Tab23[xi][yi][zl-1] = (y[zl-1] - y[zl-2])/zdist - (spli.y2[zl-1] - spli.y2[zl-2])*zdist/6 + spli.y2[zl-1]*zdist/2;
				}
			}

			for (int xi = 0; xi < xl; xi++){
				for (int yi = 0; yi < yl; yi++){
					for (int zi = 0; zi < zl; zi++){
						x[zi] = z_mi + zi*zdist;
						y[zi] = Tab12[xi][yi][zi];
					}
					Spline_interp spli(x,y); // splineinterpolate d2F/dxdy in z direction
					for (int zi = 0; zi < zl-1; zi++)
						Tab123[xi][yi][zi] = (y[zi+1] - y[zi])/zdist - (spli.y2[zi+1] - spli.y2[zi])*zdist/6 - spli.y2[zi]*zdist/2; // get cross derivatives d3F/dxdydz from spline coefficients
					Tab123[xi][yi][zl-1] = (y[zl-1] - y[zl-2])/zdist - (spli.y2[zl-1] - spli.y2[zl-2])*zdist/6 + spli.y2[zl-1]*zdist/2;
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
		void PreInterpol(NRMat3d<double[64]> *&coeff, Mat3DDoub_I &Tab){
			Mat3DDoub Tab1(xl,yl,zl), Tab2(xl,yl,zl), Tab3(xl,yl,zl);
			Mat3DDoub Tab12(xl,yl,zl), Tab13(xl,yl,zl), Tab23(xl,yl,zl), Tab123(xl,yl,zl);
			CalcDerivs(Tab,Tab1,Tab2,Tab3,Tab12,Tab13,Tab23,Tab123);

			// allocating space for the preinterpolation
			// The B*c are 5D arrays with xl x yl x zl x 4 x 4 fields
			coeff = new NRMat3d<double[64]>(xl-1,yl-1,zl-1);

			int indx, indy, indz;
			double yyy[8], yyy1[8], yyy2[8], yyy3[8], yyy12[8], yyy13[8], yyy23[8], yyy123[8]; //rectangle with values at the 4 corners for interpolation

			for (indx=0; indx < xl-1; indx++){
				for (indy=0; indy < yl-1; indy++){
					for (indz=0; indz < zl-1; indz++){
						// fill cube with values and derivatives
						// order: see Lekien, Marsden: "Tricubic interpolation in three dimensions"
						yyy[0] = Tab[indx][indy][indz];
						yyy[1] = Tab[indx+1][indy][indz];
						yyy[2] = Tab[indx][indy+1][indz];
						yyy[3] = Tab[indx+1][indy+1][indz];
						yyy[4] = Tab[indx][indy][indz+1];
						yyy[5] = Tab[indx+1][indy][indz+1];
						yyy[6] = Tab[indx][indy+1][indz+1];
						yyy[7] = Tab[indx+1][indy+1][indz+1];

						yyy1[0] = Tab1[indx][indy][indz]*xdist;
						yyy1[1] = Tab1[indx+1][indy][indz]*xdist;
						yyy1[2] = Tab1[indx][indy+1][indz]*xdist;
						yyy1[3] = Tab1[indx+1][indy+1][indz]*xdist;
						yyy1[4] = Tab1[indx][indy][indz+1]*xdist;
						yyy1[5] = Tab1[indx+1][indy][indz+1]*xdist;
						yyy1[6] = Tab1[indx][indy+1][indz+1]*xdist;
						yyy1[7] = Tab1[indx+1][indy+1][indz+1]*xdist;

						yyy2[0] = Tab2[indx][indy][indz]*ydist;
						yyy2[1] = Tab2[indx+1][indy][indz]*ydist;
						yyy2[2] = Tab2[indx][indy+1][indz]*ydist;
						yyy2[3] = Tab2[indx+1][indy+1][indz]*ydist;
						yyy2[4] = Tab2[indx][indy][indz+1]*ydist;
						yyy2[5] = Tab2[indx+1][indy][indz+1]*ydist;
						yyy2[6] = Tab2[indx][indy+1][indz+1]*ydist;
						yyy2[7] = Tab2[indx+1][indy+1][indz+1]*ydist;

						yyy3[0] = Tab3[indx][indy][indz]*zdist;
						yyy3[1] = Tab3[indx+1][indy][indz]*zdist;
						yyy3[2] = Tab3[indx][indy+1][indz]*zdist;
						yyy3[3] = Tab3[indx+1][indy+1][indz]*zdist;
						yyy3[4] = Tab3[indx][indy][indz+1]*zdist;
						yyy3[5] = Tab3[indx+1][indy][indz+1]*zdist;
						yyy3[6] = Tab3[indx][indy+1][indz+1]*zdist;
						yyy3[7] = Tab3[indx+1][indy+1][indz+1]*zdist;

						yyy12[0] = Tab12[indx][indy][indz]*xdist*ydist;
						yyy12[1] = Tab12[indx+1][indy][indz]*xdist*ydist;
						yyy12[2] = Tab12[indx][indy+1][indz]*xdist*ydist;
						yyy12[3] = Tab12[indx+1][indy+1][indz]*xdist*ydist;
						yyy12[4] = Tab12[indx][indy][indz+1]*xdist*ydist;
						yyy12[5] = Tab12[indx+1][indy][indz+1]*xdist*ydist;
						yyy12[6] = Tab12[indx][indy+1][indz+1]*xdist*ydist;
						yyy12[7] = Tab12[indx+1][indy+1][indz+1]*xdist*ydist;

						yyy13[0] = Tab13[indx][indy][indz]*xdist*zdist;
						yyy13[1] = Tab13[indx+1][indy][indz]*xdist*zdist;
						yyy13[2] = Tab13[indx][indy+1][indz]*xdist*zdist;
						yyy13[3] = Tab13[indx+1][indy+1][indz]*xdist*zdist;
						yyy13[4] = Tab13[indx][indy][indz+1]*xdist*zdist;
						yyy13[5] = Tab13[indx+1][indy][indz+1]*xdist*zdist;
						yyy13[6] = Tab13[indx][indy+1][indz+1]*xdist*zdist;
						yyy13[7] = Tab13[indx+1][indy+1][indz+1]*xdist*zdist;

						yyy23[0] = Tab23[indx][indy][indz]*ydist*zdist;
						yyy23[1] = Tab23[indx+1][indy][indz]*ydist*zdist;
						yyy23[2] = Tab23[indx][indy+1][indz]*ydist*zdist;
						yyy23[3] = Tab23[indx+1][indy+1][indz]*ydist*zdist;
						yyy23[4] = Tab23[indx][indy][indz+1]*ydist*zdist;
						yyy23[5] = Tab23[indx+1][indy][indz+1]*ydist*zdist;
						yyy23[6] = Tab23[indx][indy+1][indz+1]*ydist*zdist;
						yyy23[7] = Tab23[indx+1][indy+1][indz+1]*ydist*zdist;

						yyy123[0] = Tab123[indx][indy][indz]*xdist*ydist*zdist;
						yyy123[1] = Tab123[indx+1][indy][indz]*xdist*ydist*zdist;
						yyy123[2] = Tab123[indx][indy+1][indz]*xdist*ydist*zdist;
						yyy123[3] = Tab123[indx+1][indy+1][indz]*xdist*ydist*zdist;
						yyy123[4] = Tab123[indx][indy][indz+1]*xdist*ydist*zdist;
						yyy123[5] = Tab123[indx+1][indy][indz+1]*xdist*ydist*zdist;
						yyy123[6] = Tab123[indx][indy+1][indz+1]*xdist*ydist*zdist;
						yyy123[7] = Tab123[indx+1][indy+1][indz+1]*xdist*ydist*zdist;

						// determine coefficients of interpolation
						tricubic_get_coeff((*coeff)[indx][indy][indz], yyy, yyy1, yyy2, yyy3, yyy12, yyy13, yyy23, yyy123);
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
		long double BFieldScale(long double t){
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
		TabField3(const char *tabfile, long double Bscale, long double Escale,
				long double aNullFieldTime, long double aRampUpTime, long double aFullFieldTime, long double aRampDownTime){
			NullFieldTime = aNullFieldTime;
			RampUpTime = aRampUpTime;
			FullFieldTime = aFullFieldTime;
			RampDownTime = aRampDownTime;
			Bxc = Byc = Bzc = Vc = NULL;

			Mat3DDoub *BxTab = NULL, *ByTab = NULL, *BzTab = NULL;	// Bx/By/Bz values
			Mat3DDoub *VTab = NULL; // potential values
			ReadTabFile(tabfile,Bscale,Escale,BxTab,ByTab,BzTab,VTab); // open tabfile and read values into arrays

			CheckTab(BxTab,ByTab,BzTab,VTab); // print some info

			printf("Starting Preinterpolation ... ");
			float size = 0;
			if (BxTab){
				printf("Bx ... ");
				fflush(stdout);
				PreInterpol(Bxc,*BxTab); // precalculate interpolation coefficients for B field
				size += float(Bxc->dim1()*Bxc->dim2()*Bxc->dim3()*64*sizeof(double)/1024/1024);
				delete BxTab;
			}
			if (ByTab){
				printf("By ... ");
				fflush(stdout);
				PreInterpol(Byc,*ByTab);
				size += float(Byc->dim1()*Byc->dim2()*Byc->dim3()*64*sizeof(double)/1024/1024);
				delete ByTab;
			}
			if (BzTab){
				printf("Bz ... ");
				fflush(stdout);
				PreInterpol(Bzc,*BzTab);
				size += float(Bzc->dim1()*Bzc->dim2()*Bzc->dim3()*64*sizeof(double)/1024/1024);
				delete BzTab;
			}
			if (VTab){
				printf("V ... ");
				fflush(stdout);
				PreInterpol(Vc,*VTab);
				size += float(Vc->dim1()*Vc->dim2()*Vc->dim3()*64*sizeof(double)/1024/1024);
				delete VTab;
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
			int indx = (int)((x - x_mi)/xdist);
			int indy = (int)((y - y_mi)/ydist);
			int indz = (int)((z - z_mi)/zdist);
			if (Bscale != 0 && indx >= 0 && indx < xl - 1 && indy >= 0 && indy < yl - 1 && indz >= 0 && indz < zl - 1){
				// scale coordinates to unit cube
				x = (x - x_mi - indx*xdist)/xdist;
				y = (y - y_mi - indy*ydist)/ydist;
				z = (z - z_mi - indz*zdist)/zdist;
				// tricubic interpolation
				if (Bxc){
					B[0][0] += Bscale*tricubic_eval((*Bxc)[indx][indy][indz], x, y, z);
					B[0][1] += Bscale*tricubic_eval((*Bxc)[indx][indy][indz], x, y, z, 1, 0, 0)/xdist;
					B[0][2] += Bscale*tricubic_eval((*Bxc)[indx][indy][indz], x, y, z, 0, 1, 0)/ydist;
					B[0][3] += Bscale*tricubic_eval((*Bxc)[indx][indy][indz], x, y, z, 0, 0, 1)/zdist;
				}
				if (Byc){
					B[1][0] += Bscale*tricubic_eval((*Byc)[indx][indy][indz], x, y, z);
					B[1][1] += Bscale*tricubic_eval((*Byc)[indx][indy][indz], x, y, z, 1, 0, 0)/xdist;
					B[1][2] += Bscale*tricubic_eval((*Byc)[indx][indy][indz], x, y, z, 0, 1, 0)/ydist;
					B[1][3] += Bscale*tricubic_eval((*Byc)[indx][indy][indz], x, y, z, 0, 0, 1)/zdist;
				}
				if (Bzc){
					B[2][0] += Bscale*tricubic_eval((*Bzc)[indx][indy][indz], x, y, z);
					B[2][1] += Bscale*tricubic_eval((*Bzc)[indx][indy][indz], x, y, z, 1, 0, 0)/xdist;
					B[2][2] += Bscale*tricubic_eval((*Bzc)[indx][indy][indz], x, y, z, 0, 1, 0)/ydist;
					B[2][3] += Bscale*tricubic_eval((*Bzc)[indx][indy][indz], x, y, z, 0, 0, 1)/zdist;
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
			if (Vc &&
				(x - x_mi)/xdist > 0 && (x - x_mi - xl*xdist)/xdist < 0 &&
				(y - y_mi)/ydist > 0 && (y - y_mi - yl*ydist)/ydist < 0 &&
				(z - z_mi)/zdist > 0 && (z - z_mi - zl*zdist)/zdist < 0){
				// get coordinate index
				int indx = (int)((x - x_mi)/xdist);
				int indy = (int)((y - y_mi)/ydist);
				int indz = (int)((z - z_mi)/zdist);
				// scale coordinates to unit cube
				x = (x - x_mi - indx*xdist)/xdist;
				y = (y - y_mi - indy*ydist)/ydist;
				z = (z - z_mi - indz*zdist)/zdist;
				// tricubic interpolation
				V += tricubic_eval((*Vc)[indx][indy][indz], x, y, z);
				Ei[0] += -tricubic_eval((*Vc)[indx][indy][indz], x, y, z, 1, 0, 0)/xdist;
				Ei[1] += -tricubic_eval((*Vc)[indx][indy][indz], x, y, z, 0, 1, 0)/ydist;
				Ei[2] += -tricubic_eval((*Vc)[indx][indy][indz], x, y, z, 0, 0, 1)/zdist;
			}
		};
};


#endif // FIELD_3D_H_
