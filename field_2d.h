/**
 * \file
 * Bicubic interpolation of axisymmetric field tables.
 */

#ifndef FIELD_2D_H_
#define FIELD_2D_H_

#include "field.h"

#include "../interp2d/interp2d.h"

/**
 * Translate VECTOR (not point) components from cylindrical to cartesian coordinates.
 *
 * @param v_r Radial component of vector
 * @param v_phi Azimuthal component of vector
 * @param phi Azimuth of vector origin
 * @param v_x Returns x component of vector
 * @param v_y Returns y component of vector
 */
void CylToCart(double v_r, double v_phi, double phi, double &v_x, double &v_y){
	v_x = v_r*cos(phi) - v_phi*sin(phi);
	v_y = v_r*sin(phi) + v_phi*cos(phi);
}

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
		vector<double> rind, zind, BrTab, BphiTab, BzTab, ErTab, EphiTab, EzTab, VTab;
		interp2d *Brc, *Bphic, *Bzc, *Erc, *Ephic, *Ezc, *Vc;
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
		void ReadTabFile(const char *tabfile, double Bscale, double Escale){
			ifstream FIN(tabfile, ifstream::in);
			if (!FIN.is_open()){
				printf("\nCould not open %s!\n",tabfile);
				exit(-1);
			}
			printf("\nReading %s ",tabfile);
			int intval;
			string line;
			FIN >> m >> intval >> n;
			rind.resize(m);
			zind.resize(n);

			getline(FIN,line);
			getline(FIN,line);
			getline(FIN,line);
			bool skipy = true;
			if (line.substr(0,12) == " 2 Y [LENGU]")  skipy = false;
			getline(FIN,line);
			if (!skipy) getline(FIN,line);

			if (line.find("RBX") != string::npos){
				BrTab.resize(m*n);
				getline(FIN,line);
			}
			if (line.find("RBY") != string::npos){
				BphiTab.resize(m*n);
				getline(FIN,line);
			}
			if (line.find("RBZ") != string::npos){
				BzTab.resize(m*n);
				getline(FIN,line);
			}

			if (line.find("EX") != string::npos){
				ErTab.resize(m*n);
				getline(FIN,line);
			}
			if (line.find("EY") != string::npos){
				EphiTab.resize(m*n);
				getline(FIN,line);
			}
			if (line.find("EZ") != string::npos){
				EzTab.resize(m*n);
				getline(FIN,line);
			}

			if (line.find("RV") != string::npos){	// file contains potential?
				VTab.resize(m*n);
				getline(FIN,line);
			}

			if (!FIN || line.substr(0,2) != " 0"){
				printf("%s not found or corrupt! Exiting...\n",tabfile);
				exit(-1);
			}

			int ri = 0,zi = -1, perc = 0;
			double r, z, val;
			while (FIN.good()){
				FIN >> r;
				if (!skipy) FIN >> val;
				FIN >> z;
				if (!FIN) break;
				r *= lengthconv;
				z *= lengthconv;
				if (zi >= 0 && z < zind[zi]){
					ri++;
					zi = 0;
				}
				else zi++;

				// status if read is displayed
				PrintPercent((float)zi*ri/(n*m), perc);

				rind[ri] = r;
				zind[zi] = z;
				int i2 = INDEX_2D(ri, zi, m, n);
				if (BrTab.size() > 0){
					FIN >> val;
					BrTab[i2] = val*Bconv*Bscale;
				}
				if (BphiTab.size() > 0){
					FIN >> val;
					BphiTab[i2] = val*Bconv*Bscale;
				}
				if (BzTab.size() > 0){
					FIN >> val;
					BzTab[i2] = val*Bconv*Bscale;
				}
				if (ErTab.size() > 0){
					FIN >> val;
					ErTab[i2] = val*Econv*Escale;
				}
				if (EphiTab.size() > 0){
					FIN >> val;
					EphiTab[i2] = val*Econv*Escale;
				}
				if (EzTab.size() > 0){
					FIN >> val;
					EzTab[i2] = val*Econv*Escale;
				}
				if (VTab.size() > 0){
					FIN >> val;
					VTab[i2] = val*Escale;
				}
				FIN >> ws;
			}

			if (Bscale == 0){
				BrTab.clear();
				BphiTab.clear();
				BzTab.clear();
			}
			if (Escale == 0){
				ErTab.clear();
				EphiTab.clear();
				EzTab.clear();
				VTab.clear();
			}

			printf("\n");
			if (ri+1 != (int)m || zi+1 != (int)n){
				printf("The header says the size is %u by %u, actually it is %i by %i! Exiting...\n", m, n, ri+1, zi+1);
				exit(-1);
			}
			FIN.close();
		};


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
		void CheckTab(){
			//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB
			printf("The arrays are %u by %u.\n",m,n);
			printf("The r values go from %G to %G\n", rind[0], rind[m-1]);
			printf("The z values go from %G to %G.\n", zind[0], zind[n-1]);
			printf("rdist = %G zdist = %G\n", rind[1] - rind[0], zind[1] - zind[0]);

			double Babsmax = 0, Babsmin = 9e99, Babs;
			double Vmax = 0, Vmin = 9e99;
			for (int j=0; j < m; j++)
			{
				for (int k=0; k < n; k++)
				{
					Babs = 0;
					int i2 = INDEX_2D(j,k,m,n);
					if (BrTab.size() > 0) Babs += BrTab[i2]*BrTab[i2];
					if (BphiTab.size() > 0) Babs += BphiTab[i2]*BphiTab[i2];
					if (BzTab.size() > 0) Babs += BzTab[i2]*BzTab[i2];
					Babsmax = max(sqrt(Babs),Babsmax);
					Babsmin = min(sqrt(Babs),Babsmin);
					if (VTab.size() > 0)
					{
						Vmax = max(VTab[i2],Vmax);
						Vmin = min(VTab[i2],Vmin);
					}

				}
			}

			printf("The input table file has values of |B| from %G T to %G T and values of V from %G V to %G V\n",Babsmin,Babsmax,Vmin,Vmax);
		};


		/**
		 * Get magnetic field scale factor for a specific time.
		 *
		 * Determined by TabField::NullFieldTime, TabField::RampUpTime, TabField::FullFieldTime, TabField::RampDownTime
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
				//return 0.5 - 0.5*cos(pi*(t - NullFieldTime)/RampUpTime);

				// linear ramp
				return (t - NullFieldTime)/RampUpTime;
			}
			else if (t >= NullFieldTime + RampUpTime && t < NullFieldTime + RampUpTime + FullFieldTime)
				return 1;
			else if (t >= NullFieldTime + RampUpTime + FullFieldTime && t < NullFieldTime + RampUpTime + FullFieldTime + RampDownTime){
				// ramping down field smoothly with cosine
				//return 0.5 + 0.5*cos(pi*(t - (RampUpTime + NullFieldTime + FullFieldTime)) / RampDownTime);

				// linear ramp
				return (1 - (t - RampUpTime - NullFieldTime - FullFieldTime)/RampDownTime);
			}
			else return 0;
		};
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
				double aNullFieldTime, double aRampUpTime, double aFullFieldTime, double aRampDownTime){
			NullFieldTime = aNullFieldTime;
			RampUpTime = aRampUpTime;
			FullFieldTime = aFullFieldTime;
			RampDownTime = aRampDownTime;

			ReadTabFile(tabfile,Bscale,Escale); // open tabfile and read values into arrays

			CheckTab(); // print some info

			printf("Starting Preinterpolation ... ");
			if (BrTab.size() > 0){
				cout << "Br ... ";
				cout.flush();
				Brc = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Brc, &rind[0], &zind[0], &BrTab[0], m, n);
			}
			if (BphiTab.size() > 0){
				cout << "Bhi ... ";
				cout.flush();
				Bphic = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Bphic, &rind[0], &zind[0], &BphiTab[0], m, n);
			}
			if (BzTab.size() > 0){
				cout << "Bz ... ";
				cout.flush();
				Bzc = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Bzc, &rind[0], &zind[0], &BzTab[0], m, n);
			}
			if (VTab.size() > 0){
				cout << "V ... ";
				cout.flush();
				Vc = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Vc, &rind[0], &zind[0], &VTab[0], m, n);
				ErTab.clear(); // if there is a potential, we don't need the E-vector
				EphiTab.clear();
				EzTab.clear();
			}
			if (ErTab.size() > 0){
				cout << "Er ... ";
				cout.flush();
				Erc = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Erc, &rind[0], &zind[0], &ErTab[0], m, n);
			}
			if (EphiTab.size() > 0){
				cout << "Ephi ... ";
				cout.flush();
				Ephic = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Ephic, &rind[0], &zind[0], &EphiTab[0], m, n);
			}
			if (EzTab.size() > 0){
				cout << "Ez ... ";
				cout.flush();
				Ezc = interp2d_alloc(interp2d_bicubic, m, n);
				interp2d_init(Ezc, &rind[0], &zind[0], &EzTab[0], m, n);
			}
			printf("Done (%.f MB)\n\n", 4*float(BrTab.size() + BphiTab.size() + BzTab.size()
										 + ErTab.size() + EzTab.size() + VTab.size())
										*sizeof(double)/1024/1024);
		};

		~TabField(){
			if (BrTab.size() > 0)
				interp2d_free(Brc);
			if (BphiTab.size() > 0)
				interp2d_free(Bphic);
			if (BzTab.size() > 0)
				interp2d_free(Bzc);
			if (ErTab.size() > 0)
				interp2d_free(Erc);
			if (EphiTab.size() > 0)
				interp2d_free(Ephic);
			if (EzTab.size() > 0)
				interp2d_free(Ezc);
			if (VTab.size() > 0)
				interp2d_free(Vc);
		}

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
		void BField(long double x, long double y, long double z, long double t, long double B[4][4]){
			double r = sqrt(x*x+y*y);
			double Bscale = BFieldScale(t);
			if (Bscale != 0 && r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
				// bicubic interpolation
				double Br = 0, dBrdr = 0, dBrdz = 0, Bphi = 0, dBphidr = 0, dBphidz = 0, dBzdr = 0;
				double Bx = 0, By = 0, Bz = 0, dBxdz = 0, dBydz = 0, dBzdz = 0;
				double phi = atan2(y,x);
				if (BrTab.size() > 0){
					Br = interp2d_eval(Brc, &rind[0], &zind[0], &BrTab[0], r, z, NULL, NULL);
					dBrdr = interp2d_eval_deriv_x(Brc, &rind[0], &zind[0], &BrTab[0], r, z, NULL, NULL);
					dBrdz = interp2d_eval_deriv_y(Brc, &rind[0], &zind[0], &BrTab[0], r, z, NULL, NULL);
				}
				if (BphiTab.size() > 0){
					Bphi = interp2d_eval(Bphic, &rind[0], &zind[0], &BphiTab[0], r, z, NULL, NULL);
					dBphidr = interp2d_eval_deriv_x(Bphic, &rind[0], &zind[0], &BphiTab[0], r, z, NULL, NULL);
					dBphidz = interp2d_eval_deriv_y(Bphic, &rind[0], &zind[0], &BphiTab[0], r, z, NULL, NULL);
				}
				CylToCart(Br,Bphi,phi,Bx,By);
				B[0][0] += Bx*Bscale;
				B[1][0] += By*Bscale;
				if (r > 0){
					B[0][1] += Bscale*(dBrdr*cos(phi)*cos(phi) - dBphidr*cos(phi)*sin(phi) + (Br*sin(phi)*sin(phi) + Bphi*cos(phi)*sin(phi))/r);
					B[0][2] += Bscale*(dBrdr*cos(phi)*sin(phi) - dBphidr*sin(phi)*sin(phi) - (Br*cos(phi)*sin(phi) + Bphi*cos(phi)*cos(phi))/r);
					B[1][1] += Bscale*(dBrdr*cos(phi)*sin(phi) + dBphidr*cos(phi)*cos(phi) - (Br*cos(phi)*sin(phi) - Bphi*sin(phi)*sin(phi))/r);
					B[1][2] += Bscale*(dBrdr*sin(phi)*sin(phi) + dBphidr*cos(phi)*sin(phi) + (Br*cos(phi)*cos(phi) - Bphi*cos(phi)*sin(phi))/r);
				}
				CylToCart(dBrdz,dBphidz,phi,dBxdz,dBydz);
				B[0][3] += dBxdz*Bscale;
				B[1][3] += dBydz*Bscale;
				if (BzTab.size() > 0){
					Bz = interp2d_eval(Bzc, &rind[0], &zind[0], &BzTab[0], r, z, NULL, NULL);
					dBzdr = interp2d_eval_deriv_x(Bzc, &rind[0], &zind[0], &BzTab[0], r, z, NULL, NULL);
					dBzdz = interp2d_eval_deriv_y(Bzc, &rind[0], &zind[0], &BzTab[0], r, z, NULL, NULL);
				}
				B[2][0] += Bz*Bscale;
				B[2][1] += dBzdr*cos(phi)*Bscale;
				B[2][2] += dBzdr*sin(phi)*Bscale;
				B[2][3] += dBzdz*Bscale;
			}
		};


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
		void EField(long double x, long double y, long double z, long double t, long double &V, long double Ei[3]){
			double Vloc, dVdrj[3];
			if (VTab.size() > 0){ // prefer E-field from potential over pure E-field interpolation
				double r = sqrt(x*x+y*y);
				if (r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
					// bicubic interpolation
					Vloc = interp2d_eval(Vc, &rind[0], &zind[0], &VTab[0], r, z, NULL, NULL);
					dVdrj[0] = interp2d_eval_deriv_x(Vc, &rind[0], &zind[0], &VTab[0], r, z, NULL, NULL);
					dVdrj[2] = interp2d_eval_deriv_y(Vc, &rind[0], &zind[0], &VTab[0], r, z, NULL, NULL);
					double phi = atan2(y,x);
					V += Vloc;
					Ei[0] += -dVdrj[0]*cos(phi);
					Ei[1] += -dVdrj[0]*sin(phi);
					Ei[2] += -dVdrj[2];
				}
			}
			else if (ErTab.size() > 0 || EphiTab.size() > 0 || EzTab.size() > 0){
				double r = sqrt(x*x+y*y);
				if (r >= rind[0] && r <= rind[m-1] && z >= zind[0] && z <= zind[n-1]){
					// bicubic interpolation
					double Er = 0, Ephi = 0, Ex = 0, Ey = 0, Ez = 0;
					double phi = atan2(y,x);
					if (ErTab.size() > 0)
						Er = interp2d_eval(Erc, &rind[0], &zind[0], &ErTab[0], r, z, NULL, NULL);
					if (EphiTab.size() > 0)
						Ephi = interp2d_eval(Ephic, &rind[0], &zind[0], &EphiTab[0], r, z, NULL, NULL);
					if (EzTab.size() > 0)
						Ez = interp2d_eval(Ezc, &rind[0], &zind[0], &EzTab[0], r, z, NULL, NULL);
					CylToCart(Er,Ephi,phi,Ex,Ey);
					Ei[0] += Ex;
					Ei[1] += Ey;
					Ei[2] += Ez;
				}
			}
		};
};


#endif // FIELD_2D_H_
