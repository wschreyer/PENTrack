/**
 * \file
 * Manage all the fields.
 */

#ifndef FIELDS_H_
#define FIELDS_H_

#include "field.h"
#include "field_2d.h"
#include "field_3d.h"
#include "conductor.h"

/**
 * Contains list of all fields (2D/3D-maps, conductors, ...).
 */
struct TFieldManager{
	public:
		vector<TField*> fields; ///< list of fields
		int FieldOscillation; ///< If =1 field oscillation is turned on
		long double OscillationFraction; ///< Field oscillation amplitude
		long double OscillationFrequency; ///< Field oscillation frequency
		
		/**
		 * Constructor.
		 *
		 * Reads [FIELDS] section of configuration file and loads all field maps/conductors given there
		 *
		 * @param conf TConfig map containing field options
		 * @param aFieldOscillation Turn on field oscillations
		 * @param aOscillationFraction Amplitude of field oscillation
		 * @param aOscillationFrequency Frequency of field oscillation
		 */
		TFieldManager(TConfig &conf, int aFieldOscillation = 0, long double aOscillationFraction = 0, long double aOscillationFrequency = 0):
					FieldOscillation(aFieldOscillation), OscillationFraction(aOscillationFraction), OscillationFrequency(aOscillationFrequency){
			for (map<string, string>::iterator i = conf["FIELDS"].begin(); i != conf["FIELDS"].end(); i++){
				string type;
				string ft;
				long double Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime;
				long double Ibar, p1, p2, p3, p4, p5, p6;
				TField *f = NULL;
				istringstream ss(i->second);

				if (i->first == "2Dtable"){
					ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime;
					if (ss)
						f = new TabField(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime);
				}
				else if (i->first == "3Dtable"){
					ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime;
					if (ss)
						f = new TabField3(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime);
				}
				else if (i->first == "InfiniteWireZ"){
					ss >> Ibar >> p1 >> p2;
					if (ss)
						f = new TInfiniteWireZ(p1, p2, Ibar);
				}
				else if (i->first == "InfiniteWireZCenter"){
					ss >> Ibar;
					if (ss)
						f = new TInfiniteWireZCenter(Ibar);
				}
				else if (i->first == "FiniteWire"){
					ss >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6;
					if (ss)
						f = new TFiniteWire(p1, p2, p3, p4, p5, p6, Ibar);
				}
				else if (i->first == "FiniteWireX"){
					ss >> Ibar >> p1 >> p2 >> p3;
					if (ss)
						f = new TFiniteWireX(p1, p2, p3, Ibar);
				}
				else if (i->first == "FiniteWireY"){
					ss >> Ibar >> p1 >> p2 >> p3;
					if (ss)
						f = new TFiniteWireY(p1, p2, p3, Ibar);
				}
				else if (i->first == "FiniteWireZ"){
					ss >> Ibar >> p1 >> p2 >> p3 >> p4;
					if (ss)
						f = new TFiniteWireZ(p1, p2, p3, p4, Ibar);
				}
				else if (i->first == "FiniteWireZCenter"){
					ss >> Ibar >> p1 >> p2;
					if (ss)
						f = new TFiniteWireZCenter(p1, p2, Ibar);
				}
				else if (i->first == "FullRacetrack"){
					ss >> Ibar >> p1 >> p2 >> p3;
					if (ss)
						f = new TFullRacetrack(p1, p2, p3, Ibar);
				}

				if (f)
					fields.push_back(f);
				else{
					cout << "\nCould not load field """ << i->first << """! Did you enter invalid parameters?\n";
					exit(-1);
				}
			}
		};
		
		/**
		 * Destructor, delete all fields.
		 */
		~TFieldManager(){
			for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++)
				delete (*i);
		};
	
		/**
		 * Calculate magnetic field at a given position and time.
		 *
		 * Chooses the right map for this position and adds racetrack fields.
		 * field matrix B:
		 *	Bx,		dBxdx,	dBxdy,	dBxdz;
		 *	By,		dBydx,	dBydy,	dBydz;
		 *	Bz,		dBzdx,	dBzdy,	dBzdz;
		 *	Babs,	dBdx,	dBdy,	dBdz;
		 *
		 * @param x Cartesian x coordinate
		 * @param y Cartesian y coordinate
		 * @param z Cartesian z coordinate
		 * @param t Time
		 * @param B Returns magnetic field component matrix
		 */
		void BField (long double x, long double y, long double z, long double t, long double B[4][4]){      //B-Feld am Ort des Teilchens berechnen
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					B[i][j] = 0;

			long double BFeldSkal = BFieldScale(t);
			if (BFeldSkal != 0){
				for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++){
					(*i)->BField(x, y, z, t, B);
				}

				if (BFeldSkal != 1)
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < 4; j++)
							B[i][j] *= BFeldSkal;

				B[3][0] = sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]); // absolute value of B-Vector
				if (B[3][0]>1e-31)
				{
					B[3][1] = (B[0][0]*B[0][1] + B[1][0]*B[1][1] + B[2][0]*B[2][1])/B[3][0]; // derivatives of absolute value
					B[3][2] = (B[0][0]*B[0][2] + B[1][0]*B[1][2] + B[2][0]*B[2][2])/B[3][0];
					B[3][3] = (B[0][0]*B[0][3] + B[1][0]*B[1][3] + B[2][0]*B[2][3])/B[3][0];
				}
			}
		};
		
		/**
		 * Calculate electric field and potential at a given position.
		 *
		 * Chooses the right map for this position and returns interpolated values.
		 *
		 * @param x Cartesian x coordinate
		 * @param y Cartesian y coordinate
		 * @param z Cartesian z coordinate
		 * @param t Time
		 * @param V Return electric potential (!=0 only if a map with potential was loaded)
		 * @param Ei Returns electric field vector
		 */
		void EField(long double x, long double y, long double z, long double t, long double &V, long double Ei[3]){
			Ei[0] = Ei[1] = Ei[2] = V = 0;
			for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++){
				(*i)->EField(x, y, z, t, V, Ei);
			}
		};
	
	private:

		/**
		 * Scale magnetic field due to field oscillation.
		 *
		 * @param t Time
		 */
		long double BFieldScale(long double t){		
			if (FieldOscillation==1)
				return 1 + OscillationFraction*sin(OscillationFrequency*2*pi*t);
			return 1;
		};
		
};

#endif // FIELDS_H_
