#ifndef FIELDS_H_
#define FIELDS_H_


#include <vector>
#include <string>

#include "2dinterpolfeld.h"
#include "racetrack.h"

using namespace std;

struct TField{
	public:
		vector<TabField*> fields;
		vector<TRacetrack*> racetracks;
		// incorporate B-fieldoszillations into the code
		int FieldOscillation;        // turn field oscillation on if 1
		long double OscillationFraction, OscillationFrequency;    // Frequency in Hz
		
		TField(const char *infilename, int aFieldOscillation = 0, long double aOscillationFraction = 0, long double aOscillationFrequency = 0){
			FieldOscillation = aFieldOscillation;
			OscillationFraction = aOscillationFraction;
			OscillationFrequency = aOscillationFrequency;
			ifstream infile(infilename);
			string line;
			char c;
			while (infile.good()){ // read field table filenames
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if ((infile.peek() == '[') && getline(infile,line).good()){	// parse infile for section header
					if (line.compare(0,8,"[FIELDS]") == 0)
						LoadFieldsSection(infile);
					else if (line.compare(0,12,"[RACETRACKS]") == 0)
						LoadRacetrackSection(infile);
					else getline(infile,line);
				}
				else getline(infile,line);
			}
		};
		
		~TField(){
			for (vector<TabField*>::iterator i = fields.begin(); i != fields.end(); i++)
				delete (*i);	
			for (vector<TRacetrack*>::iterator i = racetracks.begin(); i != racetracks.end(); i++)
				delete (*i);	
		};
	
		// fill B field matrix with values
		/* field matrix:
			Bx		dBxdx	dBxdy	dBxdz
			By		dBydx	dBydy	dBydz
			Bz		dBzdx	dBzdy	dBzdz
			Babs	dBdx	dBdy	dBdz
		*/
		void BFeld (long double x, long double y, long double z, long double t, long double B[4][4]){      //B-Feld am Ort des Teilchens berechnen
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					B[i][j] = 0;

			if (fields.size() > 0){
				long double BFeldSkal = BFieldScale(t);
				for (vector<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){
					if ((*i)->BInterpol(t, x, y, z, B))
						break;
/*					else if (++i == fields.end())
						Log("\nThe particle has left fieldval boundaries: x=%LG, y=%LG z=%LG!", x, y, z);
*/				}
				if (BFeldSkal != 1)
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < 4; j++)
							B[i][j] *= BFeldSkal;
			}
			for (vector<TRacetrack*>::iterator i = racetracks.begin(); i != racetracks.end(); i++)
				(*i)->BFeld(x,y,z,B);

			B[3][0] = sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]); // absolute value of B-Vector
			if (B[3][0]>1e-31)
			{
				B[3][1] = (B[0][0]*B[0][1] + B[1][0]*B[1][1] + B[2][0]*B[2][1])/B[3][0]; // derivatives of absolute value
				B[3][2] = (B[0][0]*B[0][2] + B[1][0]*B[1][2] + B[2][0]*B[2][2])/B[3][0];
				B[3][3] = (B[0][0]*B[0][3] + B[1][0]*B[1][3] + B[2][0]*B[2][3])/B[3][0];
			}
			
		};
		
		// fill E field-vector with values
		void EFeld(long double x, long double y, long double z, long double &V, long double Ei[3]){
			Ei[0] = Ei[1] = Ei[2] = V = 0;
			for (vector<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){
				if ((*i)->EInterpol(x, y, z, V, Ei))
					break;
	/*			else if (++i == fields.end())
					Log("\nThe particle has left fieldval boundaries: x=%LG, y=%LG z=%LG!", x, y, z);
	*/		}
		};
	
	private:
		void LoadFieldsSection(ifstream &infile){
			char c;
			string line;
			cout << "Loading table fields..." << endl;
			do{	// parse table field list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue;	// skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				string ft;
				long double Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime;
				infile >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime;
				TabField *tf;
				tf = new TabField(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime);
				fields.push_back(tf);
			}while(infile.good() && getline(infile,line).good());
		};
		
		void LoadRacetrackSection(ifstream &infile){
			char c;
			string line;
			cout << "Loading racetracks..." << endl;
			do{	// parse STLfile list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue;	// skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				string type;
				long double Ibar, p1, p2, p3, p4, p5, p6;
				TRacetrack *rt = NULL;
				infile >> type;
				if (type == "InfiniteWireZ"){
					infile >> Ibar >> p1 >> p2;
					rt = new TInfiniteWireZ(p1, p2, Ibar);
				}
				else if (type == "InfiniteWireZCenter"){
					infile >> Ibar;
					rt = new TInfiniteWireZCenter(Ibar);
				}
				else if (type == "FiniteWire"){
					infile >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6;
					rt = new TFiniteWire(p1, p2, p3, p4, p5, p6, Ibar);
				}
				else if (type == "FiniteWireX"){
					infile >> Ibar >> p1 >> p2 >> p3;
					rt = new TFiniteWireX(p1, p2, p3, Ibar);
				}
				else if (type == "FiniteWireY"){
					infile >> Ibar >> p1 >> p2 >> p3;
					rt = new TFiniteWireY(p1, p2, p3, Ibar);
				}
				else if (type == "FiniteWireZ"){
					infile >> Ibar >> p1 >> p2 >> p3 >> p4;
					rt = new TFiniteWireZ(p1, p2, p3, p4, Ibar);
				}
				else if (type == "FiniteWireZCenter"){
					infile >> Ibar >> p1 >> p2;
					rt = new TFiniteWireZCenter(p1, p2, Ibar);
				}
				else if (type == "FullRacetrack"){
					infile >> Ibar >> p1 >> p2 >> p3;
					rt = new TFullRacetrack(p1, p2, p3, Ibar);
				}
				if (rt)
					racetracks.push_back(rt);
			}while(infile.good() && getline(infile,line).good());
		};
	
		//B-Feld nach Wunsch skalieren, BFeldSkal wird an alle B Komp und deren Ableitungen multipliziert
		long double BFieldScale(long double t){		
			if (FieldOscillation==1)
				return 1 + OscillationFraction*sin(OscillationFrequency*2*pi*t);
			return 1;
		};
		
};

#endif // FIELDS_H_
