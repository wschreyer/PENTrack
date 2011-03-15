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
		int bfeldwahl;
		long double BFeldSkalGlobal, EFeldSkal;
		// incorporate B-fieldoszillations into the code
		int FieldOscillation;        // turn field oscillation on if 1
		long double OscillationFraction, OscillationFrequency;    // Frequency in Hz
		
		TField(const char *infilename, int abfeldwahl, long double aBFeldSkalGlobal, long double aEFeldSkal,
				int aFieldOscillation = 0, long double aOscillationFraction = 0, long double aOscillationFrequency = 0){
			bfeldwahl = abfeldwahl;
			BFeldSkalGlobal = aBFeldSkalGlobal;
			EFeldSkal = aEFeldSkal;	
			FieldOscillation = aFieldOscillation;
			OscillationFraction = aOscillationFraction;
			OscillationFrequency = aOscillationFrequency;
			if ((bfeldwahl == 0 || bfeldwahl == 2) && (BFeldSkalGlobal != 0 || EFeldSkal != 0))
			{
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
			}
			else if(bfeldwahl==4)
			{
	//			ReadMagnets(("in/coils.cond").c_str());
			
				printf("\n \n Test of integration\n");
				//	long double TestInt;
				long double B[4][4];
				for (int a = 0;a<1;a++)
				{
				
					BFeld(0.3,0,0.1, 500.0, B);
					printf("T\n");
				}
				printf("Br = %.17LG \n",B[0][0]);
				printf("dBrdr = %.17LG \n",B[0][1]);
				printf("dBrdz = %.17LG \n",B[0][3]);
				printf("Bz = %.17LG \n",B[2][0]);
				printf("dBzdr = %.17LG \n",B[2][1]);
				printf("dBzdz = %.17LG \n",B[2][3]);	
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
				
			if (bfeldwahl != 1){
				long double BFeldSkal = BFieldScale(t);
				switch (bfeldwahl)
				{				
					case 0:	if (BFeldSkal != 0){
								for (vector<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){
									if ((*i)->BInterpol(x, y, z, B))
										break;
		/*							else if (++i == fields.end())
										Log("\nThe particle has left fieldval boundaries: x=%LG, y=%LG z=%LG!", x, y, z);
		*/						}
								if (BFeldSkal != 1)
									for (int i = 0; i < 3; i++)
										for (int j = 0; j < 4; j++)
											B[i][j] *= BFeldSkal;												
							}
							RacetrackField(x,y,z,B);
							break;
							
					case 3: Banalytic(x, y, z, t, B);
							break;
							
					case 4: //BForbes(x, y, z, t, Bi, dBidrj);
							RacetrackField(x,y,z, B);
							break;
					default: printf("\nbfeldwahl = %d unknown!\n",bfeldwahl);
							exit(-1);
				}   
				B[3][0] = sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]); // absolute value of B-Vector
				if (B[3][0]>1e-31)
				{			
					B[3][1] = (B[0][0]*B[0][1] + B[1][0]*B[1][1] + B[2][0]*B[2][1])/B[3][0]; // derivatives of absolute value
					B[3][2] = (B[0][0]*B[0][2] + B[1][0]*B[1][2] + B[2][0]*B[2][2])/B[3][0];
					B[3][3] = (B[0][0]*B[0][3] + B[1][0]*B[1][3] + B[2][0]*B[2][3])/B[3][0];
				}
			}
			
		};
		
		// fill E field-vector with values
		void EFeld(long double x, long double y, long double z, long double &V, long double Ei[3]){
			Ei[0] = Ei[1] = Ei[2] = V = 0;
			if (EFeldSkal != 0 && (bfeldwahl == 0 || bfeldwahl == 2)){
				for (vector<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){
					if ((*i)->EInterpol(x, y, z, V, Ei))
						break;
		/*			else if (++i == fields.end())
						Log("\nThe particle has left fieldval boundaries: x=%LG, y=%LG z=%LG!", x, y, z);
		*/		}
				V *= EFeldSkal;
				Ei[0] *= EFeldSkal;
				Ei[1] *= EFeldSkal;
				Ei[2] *= EFeldSkal;												
		    }
			return;
		};
	
	private:
		void LoadFieldsSection(ifstream &infile){
			char c;
			string line;
			do{	// parse table field list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue;	// skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				string ft;
				infile >> ft;
				TabField *tf;
				tf = new TabField(ft.c_str());
				fields.push_back(tf);
			}while(infile.good() && getline(infile,line).good());
		};
		
		void LoadRacetrackSection(ifstream &infile){
			char c;
			string line;
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
			int phase = ExpPhase(t);
			long double result = 0;
			switch (phase){
				case FILLING_PHASE:
				case CLEANING_PHASE:
				case COUNTING_PHASE: 
					return 0;
				case RAMPUP_PHASE:
					// ramping up field smoothly with cosine
//					result = (0.5 - 0.5*cos(pi*(t - CleaningTime - FillingTime)/RampUpTime)) * BFeldSkalGlobal;

					// linear ramp
					result = (t - CleaningTime - FillingTime)/RampUpTime * BFeldSkalGlobal;
					break;
				case FULLFIELD_PHASE:
					result = BFeldSkalGlobal;
					break;
				case RAMPDOWN_PHASE:
					// ramping down field smoothly with cosine
					//result = (0.5 + 0.5*cos(pi*(t - (RampUpTime + CleaningTime + FillingTime + FullFieldTime)) / RampDownTime)) * BFeldSkalGlobal;

					// linear ramp
					result = (1 - (t - RampUpTime - CleaningTime - FillingTime - FullFieldTime)/RampDownTime) * BFeldSkalGlobal;
					break;
				default: result = BFeldSkalGlobal;
			}
			if (FieldOscillation==1)
				result *= 1 + OscillationFraction*sin(OscillationFrequency*2*pi*t);		
			return result;
		};
		
		// produce linearly rising B-field in z-direction and linearly changing in time
		void Banalytic(long double x,long double y,long double z, long double t, long double B[4][4]){
			long double dBx = 5;
			long double x0 = 1.1;
			long double alpha = 0.5;
			long double beta = 0.025;
			
			
			// transform Bfield to higher z-values
			long double offset = -0.0;
			z = z+offset;
			
			B[0][0] += (-dBx*(x0*z-0.5*z*z)+6)*(alpha+beta*t);
		    B[0][3] += -dBx*(alpha+beta*t)*(x0-z);
		};	
		
		void RacetrackField(long double x, long double y, long double z, long double B[4][4]){
			for (vector<TRacetrack*>::iterator i = racetracks.begin(); i != racetracks.end(); i++)
				(*i)->BFeld(x,y,z,B);
		};

};

#endif // FIELDS_H_
