#ifndef GLOBALS_H_
#define GLOBALS_H_

#define KENNZAHL_UNKNOWN 0
#define KENNZAHL_NOT_FINISH 1
#define KENNZAHL_HIT_BOUNDARIES 2
#define KENNZAHL_NRERROR 3
#define KENNZAHL_DECAYED 4
#define KENNZAHL_INITIAL_NOT_FOUND 5

#define NEUTRON 1
#define PROTON 2
#define BF_ONLY 3
#define BF_CUT 4
#define ELECTRON 6

#define POLARISATION_GOOD 1
#define POLARISATION_BAD 2
#define POLARISATION_NONE 3

#define OUTPUT_EVERYTHING 1
#define OUTPUT_ENDPOINTS 2
#define OUTPUT_NOTHING 5
#define OUTPUT_EVERYTHINGandSPIN 3
#define OUTPUT_ENDPOINTSandSPIN 4


#include <string>

using namespace std;

int Log(const char* format, ...); // works like printf, prints to logfile and, if jobnumber == 0, to stdout
void percent(long double x, long double x1, long double len, int &perc);
void hunt(long double xx[], int n, long double x, int *jlo);

extern FILE *LOGSCR;

// physical constants
const long double pi = 3.14159265359;
const long double ele_e = 1.602176487E-19;      //elementary charge in SI
const long double gravconst = 9.80665;	// g
const long double conv = pi/180.;
const long double mu0 = 4*pi*1e-7;      //permeability,
const long double m_n = 1.674927211E-27/ele_e;	//neutron mass (eV/c^2)
const long double m_p = 1.672621637E-27/ele_e;	//proton mass (eV/c^2)
const long double m_e = 9.10938215e-31/ele_e;	//electron mass
const long double c_0 = 299792458;	// lightspeed
const long double hquer = 1.05457266e-34;	// planck constant (Js)
const long double mu_nSI = 0.96623641e-26;	// Neutron magn Mom (in J/T)
extern long double tau;	// neutron lifetime

//set the duration of the experiment
extern long double FillingTime;										//Filling in UCN
extern long double CleaningTime;                        // cleaning without field
extern long double RampUpTime;                          // ramping up coils
extern long double FullFieldTime;                       // storing in full field
extern long double RampDownTime;                        // ramping down coils
extern long double EmptyingTime;                        // emptying without field
extern long double StorageTime;                     // time when ramping down shall start, if xend > storage time, let neutron decay
extern int ExpPhase;  															// current experiment phase

extern int jobnumber;
extern string inpath, outpath; 
extern int ausgabewunsch; // user choice for Ausgabewunsch


#endif /*GLOBALS_H_*/
