#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <set>
#include <vector>

#include "particle.h"

using namespace std;

struct material{
	string name;
	long double FermiReal, FermiImag, DiffProb;
};

struct solid{
	string name;
	material mat;
	unsigned kennz;
	set<int> ignoretimes;
};


void LoadGeometry(const char *geometryin);
void RandomPointInSourceVolume(long double &r, long double &phi, long double &z, long double &alpha, long double &gamma);
bool ReflectCheck(TParticle *particle, long double x1, long double &h, long double *y1, long double *y2, int &itercount); 
void AbsorbCheck(TParticle *particle, long double *y1, long double *y2, long double H);


extern int reflektlog;
extern int reflekt, diffuse;
extern vector<int> kennz_counter[3];
extern vector<solid> solids;	// dynamic array to associate solids with material/kennzahl/name
extern long double Emin_n;

#endif /*GEOMETRY_H_*/
