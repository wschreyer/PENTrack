#include "particle.h"

void LoadMCGenerator(const char *infile);
void MCStartwerte(TParticle *particle);
void MCZerfallsstartwerte(TParticle *neutron, TParticle *proton, TParticle *electron);

long double UniformDist(long double min, long double max);
long double CosDist(long double min, long double max);
long double SquareDist(long double min, long double max);
void IsotropicDist(long double &alpha, long double &gamma);
long double ProtonSpectrum();
long double ElectronSpectrum();

extern unsigned long int monthinmilliseconds; // RandomSeed
extern int decay, polarisation;	// decay on or off?, polarisation
extern long double decayoffset;
