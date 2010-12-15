#include "nr/nr3.h"

// functions to save and print neutron distribution
void prepndist();
void fillndist(long double x1, VecDoub_I &y1, long double x2, VecDoub_I &y2);
void outndist(const char *ndistfile);

extern int neutdist;
