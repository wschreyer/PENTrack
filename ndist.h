#ifndef NDIST_H_
#define NDIST_H_

extern int neutdist;

// prepare neutron density histogram matrix
void prepndist();

/// fill neutron density histogram
void fillndist(double x1, double y1[6], double x2, double y2[6]);

// print neutron density distribution in ndist.out
void outndist(const char *ndistfile);

#endif
