
// root stuff
//#include "TROOT.h"
//#include "TNtupleD.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TNetFile.h"
//#include "TRandom.h"
//#include "TBranch.h"
//#include "TClonesArray.h"
//#include "TStopwatch.h"
//#include "TMath.h"
#include "TF1.h"
//#include "TF2.h"


extern int ReadMagnets(void);

extern long double BForbes(long double r, long double phi, long double z, long double t);
extern long double trapzd(long double (*func)(long double), long double a, long double b, int n);
extern long double qromb(long double (*func)(long double), long double a, long double b);
extern long double BrINT(long double beta, long double sign1, long double sign2, long double r, long double phi, long double z, long double a, long double b, long double R_0);
extern long double OneCoil(long double r, long double phi, long double z, long double a, long double b, long double R_0, long double J_0, long double zoffset);
extern long double OneCoilRoot(long double r, long double phi, long double z, long double a, long double b, long double R_0, long double J_0, long double zoffset);
extern long double BrFoInt(long double beta);
extern long double dfridr(long double (*func)(long double), long double x, long double h, long double *err);
extern long double BzFoInt1(long double beta);
extern long double BzFoInt2(long double beta);
extern long double BzFoInt3(long double beta);
extern long double BzFoInt4(long double beta);
extern Double_t BrFoIntRoot(Double_t *beta, Double_t *par);
extern Double_t BrRoot(Double_t *x, Double_t *par);
extern Double_t BrRootR(Double_t *x, Double_t *par);
extern Double_t BrRootZ(Double_t *x, Double_t *par);
extern Double_t BzRoot(Double_t *x, Double_t *par);
extern Double_t BzRootR(Double_t *x, Double_t *par);
extern Double_t BzRootZ(Double_t *x, Double_t *par);
extern Double_t BzFoInt1Root(Double_t *beta, Double_t *par);
extern Double_t BzFoInt2Root(Double_t *beta, Double_t *par);
extern Double_t BzFoInt3Root(Double_t *beta, Double_t *par);
extern Double_t BzFoInt4Root(Double_t *beta, Double_t *par);
