#ifndef DINTERPOLFELD_H_
#define DINTERPOLFELD_H_

#include <vector>

#include "nr/nr3.h"

using namespace std;

// class for field interpolation, create one for every table file you want to use
class TabField{
	private:
		int m,n;	// size of the arrays
		long double rdist, zdist, r_mi, z_mi;	// distance between coordinates
		NRmatrix<Doub[4][4]> Brc, Bphic, Bzc, Erc, Ephic, Ezc, Vc;
		void ReadTabFile(const char *tabfile, MatDoub_O &BrTab, MatDoub_O &BphiTab, MatDoub_O &BzTab, MatDoub_O &ErTab, MatDoub_O &EphiTab, MatDoub_O &EzTab, MatDoub_O &VTab);
		void CalcDerivs(MatDoub_I &Tab, MatDoub_O &Tab1, MatDoub_O &Tab2, MatDoub_O &Tab12);
		void CheckTab(MatDoub_I &BrTab, MatDoub_I &BphiTab, MatDoub_I &BzTab, MatDoub_I &ErTab, MatDoub_I &EphiTab, MatDoub_I &EzTab, MatDoub_I &VTab);
		void PreInterpol(NRmatrix<Doub[4][4]> &coeff, MatDoub_I &Tab);		
	public:
		TabField(const char *tabfile);
		bool BInterpol(long double x, long double y, long double z, long double B[4][4]);
		bool EInterpol(long double x, long double y, long double z, long double &V, long double Ei[3]);
};


#endif // DINTERPOLFELD_H_
