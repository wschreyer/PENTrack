// class for field interpolation, create one for every table file you want to use
class TabField{
	private:
		int n,m;	// size of the arrays
		long double *rind, *zind, rdist, zdist;	// coordinate arrays, distance between coordinates
		long double conv_rA, conv_rB, conv_zA, conv_zB;	// conversion factors: array index <-> coordinate
		long double ****Brc, ****Bzc, ****Erc, ****Ezc; // 4D-arrays of interpolation coefficients
		void ReadTabFile(const char *tabfile, long double **&BrTab, long double **&BzTab, long double **&ErTab, long double **&EzTab);
		void CalcDeriv4th(long double **Tab, long double **Tab1, long double **Tab2,long double **Tab12);
		void CheckTab(long double **BrTab, long double **BzTab, long double **ErTab, long double **EzTab);
		void PreInterpol(long double ****&Brc, long double ****&Bzc, long double **BrTab, long double **BzTab);		
	public:
		TabField(const char *tabfile);
		~TabField();
		bool BInterpol(long double r, long double z, long double Bi[3], long double dBidrj[3][3]);
		bool EInterpol(long double r, long double z, long double Ei[3], long double dEidrj[3][3]);
};

