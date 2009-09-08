class TabField{
	private:
		int n,m;
		long double *rind, *zind, rdist, zdist;
		long double conv_rA, conv_rB, conv_zA, conv_zB;
		long double ****Brc, ****Bzc, ****Erc, ****Ezc; 
		void ReadTabFile(const char *tabfile, long double **&BrTab, long double **&BzTab, long double **&ErTab, long double **&EzTab);
		void CalcDeriv4th(long double **BrTab, long double **BzTab, long double **BrTab1, long double **BrTab2,long double **BrTab12,
																	long double **BzTab1, long double **BzTab2, long double **BzTab12);
		void CheckTab(long double **BrTab, long double **BzTab, long double **ErTab, long double **EzTab);
		void PreInterpol(long double ****&Brc, long double ****&Bzc, long double **BrTab, long double **BzTab);		
	public:
		TabField(const char *tabfile);
		~TabField();
		bool BInterpol(long double r, long double z);
		bool EInterpol(long double r, long double z);
};


void bcucof(long double y[], long double y1[], long double y2[], long double y12[], long double d1, long double d2,long double **c);
void bcuint_new(int indr, int indz,long double ****c, long double d1, long double d2, long double x1l, long double x1u, long double x2l, long double x2u, long double x1, long double x2, long double *ansy, long double *ansy1, long double *ansy2);

