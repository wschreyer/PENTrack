extern void PrepareBField(); // Read fieldval tab and preinterpolate or read coils.cond
extern void PrepIntpol(int k);
extern void BInterpol(long double r_n, long double phi, long double z_n);
extern void EInterpol(long double r_n, long double phi, long double z_n);
extern long double Bnull();
extern long double Banalytic(long double r,long double phi,long double z, long double t);
extern void Preinterpol(int k);

void GetDim(int *m, int *n);
void CalcDeriv4th(long double rind[], long double zind[], long double **BrTab, long double **BzTab, long double **BphiTab, long double **BrTab1, long double **BrTab2, long double **BrTab12, long double **BzTab1, long double **BzTab2, long double **BzTab12, long double **BphiTab1, long double **BphiTab2,long double **BphiTab12,int m, int n);
int readWert(long double rind[], long double zind[], long double **BrTab,long double **BphiTab, long double **BzTab,long double **ErTab, long double **EphiTab,long double **EzTab, int *m, int *n);
void hunt(long double xx[], int n, long double x, int *jlo);
void bcucof(long double y[], long double y1[], long double y2[], long double y12[], long double d1, long double d2,long double **c);
void bcuint_new(int indr, int indz,long double ****c, long double d1, long double d2, long double x1l, long double x1u, long double x2l, long double x2u, long double x1, long double x2, long double *ansy, long double *ansy1, long double *ansy2);
void Racetrack(long double r,long double phi,long double z,  long double lx, long double ly, long double I_rt);
void CenterCurrent(long double r,long double phi,long double z,  long double I_rt);
void StraightWireField(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z);
void BarRaceTrack(long double r_current, long double phi_current, long double z_current, long double I_bar);

extern long double Emin_n;
