extern void ReflektDiffuse(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t);
extern void ReflektEckig(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t);
extern int Absorber(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t);
extern long double Transmission(long double Er, long double Mf, long double Pf);
extern void StatAbsorbtion(long double v, long double FermiReal, long double FermiImag);
extern int ReflexTop(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,
										long double ProbDiffuse, long double FermiReal, long double FermiImag);
extern int ReflexBottom(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,
										long double ProbDiffuse, long double FermiReal, long double FermiImag);
extern int ReflexOut(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,
										long double ProbDiffuse, long double FermiReal, long double FermiImag);
extern int ReflexIn(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,
										long double ProbDiffuse, long double FermiReal, long double FermiImag);
extern int ReflexGeneral(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,
										long double ProbDiffuse, long double FermiReal, long double FermiImag,long double rNorm,long double zNorm,long double phiNorm);