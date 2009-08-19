extern void LoadGeometry();
extern void RandomPointInSourceVolume(long double &r, long double &phi, long double &z);
extern bool InSourceVolume(long double r, long double phi, long double z);
extern short ReflectCheck(long double x1, long double *y1, long double &x2, long double *y2); 
extern void IncrementCodes(int kennz);
extern void OutputCodes(int iMC);
