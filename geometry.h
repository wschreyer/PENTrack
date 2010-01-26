extern void LoadGeometry();
extern void RandomPointInSourceVolume(long double &r, long double &phi, long double &z, long double &alpha, long double &gamma);
extern bool InSourceVolume(long double r, long double phi, long double z);
extern short ReflectCheck(long double x1, long double *y1, long double &x2, long double *y2, int &itercount); 
extern void IncrementCodes(int kennz);
extern void OutputCodes(int iMC);
extern void Snapshooter(long double x2, long double *ystart, long double H);

extern int reflektlog;

void RotateVector(long double v[3], long double n[3]); // rotate coordinate system so, that new z-axis lies on NORMALIZED vector n
