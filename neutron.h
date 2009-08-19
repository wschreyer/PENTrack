extern void BFeld (long double rloc, long double philoc, long double zloc, long double t);
extern void EFeld(long double rloc, long double philoc, long double zloc);
extern void CylKartCoord(long double Wr, long double Wphi, long double Wz, long double phi, long double *Wx0, long double *Wy0, long double *Wz0);
extern void KartCylCoord(long double Wx, long double Wy, long double Wz, long double phi, long double *Wr0, long double *Wphi0, long double *Wz0);
extern void SwitchField(long double t);
extern void PrintConfig(void);
extern long double AbsValueCart(long double x, long double y, long double z);
