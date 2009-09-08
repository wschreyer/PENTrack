extern void PrepareFields(); // Read fieldval tab and preinterpolate or read coils.cond
extern void BFeld (long double rloc, long double philoc, long double zloc, long double t);
extern void EFeld(long double rloc, long double philoc, long double zloc);
extern void FreeFields();

void SwitchField(long double t);
long double Bnull();
long double Banalytic(long double r,long double phi,long double z, long double t);

