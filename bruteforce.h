long double BruteForceIntegration(long double x, long double *y, long double Br, long double Bphi, long double Bz, long double Bws); // integrate spin flip probability

// variables for BruteForce Bloch integration BEGIN
extern long double BFTargetB;     // smallest value of Babs during step, Babs < BFTargetB => integrate,

