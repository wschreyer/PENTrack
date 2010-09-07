#include <list>
#include <string>

using namespace std;

void PrepareFields(list<string> &fieldvaltab); // Read fieldval tab and preinterpolate or read coils.cond
void BFeld (long double rloc, long double philoc, long double zloc, long double t);
void EFeld(long double rloc, long double philoc, long double zloc);
void FreeFields();

void PrintBFieldCut(const char *outfile);	// print cut through BField to file
void PrintBField(const char *outfile);	// print Bfield to file and investigate ramp heating

extern int bfeldwahl, Feldcount;
extern long double BFeldSkalGlobal, BFeldSkal, EFeldSkal;

// incorporate B-fieldoszillations into the code
extern int FieldOscillation;        // turn field oscillation on if 1
extern long double OscillationFraction, OscillationFrequency;

extern long double dBrdr, dBrdz, dBzdr, dBzdz, Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
extern long double dBphidr, dBphidz, dBrdphi, dBzdphi, dBphidphi;
extern long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field

extern long double BCutPlanePoint[3], BCutPlaneNormalAlpha, BCutPlaneNormalGamma, BCutPlaneSampleDist; // plane for cut through BField
extern int BCutPlaneSampleCount;
