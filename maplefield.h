/*****************************************
     Bwso erstellt nicht mehr folgende Matrix ar
               [1]       [2]         [3]
      [1]      Br        Bphi        Bz
      [2]      dBr/dr    dBphi/dr    dBz/dr
      [3]      dBr/dphi  dBphi/dphi  dBz/dphi
      [4]      dBr/dz    dBphi/dz    dBz/dz
***********************************************/


long double mapleBrp(long double r,long double z);
long double mapleBzp(long double r,long double z);
long double mapleBrdr(long double r,long double z);
long double mapleBrdz(long double r,long double z);
long double mapleBzdr(long double r,long double z);
long double mapleBzdz(long double r,long double z);
extern void Bwsomaple(long double r,long double phi,long double z, long double Ibar, long double Bi[3], long double dBidrj[3][3]); 


