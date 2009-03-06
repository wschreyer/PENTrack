/*****************************************
     Bwso erstellt nicht mehr folgende Matrix ar
               [1]       [2]         [3]
      [1]      Br        Bphi        Bz
      [2]      dBr/dr    dBphi/dr    dBz/dr
      [3]      dBr/dphi  dBphi/dphi  dBz/dphi
      [4]      dBr/dz    dBphi/dz    dBz/dz
***********************************************/


extern long double mapleBrp(long double r,long double z);
extern long double mapleBzp(long double r,long double z);
extern long double mapleBrdr(long double r,long double z);
extern long double mapleBrdz(long double r,long double z);
extern long double mapleBzdr(long double r,long double z);
extern long double mapleBzdz(long double r,long double z);
extern void Bwsomaple(long double r,long double phi,long double z); 


