#include <cmath>
#include "globals.h"
#include "maplefield.h"

/*****************************************
     Bwso erstellt nicht mehr folgende Matrix ar
               [1]       [2]         [3]
      [1]      Br        Bphi        Bz
      [2]      dBr/dr    dBphi/dr    dBz/dr
      [3]      dBr/dphi  dBphi/dphi  dBz/dphi
      [4]      dBr/dz    dBphi/dz    dBz/dz
***********************************************/

long double mapleBrp(long double r,long double z)
{
  {
    return(-0.39213E1*expl(-0.56331E2*r+0.67597E1)*cosl(0.5236E2*z-0.78052)+
0.39213E1*expl(-0.3688E2+0.73761E2*r)*cosl(0.5236E2*z-0.80105)+0.39213E1*cosl(
0.5236E2*r+0.11165E1)*expl(-0.4618E2*z));
  }
}



long double mapleBzp(long double r,long double z)
{
  {   
    return((0.74891138E-1*expl(-0.56331E2*r+0.67597E1)*sinl(0.5236E2*z-0.78052)
-0.74891138E-1*expl(-0.3688E2+0.73761E2*r)*sinl(0.5236E2*z-0.80105)+0.84913382E-1
*cosl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)+r*(-0.42186927E1*expl(-0.56331E2*r+
0.67597E1)*sinl(0.5236E2*z-0.78052)-0.55240452E1*expl(-0.3688E2+0.73761E2*r)*sinl(
0.5236E2*z-0.80105)-0.44460647E1*sinl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)))/r);
  }
}


long double mapleBrdr(long double r,long double z)
{
  {
    return(0.22089075E3*expl(-0.56331E2*r+0.67597E1)*cosl(0.5236E2*z-0.78052)+
0.28923901E3*expl(-0.3688E2+0.73761E2*r)*cosl(0.5236E2*z-0.80105)-0.20531927E3*
sinl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z));
  }
}



long double mapleBrdz(long double r,long double z)
{
  {
   return(0.20531927E3*expl(-0.56331E2*r+0.67597E1)*sinl(0.5236E2*z-0.78052)
-0.20531927E3*expl(-0.3688E2+0.73761E2*r)*sinl(0.5236E2*z-0.80105)-0.18108563E3*
cosl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z));
  }
}



long double mapleBzdr(long double r,long double z)
{
  {
    return((-0.84373854E1*expl(-0.56331E2*r+0.67597E1)*sinl(0.5236E2*z-0.78052)
-0.1104809E2*expl(-0.3688E2+0.73761E2*r)*sinl(0.5236E2*z-0.80105)-0.88921294E1*
sinl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)+r*(0.23764318E3*expl(-0.56331E2*r+
0.67597E1)*sinl(0.5236E2*z-0.78052)-0.4074591E3*expl(-0.3688E2+0.73761E2*r)*sinl(
0.5236E2*z-0.80105)-0.23279595E3*cosl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)))/r
-(0.74891138E-1*expl(-0.56331E2*r+0.67597E1)*sinl(0.5236E2*z-0.78052)
-0.74891138E-1*expl(-0.3688E2+0.73761E2*r)*sinl(0.5236E2*z-0.80105)+0.84913382E-1
*cosl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)+r*(-0.42186927E1*expl(-0.56331E2*r+
0.67597E1)*sinl(0.5236E2*z-0.78052)-0.55240452E1*expl(-0.3688E2+0.73761E2*r)*sinl(
0.5236E2*z-0.80105)-0.44460647E1*sinl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)))/(
r*r));
  }
}



long double mapleBzdz(long double r,long double z)
{
  {
     return((0.39213E1*expl(-0.56331E2*r+0.67597E1)*cosl(0.5236E2*z-0.78052)
-0.39213E1*expl(-0.3688E2+0.73761E2*r)*cosl(0.5236E2*z-0.80105)-0.39213E1*cosl(
0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)+r*(-0.22089075E3*expl(-0.56331E2*r+
0.67597E1)*cosl(0.5236E2*z-0.78052)-0.28923901E3*expl(-0.3688E2+0.73761E2*r)*cosl(
0.5236E2*z-0.80105)+0.20531927E3*sinl(0.5236E2*r+0.11165E1)*expl(-0.4618E2*z)))/r);
  }
}






void Bwsomaple(long double r,long double phi,long double z, long double Ibar, long double Bi[3], long double dBidrj[3][3])
{
	Bi[0] += mapleBrp(r,z);
	Bi[2] += mapleBzp(r,z);
	dBidrj[0][0] += mapleBrdr(r,z);
	dBidrj[0][2] += mapleBrdz(r,z);
	dBidrj[2][0] += mapleBzdr(r,z);
	dBidrj[2][2] += mapleBzdz(r,z);

  	Bi[1] += Ibar*mu0/(2*pi*r);
  	dBidrj[1][0] += -Bi[1]/r;

  	return;
}

