#include <cmath>

#include "adiabacity.h"
#include "globals.h"

// Wahrscheinlichkeit f�r keinen Spinflip nach Rabi (im Matora Paper)
long double rabiplus(double Br,double Bz,double dBrdr,double dBrdz,double dBzdr,double dBzdz,double vr_n,double vz_n,double t)
{
        long double t3,t5,t6,t9,t11,t12,t13,t14,t16,t17,t22,t23,t27,t33,t34,rabiplus;
        long double gamm = 1.83247188e+8;

        t3 = dBrdr*vr_n+dBrdz*vz_n;
        t5 = Br+t3*t;
      t6 = t5*t5;
      t9 = dBzdr*vr_n+dBzdz*vz_n;
      t11 = Bz+t9*t;
      t12 = t11*t11;
      t13 = t6+t12;
      t14 = pow(t13,0.5l);
      t16 = gamm*gamm;
      t17 = pow(t13,0.1E1l);
      t22 = pow(t3*t11-t9*t5,2.0l);
      t23 = t13*t13;
      t27 = sqrt(t16*t17+t22/t23);
      t33 = cos(0.5l*t27*t);
      t34 = (gamm*t14+t27)/t27*t33;

      rabiplus = t34 * t34;

   return rabiplus;
      }

// Wahrscheinlichkeit f�r einen Spinflip nach Rabi (im Matora Paper)
long double rabimin(double Br,double Bz,double dBrdr,double dBrdz,double dBzdr,double dBzdz,double vr_n,double vz_n,double t)
{
        long double t3,t6,t8,t11,t13,t14,t15,t16,t19,t20,t22,t23,t27,t31,t33,rabimin;
        long double gamm = 1.83247188e+8;


        t3 = dBrdr*vr_n+dBrdz*vz_n;
        t6 = dBzdr*vr_n+dBzdz*vz_n;
      t8 = Bz+t6*t;
      t11 = Br+t3*t;
      t13 = t3*t8-t6*t11;
      t14 = t11*t11;
      t15 = t8*t8;
      t16 = t14+t15;
      t19 = gamm*gamm;
      t20 = pow(t16,0.1E1l);
      t22 = t13*t13;
      t23 = t16*t16;
      t27 = sqrt(t19*t20+t22/t23);
      t31 = sin(0.5l*t27*t);
      t33 = t13/t16/t27*t31;


      rabimin = t33 * t33;

	return rabimin;
}

// Wahrscheinlichkeit f�r einen Spinflip nach Vladimirsky (Soviet Physics JETP, Volume 12, Number 4, April 1961)
long double vladimirsky(long double r, long double Br,long double Bphi, long double Bz,long double dBrdr, long double dBrdphi, long double dBrdz, long double dBphidr, long double dBphidphi, long double dBphidz, long double dBzdr, long double dBzdphi, long double dBzdz, long double Bws, long double vr_n, long double phidot, long double vz_n){
	long double vabs, dBdt_par, dBdt_perp, dBdt_r, dBdt_phi, dBdt_z, dBdt_square, W;

    dBdt_r =   dBrdr  *vr_n + (dBrdphi-Bphi)*phidot + dBrdz  *vz_n;
    dBdt_phi = dBphidr*vr_n + (dBphidphi+Br)*phidot + dBphidz*vz_n;
	dBdt_z =   dBzdr  *vr_n + (dBzdphi)     *phidot + dBzdz  *vz_n;
	dBdt_square = dBdt_r*dBdt_r+dBdt_phi*dBdt_phi+dBdt_z*dBdt_z;

	vabs = sqrt(vr_n*vr_n + phidot*phidot*r*r + vz_n*vz_n);
	// component of dBdt parallel to B
	dBdt_par = (1.0/vabs) * (dBdt_r*vr_n + dBdt_phi*r*phidot + dBdt_z*vz_n);
	// component of dBdt perpendicular to B
	dBdt_perp = sqrtl(dBdt_square-dBdt_par*dBdt_par);

	// spin flip probability according to Vladimirsky
	W = exp(-pi*mu_nSI*Bws*Bws/(hquer*dBdt_perp));
	
	if (W>1){
		printf("Schei�e!!!\n");
	}

    return W;
}

long double thumbrule(long double Br,long double Bphi, long double Bz,long double dBrdr, long double dBrdphi, long double dBrdz, long double dBphidr, long double dBphidphi, long double dBphidz, long double dBzdr, long double dBzdphi, long double dBzdz, long double Bws, long double vr_n, long double phidot, long double vz_n){

	long double dBdt, dBdt_r, dBdt_phi, dBdt_z;

    dBdt_r =   dBrdr  *vr_n + (dBrdphi-Bphi)*phidot + dBrdz  *vz_n;
    dBdt_phi = dBphidr*vr_n + (dBphidphi+Br)*phidot + dBphidz*vz_n;
    dBdt_z =   dBzdr  *vr_n + (dBzdphi)     *phidot + dBzdz  *vz_n;

    dBdt = sqrt(dBdt_r*dBdt_r+dBdt_phi*dBdt_phi+dBdt_z*dBdt_z);

	// Adiabacity mit Daumenformel

    if(Bws!=0) return (dBdt*hquer)/(2*mu_nSI*Bws*Bws);
    else return 1e31;
}
