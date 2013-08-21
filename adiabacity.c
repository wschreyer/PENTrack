#include <cstdio>
#include <cmath>

#include "adiabacity.h"
#include "globals.h"

// Probability for spin flip according to Rabi (in Matora Paper)
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

// Probability for spin flip according to Rabi (im Matora Paper)
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

// Probability for spin flip according to Vladimirsky (Soviet Physics JETP, Volume 12, Number 4, April 1961)
long double vladimirsky(long double Bx,long double By, long double Bz,long double dBxdx, long double dBxdy, long double dBxdz, long double dBydx, long double dBydy, long double dBydz, long double dBzdx, long double dBzdy, long double dBzdz, long double Bws, long double vx, long double vy, long double vz){
	long double vabs, dBdt_par, dBdt_perp, dBdt_x, dBdt_y, dBdt_z, dBdt_square, W;

    dBdt_x = dBxdx*vx + dBxdy*vy + dBxdz*vz;
    dBdt_y = dBydx*vx + dBydy*vy + dBydz*vz;
	dBdt_z = dBzdx*vx + dBzdy*vy + dBzdz*vz;
	dBdt_square = dBdt_x*dBdt_x+dBdt_y*dBdt_y+dBdt_z*dBdt_z;

	vabs = sqrt(vx*vx + vy*vy + vz*vz);
	// component of dBdt parallel to B
	dBdt_par = (1.0/vabs) * (dBdt_x*vx + dBdt_y*vy + dBdt_z*vz);
	// component of dBdt perpendicular to B
	dBdt_perp = sqrtl(dBdt_square-dBdt_par*dBdt_par);

	// spin flip probability according to Vladimirsky
	W = exp(pi*mu_nSI*Bws*Bws/(hbar*dBdt_perp));
	
	if (W>1){
		printf("Scheisse!!!\n");
	}

    return W;
}

long double thumbrule(long double Bx, long double By, long double Bz,long double dBxdx, long double dBxdy, long double dBxdz, long double dBydx, long double dBydy, long double dBydz, long double dBzdx, long double dBzdy, long double dBzdz, long double Bws, long double vx, long double vy, long double vz){

	long double dBdt, dBdt_x, dBdt_y, dBdt_z;

    dBdt_x = dBxdx*vx + dBxdy*vy + dBxdz*vz;
    dBdt_y = dBydx*vx + dBydy*vy + dBydz*vz;
    dBdt_z = dBzdx*vx + dBzdy*vy + dBzdz*vz;

    dBdt = sqrt(dBdt_x*dBdt_x+dBdt_y*dBdt_y+dBdt_z*dBdt_z);

	// Adiabacity mit Daumenformel

    if(Bws!=0) return (dBdt*hbar)/(-2*mu_nSI*Bws*Bws);
    else return 1e31;
}
