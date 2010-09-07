#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <cmath>

#include "globals.h"
#include "fields.h"

class TParticle{
	public:
		int kennz, particlenumber, protneut;
		long double xstart, xend, dt; // start time, max simulation time, actual simulation time
		long double Hstart, Hend, Hmax, NeutEnergie;
		long double rstart, rend, phistart, phiend, zstart, zend, alphastart, alphaend, gammastart, gammaend, vstart, vend;
		long double h1, hmax, trajlength, BahnPointSaveTime;
		long double vlad, vladmax, vladtotal, frac, thumbmax, BFsurvprob, decayerror;
		int polarisation;
		int nrefl, NSF, nok, nbad;

		TParticle(int aprotneut, int number){
			kennz = KENNZAHL_UNKNOWN; particlenumber = number;
			Hmax = 0; BFsurvprob = 1; trajlength = decayerror = vladtotal = vladmax = thumbmax = 0;
			nrefl = NSF = nok = nbad = 0;
			protneut = aprotneut;
			switch (protneut){
				case NEUTRON: h1 = 5e-5; BahnPointSaveTime = 5e-7; break;
				case PROTON: h1 = 1e-8; BahnPointSaveTime = 1e-8; break;
				case ELECTRON: h1 = 2e-12; BahnPointSaveTime = 1e-11; break;
				default: Log("Particle type %i unknown!",protneut); exit(-1);
			}
		};

		long double Energy(long double x, long double *y) const{
			long double v;
			switch (protneut){
				case NEUTRON:
					BFeld(y[1], y[5], y[3], x); 
					v = sqrt(y[2]*y[2] + y[4]*y[4] + y[1]*y[6]*y[1]*y[6]);
					return m_n*gravconst*y[3] + 0.5*m_n*v*v - polarisation*mu_nSI/ele_e*Bws;
				case PROTON:
					v = sqrt(y[2]*y[2] + y[4]*y[4] + y[1]*y[6]*y[1]*y[6]);
					return 0.5*m_p*v*v;
				case ELECTRON:
					v = sqrt(y[2]*y[2] + y[4]*y[4] + y[1]*y[6]*y[1]*y[6]);
					return c_0*c_0  * m_e *  (1/sqrt(1 - v*v/(c_0*c_0)) - 1);
				default: return 0;
			}
		};

		void Fillystart(long double *y){
			y[1] = rstart;
			y[3] = zstart;
			y[5] = phistart;
			long double Ekin;
			switch (protneut){
				case NEUTRON:
					BFeld(rstart, phistart, zstart, xstart);
					Ekin = Hstart - m_n*gravconst*zstart + polarisation*mu_nSI/ele_e*Bws;
					if (Ekin > 0)
						vstart = sqrt(2/m_n*Ekin);
					else
						vstart = 0;
					break;
				case PROTON:
					vstart = sqrt(2/m_p*Hstart);
					break;
				case ELECTRON:
					long double gammarel = Hstart/m_e/c_0/c_0 + 1;
					vstart = c_0 * sqrtl(1-(1/(gammarel*gammarel)));     // relatistic velocity from E_kin in eV
					break;
			}
			y[2] = vstart*sin(gammastart)*cos(alphastart - phistart);
			y[4] = vstart*cos(gammastart);
			if (rstart != 0) y[6] = vstart*sin(gammastart)*sin(alphastart - phistart)/rstart;
			else y[6] = 0;
		};

		void derivs(long double x, long double *y, long double *dydx) const{
			BFeld(y[1],y[5],y[3], x);
			dydx[1]= y[2];
			dydx[3]= y[4];
			dydx[5]= y[6];
			switch (protneut){
				case NEUTRON:
					dydx[2]= y[1]*(y[6]*y[6]) + polarisation*mu_nSI/ele_e/m_n*dBdr;
					dydx[4]= polarisation*mu_nSI/ele_e/m_n*dBdz - gravconst; // [m/s^2]
					dydx[6]= -2*y[2]*y[6]/y[1] + polarisation*mu_nSI/ele_e/m_n/y[1]*dBdphi/y[1];
					break;
				case PROTON:
					EFeld(y[1],y[5],y[3]);
					dydx[2]= y[1]*(y[6]*y[6]) + (Bz*y[1]*y[6] - Bphi*y[4] + Er)/m_p;
					dydx[4]= (Ez + Bphi*y[2] - y[1]*Br*y[6])/m_p - gravconst; 
					dydx[6]= -2*y[2]*y[6]/y[1] + (Ephi + Br*y[4] - Bz*y[2])/y[1]/m_p;
					break;
				case ELECTRON:
					EFeld(y[1],y[5],y[3]);
					long double Qm0 = -1.0/m_e*sqrtl(1 - (y[2]*y[2] + y[1]*y[1]*y[6]*y[6] + y[4]*y[4])/(c_0*c_0));
					dydx[2]= y[1]*(y[6]*y[6]) + Qm0*(Bz*y[1]*y[6] - Bphi*y[4] + Er);
					dydx[4]= Qm0*(Ez + Bphi*y[2]-y[1]*Br*y[6]) - gravconst;
					dydx[6]= -2*y[2]*y[6]/y[1] + Qm0*(Ephi + Br*y[4] - Bz*y[2])/y[1];
					break;
			}
		};
};


#endif /*PARTICLE_H_*/
