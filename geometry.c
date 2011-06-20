#include <complex>

#include "globals.h"
#include "geometry.h"

solid defaultsolid = {"Restgas",{"Restgas",0,0,0},0}; // solid in which particle start

// transmission through a wall with Fermi potential  Mf (real) and PF (im) for a UCN with energy Er perpendicular to the wall (all in neV)
long double Transmission(const long double Er, const long double Mf, const long double Pf){
	complex<long double> V(Mf,Pf), sEV = sqrt(Er - V);
	long double sE = sqrt(Er);
	return 1 - norm((sE - sEV)/(sE + sEV));
}

// absorption probability of neutron flying distance l (in m) through Fermi potential Mf + i*Pf (all in neV)
long double Absorption(const long double E, const long double Mf, const long double Pf, const long double l){
	complex<long double> V(Mf,Pf);
	complex<long double> i(0,1);
	complex<long double> k = sqrt(2*m_n*1e-9*(E - V)); // complex wave vector
	complex<long double> wavefunc = exp(i*k*l/(hquer/ele_e));
	return 1 - (real(wavefunc)*real(wavefunc) + imag(wavefunc)*imag(wavefunc));
}
	
