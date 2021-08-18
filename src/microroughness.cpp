#include "microroughness.h"

#include <complex>

#include "optimization.h"
#include "specialfunctions.h"
#include <boost/numeric/odeint.hpp>

#include "globals.h"

using namespace std;

namespace MR{

bool MRValid(const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength){
	if (RMSroughness == 0 || correlationLength == 0)
		return false;
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; // velocity squared
	double E = 0.5*m_n*v2; // kinetic energy
	double ki = sqrt(2*m_n*E)*ele_e/hbar; // wave number in first solid
	if (Estep >= 0){
		double kc = sqrt(2*m_n*Estep)*ele_e/hbar; // critical wave number of potential wall
		if (2*RMSroughness*ki < 2 && 2*RMSroughness*kc < 2)
			return true;
	}
	else{
		double kt = sqrt(2*m_n*(E - Estep))*ele_e/hbar; // wave number in second solid
		if (2*RMSroughness*ki < 2 && 2*RMSroughness*kt < 2)
			return true;
	}
	std::cout << "MR model not applicable. Falling back to specular reflection!\n";
	return false;
}

double MRDist(const bool transmit, const bool integral, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength, const double theta, const double phi){
	if (theta < 0 || theta > pi/2)
		return 0;
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; // velocity squared
	double vnormal = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]; // velocity projected onto surface normal
	double E = 0.5*m_n*v2; // kinetic energy
	if (transmit && E <= Estep) // if particle is transmitted: check if energy is higher than potential wall
		return 0;

	double ki = sqrt(2*m_n*E)*ele_e/hbar; // wave number in first solid
	std::complex<double> kc = std::sqrt(2*(double)m_n*std::complex<double>(Estep, 0.))*(double)ele_e/(double)hbar; // critical wave number of potential wall
	double kt = 0;

	double costheta_i = abs(vnormal/sqrt(v2)); // cosine of angle between normal and incoming velocity vector
	std::complex<double> Si = 2*costheta_i/(costheta_i + sqrt(costheta_i*costheta_i - kc*kc/ki/ki)); // specularly transmitted amplitude
	std::complex<double> So;
	if (transmit){
		kt = sqrt(2*m_n*(E - Estep))*ele_e/hbar; // wave number in second solid
		if (-real(kc*kc)/kt/kt > cos(theta)*cos(theta)) // this can happen if Estep < 0
			return 0;
		So = 2*cos(theta)/(cos(theta) + sqrt(cos(theta)*cos(theta) + kc*kc/kt/kt)); // diffusely transmitted amplitude
	}
	else
		So = 2*cos(theta)/(cos(theta) + sqrt(cos(theta)*cos(theta) - kc*kc/ki/ki)); // diffusely reflected amplitude

	double b = RMSroughness, w = correlationLength;

	double theta_i = acos(costheta_i);
	double Fmu = 0; // fourier transform of roughness correlation function
	if (!integral && !transmit)
		Fmu = exp(-w*w/2*ki*ki*(sin(theta_i)*sin(theta_i) +       sin(theta)*sin(theta) -                                 2*sin(theta_i)*sin(theta)*cos(phi)));
	else if (!integral && transmit)
		Fmu = exp(-w*w/2*(ki*ki*sin(theta_i)*sin(theta_i) + kt*kt*sin(theta)*sin(theta) -                           2*ki*kt*sin(theta_i)*sin(theta)*cos(phi)));
	else if (integral && !transmit) // if integral is set, precalculate phi-integral using modified Bessel function of first kind
		Fmu = exp(-w*w/2*ki*ki*(sin(theta_i)*sin(theta_i) +       sin(theta)*sin(theta))) * 2*pi*alglib::besseli0(w*w*ki*ki*sin(theta_i)*sin(theta));
	else if (integral && transmit)
		Fmu = exp(-w*w/2*(ki*ki*sin(theta_i)*sin(theta_i) + kt*kt*sin(theta)*sin(theta))) * 2*pi*alglib::besseli0(w*w*ki*kt*sin(theta_i)*sin(theta));
	double factor = real(kc*kc*kc*kc)*b*b*w*w/8/pi/costheta_i;
	if (transmit)
		factor *= kt/ki;
	if (integral)
		factor *= sin(theta); // integral over sin(theta) dtheta dphi
	double result = factor*norm(Si)*norm(So)*Fmu;
	//	cout << integral << transmit << ' ' << theta_i << ' ' << theta << ' ' << phi << ' ' << ki << ' ' << kt << ' ' << kc << ' ' << factor << ' ' << Si << ' ' << So << ' ' << Fmu << ' ' << result << '\n';
	return result; // return probability
}

/**
 * struct containing all collision parameters, used to pass them to the wrapper functions below
 */
struct TMRParams{
	const bool transmit; ///< True, if the scattering distribution for transmission through the surface should be returned
	const bool integral; ///< True, if the phi-integrated scattering distribution should be returned
	const double *v; ///< Icident velocity vector
	const double *normal; ///< Normal vector
	const double Estep; ///< Potential step
	const double RMSroughness; ///< RMS roughness of surface
	const double correlationLength; ///< Correlation length of surface
};

double MRProb(const bool transmit, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength){
	vector<double> total(1, 0);
	auto integrand = [transmit, v, normal, Estep, RMSroughness, correlationLength](const vector<double> &dummy, std::vector<double> &result, const double theta){ // use lambda expression to define local function that has the proper parameters for odeint
		result[0] = MRDist(transmit, true, v, normal, Estep, RMSroughness, correlationLength, theta, 0);
	};
	boost::numeric::odeint::integrate(integrand, total, 0.0, (double)pi/2, 0.01); // integrate "differential equation" dP/dtheta = MRdist(theta) from 0 to pi/2
	return total[0];
}

/**
 * Wrapper function of MRDist, is passed to numerical-optimization routine from alglib library
 * 
 * @param x Array of function parameters, in this case only one: polar angle
 * @param y Returns MRDist value at theta=x[0] and phi=0
 * @param ptr Pointer data, used to pass TMRParams struct
 */
void NegMRDist(const alglib::real_1d_array &x, double &y, void *ptr){
	TMRParams *p = static_cast<TMRParams*>(ptr);
	y = -MRDist(p->transmit, p->integral, p->v, p->normal, p->Estep, p->RMSroughness, p->correlationLength, x[0], 0);
}


double MRDistMax(const bool transmit, const double v[3], const double normal[3], const double Estep, const double RMSroughness, const double correlationLength){
	alglib::mincgstate s;
	alglib::real_1d_array theta = "[0.7853981635]";
	alglib::mincgcreatef(1, theta, 1e-10, s);
	alglib::mincgsetcond(s, 0, 0, 0, 0);
	TMRParams p = {transmit, false, v, normal, Estep, RMSroughness, correlationLength};
	alglib::mincgoptimize(s, NegMRDist, NULL, &p);
	alglib::mincgreport r;
	alglib::mincgresults(s, theta, r);
//	cout << theta[0] << " min: " << r.terminationtype << " in " << r.iterationscount << '\n';
//	if (theta[0] < 0 || theta[0] > pi/2)
//		throw std::runtime_error("theta out of bounds");
//	if (r.terminationtype < 0)
//		throw std::runtime_error("microroughness maximization did not converge properly");
	return MRDist(transmit, false, v, normal, Estep, RMSroughness, correlationLength, theta[0], 0);
}

}
