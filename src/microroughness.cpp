#include <complex>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

#include "microroughness.h"

bool TMicroRoughness::MRValid(const double v[3], const double normal[3], solid *leaving, solid *entering){
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; // velocity squared
	double vnormal = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]; // velocity projected onto surface normal
	double E = 0.5*m_n*v2; // kinetic energy
	double Estep = entering->mat.FermiReal*1e-9 - leaving->mat.FermiReal*1e-9;
	double ki = sqrt(2*m_n*E)*ele_e/hbar; // wave number in first solid
	material *mat;
	if (vnormal > 0)
		mat = &leaving->mat;
	else
		mat = &entering->mat;
	if (Estep >= 0){
		double kc = sqrt(2*m_n*Estep)*ele_e/hbar; // critical wave number of potential wall
		if (2*mat->RMSRoughness*ki < 1 && 2*mat->RMSRoughness*kc < 1)
			return true;
	}
	else{
		double kt = sqrt(2*m_n*(E - Estep))*ele_e/hbar; // wave number in second solid
		if (2*mat->RMSRoughness*ki < 1 && 2*mat->RMSRoughness*kt < 1)
			return true;
	}
	std::cout << "MR model not applicable on material " << mat->name << ". Falling back to Lambert model!\n";
	return false;
}

double TMicroRoughness::MRDist(bool transmit, bool integral, const double v[3], const double normal[3], solid *leaving, solid *entering, double theta, double phi){
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; // velocity squared
	double vnormal = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]; // velocity projected onto surface normal
	double E = 0.5*m_n*v2; // kinetic energy
	double Estep = entering->mat.FermiReal*1e-9 - leaving->mat.FermiReal*1e-9; // potential wall
	double costheta_i = abs(vnormal/sqrt(v2)); // cosine of angle between normal and incoming velocity vector
	if (transmit && E*costheta_i*costheta_i <= Estep) // if particle is transmitted: check if energy is higher than potential wall
		return 0;

	double ki = sqrt(2*m_n*E)*ele_e/hbar; // wave number in first solid
	std::complex<double> kc = std::sqrt(2*(double)m_n*std::complex<double>(Estep, 0.))*(double)ele_e/(double)hbar; // critical wave number of potential wall
	double kt = 0;

	std::complex<double> Si = 2*costheta_i/(costheta_i + sqrt(costheta_i*costheta_i - kc*kc/ki/ki)); // specularly transmitted amplitude
	std::complex<double> So;
	if (transmit){
		kt = sqrt(2*m_n*(E - Estep))*ele_e/hbar; // wave number in second solid
		So = 2*cos(theta)/(cos(theta) + sqrt(cos(theta)*cos(theta) + kc*kc/kt/kt)); // diffusely transmitted amplitude
	}
	else
		So = 2*cos(theta)/(cos(theta) + sqrt(cos(theta)*cos(theta) - kc*kc/ki/ki)); // diffusely reflected amplitude

	double b, w;
	if (vnormal > 0){
		b = leaving->mat.RMSRoughness;
		w = leaving->mat.CorrelLength;
	}
	else{
		b = entering->mat.RMSRoughness;
		w = entering->mat.CorrelLength;
	}

	double theta_i = acos(costheta_i);
	double Fmu = 0; // fourier transform of roughness correlation function
	if (!integral && !transmit)
		Fmu = exp(-w*w/2*ki*ki*(sin(theta_i)*sin(theta_i) + sin(theta)*sin(theta) - 2*sin(theta_i)*sin(theta)*cos(phi)));
	else if (!integral && transmit)
		Fmu = exp(-w*w/2*(ki*ki*sin(theta_i)*sin(theta_i) + kt*kt*sin(theta)*sin(theta) - 2*ki*kt*sin(theta_i)*sin(theta)*cos(phi)));
	else if (integral && !transmit) // if integral is set, precalculate phi-integral using modified Bessel function of first kind
		Fmu = exp(-w*w/2*ki*ki*(sin(theta_i)*sin(theta_i) + sin(theta)*sin(theta))) * 2*pi*boost::math::cyl_bessel_i(0, w*w*ki*ki*sin(theta_i)*sin(theta));
	else if (integral && transmit)
		Fmu = exp(-w*w/2*(ki*ki*sin(theta_i)*sin(theta_i) + kt*kt*sin(theta)*sin(theta))) * 2*pi*boost::math::cyl_bessel_i(0, w*w*ki*kt*sin(theta_i)*sin(theta));
	double factor = norm(kc)*norm(kc)*b*b*w*w/8/pi/costheta_i;
	if (transmit)
		factor *= kt/ki;
	if (integral)
		factor *= sin(theta); // integral over sin(theta) dtheta dphi
	double result = factor*norm(Si)*norm(So)*Fmu;
	//	cout << integral << transmit << ' ' << theta_i << ' ' << theta << ' ' << phi << ' ' << ki << ' ' << kt << ' ' << kc << ' ' << factor << ' ' << Si << ' ' << So << ' ' << Fmu << ' ' << result << '\n';
	return result; // return probability
}

void TMicroRoughness::operator()(std::vector<double> state, std::vector<double> &derivs, double x){
	derivs[0] = MRDist(ftransmit, fintegral, fv, fnormal, fleaving, fentering, x, 0);
}

void TMicroRoughness::MRDist(double theta, double xminusa, double bminusx, double &y, void *params){
	TMicroRoughness *p = (TMicroRoughness*)params;
	y = p->MRDist(p->ftransmit, p->fintegral, p->fv, p->fnormal, p->fleaving, p->fentering, theta, 0);
}

void TMicroRoughness::NegMRDist(const alglib::real_1d_array &x, double &f, void *params){
	TMicroRoughness *p = (TMicroRoughness*)params;
	f = -p->MRDist(p->ftransmit, p->fintegral, p->fv, p->fnormal, p->fleaving, p->fentering, x[0], 0);
}

double TMicroRoughness::MRProb(bool transmit, const double v[3], const double normal[3], solid *leaving, solid *entering){
	ftransmit = transmit;
	fv = v;
	fnormal = normal;
	fleaving = leaving;
	fentering = entering;

	double prob;
	alglib::autogkstate s;
	alglib::autogksmooth(0, pi/2, s);
	alglib::autogkintegrate(s, &MRDist, this);
	alglib::autogkreport r;
	alglib::autogkresults(s, prob, r);
	return prob;

}

double TMicroRoughness::MRDistMax(bool transmit, const double v[3], const double normal[3], solid *leaving, solid *entering){
	alglib::minbleicstate s;
	alglib::real_1d_array theta = "[0.7853981635]", bndl = "[0]", bndu = "[1.570796327]";
	alglib::minbleiccreatef(1, theta, 1e-6, s);
	alglib::minbleicsetbc(s, bndl, bndu);
	alglib::minbleicsetcond(s, 1e-6, 0, 0, 0);
	alglib::minbleicoptimize(s, &NegMRDist, NULL, this);
	alglib::minbleicreport r;
	alglib::minbleicresults(s, theta, r);
//	cout << "min: " << r.terminationtype << " in " << r.iterationscount << '\n';
	return MRDist(transmit, false, v, normal, leaving, entering, theta[0], 0);
}
