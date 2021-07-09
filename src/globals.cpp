#include "globals.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include <CGAL/Simple_cartesian.h>

#include <boost/format.hpp>

std::atomic<bool> quit(false);

const long double pi = 3.1415926535897932384626L; ///< Pi
const long double ele_e = 1.602176487E-19L; ///< elementary charge [C]
const long double gravconst = 9.80665L; ///< g [m/s]
const long double boltzconst = 1.38064852E-23L/ele_e; /// Boltzmann's constant [ev/K] from http://physics.nist.gov/cgi-bin/cuu/Value?sigma
const long double avogadroconst = 6.022140857E23L; /// Avogadro's constant [ mol^-1 ] from http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avagadro%27s+number
const long double conv = pi/180.L; ///< deg to rad conversion factor
const long double mu0 = 4*pi*1e-7L; ///< magnetic permeability [Vs/Am]
const long double m_n = 1.674927211E-27L/ele_e; ///< neutron mass [eV/c^2]
const long double m_p = 1.672621637E-27L/ele_e; ///< proton mass [eV/c^2]
const long double m_e = 9.10938215e-31L/ele_e; ///< electron mass [eV/c^2]
const long double m_hg = 198.96828064/(1000*avogadroconst)/ele_e; ///< mercury-199 mass [ev/c^2] (http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Hg)
const long double m_xe = 128.9047808611/(1000*avogadroconst)/ele_e; ///< xenon-129 mass [ev/c^2] (http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Xe)
const long double c_0 = 299792458.L; ///< light speed [m/s]
const long double hbar = 1.05457266e-34L; ///< planck constant [Js]
const long double mu_nSI = -0.96623641e-26L;	///< Neutron magnetic moment [J/T]
const long double mu_hgSI = 2.555118e-27L; ///< mercury-199 magnetic moment [J/T]
const long double mu_xeSI = -3.392939e-27L; ///< Xenon-129 avogadroconstmagnetic moment [J/T]
const long double gamma_n = -1.83247185e8L; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]
const long double gamma_hg = 4.76901003e7L; ///< from: http://www.sciencedirect.com/science/article/pii/S0370269314007692 [ 1/Ts ]
const long double gamma_xe = -7.399707336e7L; ///< from: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio [ 1/Ts ]

long long int jobnumber = 0; ///< job number, read from command line paramters, used for parallel calculations
boost::filesystem::path configpath = boost::filesystem::current_path() / "in/config.in"; ///< path to configuration files, read from command line paramters
boost::filesystem::path outpath = boost::filesystem::current_path() / "out/"; ///< path where the log file should be saved to, read from command line parameters


// rotate vector into new coordinate whose z axis is parallel to n and whose x axis is defined by the projection of x onto the plane defined by n (active transformation)
void RotateVector(double v[3], const double n[3], const double x[3])
{
	typedef CGAL::Simple_cartesian<double> K;
	CGAL::Vector_3<K> vv(v[0], v[1], v[2]), xv, zv(n[0], n[1], n[2]);
	CGAL::Plane_3<K> plane(CGAL::ORIGIN, zv);
	if (x == NULL) // if no x is given, choose random vector on plane as x
		xv = plane.base1();
	else{
		CGAL::Point_3<K> xp(x[0], x[1], x[2]);
		xp = plane.projection(xp); // project x onto plane defined by normal n
		xv = xp - CGAL::ORIGIN;
		if (xv.squared_length() < 1e-30)
			xv = plane.base1(); // if x is parallel to n choose some random vector on plane as x axis
	}
	xv = xv/sqrt(xv.squared_length()); // build orthonormal basis of new coordinate system
	zv = zv/sqrt(zv.squared_length());
	CGAL::Vector_3<K> yv = CGAL::cross_product(zv, xv); // new y-axis corresponds to cross product of normal and velocity
	CGAL::Aff_transformation_3<K> rotmatrix(xv.x(), yv.x(), zv.x(), xv.y(), yv.y(), zv.y(), xv.z(), yv.z(), zv.z());
//	std::cout << "(" << xv << "," << yv << "," << zv << ") + " << vv << " = ";
	vv = vv.transform(rotmatrix); // transform v into new coordinate system
//	std::cout << vv << std::endl;
	v[0] = vv.x();
	v[1] = vv.y();
	v[2] = vv.z();
}

//======== Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta ======================================================
void BOOST(const std::vector<double> &beta, std::vector<double> &p){
   //Boost this Lorentz vector (copy&paste from ROOT)
   double b2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
   double gamma = 1.0 / sqrt(1.0 - b2);
   double bp = beta[0]*p[1] + beta[1]*p[2] + beta[2]*p[3];
   double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

   p[1] = (p[1] + gamma2*bp*beta[0] + gamma*beta[0]*p[0]);
   p[2] = (p[2] + gamma2*bp*beta[1] + gamma*beta[1]*p[0]);
   p[3] = (p[3] + gamma2*bp*beta[2] + gamma*beta[2]*p[0]);
   p[0] = (gamma*(p[0] + bp));
}
//======== end of BOOST ====================================================================================================

// energy distribution of protons (0 < E < 750 eV)
// proton recoil spectrum from "Diplomarbeit M. Simson"
// result always < 1!
double ProtonBetaSpectrum(const double E){
	double DeltaM = m_n - m_p;
	double Xi = m_e / DeltaM;
	double Sigma = 1 - 2*E*m_n / pow(DeltaM, 2) / pow(c_0, 2);
	double g1 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 				4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double g2 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 	4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double littlea = -0.1017;
	return 0.5*(g1 + littlea * g2);
}


// energy distribution of electrons (0 < E < 782 keV)
// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
// result always < 1!
double ElectronBetaSpectrum(const double E){
	double Qvalue = 0.782; //[MeV]
	return 8.2*sqrt(E*1e-6*E*1e-6 + 2*E*1e-6*m_e*c_0*c_0*1e-6) * pow(Qvalue - E*1e-6, 2) * (E*1e-6 + m_e*c_0*c_0*1e-6);
}

// energy distribution of comagnetometer gases using Maxwell-Boltzmann distribution
// from en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
// result always < 1!
double MaxwellBoltzSpectrum (const double T, const double E) {
	double kT = boltzconst*T;  
	return 2*sqrt(E/pi)/sqrt(kT*kT*kT)*exp(-E/kT) / (2*sqrt(0.5/pi)/kT*exp(-0.5)); // return distribution divided by its maximum at E/kt = 0.5
}


std::string ResolveFormula(const std::string &formulaName, const std::map<std::string, std::string> &formulas){
	auto i = formulas.find(formulaName);
	if (i == formulas.end()){
		return formulaName;
	}
	else{
		return i->second;
	}
}