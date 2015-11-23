#include <iostream>
#include <fstream>
#include <cmath>

#include <CGAL/Simple_cartesian.h>

#include "globals.h"

long long int jobnumber = 0; ///< job number, read from command line paramters, used for parallel calculations
std::string inpath = "."; ///< path to configuration files, read from command line paramters
std::string outpath = "."; ///< path where the log file should be saved to, read from command line parameters

// print progress in percent
void PrintPercent(double percentage, int &lastprint){
	// write status to console
	// one point per 2 percent of endtime
	while (lastprint < percentage*100){
		lastprint += 2;
		if (lastprint % 10 == 0)
			std::cout << lastprint << "%";
		else
			std::cout << ".";
//		std::cout.flush();
	}
}

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
void BOOST(double beta[3], double p[4]){
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

// energy distribution of protons (0 < 'Energie' < 750 eV)
// proton recoil spectrum from "Diplomarbeit M. Simson"
// result always < 1!
double ProtonBetaSpectrum(double E){
	double DeltaM = m_n - m_p;
	double Xi = m_e / DeltaM;
	double Sigma = 1 - 2*E*m_n / pow(DeltaM, 2) / pow(c_0, 2);
	double g1 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 				4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double g2 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 	4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double littlea = -0.1017;
	return 0.5*(g1 + littlea * g2);
}


// energy distribution of electrons (0 < 'Energie' < 782 keV)
// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
// result always < 1!
double ElectronBetaSpectrum(double E){
	double Qvalue = 0.782; //[MeV]
	return 8.2*sqrt(E*1e-6*E*1e-6 + 2*E*1e-6*m_e*c_0*c_0*1e-6) * pow(Qvalue - E*1e-6, 2) * (E*1e-6 + m_e*c_0*c_0*1e-6);
}

// energy distribution of comagnetometer gases using Maxwell-Boltzmann distribution
// from "http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html"
// result always < 1!
double MaxwellBoltzSpectrum (double T, double E) {
	double kT = boltzconst*T;  
	return 2*sqrt(E/pi)*sqrt(kT*kT*kT)*exp(-E/kT);
}

typedef std::map<std::string, std::map<std::string, std::string> > TConfig;

//read variables from *.in file into map
void ReadInFile(const char *inpath, TConfig &vars){
	std::ifstream infile(inpath);
	char c;
	std::string rest,section,key;
	while (infile && (infile >> std::ws) && (c = infile.peek())){
		if (c == '[' && infile.ignore()){
			if (infile.peek() == '/'){
				section = "";
			}
			else{
				std::getline(infile, section, ']');
//				std::cout << "\nsection: " << section.c_str() << '\n';
			}
			std::getline(infile,rest);
		}
		else if (c == '#')
			std::getline(infile,rest);
		else if (section != ""){
			infile >> key;
			std::getline(infile,rest);
			if (infile){
				std::string::size_type l = rest.find('#');
				if (l == std::string::npos)
					vars[section][key] = rest;
				else
					vars[section][key] = rest.substr(0,l);
//				std::cout << key << " " << vars[section][key] << '\n';
			}
		}
		else
			std::getline(infile,rest);
	}
}
