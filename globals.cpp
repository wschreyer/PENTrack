#include <iostream>
#include <fstream>
#include <cmath>

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

// rotate vector into new coordinate sys whose z-axis lies on NORMALIZED vector n (active transformation)
void RotateVector(double v[3], const double n[3])
{
	double cosalpha = n[2], sinalpha = sqrt(1 - cosalpha*cosalpha);	// rotation angle (angle between z and n)
	if (sinalpha > 1e-30){ // when normal not parallel to z-axis rotate new velocity into the coordinate system where the normal is the z-axis
		double a[2] = {-n[1]/sinalpha, n[0]/sinalpha};	// rotation axis (z cross n), a[2] = 0
		double vtemp[3] = {v[0],v[1],v[2]};
		// rotate velocity vector
		v[0] = (cosalpha + a[0]*a[0]*(1 - cosalpha))*	vtemp[0] +  a[0]*a[1]*(1 - cosalpha)*				vtemp[1] + a[1]*sinalpha*	vtemp[2];
		v[1] =  a[1]*a[0]*(1 - cosalpha)*				vtemp[0] + (cosalpha + a[1]*a[1]*(1 - cosalpha))*	vtemp[1] - a[0]*sinalpha*	vtemp[2];
		v[2] = -a[1]*sinalpha*							vtemp[0] +  a[0]*sinalpha*							vtemp[1] + cosalpha*		vtemp[2];
	}
	else if (cosalpha < 0){
		v[0] = -v[0];
		v[1] = -v[1];
		v[2] = -v[2];
	}
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
