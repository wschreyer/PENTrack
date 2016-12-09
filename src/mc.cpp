#include "mc.h"

#include <chrono>

#include "exprtk.hpp"

#include "globals.h"

const int PIECEWISE_LINEAR_DIST_INTERVALS = 1000;

TMCGenerator::TMCGenerator(const uint64_t aseed){
	if (aseed == 0){
		// get high-resolution timestamp to generate seed
		using namespace std::chrono;
		seed = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
	}
	else
		seed = aseed;
	std::cout << "Random Seed: " << seed << "\n\n";
	rangen.seed(seed);

	ProtonBetaSpectrumDist = ParseDist(&ProtonBetaSpectrum, 0., 750.); // create piecewise linear distribution from predefined energy spectrum of protons from free-neutron decay
	ElectronBetaSpectrumDist = ParseDist(&ElectronBetaSpectrum, 0., 782000.); // create piecewise linear distribution from predefined energy spectrum of protons from free-neutron decay
}


double TMCGenerator::UniformDist(const double min, const double max) const{
	if (min == max)
		return min;
	std::uniform_real_distribution<double> unidist(min, max);
	return unidist(rangen);
}

double TMCGenerator::SinDist(const double min, const double max) const{
	return acos(cos(min) - UniformDist(0,1) * (cos(min) - cos(max)));
}

double TMCGenerator::SinCosDist(const double min, const double max) const{
	return acos(sqrt(UniformDist(0,1)*(cos(max)*cos(max) - cos(min)*cos(min)) + cos(min)*cos(min)));
}

double TMCGenerator::SquareDist(const double min, const double max) const{
	return pow(UniformDist(0,1)*(pow(max,3) - pow(min,3)) + pow(min,3),1.0L/3.0);
}

double TMCGenerator::LinearDist(const double min, const double max) const{
	return sqrt(UniformDist(0,1)*(max*max - min*min) + min*min);
}

double TMCGenerator::SqrtDist(const double min, const double max) const{
	return pow((pow(max, 1.5) - pow(min, 1.5))*UniformDist(0,1) + pow(min, 1.5), 2.0/3.0);
}

double TMCGenerator::NormalDist(const double mean, const double sigma) const{
	std::normal_distribution<double> normaldist(mean, sigma);
	return normaldist(rangen);
}

double TMCGenerator::ExpDist(const double exponent) const{
	std::exponential_distribution<double> expdist(exponent);
	return expdist(rangen);
}

void TMCGenerator::IsotropicDist(double &phi, double &theta) const{
	phi = UniformDist(0,2*pi);
	theta = SinDist(0,pi);
}


//double TMCGenerator::NeutronSpectrum() const{
//		return SqrtDist(0, 300e-9);

/*
	//neutron energy spectrum for PENeLOPE (storage only) 180cm above source and 10cm absorber
	long double x,y;
	for(;;){
		x = UniformDist(0, 130);
		y = UniformDist(0,500);
		if ((x <= 53.2 && y <= (0.0817*x*x + 5.106*x - 4))
			|| (x > 53.2 && x <= 71.4 && y <= 487)
			|| (x > 71.4 && y <= (1.727e-4*x*x*x*x - 0.07839*x*x*x + 13.2296*x*x - 985.9*x + 27485)))
			return x*1e-9;
	}
*/
/*
	//neutron energy spectrum for PENeLOPE (storage+buffer) 180cm above source and 10cm absorber
	long double x,y;
	for(;;){
		x = UniformDist(-18, 130);
		y = UniformDist(0,585);
		if ((x <= 56 && y <= (8*x + 140))
			|| (x > 56 && x <= 70 && y <= -3.87*x + 802)
			|| (x > 70 && y <= (-3.719e-5*x*x*x*x + 7.585e-3*x*x*x + 0.1614*x*x - 113.3*x + 5.951e3)))
			return x*1e-9;
	}
*/
/*
	//low-field-seeker spectrum after ramping (storage+buffer) 180cm above source and 10cm absorber (outer&inner)
	long double x,y;
	for(;;){
		x = UniformDist(-12, 105);
		y = UniformDist(0,1900);
		if ((x <= 75 && y <= 22*x + 246)
			|| (x > 75 && x <= 87 && y <= 1900)
			|| (x > 87 && y <= -116.4*x + 12092))
			return x*1e-9;
	}
*/
/*
	//low-field-seeker spectrum after ramping (storage only) 180cm above source and 10cm absorber (outer&inner)
	long double x,y;
	for(;;){
		x = UniformDist(11,107);
		y = 6.74339e-11*pow(x,8) - 2.56528e-8*pow(x,7) + 3.89341e-6*pow(x,6) - 0.000301768*pow(x,5) + 0.0126685*pow(x,4) - 0.284483*pow(x,3) + 3.54352*pow(x,2) - 23.9924*x + 72.2083;
		// polynom 8th order
		if (UniformDist(0,1060)  < y)
			return x*1e-9;
	}
*/
/*
	//high-field-seeker spectrum after ramping (storage+buffer) 180cm above source and 10cm absorber (outer&inner)
	long double x,y;
	for(;;){
		x = UniformDist(-200,50);
		y = UniformDist(0,751);
		// two gaussians
		if ((x < -35 && y < 561*exp(-(x + 39.30)*(x + 39.30)/2/43.50/43.50))
		|| (x >= -35 && y < 751*exp(-(x + 17.60)*(x + 17.60)/2/19.69/19.69)))
			return x*1e-9;
	}
*/
/*
	//spectrum entering buffervolume from guide
	long double x,y;
	for(;;){
		x = UniformDist(0,150);
		y = -3.142e-10*pow(x,7) + 9.549e-8*pow(x,6) - 6.808e-6*pow(x,5) - 0.000414*pow(x,4) + 0.0635*pow(x,3) - 2.672*pow(x,2) + 131.12*x - 23.6645;
		if (UniformDist(0,5120) < y)
			return x*1e-9;
	}
*/
/*
	//spectrum leaving horizontal guide
	long double x,y;
	for(;;){
		x = UniformDist(100,300);
		y = 1.96616e-6*pow(x,5) - 0.00204264*pow(x,4) + 0.834378*pow(x,3) - 167.958*pow(x,2) + 16674.8*x - 639317;
		if (UniformDist(0,14000) < y)
			return x*1e-9;
	}
*/
/*
	//lfs spectrum after ramping with 0.5B and absorber 34cm above storage bottom
	long double x,y;
	for(;;){
		x = UniformDist(10,58);
		y = UniformDist(0,180);
		if ((x <= 49 && y < 4.86222*x - 58.4666)
		|| (x > 49 && y < exp(-0.616279*x + 35.3488)))
			return x*1e-9;
	}
*/
/*
	//neutron energy spectrum by Gerd Petzoldt
	long double x,y,cutoff = 900,decay = 0.05;
	for(;;){
		x = UniformDist(nini.EnergieS, nini.EnergieE);
		y = UniformDist(0,cutoff);
		if (x <= cutoff && y <= x) return (x+100)*1e-9; // linear spectrum until cutoff
		else if (x > cutoff && y < cutoff*exp(-decay*(x - cutoff))) return (x+100)*1e-9; // epxonential decay after cutoff
	}
*/
	/*
	// neutron energy distribution for AbEx@ILL
	long double p1 = -77.6138, p2 = -1.34704, p3 = 0.00739579, p4 = 0.00012494, p5 = -1.88103e-6 , p6 = 8.52798e-9;
	do
	{	cout << "above distribution... dicing on..." << '\n';
		x = mt_get_double(v_mt_state) * 300e-9;
		WktTMP = mt_get_double(v_mt_state)*140;
		cout << "AbEx energy dist: E = " << x << " Wkt = " << WktTMP << '\n';
	}while(WktTMP > (-1*(p1 + 2*p2 * x+ 3*p3 * pow(x,2) + 4*p4 * pow(x,3) + 5*p5 * pow(x,4) + 6*p6 * pow(x,5))));

	cout << "below! Take this value E = " << x <<  '\n';
	Energie = NeutEnergie = x*1e-9;
	// END AbEx@ILL
	*/

//	return 0;
//}

double TMCGenerator::DicePolarisation(const double polarisation) const{
	if(UniformDist(0,1) < 0.5*(1. + polarisation))
		return 1;
	else
		return -1;
}

void TMCGenerator::NeutronDecay(const double v_n[3], double &E_p, double &E_e, double &phi_p, double &phi_e, double &theta_p, double &theta_e, double &pol_p, double &pol_e) const
{
	double m_nue = 1 / pow(c_0, 2); // [eV/c^2]

 	double pabs, beta1, beta2, delta1, delta2; // 3-momentums (abs+directions)
 	double p[4], e[4];

//-------- Step 1 ----------------------------------------------------------------------------------------------------------
/*
 * p   = (E/c, p_x, p_y, p_z)
 * E   = E_kin + m*c^2
 * p^2 = (E/c)^2 - (m*c)^2
 *
 */
	// energy of proton
 	double Ekin = ProtonBetaSpectrumDist(rangen);
	p[0] = Ekin/c_0 + m_p*c_0;
	// momentum norm of proton
	pabs = sqrt(p[0]*p[0] - m_p*m_p*c_0*c_0);

//-------- Step 2 ----------------------------------------------------------------------------------------------------------
	// isotropic emission characteristics (in the rest frame of the neutron)
	IsotropicDist(beta1,delta1);
	// 3-momentum of the proton
	p[1] = pabs * sin(delta1) * cos(beta1);
	p[2] = pabs * sin(delta1) * sin(beta1);
	p[3] = pabs * cos(delta1);

//-------- Step 3 ----------------------------------------------------------------------------------------------------------
	// calculate intermediate virtual state via 4-momentum conservation
	double virt[4] = {(double)m_n*(double)c_0 - p[0], -p[1], -p[2], -p[3]}; // 4-momentum of intermediate virtual state (n - p)
	double m2_virt = virt[0]*virt[0] - virt[1]*virt[1] - virt[2]*virt[2] - virt[3]*virt[3]; // squared mass of virtual state


//-------- Step 4 ----------------------------------------------------------------------------------------------------------

	// energy of electron from two-body-decay of virtual state
	e[0] = (m2_virt + m_e*m_e*c_0*c_0 - m_nue*m_nue*c_0*c_0)/2/sqrt(m2_virt);
	pabs = sqrt(e[0]*e[0] - m_e*m_e*c_0*c_0);

//-------- Step 5 ----------------------------------------------------------------------------------------------------------
	// isotropic emission characteristics (in the rest frame of the virtual state)
	IsotropicDist(beta2, delta2);
	// 3-momentum of the electron
	e[1] = pabs * sin(delta2) * cos(beta2);
	e[2] = pabs * sin(delta2) * sin(beta2);
	e[3] = pabs * cos(delta2);

//-------- Step 6 ----------------------------------------------------------------------------------------------------------
	// boost e into moving frame of virtual state
	double beta[3] = {virt[1]/virt[0], virt[2]/virt[0], virt[3]/virt[0]};
	BOOST(beta,e);

//-------- Step 7 ----------------------------------------------------------------------------------------------------------
	// get 4-momentum of neutrino via 4-momentum conservation
	double nue[4] = {virt[0] - e[0], virt[1] - e[1], virt[2] - e[2], virt[3] - e[3]};

//-------- Step 8 ----------------------------------------------------------------------------------------------------------
	// boost p,e,nu into moving frame of neutron
	beta[0] = v_n[0]/c_0;
	beta[1] = v_n[1]/c_0;
	beta[2] = v_n[2]/c_0;
	BOOST(beta,e);
	BOOST(beta,p);
	BOOST(beta,nue);

//-------- Step 9 ----------------------------------------------------------------------------------------------------------
	// use neutrino mass for check
//		long double decayerror = m_nue*c_0*c_0 - sqrt(nue[0]*nue[0] - nue[1]*nue[1] - nue[2]*nue[2] - nue[3]*nue[3])*c_0;
//		printf("\n   +++ decay error : %LG eV +++\n", decayerror);


//-------- Finished --------------------------------------------------------------------------------------------------------
	double beta_p2 = (p[1]*p[1] + p[2]*p[2] + p[3]*p[3])/p[0]/p[0];
	E_p = m_p*c_0*c_0*beta_p2*(0.5 + beta_p2*(3./8. + beta_p2*(5./16.))); // use series expansion for low energy proton calculations
	E_e = e[0]*c_0 - m_e*c_0*c_0; // use fully relativistic formula for high energy electron
	phi_p = atan2(p[2], p[1]);
	phi_e = atan2(e[2], e[1]);
	theta_p = acos(p[3]/sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
	theta_e = acos(e[3]/sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]));
	pol_p = UniformDist(0, 1) < 0.5 ? 1 : -1;
	pol_e = UniformDist(0, 1) < 0.5 ? 1 : -1;
}


void TMCGenerator::tofDist(double &Ekin, double &phi, double &theta) const{
	long double v_x,v_y, v_tot, xz_ang, v_yaxis, v_xaxis, v_zaxis;
	v_tot = sqrt(2*Ekin/m_n);
//	int i =0;
	for(;;){
//		i++;
		v_x = UniformDist(0, v_tot);
		v_y = (1/(2.270*v_x + 0.0122*pow(v_x,2)))*exp(-pow((log(2.270/v_x +0.0122) +1.4137),2)/(2*pow(0.3420,2)));

		if(UniformDist(0,0.04570423)<v_y){
			v_yaxis = v_x;
			xz_ang = UniformDist(0,2*pi);
			v_xaxis = sqrt(pow(v_tot,2)-pow(v_yaxis,2))*cos(xz_ang);
			v_zaxis = sqrt(pow(v_tot,2)-pow(v_yaxis,2))*sin(xz_ang);
			phi = atan2(v_yaxis,v_xaxis);
			theta = acos(v_zaxis/v_tot);
			//check x component
			if(v_xaxis - v_tot*cos(phi)*sin(theta) > 0.01){
				std::cout<< "There was an error in the calculation of v_x \n";
//					sleep(1);
			}else if(v_yaxis - v_tot*sin(phi)*sin(theta) > 0.01){
				std::cout<< "There was an error in the calculation of  v_y \n";
//					sleep(1);
			}else if(v_zaxis - v_tot*cos(theta) > 0.01){
				std::cout<< "There was an error in the calculation of v_z \n";
//					sleep(1);
			}
			return;
		}
/*		if(i > 100){ //show that there is a problem by setting phi and theta to 100
			phi = i;
			theta = i;
			return;
		}
*/

	}
}



std::piecewise_linear_distribution<double> TMCGenerator::ParseDist(const std::string &func, const double range_min, const double range_max) const{
	if (range_min > range_max)
		throw std::runtime_error("Invalid range used to parse distribution");
	std::vector<double> range, dist;
	double x;
	exprtk::symbol_table<double> symbol_table;
	symbol_table.add_variable("x", x);
	symbol_table.add_function("ProtonBetaSpectrum", ProtonBetaSpectrum);
	symbol_table.add_function("ElectronBetaSpectrum", ElectronBetaSpectrum);
	symbol_table.add_function("MaxwellBoltzSpectrum", MaxwellBoltzSpectrum);
	symbol_table.add_constants();
	
	exprtk::expression<double> expression;
	expression.register_symbol_table(symbol_table);

	exprtk::parser<double> parser;
	parser.compile(func, expression);

	for (x = range_min; x <= range_max; x += (range_max - range_min)/PIECEWISE_LINEAR_DIST_INTERVALS){ // go through parameter range
		range.push_back(x); // create list of parameter values
		dist.push_back(expression.value()); // create list of function values
		if (range_min == range_max){ // if min and max of parameter range are equal
			x = std::nextafter(range_max, std::numeric_limits<double>::max()); // add another parameter value slightly larger than min, else piecewise linear distribution defaults to range 0..1
			range.push_back(x);
			dist.push_back(expression.value());
			break; // stop list creation
		}
	}

	return std::piecewise_linear_distribution<double>(range.begin(), range.end(), dist.begin()); // create piecewise linear distribution from parameter- and distribution-value lists
}

std::piecewise_linear_distribution<double> TMCGenerator::ParseDist(double (*func)(double), const double range_min, const double range_max) const{
	if (range_min > range_max)
		throw std::runtime_error("Invalid range used to parse distribution");
	std::vector<double> range, dist;
	for (double x = range_min; x <= range_max; x += (range_max - range_min)/PIECEWISE_LINEAR_DIST_INTERVALS){ // go through parameter range
		range.push_back(x); // create list of parameter values
		dist.push_back((*func)(x)); // create list of function values
		if (range_min == range_max){ // if min and max of parameter range are equal
			x = std::nextafter(range_max, std::numeric_limits<double>::max()); // add another parameter value slightly larger than min, else piecewise linear distribution defaults to range 0..1
			range.push_back(x);
			dist.push_back((*func)(x));
			break; // stop list creation
		}
	}

	return std::piecewise_linear_distribution<double>(range.begin(), range.end(), dist.begin()); // create piecewise linear distribution from parameter- and distribution-value lists
}
