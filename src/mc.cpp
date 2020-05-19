#include "mc.h"

#include <chrono>
#include <iostream>
#include "exprtk.hpp"

#include "globals.h"

const int PIECEWISE_LINEAR_DIST_INTERVALS = 1000;

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
/*
void TMCGenerator::tofDist(double &Ekin, double &phi, double &theta) const{
	long double v_x,v_y, v_tot, xz_ang, v_yaxis, v_xaxis, v_zaxis;
	v_tot = sqrt(2*Ekin/m_n);
//	int i =0;
	for(;;){
//		i++;
		std::uniform_real_distribution<double> vdist(0, 25);
		v_x = vdist(rangen);
		v_y = (1/(2.270*v_x + 0.0122*pow(v_x,2)))*exp(-pow((log(2.270/v_x +0.0122) +1.4137),2)/(2*pow(0.3420,2)));

		std::uniform_real_distribution<double> vydist(0., 0.04570423);
		if(vydist(rangen)<v_y){
			v_yaxis = v_x;
			std::uniform_real_distribution<double> phidist(0, 2.*pi);
			xz_ang = phidist(rangen);
			std::sincos_distribution<double> sincosdist(0, acos(v_yaxis/25));
			theta = sincosdist(rangen);
			v_tot = v_yaxis/cos(theta);
			v_xaxis = sqrt(pow(v_tot,2)-pow(v_yaxis,2))*cos(xz_ang);
			v_zaxis = sqrt(pow(v_tot,2)-pow(v_yaxis,2))*sin(xz_ang);
			phi = atan2(v_yaxis,v_xaxis);
			theta = acos(v_zaxis/v_tot);
			Ekin = 0.5*m_n*v_tot*v_tot;
			//check x component
			if(v_xaxis - v_tot*cos(phi)*sin(theta) > 0.01){
				throw std::runtime_error("There was an error in the calculation of v_x");
//					sleep(1);
			}else if(v_yaxis - v_tot*sin(phi)*sin(theta) > 0.01){
				throw std::runtime_error("There was an error in the calculation of  v_y");
//					sleep(1);
			}else if(v_zaxis - v_tot*cos(theta) > 0.01){
				throw std::runtime_error("There was an error in the calculation of v_z");
//					sleep(1);
			}
			return;
		}
		if(i > 100){ //show that there is a problem by setting phi and theta to 100
			phi = i;
			theta = i;
			return;
		}


	}
}

*/

template<typename UnaryFunction>
std::piecewise_linear_distribution<double> parse_distribution(UnaryFunction f, const double range_min, const double range_max){
	if (range_min > range_max)
		throw std::runtime_error("Invalid range used to parse distribution");
	double rmin = range_min;
	double rmax = range_max;
	int nw = PIECEWISE_LINEAR_DIST_INTERVALS;
	if (range_min == range_max){
		rmax = std::nextafter(range_max, std::numeric_limits<double>::max()); // use slightly larger double value for range_max, else distribution will default to range 0..1
		nw = 1;
	}

	return std::piecewise_linear_distribution<double>(nw, rmin, rmax, f);
}


std::piecewise_linear_distribution<double> parse_distribution(const std::string &func, const double range_min, const double range_max){
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
	if (not parser.compile(func, expression)){
	    throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing formula '" + func + "': " + parser.get_error(0).diagnostic);
	}
	return parse_distribution(
			[&x, &expression](const double px){
				x = px;
				return expression.value();
			},
			range_min,
			range_max
	);
}

std::piecewise_linear_distribution<double> proton_beta_distribution = parse_distribution(ProtonBetaSpectrum, 0., 750.);
std::piecewise_linear_distribution<double> electron_beta_distribution = parse_distribution(ElectronBetaSpectrum, 0., 782000.);


