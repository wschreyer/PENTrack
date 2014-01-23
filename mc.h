/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <string>
#include <map>

#include "nr/ran.h"
#include "muParser.h"
#include "globals.h"

/**
 * For each section in particle.in such a struct is created containing all user options
 */
struct TParticleConfig{
	double tau; ///< lifetime
	double tmax; ///< max. simulation time
	double lmax; ///< max. trajectory length
	int polarization; ///< initial polarization
	double Emin; ///< min. initial energy
	double Emax; ///< max. initial energy
	mu::Parser spectrum; ///< Parsed energy spectrum given by user
	double phi_v_min; ///< Parsed minimum for initial azimuthal angle of velocity given by user
	double phi_v_max; ///< Parsed maximum for initial azimuthal angle of velocity given by user
	mu::Parser phi_v; ///< Parsed initial azimuthal angle distribution of velocity given by user
	double theta_v_min; ///< Parsed minimum for initial polarl angle of velocity given by user
	double theta_v_max; ///< Parsed maximum for initial polar angle of velocity given by user
	mu::Parser theta_v; ///< Parsed initial polar angle distribution of velocity given by user
};

/**
 * Class to generate random numbers in several distributions.
 */
class TMCGenerator{
private:
	Ran *rangen; ///< random number generator
	double xvar; ///< x variable for formula parser
	map<string, TParticleConfig> pconfigs;

public:
	Ullong seed; ///< initial random seed

	/**
	 * Constructor.
	 *
	 * Create random seed and read infile.
	 *
	 * @param infile Path to configuration file
	 */
	TMCGenerator(const char *infile){
		// get high resolution timestamp to generate seed
		timespec highrestime;
		clock_gettime(CLOCK_REALTIME, &highrestime);
		seed = (Ullong)highrestime.tv_sec * (Ullong)1000000000 + (Ullong)highrestime.tv_nsec;
		printf("Random Seed: %Lu\n\n", seed);
		rangen = new Ran(seed);

		TConfig invars; ///< contains variables from *.in file
		ReadInFile(infile, invars);
		for (TConfig::iterator i = invars.begin(); i != invars.end(); i++){
			if (i->first != "all")
				i->second = invars["all"]; // default particle specific settings to "all"-settings
		}
		ReadInFile(infile, invars); // read particle.in again to overwrite defaults for specific particle settings

		for (TConfig::iterator i = invars.begin(); i != invars.end(); i++){
			istringstream(i->second["Emax"]) >> pconfigs[i->first].Emax;
			istringstream(i->second["Emin"]) >> pconfigs[i->first].Emin;
			istringstream(i->second["lmax"]) >> pconfigs[i->first].lmax;
			istringstream(i->second["polarization"]) >> pconfigs[i->first].polarization;
			istringstream(i->second["tau"]) >> pconfigs[i->first].tau;
			istringstream(i->second["tmax"]) >> pconfigs[i->first].tmax;
			istringstream(i->second["phi_v_min"]) >> pconfigs[i->first].phi_v_min;
			istringstream(i->second["phi_v_max"]) >> pconfigs[i->first].phi_v_max;
			istringstream(i->second["theta_v_min"]) >> pconfigs[i->first].theta_v_min;
			istringstream(i->second["theta_v_max"]) >> pconfigs[i->first].theta_v_max;
			try{
				pconfigs[i->first].spectrum.DefineVar("x", &xvar);
				pconfigs[i->first].spectrum.SetExpr(i->second["spectrum"]);
				pconfigs[i->first].spectrum.DefineFun("ProtonBetaSpectrum", &ProtonBetaSpectrum);
				pconfigs[i->first].spectrum.DefineFun("ElectronBetaSpectrum", &ElectronBetaSpectrum);
				pconfigs[i->first].phi_v.DefineVar("x", &xvar);
				pconfigs[i->first].phi_v.SetExpr(i->second["phi_v"]);
				pconfigs[i->first].theta_v.DefineVar("x", &xvar);
				pconfigs[i->first].theta_v.SetExpr(i->second["theta_v"]);
			}
			catch (mu::Parser::exception_type &exc){
				cout << exc.GetMsg();
				exit(-1);
			}
		}
	};

	/**
	 * Destructor.
	 *
	 * Delete random number generator.
	 */
	~TMCGenerator(){
		delete rangen;
	};

	/// return uniformly distributed random number in [min..max]
	long double UniformDist(long double min, long double max){
		return rangen->doub()*(max - min) + min;
	};
	

	/// return sine distributed random number in [min..max] (in rad!)
	long double SinDist(long double min, long double max){
		return acos(cos(min) - UniformDist(0,1) * (cos(min) - cos(max)));
	};
	

	/// return sin(x)*cos(x) distributed random number in [min..max] (0 <= min/max < pi/2!)
	long double SinCosDist(long double min, long double max){
		return acos(sqrt(UniformDist(0,1)*(cos(max)*cos(max) - cos(min)*cos(min)) + cos(min)*cos(min)));
	};
	

	/// return x^2 distributed random number in [min..max]
	long double SquareDist(long double min, long double max){
		return pow(UniformDist(0,1)*(pow(max,3) - pow(min,3)) + pow(min,3),1.0L/3.0);
	};
	

	/// return linearly distributed random number in [min..max]
	long double LinearDist(long double min, long double max){
		return sqrt(UniformDist(0,1)*(max*max - min*min) + min*min);
	};


	/// return sqrt(x) distributed random number in [min..max]
	long double SqrtDist(long double min, long double max){
		return pow((pow(max, 1.5L) - pow(min, 1.5L))*UniformDist(0,1) + pow(min, 1.5L), 2.0L/3.0);
	};
	

	/**
	 * Create isotropically distributed 3D angles.
	 *
	 * @param phi Azimuth
	 * @param theta Polar angle
	 */
	void IsotropicDist(long double &phi, long double &theta){
		phi = UniformDist(0,2*pi);
		theta = SinDist(0,pi);
	};
	

	/// energy distribution of UCNs
	long double NeutronSpectrum(){
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

		return 0;
	};

	/**
	 * Energy distribution for each particle type
	 */
	double Spectrum(const string &particlename){
		double y;
		TParticleConfig *pconfig = &pconfigs[particlename];
		for (;;){
			xvar = UniformDist(pconfig->Emin, pconfig->Emax);
			try{
				y = pconfig->spectrum.Eval();
			}
			catch(mu::Parser::exception_type &exc){
				cout << exc.GetMsg();
				exit(-1);
			}
			if (UniformDist(0,1) < y)
				return xvar;
		}
		return 0;
	};


	/**
	 * Angular velocity distribution for each particle type
	 * 
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 * @param phi_v Returns velocity azimuth
	 * @param theta_v Returns velocity polar angle
	 */
	void AngularDist(const string &particlename, double &phi_v, double &theta_v){
		double y;
		TParticleConfig *pconfig = &pconfigs[particlename];
		for (;;){
			xvar = UniformDist(pconfig->phi_v_min, pconfig->phi_v_max);
			try{
				y = pconfig->phi_v.Eval();
			}
			catch(mu::Parser::exception_type &exc){
				cout << exc.GetMsg();
				exit(-1);
			}
			if (UniformDist(0,1) < y)
				phi_v = xvar;
		}
		for (;;){
			xvar = UniformDist(pconfig->theta_v_min, pconfig->theta_v_max);
			try{
				y = pconfig->theta_v.Eval();
			}
			catch(mu::Parser::exception_type &exc){
				cout << exc.GetMsg();
				exit(-1);
			}
			if (UniformDist(0,1) < y)
				theta_v = xvar;
		}
	}


	/**
	 * Lifetime of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns lifetimes using an exponentially decaying or flat distribution, depending on user choice in particle.in
	 */
	long double LifeTime(const string &particlename){
		double tau = pconfigs[particlename].tau;
		if (tau != 0)
			return -tau * log(UniformDist(0,1));
		else
			return pconfigs[particlename].tmax;
	};
	

	/**
	 * Max. trajectory length of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns max. trajectory length, depending on user choice in particle.in
	 */
	long double MaxTrajLength(const string &particlename){
		return pconfigs[particlename].lmax;
	};
	

	/**
	 * Initial polarisation of different particles
	 *
	 * @param particlename Name of the particle (chooses corresponding section in particle.in)
	 *
	 * @return Returns deiced of fixed polarisation (-1, 0, 1), depending on user choice in particle.in
	 */
	int DicePolarisation(const string &particlename){
		int p = pconfigs[particlename].polarization;
		if (p == 0){
			if(UniformDist(0,1) < 0.5)
				return -1;
			else 
				return 1;
		}
		else
			return p;
	};


	/**
	 * Simulate neutron beta decay.
	 *
	 * Calculates velocities of neutron decay products proton and electron.
	 *
	 * Reaction(s):
	 *
	 *     n0  ->  p+  +  R-
	 *
	 *     R-  ->  e-  +  nue
	 * 
	 * Procedure:
	 * (1) Dice the energy of the decay proton according to #ProtonBetaSpectrum in the rest frame of the neutron and
	 *     calculate the proton momentum via the energy momentum relation.
	 * (2) Dice isotropic orientation of the decay proton.
	 * (3) Calculate 4-momentum of rest R- via 4-momentum conservation.
	 * (4) Get fixed electron energy from two body decay of R-.
	 * (5) Dice isotropic electron orientation in the rest frame of R-.
	 * (6) Lorentz boost electron 4-momentum into moving frame of R-.
	 * (7) Calculate neutrino 4-momentum via 4-momentum conservation.
	 * (8) Boost all 4-momentums into moving neutron frame.
	 * 
	 * Cross-check:
	 * (9) Print neutrino 4-momentum invariant mass (4-momentum square, should be zero).
	 * 
	 * @param v_n Velocity of decayed neutron
	 * @param v_p Returns proton velocity
	 * @param v_e Returns electron velocity
	 */
	void NeutronDecay(long double v_n[3], long double v_p[3], long double v_e[3])
	{
		long double m_nue = 1 / pow(c_0, 2); // [eV/c^2]

	 	long double pabs, beta1, beta2, delta1, delta2; // 3-momentums (abs+directions)
	 	long double p[4], e[4];
	 	
	//-------- Step 1 ----------------------------------------------------------------------------------------------------------
	/*
	 * E   = E_kin + m*c^2
	 * p^2 = (E/c)^2 - (m*c)^2
	 *
	 */
		// energy of proton
	 	for (;;){
	 		long double E = UniformDist(0,751);
	 		if (UniformDist(0,1) < ProtonBetaSpectrum(E)){
	 			p[0] = E/c_0 + m_p*c_0;
	 			break;
	 		}
	 	}
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
		long double virt[4] = {m_n*c_0 - p[0], -p[1], -p[2], -p[3]}; // 4-momentum of intermediate virtual state (n - p)
		long double m2_virt = virt[0]*virt[0] - virt[1]*virt[1] - virt[2]*virt[2] - virt[3]*virt[3]; // squared mass of virtual state 		
			
	
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
		long double beta[3] = {virt[1]/virt[0], virt[2]/virt[0], virt[3]/virt[0]};
		BOOST(beta,e);
			
	//-------- Step 7 ----------------------------------------------------------------------------------------------------------
		// get 4-momentum of neutrino via 4-momentum conservation
		long double nue[4] = {virt[0] - e[0], virt[1] - e[1], virt[2] - e[2], virt[3] - e[3]};
	
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
		v_p[0] = p[1]/p[0]*c_0;
		v_p[1] = p[2]/p[0]*c_0;
		v_p[2] = p[3]/p[0]*c_0;
		v_e[0] = e[1]/e[0]*c_0;
		v_e[1] = e[2]/e[0]*c_0;
		v_e[2] = e[3]/e[0]*c_0;
	};
};

#endif /*MC_H_*/
