/**
 * \file
 * All about random numbers.
 */

#ifndef MC_H_
#define MC_H_

#include <cstdlib>
#include <sys/time.h>
#include <string>

#include "mersenne/mt.h"
#include "globals.h"

using namespace std;


/**
 * Class to generate random numbers in several distributions.
 */
struct TMCGenerator{
	/// struct to store limits for initial conditions read from all3inone.in
	struct initial{
		long double EnergieS; ///< minimum of energy spectrum
		long double EnergieE; ///< maximum of energy spectrum
		long double alphas; ///< minimum of starting angle between position vector and velocity vector projected onto xy-plane
		long double alphae; ///< maximum of starting angle between position vector and velocity vector projected onto xy-plane
		long double gammas; ///< minimum of starting angle between z axis and velocity vector
		long double gammae; ///< maximum of starting angle between z axis and velocity vector
		long double trajend; ///< max trajectory length
		long double xend; ///< max. lifetime
	};
	initial nini; ///< Initial condition limits for neutrons
	initial pini; ///< Initial condition limits for protons
	initial eini; ///< Initial condition limits for electrons
	unsigned long int monthinmilliseconds; ///< Random seed, derived from program start time
	int decay; ///< decay on (1/2) or off (0)
	int polarisation; ///< Polarisation parallel (1), antiparallel (2) or random (4)
	long double tau_n; ///< neutron lifetime
	long double decayoffset; ///< start decay timer after this time

	mt_state_t v_mt_state; ///< mersenne twister random number generator
	

	/**
	 * Constructor.
	 *
	 * Create random seed and read all3inone.in.
	 *
	 * @param infile Path to configuration file
	 * @param apolarisation Polarisation parallel (1), antiparallel (2) or random (4)
	 * @param adecay decay on (1/2) or off (0)
	 * @param adecayoffset Start decay timer after this time
	 * @param NeutLifetime Neutron lifetime
	 */
	TMCGenerator(const char *infile, int apolarisation, int adecay, long double adecayoffset, long double NeutLifetime){
		decay = adecay;
		polarisation = apolarisation;
		decayoffset = adecayoffset;
		tau_n = NeutLifetime;
		time_t mytime;
		tm *monthday;
		timeval daysec;
	
		// for random numbers we need to set an initial seed
		mytime = time(NULL);
		
		monthday = localtime(&mytime);
		monthinmilliseconds = (unsigned long int)(monthday->tm_mday)*24*60*60*1000; // add day in ms
		monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_hour)*60*60*1000; // add hour in ms
		monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_min)*60*1000; // add minute in ms 
		monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_sec)*1000;  // add second in ms	
		gettimeofday(&daysec, NULL);
		monthinmilliseconds = monthinmilliseconds + daysec.tv_usec/1000; // add milliseconds

		printf("Random Seed: %lu\n\n", monthinmilliseconds);
	
		mt_set (&v_mt_state, monthinmilliseconds);
		
		int i = 0, ncont;
	    string cline;
		string path;
	
		// creating 'inistream' to all3inone.in
		ifstream inistream(infile);
		if(!inistream.is_open())
		{	
			printf("Can't open %s\n", infile);
			exit(-1);
		}
		
		// looping  over the lines of all3inone.in and reading them in
		for(i = 1; i < 31; i++)
		{
			getline(inistream, cline);
			ncont = 42;
			switch (i)
			{
				case  2:	ncont = sscanf(cline.c_str(), "%LG %LG ", &nini.EnergieS, &nini.EnergieE);
							break;
				case  3:	ncont = sscanf(cline.c_str(), "%LG %LG ", &nini.alphas, &nini.alphae);
							nini.alphas *= conv;
							nini.alphae *= conv;
							break;
				case  4:	ncont = sscanf(cline.c_str(), "%LG %LG ", &nini.gammas, &nini.gammae);
							nini.gammas *= conv;
							nini.gammae *= conv;
							break;
				case 5:		ncont = sscanf(cline.c_str(), "%LG %LG ", &nini.trajend, &nini.xend);
							break;
				case 7:		ncont = sscanf(cline.c_str(), "%LG %LG ", &pini.EnergieS, &pini.EnergieE);
							break;
				case 8:		ncont = sscanf(cline.c_str(), "%LG %LG ", &pini.alphas, &pini.alphae);
							pini.alphas *= conv;
							pini.alphae *= conv;
							break;
				case 9:		ncont = sscanf(cline.c_str(), "%LG %LG ", &pini.gammas, &pini.gammae);
							pini.gammas *= conv;
							pini.gammae *= conv;
							break;
				case 10:	ncont = sscanf(cline.c_str(), "%LG %LG ", &pini.trajend, &pini.xend);
							break;
				case 12:	ncont = sscanf(cline.c_str(), "%LG %LG ", &eini.EnergieS, &eini.EnergieE);
							break;
				case 13:	ncont = sscanf(cline.c_str(), "%LG %LG ", &eini.alphas, &eini.alphae);
							eini.alphas *= conv;
							eini.alphae *= conv;
							break;
				case 14:	ncont = sscanf(cline.c_str(), "%LG %LG ", &eini.gammas, &eini.gammae);
							eini.gammas *= conv;
							eini.gammae *= conv;
							break;
				case 15:	ncont = sscanf(cline.c_str(), "%LG %LG ", &eini.trajend, &eini.xend);
							break;
			}
			if(ncont < 1) printf("an error occured while reading the %i. item in line %i of all3inone.in", ncont, i);
	    }
	    
		// starting value (S) > final value (E)?? EXIT!!
		for(int i = 1; i < 4; i++)
		{	struct initial particleini;
			switch(i)
			{	case 1:	particleini = nini;
						break;	
				case 2:	particleini = pini;
						break;
				case 3:	particleini = eini;
						break;
			}
			if((particleini.EnergieS > particleini.EnergieE) || (particleini.alphas > particleini.alphae) || (particleini.gammas > particleini.gammae))
			{	printf("\n\n\nERROR: Two or more initial values are inconsistent.\n"
				             "       Check ALL starting and final values in all3inone.in!\n\n"
				             "EXIT!!\n");
				exit(-1);
			}
		}
	};


	/// return uniformly distributed random number in [min..max]
	long double UniformDist(long double min, long double max){
		return mt_get_double(&v_mt_state)*(max - min) + min;
	};
	

	/// return sine distributed random number in [min..max] (in rad!)
	long double SinDist(long double min, long double max){
		return acos(cos(min) - mt_get_double(&v_mt_state) * (cos(min) - cos(max)));	
	};
	

	/// return sin(x)*cos(x) distributed random number in [min..max] (0 <= min/max < pi/2!)
	long double SinCosDist(long double min, long double max){
		return acos(sqrt(mt_get_double(&v_mt_state)*(cos(max)*cos(max) - cos(min)*cos(min)) + cos(min)*cos(min)));
	};
	

	/// return x^2 distributed random number in [min..max]
	long double SquareDist(long double min, long double max){
		return powl(mt_get_double(&v_mt_state) * (pow(max,3) - pow(min,3)) + pow(min,3),1.0/3.0);
	};
	

	/// return linearly distributed random number in [min..max]
	long double LinearDist(long double min, long double max){
		return sqrt(mt_get_double(&v_mt_state) * (max*max - min*min) + min*min);
	};


	/// return sqrt(x) distributed random number in [min..max]
	long double SqrtDist(long double min, long double max){
		return powl((powl(max, 1.5) - powl(min, 1.5)) * mt_get_double(&v_mt_state) + powl(min, 1.5), 2.0/3.0);
	};
	

	/**
	 * Create isotropically distributed 3D angles.
	 *
	 * @param alpha Angle between vector and x-axis projected onto xy-plane
	 * @param gamma Angle between vector and z-axis
	 */
	void IsotropicDist(long double &alpha, long double &gamma){
		alpha = UniformDist(0,2*pi);
		gamma = SinDist(0,pi);
	};
	

	/// energy distribution of UCNs (usually in [nini::EnergieS..nini::EnergieE])
	long double NeutronSpectrum(){
		return SqrtDist(nini.EnergieS*1e-9, nini.EnergieE*1e-9);

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
			x = EnergieS*1e9 + mt_get_double(v_mt_state) * (EnergieE*1e9 - EnergieS*1e9);
			WktTMP = mt_get_double(v_mt_state)*140;
			cout << "AbEx energy dist: E = " << x << " Wkt = " << WktTMP << '\n';
		}while(WktTMP > (-1*(p1 + 2*p2 * x+ 3*p3 * pow(x,2) + 4*p4 * pow(x,3) + 5*p5 * pow(x,4) + 6*p6 * pow(x,5))));
	
		cout << "below! Take this value E = " << x <<  '\n';
		Energie = NeutEnergie = x*1e-9;
		// END AbEx@ILL
		*/
	};
	

	/**
	 * Energy distribution of protons [0..750 eV]
	 *
	 * Proton recoil spectrum from neutron beta decay from "Diplomarbeit M. Simson"
	 */
	long double ProtonSpectrum(){
		long double Wahrschl, WktTMP, Energie;
		// dice as long as point is above curve
		do
		{	//cout << "above distribution... dicing on..." << '\n';
			Energie = UniformDist(0, 751); // constant distribution between 0 and 751 eV // [eV]
			WktTMP = UniformDist(0, 2);
			
			long double DeltaM = m_n - m_p;
			long double Xi = m_e / DeltaM;
			long double Sigma = 1 - 2*Energie * m_n / powl(DeltaM, 2) / powl(c_0, 2);
			long double g1 = powl(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
			long double g2 = powl(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
			long double littlea = -0.1017;
			Wahrschl = g1 + littlea * g2;
			//cout << "Proton energy dist: E = " << Energie << " decay rate: " << Wahrschl << " diced value = " << WktTMP << '\n';
	
		}while(WktTMP > Wahrschl);
		return Energie;
	};
	

	/**
	 * Energy distribution of electrons [0..782 keV]
	 *
	 * Electron energy spectrum from neutron beta decay from http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html
	 */
	long double ElectronSpectrum(){
		long double Elvert, Energie, WktTMP;
		long double Qvalue = 782e3;
		// dice as long as point is above curve
		do
		{
			Energie = UniformDist(0, 782e3); // constant distribution between 0 and 800 keV
			WktTMP = UniformDist(0, 0.12);
	
			Elvert = sqrtl(Energie*1e-6 * Energie*1e-6 + 2*Energie*1e-6 * m_e * c_0 * c_0*1e-6) * pow(Qvalue*1e-6 - Energie*1e-6, 2) * (Energie*1e-6 + m_e * c_0 * c_0*1e-6);			
	//		cout << "Electron energy dist: E = " << Energie << " decay rate: " << Elvert << " diced value = " << WktTMP << '\n';
	 
		}while(WktTMP > Elvert);
		return Energie;
	};
	

	/**
	 * Lifetime of different particles
	 *
	 * For neutrons it either returns nini::xend (from all3inone.in, if #decay=0) or a random exponential decay lifetime (if #decay =1 or =2 in config.in) with decay time #tau_n plus #decayoffset.
	 * For electrons and protons it always returns eini::xend and pini::xend.
	 */
	long double LifeTime(int protneut){
		switch (protneut){
			case NEUTRON:
				if (decay != 0)
					return decayoffset - tau_n * log(mt_get_double(&v_mt_state));
				else
					return nini.xend;
			case PROTON:
				return pini.xend;
			case ELECTRON:
				return eini.xend;
		}
		return 0;
	};
	

	/// max. trajectory length for different particles
	long double MaxTrajLength(int protneut){
		switch (protneut){
			case NEUTRON:
				return nini.trajend;
			case PROTON:
				return pini.trajend;
			case ELECTRON:
				return eini.trajend;
		}
		return 9e99;
	};
	

	/// get neutron polarisation, either diced (#polarisation==4) or fixed (#polarisation==1/2)
	int DicePolarisation(){
		if (polarisation == 4)
		{
			if(UniformDist(0,1) < 0.5)
				return -1;
			else 
				return 1;
		}
		else if (polarisation == POLARISATION_GOOD)
			return 1;
		else if (polarisation == POLARISATION_BAD)
			return -1;
		return 0;	
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
	 * (1) Dice the energy of the decay proton according to #ProtonSpectrum in the rest frame of the neutron and
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
		long double m_nue = 1 / powl(c_0, 2); // [eV/c^2]

	 	long double pabs, beta1, beta2, delta1, delta2; // 3-momentums (abs+directions)
	 	long double p[4], e[4];
	 	
	//-------- Step 1 ----------------------------------------------------------------------------------------------------------
	/*
	 * E   = E_kin + m*c^2
	 * p^2 = (E/c)^2 - (m*c)^2
	 *
	 */
		// energy of proton
		p[0] = ProtonSpectrum() / c_0 + m_p * c_0;
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
		long double decayerror = m_nue*c_0*c_0 - sqrt(nue[0]*nue[0] - nue[1]*nue[1] - nue[2]*nue[2] - nue[3]*nue[3])*c_0;
		printf("\n   +++ decay error : %LG eV +++\n", decayerror);
	
	
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
