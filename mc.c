#include <iostream>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <string>

#include "mc.h"
#include "globals.h"
#include "mersenne/mt.h"
#include "particle.h"
#include "geometry.h"

using namespace std;

struct initial												// record for initial values of all 3 particle types (for more details see "all3inone.in")
{	long double EnergieS, EnergieE;				// E
	long double alphas, alphae;						// alpha
	long double gammas, gammae;						// gamma
	long double hmax, xend;									// delx | xend
} nini, pini, eini;						// one 'initial' for each particle type

void Startbed(const char *infile); // read initial values from all3inone.in
void BOOST(long double beta[3], long double v[4]);


unsigned long int monthinmilliseconds; // RandomSeed
int decay = 1, polarisation = 0; //decay on or off?
long double decayoffset=0;

mt_state_t v_mt_state; //mersenne twister state var

// load random number generator and limits from in-file
void LoadMCGenerator(const char *infile){
	time_t mytime;
	tm *monthday;
	timeval daysec;
	ldiv_t divresult; 

	// for random numbers we need a statevar + we need to set an initial seed
	mytime = time(NULL);
	
	monthday = localtime(&mytime);
	monthinmilliseconds = (unsigned long int)(monthday->tm_mday)*24*60*60*1000; // add day in ms
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_hour)*60*60*1000; // add hour in ms
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_min)*60*1000; // add minute in ms 
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(monthday->tm_sec)*1000;  // add second in ms	
	gettimeofday(&daysec, 0);
	divresult = div((long int)(daysec.tv_usec), (long int)(1000));
	monthinmilliseconds = monthinmilliseconds + (unsigned long int)(divresult.quot); // add milliseconds

	mt_set (&v_mt_state, monthinmilliseconds);
	
	Startbed(infile);	
}


//======== read in of the initial values (out of all3inone.in) and the dimensions (out of dimensions.in) ========
void Startbed(const char *infile)
{	int i = 0, ncont;
    char cline[200];    
	string path;

	// creating 'inistream' to all3inone.in
	FILE *inistream = fopen(infile, "r");
	if(inistream == NULL)
	{	cout << "Can't open " << infile << "\n";
		exit(-1);
	}
	
	// looping  over the lines of all3inone.in and reading them in
	for(i = 1; i < 31; i++)
	{	fgets(cline, 200, inistream);
		ncont = 42;
		switch (i)
		{
			case  2:	ncont = sscanf(cline, "%LG %LG ", &nini.EnergieS, &nini.EnergieE);
						break;
			case  3:	ncont = sscanf(cline, "%LG %LG ", &nini.alphas, &nini.alphae);
						nini.alphas *= conv;
						nini.alphae *= conv;
						break;
			case  4:	ncont = sscanf(cline, "%LG %LG ", &nini.gammas, &nini.gammae);
						nini.gammas *= conv;
						nini.gammae *= conv;
						break;
			case 5:		ncont = sscanf(cline, "%LG %LG ", &nini.hmax, &nini.xend);
						break;
			case 7:		ncont = sscanf(cline, "%LG %LG ", &pini.EnergieS, &pini.EnergieE);
						break;
			case 8:	ncont = sscanf(cline, "%LG %LG ", &pini.alphas, &pini.alphae);
						pini.alphas *= conv;
						pini.alphae *= conv;
						break;
			case 9:	ncont = sscanf(cline, "%LG %LG ", &pini.gammas, &pini.gammae);
						pini.gammas *= conv;
						pini.gammae *= conv;
						break;
			case 10:	ncont = sscanf(cline, "%LG %LG ", &pini.hmax, &pini.xend);
						break;
			case 12:	ncont = sscanf(cline, "%LG %LG ", &eini.EnergieS, &eini.EnergieE);
						break;
			case 13:	ncont = sscanf(cline, "%LG %LG ", &eini.alphas, &eini.alphae);
						eini.alphas *= conv;
						eini.alphae *= conv;
						break;
			case 14:	ncont = sscanf(cline, "%LG %LG ", &eini.gammas, &eini.gammae);
						eini.gammas *= conv;
						eini.gammae *= conv;
						break;
			case 15:	ncont = sscanf(cline, "%LG %LG ", &eini.hmax, &eini.xend);
						break;
		}
		if(ncont < 1) Log("an error occourd while reading the %i. item in line %i of all3inone.in", ncont, i);
    }
    
    // closing 'inistream' to all3inone.in
 	fclose(inistream);
	
	// check if lower neutron energy boundary is possible, if not set to possible value
	if((FillingTime == 0) && (CleaningTime == 0) && (RampUpTime == 0) && ((FullFieldTime != 0) || (RampDownTime != 0)))
	{	
		// minimal necessary energy (Emin_n) > maximum input energy (nini.EnergieE)?? EXIT !!
		if(Emin_n*1e9 > nini.EnergieE)
		{	Log("\n\n\nERROR: Emin_n (= %.5LG neV) > nini.EnergieE (= %.5LG neV)\n"
			             "       Check the energy range of the neutrons and / or the B-field configuration!\n\n"
			             "EXIT!!\n", (Emin_n*1e9), nini.EnergieE);
			exit(-1);
		}
	
		nini.EnergieS = max(Emin_n*1e9, nini.EnergieS); // higher energy of the neutron
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
		{	Log("\n\n\nERROR: Two or more initial values are inconsistent.\n"
			             "       Check ALL starting and final values in all3inone.in!\n\n"
			             "EXIT!!\n");
			exit(-1);
		}
	}

	// writing parameters to screen
/*	Log("\nStart parameters:\n"
	       "  Energy (min, max): %.17LG neV/eV/keV, %.17LG neV/eV/keV\n"
	       "  Maximum runtime: %.17LG s\n"
	       "  Alpha (min, max): %.17LG degree, %.17LG degree\n"
	       "  Gamma (min, max): %.17LG degree, %.17LG degree\n"
	       "  Filling time: %LG s\n"
	       "  Cleaning time: %.17LG s\n"
	       "  Ramp up time: %.17LG s\n"
	       "  Full field time: %.17LG s\n"
	       "  RampDownTime: %.17LG s\n"
	       "  B field scaling factor: %.17LG\n"
	       "  E field scaling factor: %.17LG\n",
	       EnergieS, EnergieE,
	       xend,
	       alphas/conv, alphae/conv,
	       gammas/conv,  gammae/conv,
	       FillingTime,
	       CleaningTime,
	       RampUpTime,
	       FullFieldTime,
	       RampDownTime,
	       BFeldSkalGlobal,
	       EFeldSkal);*/
}
//======== end of Startbed ======================================================================================



// return uniformly distributed random number in [min..max]
long double UniformDist(long double min, long double max){
	return mt_get_double(&v_mt_state)*(max - min) + min;
}

// return cosine distributed random number in [min..max] (in rad!)
long double CosDist(long double min, long double max){
	return acos(cos(min) - mt_get_double(&v_mt_state) * (cos(min) - cos(max)));	
}

// return x^2 distributed random number in [min..max]
long double SquareDist(long double min, long double max){
	return powl(mt_get_double(&v_mt_state) * (pow(max,3) - pow(min,3)) + pow(min,3),1.0/3.0);
}

// return sqrt(x) distributed random number in [min..max]
long double SqrtDist(long double min, long double max){
	return powl((powl(max, 1.5) - powl(min, 1.5)) * mt_get_double(&v_mt_state) + powl(min, 1.5), 2.0/3.0);
}

// write isotropically distributed 3D-angles into alpha and gamma
// alpha = angle between x-axis and vector projected onto xy-plane
// gamma = angle between z-axis and vector
void IsotropicDist(long double &alpha, long double &gamma){
	alpha = UniformDist(0,2*pi);
	gamma = CosDist(0,pi);
}

// energy distribution of protons (0 < 'Energie' < 750 eV)			
// proton recoil spectrum from "Diplomarbeit M. Simson"
long double ProtonSpectrum(){
	long double Wahrschl, WktTMP, Energie;
	// dice as long as point is above curve
	do
	{	//cout << "above distribution... dicing on..." << endl;
		Energie = pini.EnergieS + mt_get_double(&v_mt_state) * (pini.EnergieE - pini.EnergieS); // constant distribution between 0 and 800 eV // [eV]
		WktTMP = mt_get_double(&v_mt_state)*2;
		
		long double DeltaM = m_n - m_p;
		long double Xi = m_e / DeltaM;
		long double Sigma = 1 - 2*Energie * m_n / powl(DeltaM, 2) / powl(c_0, 2);
		long double g1 = powl(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
		long double g2 = powl(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
		long double littlea = -0.1017;
		Wahrschl = g1 + littlea * g2;
		//cout << "Proton energy dist: E = " << Energie << " decay rate: " << Wahrschl << " diced value = " << WktTMP << endl;

	}while(WktTMP > Wahrschl);
	return Energie;
}

// energy distribution of electrons (0 < 'Energie' < 782 keV)			
// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
long double ElectronSpectrum(){
	long double Elvert, Energie, WktTMP;
	long double Qvalue = 782e3;
	// dice as long as point is above curve
	do
	{
		Energie = eini.EnergieS*1e3 + mt_get_double(&v_mt_state) * (eini.EnergieE - eini.EnergieS)*1e3; // constant distribution between 0 and 800 keV
		WktTMP = mt_get_double(&v_mt_state) * 0.12;

		Elvert = sqrtl(Energie*1e-6 * Energie*1e-6 + 2*Energie*1e-6 * m_e * c_0 * c_0*1e-6) * pow(Qvalue*1e-6 - Energie*1e-6, 2) * (Energie*1e-6 + m_e * c_0 * c_0*1e-6);			
		cout << "Electron energy dist: E = " << Energie << " decay rate: " << Elvert << " diced value = " << WktTMP << endl;
 
	}while(WktTMP > Elvert);
	return Energie;
}

/*
 * With in the initial conditions random starting values (energie, position) and a random decay time will be diced
 * and crosschecked (if the energy is possible at this starting point)
 */
void MCStartwerte(TParticle *particle)
{	
	initial *ini;
	particle->xstart = 0;

	// square root neutron energy distribution
	particle->NeutEnergie = SqrtDist(nini.EnergieS*1e-9, nini.EnergieE*1e-9);

	if (particle->protneut == NEUTRON){
		ini = &nini;
		particle->Hstart = particle->NeutEnergie;
		// random start time during filling period
		if(FillingTime > 0)
			particle->xstart = FillingTime * (mt_get_double(&v_mt_state));
		// set time span till decay
		if(decay != 0)
			particle->xend = decayoffset - tau * log(mt_get_double(&v_mt_state)); // time span till neutron decay
		else particle->xend = nini.xend;
	}
	else if (particle->protneut == PROTON){
		ini = &pini;
		particle->Hstart = ProtonSpectrum();			
		particle->xend = pini.xend;
						
	}		
	else if (particle->protneut == ELECTRON){
		ini = &eini;
		particle->Hstart = ElectronSpectrum();
		particle->xend = eini.xend;
	}
	
	particle->alphastart = UniformDist(ini->alphas,ini->alphae);
	particle->gammastart = CosDist(ini->gammas,ini->gammae);
	particle->hmax = ini->hmax;
	
	if (polarisation == 4)
	{
		if(UniformDist(0,1) < 0.5)
			particle->polarisation = -1;
		else 
			particle->polarisation = 1;
	}
	else if (polarisation == POLARISATION_GOOD)
		particle->polarisation = -1;
	else if (polarisation == POLARISATION_BAD)
		particle->polarisation = 1;
	else if (polarisation == POLARISATION_NONE)
		particle->polarisation = 0;
	
	/*
	// neutron energy distribution for AbEx@ILL
	long double p1 = -77.6138, p2 = -1.34704, p3 = 0.00739579, p4 = 0.00012494, p5 = -1.88103e-6 , p6 = 8.52798e-9;
	do
	{	cout << "above distribution... dicing on..." << endl;
		x = EnergieS*1e9 + mt_get_double(v_mt_state) * (EnergieE*1e9 - EnergieS*1e9);
		WktTMP = mt_get_double(v_mt_state)*140;
		cout << "AbEx energy dist: E = " << x << " Wkt = " << WktTMP << endl; 
	}while(WktTMP > (-1*(p1 + 2*p2 * x+ 3*p3 * pow(x,2) + 4*p4 * pow(x,3) + 5*p5 * pow(x,4) + 6*p6 * pow(x,5))));

	cout << "below! Take this value E = " << x <<  endl;
	Energie = NeutEnergie = x*1e-9;
	// END AbEx@ILL
	*/
	
	Log("\nDice starting position for E_neutron = %LG neV ... (Please be patient!)\n",particle->NeutEnergie*1e9);
	// repeat dicing for starting point until particle is valid (energy is possible at this specific point) because neutrons
	// are filled in according to the energy spectrum above, the points in the trap that are reachable are determined by that
	long int nroll = 0; // counter for the number of rolled dices 
	long double crit;   // crit = initial energie - potenial energy by gravitation + potential energie by B-field  // [eV]
	while(1)
	{	
		nroll++;
		// starting point
		RandomPointInSourceVolume(particle->rstart, particle->phistart, particle->zstart, particle->alphastart, particle->gammastart);
					
		// check if neutron could possiblly reached this positon by its own energy
		BFeld(particle->rstart, particle->phistart, particle->zstart, particle->xstart);
		crit = particle->NeutEnergie - m_n*gravconst*particle->zstart + particle->polarisation*mu_nSI/ele_e*Bws; // crit = initial energie - potenial energy by gravitation + potential energie by B-field // [eV]
		//printf("Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n", NeutEnergie*1e9, m_n * gravconst * z_n*1e9, (mu_nSI / ele_e) * Bws*1e9, crit*1e9);
		//fprintf(LOGSCR,"Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie * 1e9, m_n * gravconst * z_n*1e9, (mu_nSI / ele_e) * Bws*1e9, crit*1e9);

		if (crit >= 0.0){	
			Log("... %li dice(s) rolled\n", nroll);
			break;
		}
		else if(nroll == 42000000)
		{	
			particle->kennz = KENNZAHL_INITIAL_NOT_FOUND;
			Log("... ABORT: Failed %li times to find a compatible spot!! NO particle will be simulated!!\n\n", nroll);
			return;
		}		
	}
}

//======== Determination of starting values (energie, position, orientation) of the decay products =========================

/*
 * The starting values (energie, position, orientation) of the proton and electron from the decay of the latest simulated 
 * neutron will be set. Initial parameters of the products will be concerned, but due to the Lorentz boost the energy may 
 * exceed the limits.
 * 
 * For more details see the comments below.
 * 
 */
 
void MCZerfallsstartwerte(TParticle *neutron, TParticle *proton, TParticle *electron)
{	
	long double m_nue = 1 / powl(c_0, 2); // [eV/c^2]
	
	// keep polarisation
	proton->polarisation = electron->polarisation = neutron->polarisation;
	proton->xend = pini.xend;
	proton->hmax = pini.hmax;
	electron->xend = eini.xend;
	electron->hmax = eini.hmax;
	
	// restore neutron energy 'NeutEnergie' from decayed neutron
	proton->NeutEnergie = electron->NeutEnergie = neutron->NeutEnergie; // [eV]
	
/*
 * The decay products will start at the point where the neutron decayed.
 * 
 */
 	proton->rstart = electron->rstart = neutron->rend;
 	proton->phistart = electron->phistart = neutron->phiend;
 	proton->zstart = electron->zstart = neutron->zend;
 	proton->xstart = electron->xstart = neutron->xstart + neutron->dt;

/*
 * Reaction(s): 
 *     A   ->  B   +  C
 *     n0  ->  p+  +  R-
 *     R-  ->  e-  +  nue
 * 
 * Procedure:
 * (1) Dicing the energy of the decay proton according to the characteristic spectrum in the rest frame of the neutron and
 *     calculating the proton momentum via the energy momentum relation. 
 * (2) Dicing isotrope orientation of the decay proton
 * (3) Calculating 4-momentum of R- via 4-momentum conservation
 * (4) Getting fixed electron energy from two body decay of R-
 * (5) Dicing isotrope electron orientation in the rest frame of R-
 * (6) Lorentz boosting electron 4-momentum into moving frame of R-
 * (7) Calculating neutrino 4-momentum via 4-momentum conservation
 * (8) Boosting all 4-momentums into moving neutron frame
 * 
 * Cross-check:
 * (9) Check if neutrino 4-momentum has the correct invariant mass (4-momentum square)
 * 
 */
 	long double p[4], e[4]; // 4-momentums
 	long double pabs, beta1, beta2, delta1, delta2; // 3-momentums (abs+directions)
 	
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
	// isotropic emission characteristics (in the rest frame of the neutron with the z''-axis as its track)
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
	long double vend_n = sqrt(2*neutron->Hend/m_n);
	beta[0] = vend_n/c_0*sin(neutron->gammaend)*cos(neutron->alphaend + neutron->phiend);
	beta[1] = vend_n/c_0*sin(neutron->gammaend)*sin(neutron->alphaend + neutron->phiend);
	beta[2] = vend_n/c_0*cos(neutron->gammaend);
	BOOST(beta,e);
	BOOST(beta,p);
	BOOST(beta,nue);
	
	proton->Hstart = p[0]*c_0 - m_p*c_0*c_0;
	proton->gammastart = acos(p[3]/sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
	proton->alphastart = fmod(atan2(p[2],p[1]) - proton->phistart,2*pi);
	
	electron->Hstart = e[0]*c_0 - m_e*c_0*c_0;
	electron->gammastart = acos(e[3]/sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]));
	electron->alphastart = fmod(atan2(e[2],e[1]) - electron->phistart,2*pi);
	

//-------- Step 9 ----------------------------------------------------------------------------------------------------------
	// use neutrino mass for check
	neutron->decayerror = proton->decayerror = electron->decayerror = m_nue*c_0*c_0 - sqrt(nue[0]*nue[0] - nue[1]*nue[1] - nue[2]*nue[2] - nue[3]*nue[3])*c_0;
	Log("\n   +++ decay error : %LG eV +++\n", neutron->decayerror);


//-------- Finished --------------------------------------------------------------------------------------------------------
		
	return;
}
//======== end of MCZerfallsstartwerte =====================================================================================

//======== Lorentz boost of four-vector p into frame moving in arbitrary direction with v/c = beta ======================================================
void BOOST(long double beta[3], long double p[4]){
   //Boost this Lorentz vector (copy&paste from ROOT)
   long double b2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
   long double gamma = 1.0 / sqrt(1.0 - b2);
   long double bp = beta[0]*p[1] + beta[1]*p[2] + beta[2]*p[3];
   long double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

   p[1] = (p[1] + gamma2*bp*beta[0] + gamma*beta[0]*p[0]);
   p[2] = (p[2] + gamma2*bp*beta[1] + gamma*beta[1]*p[0]);
   p[3] = (p[3] + gamma2*bp*beta[2] + gamma*beta[2]*p[0]);
   p[0] = (gamma*(p[0] + bp));}
//======== end of BOOST ====================================================================================================

