#include "main.h"
//#include "vars.h"

bool noparticle = false;

/*
 * With in the initial conditions random starting values (energie, position) and a random decay time will be diced
 * and crosschecked (if the energy is possible at this starting point)
 */

void MCStartwerte(long double delx)
{	
	if((decay.on == 2) && ((protneut == PROTON) || (protneut == ELECTRONS)))
	{	MCZerfallsstartwerte();
		return;
	}
	
	long int nroll = 0; // counter for the number of rolled dices 
	long double crit;   // crit = initial energie - potenial energy by gravitation + potential energie by B-field  // [eV]
	long double WktTMP;

	// examine correlation of fields with proton collection 
	 // BFeldSkalGlobal = mt_get_double(v_mt_state);
	 // EFeldSkal = mt_get_double(v_mt_state);
	 // EFeldSkal = 0;
	
	// dicing of spin orientation towards magnetic field
	if (polarisationsave == 4)
	{
		long double poldice = mt_get_double(v_mt_state);
		if(poldice <= 0.5)
		{
			polarisation = 1;
			hfs = -1;
			mu_n = hfs * mu_nSI / ele_e;
			mumB = mu_n / M; // [c^2/T]
		}
		else if(poldice > 0.5) 
		{
			polarisation = 2;
			hfs = 1;
			mu_n = hfs * mu_nSI / ele_e;		
			mumB = mu_n / M; // [c^2/T]
		}
	}

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


	// square root neutron energy distribution
	if ((protneut != PROTON)&&(protneut != ELECTRONS))
	{	//Energie = powl(mt_get_double(v_mt_state) * (EnergieE * EnergieE - EnergieS * EnergieS) + EnergieS * EnergieS, 0.5); // squared weighted energy
		Energie = powl((powl(nini.EnergieE*1e-9, 1.5) - powl(nini.EnergieS*1e-9, 1.5)) * mt_get_double(v_mt_state) + powl(nini.EnergieS*1e-9, 1.5), 2.0/3.0); // [eV]
		NeutEnergie = Energie; // [eV]
	}		

	// proton energy distribution
	else if (protneut == PROTON)
	{   // energy distribution of protons (0 < 'Energie' < 750 eV)			
		// proton recoil spectrum from "Diplomarbeit M. Simson"
		long double Wahrschl;
		// dice as long as point is above curve
		do
		{	//cout << "above distribution... dicing on..." << endl;
			Energie = pini.EnergieS + mt_get_double(v_mt_state) * (pini.EnergieE - pini.EnergieS); // constant distribution between 0 and 800 eV // [eV]
			WktTMP = mt_get_double(v_mt_state)*2;
			
			long double DeltaM = m_n - m_p;
			long double Xi = m_e / DeltaM;
			long double Sigma = 1 - 2*Energie * m_n / powl(DeltaM, 2) / powl(c_0, 2);
			long double g1 = powl(((Sigma - powl(Xi, 2)) / Sigma), 2) * sqrtl(1 - Sigma) * (4*(1 + powl(Xi, 2) / Sigma) - 4/3*((Sigma - powl(Xi, 2)) / Sigma * (1 - Sigma)));
			long double g2 = powl(((Sigma - powl(Xi, 2)) / Sigma), 2) * sqrtl(1 - Sigma) * (4*(1 + powl(Xi, 2) / Sigma - 2 * Sigma) - 4/3*((Sigma - powl(Xi, 2)) / Sigma * (1 - Sigma)));
			long double littlea = -0.1017;
			Wahrschl = g1 + littlea * g2;
			//cout << "Proton energy dist: E = " << Energie << " decay rate: " << Wahrschl << " diced value = " << WktTMP << endl;

		}while(WktTMP > Wahrschl);			

		// temporary for energy dependence test of proton extraction efficiency 
		//NeutEnergie =  15.0e-9 + mt_get_double(v_mt_state) * (110.0e-9 - 15.0e-9); // equal weighted neutron energy		
		// END temporary for energy dependence test of proton extraction efficiency
	
																					// add 1neV to avoid infinite initial condition dicing
		NeutEnergie  = powl((powl(nini.EnergieE*1e-9, 1.5) - powl(nini.EnergieS*1e-9 + 1e-9, 1.5)) * mt_get_double(v_mt_state) + powl(nini.EnergieS*1e-9 + 1e-9, 1.5), 2.0/3.0); // [eV]
	}
	
	// electron energy distribution
	else if (protneut == ELECTRONS)
	{	// energy distribution of electrons (0 < 'Energie' < 782 keV)			
		// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
		long double Elvert;
		long double Qvalue = 782e3;
		// dice as long as point is above curve
		do
		{
			Energie = eini.EnergieS*1e3 + mt_get_double(v_mt_state) * (eini.EnergieE - eini.EnergieS)*1e3; // constant distribution between 0 and 800 keV
			WktTMP = mt_get_double(v_mt_state) * 0.12;

			Elvert = sqrtl(Energie*1e-6 * Energie*1e-6 + 2*Energie*1e-6 * m_e * c_0 * c_0*1e-6) * powl(Qvalue*1e-6 - Energie*1e-6, 2) * (Energie*1e-6 + m_e * c_0 * c_0*1e-6);			
			cout << "Electron energy dist: E = " << Energie << " decay rate: " << Elvert << " diced value = " << WktTMP << endl;
 
		}while(WktTMP > Elvert);

																					// add 1neV to avoid infinite initial condition dicing
		NeutEnergie  = powl( (powl(nini.EnergieE*1e-9, 1.5) - powl(nini.EnergieS*1e-9 + 1e-9, 1.5)) * mt_get_double(v_mt_state) + powl(nini.EnergieS*1e-9 + 1e-9, 1.5), 2.0/3.0); // [eV]			
	}


	Log("\nDice starting position for E_neutron = %LG neV ... (Please be patient!)\n",NeutEnergie*1e9);
	// repeat dicing for starting point until particle is valid (energy is possible at this specific point) because neutrons
	// are filled in according to the energy spectrum above, the points in the trap that are reachable are determined by that
	nroll = 0;
	noparticle = false;
	while(1)
	{	nroll++;
		// starting point
		RandomPointInSourceVolume(r_n,phi_n,z_n,alpha,gammaa);
					
		// check if neutron could possiblly reached this positon by its own energy
		BFeld(r_n, phi_n, z_n, 0.0);
		crit = NeutEnergie - m_n * gravconst * z_n + mu_n * Bws; // crit = initial energie - potenial energy by gravitation + potential energie by B-field // [eV]
		//printf("Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n", NeutEnergie*1e9, m_n * gravconst * z_n*1e9, (mu_nSI / ele_e) * Bws*1e9, crit*1e9);
		//fprintf(LOGSCR,"Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie * 1e9, m_n * gravconst * z_n*1e9, (mu_nSI / ele_e) * Bws*1e9, crit*1e9);

		if (crit >= 0.0){	
			Log("... %li dice(s) rolled\n", nroll);
			break;
		}
		if(nroll == 42000000)
		{	kennz = KENNZAHL_INITIAL_NOT_FOUND;
			noparticle = true;
			Log("... ABORT: Failed %li times to find a compatible spot!! NO particle will be simulated!!\n\n", nroll);
			return;
		}		
	}

	// for checking the dependence of spin flip on inner rod current
	if((protneut == NEUTRON) && (DiceRodField == 1))
		RodFieldMultiplicator = (mt_get_double(v_mt_state));

	xstart = 0;
	// random start time from 0 to 2e-5 for protons
	if((protneut == PROTON) || (protneut == ELECTRONS))
		xstart = (mt_get_double(v_mt_state))*2e-5;

	// random start time during filling period
	if((FillingTime > 0) && (protneut == NEUTRON))
		xstart = FillingTime * (mt_get_double(v_mt_state));

	// set time span till decay
	if((decay.on) && (protneut == NEUTRON))
	{	
		xend = decayoffset - tau * log(mt_get_double(v_mt_state)); // time span till neutron decay
		//if(decay.on == 2) xend = 0.1; // ONLY FOR TESTING !! DELETE THIS LINE !!
	}
	
	return ;
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
 
void MCZerfallsstartwerte()
{	long double m_nue = 1 / powl(c_0, 2); // [eV/c^2]
	long double WktTMP;
	
	// keep polarisation
	polarisation = decay.Npolarisation;

	// keep scaling factor for the global E- and B-field	// !?
	//EFeldSkal = EFeldSkalSave;							// !?
	//BFeldSkalGlobal = BFeldSkalGlobalSave;				// !?

	// keep 'RodFieldMultiplicator' (no redicing)
	
	// restore neutron energy 'NeutEnergie' from decayed neutron
	NeutEnergie = decay.NH*1e-9; // [eV]
	
/*
 * The decay products will start at the point where the neutron decayed.
 * 
 */
 	z_n = decay.Nz;
	r_n = decay.Nr;
	phi_n = decay.Nphi;

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
 	
	if(protneut == PROTON)
	{
//-------- Step 1 ----------------------------------------------------------------------------------------------------------
		// energy distribution of protons (0 < 'Energie' < 800 eV)			
		// proton recoil spectrum from "Diplomarbeit M. Simson"
		long double Wahrschl;
		// dice as long as point is above curve
		do
		{	Energie = pini.EnergieS + mt_get_double(v_mt_state) * (pini.EnergieE - pini.EnergieS); // constant distribution between 0 and 800 eV // [eV]
			WktTMP = mt_get_double(v_mt_state)*2;			
			long double DeltaM = m_n - m_p;
			long double Xi = m_e / DeltaM;
			long double Sigma = 1 - 2 * Energie * m_n / powl(DeltaM, 2) / powl(c_0, 2);
			long double g1 = powl(((Sigma - powl(Xi, 2)) / Sigma), 2) * sqrtl(1 - Sigma) * (4*(1 + powl(Xi, 2) / Sigma) - 4/3*((Sigma - powl(Xi, 2)) / Sigma * (1 - Sigma)));
			long double g2 = powl(((Sigma - powl(Xi, 2)) / Sigma), 2) * sqrtl(1 - Sigma) * (4*(1 + powl(Xi, 2) / Sigma - 2 * Sigma) - 4/3*((Sigma - powl(Xi, 2)) / Sigma * (1 - Sigma)));
			long double littlea = -0.1017;
			Wahrschl = g1 + littlea * g2;
			//cout << "Proton energy dist: E = " << Energie << " decay rate: " << Wahrschl << " diced value = " << WktTMP << endl;
		}while(WktTMP > Wahrschl);
		
/*
 * E   = E_kin + m*c^2
 * p^2 = (E/c)^2 - (m*c)^2
 *
 */
 		// energy of proton
		p[0] = Energie / c_0 + m_p * c_0;
		// momentum norm of proton
		pabs = sqrt(p[0]*p[0] - m_p*m_p*c_0*c_0);

//-------- Step 2 ----------------------------------------------------------------------------------------------------------
		// isotropic emission characteristics (in the rest frame of the neutron with the z''-axis as its track)
		beta1 = 2*pi * (mt_get_double(v_mt_state));
		delta1 = acosl(1 - 2 * (mt_get_double(v_mt_state)));
		// 3-momentum of the proton
		p[1] = pabs * sin(delta1) * cos(beta1);
		p[2] = pabs * sin(delta1) * sin(beta1);
		p[3] = pabs * cos(delta1);
		
//-------- Step 3 ----------------------------------------------------------------------------------------------------------
		// calculate intermediate virtual state via 4-momentum conservation
 		long double virt[4] = {m_n*c_0 - p[0], -p[1], -p[2], -p[3]}; // 4-momentum of intermediate virtual state (n - p)
 		long double m2_virt = virt[0]*virt[0] - virt[1]*virt[1] - virt[2]*virt[2] - virt[3]*virt[3]; // squared mass of virtual state 		
		

//-------- Step 4 ----------------------------------------------------------------------------------------------------------
/*
		// energy distribution of electrons (0 < 'Energie' < 782 keV)			
		// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
		long double Elvert;
		long double Qvalue = (m_n - m_p - m_e - m_nue) * powl(c_0, 2);
		//                 = 782e3; // [eV]
		// dice as long as point is above curve
		do
		{	Energie = eini.EnergieS*1e3 + mt_get_double(v_mt_state) * (eini.EnergieE - eini.EnergieS)*1e3; // constant distribution between 0 and 800 keV
			WktTMP = mt_get_double(v_mt_state) * 0.12;
			Elvert = sqrtl(Energie*1e-6 * Energie*1e-6 + 2*Energie*1e-6 * m_e * c_0 * c_0*1e-6) * powl(Qvalue*1e-6 - Energie*1e-6, 2) * (Energie*1e-6 + m_e * c_0 * c_0*1e-6);			
			//cout << "Electron energy dist: E = " << Energie << " decay rate: " << Elvert << " diced value = " << WktTMP << endl; 
		}while(WktTMP > Elvert);
		
		decay.EEnergie = Energie;
*/

 		// energy of electron from two-body-decay of virtual state
		e[0] = (m2_virt + m_e*m_e*c_0*c_0 - m_nue*m_nue*c_0*c_0)/2/sqrt(m2_virt);
		pabs = sqrt(e[0]*e[0] - m_e*m_e*c_0*c_0);

//-------- Step 5 ----------------------------------------------------------------------------------------------------------
		// isotropic emission characteristics (in the rest frame of the virtual state)
		beta2 = 2*pi * (mt_get_double(v_mt_state));
		delta2 = acosl(1 - 2 * (mt_get_double(v_mt_state)));
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
		beta[0] = decay.Nv/c_0*sin(decay.Ngamma)*cos(decay.Nalpha + decay.Nphi);
		beta[1] = decay.Nv/c_0*sin(decay.Ngamma)*sin(decay.Nalpha + decay.Nphi);
		beta[2] = decay.Nv/c_0*cos(decay.Ngamma);
		BOOST(beta,e);
		BOOST(beta,p);
		BOOST(beta,nue);
		
		decay.PEnergie = p[0]*c_0 - m_p*c_0*c_0;
		decay.Pgamma = acos(p[3]/sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));
		decay.Palpha = fmod(atan2(p[2],p[1]) - decay.Nphi,2*pi);
		
		decay.EEnergie = e[0]*c_0 - m_e*c_0*c_0;
		decay.Egamma = acos(e[3]/sqrt(e[1]*e[1] + e[2]*e[2] + e[3]*e[3]));
		decay.Ealpha = fmod(atan2(e[2],e[1]) - decay.Nphi,2*pi);
		

//-------- Step 9 ----------------------------------------------------------------------------------------------------------
		// use neutrino mass for check
		decay.error = m_nue*c_0*c_0 - sqrt(nue[0]*nue[0] - nue[1]*nue[1] - nue[2]*nue[2] - nue[3]*nue[3])*c_0;
		Log("\n   +++ decay error : %LG eV +++\n", decay.error);


//-------- Finished --------------------------------------------------------------------------------------------------------
		
		Energie = decay.PEnergie;
		alpha = decay.Palpha;
		gammaa = decay.Pgamma;
	}
	else if(protneut == ELECTRONS)
	{	
		Energie = decay.EEnergie;
		alpha = decay.Ealpha;
		gammaa = decay.Egamma;
	}

	// starting time = time of neutron decay
	xstart = decay.Nx;
	
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

