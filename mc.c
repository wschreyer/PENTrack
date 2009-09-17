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


	Log("\nDice starting position ... (Please be patient!)\n");
	// repeat dicing for starting point until particle is valid (energy is possible at this specific point) because neutrons
	// are filled in according to the energy spectrum above, the points in the trap that are reachable are determined by that
	nroll = 0;
	noparticle = false;
	while(1)
	{	nroll++;
		// starting point
		RandomPointInSourceVolume(r_n,phi_n,z_n);
					
		// check if neutron could possiblly reached this positon by its own energy
		BFeld(r_n, phi_n * conv, z_n, 0.0);
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

	alpha = alphas + (mt_get_double(v_mt_state)) * (alphae - alphas); // constant angular distribution
	gammaa = acosl(cosl(gammas*conv) - mt_get_double(v_mt_state) * (cosl(gammas * conv) - cosl(gammae * conv))) / conv; // isotropic emission characteristics

	// for checking the dependence of spin flip on inner rod current
	if((protneut == NEUTRON) && (DiceRodField == 1))
		RodFieldMultiplicator = (mt_get_double(v_mt_state));

	// random start time from 0 to 2e-5 for protons
	if(protneut == (PROTON || ELECTRONS))
		xstart = (mt_get_double(v_mt_state))*2e-5;

	// random start time during filling period
	if((FillingTime > 0) && (protneut == NEUTRON))
		xstart = FillingTime * (mt_get_double(v_mt_state));

	// set time span till decay
	if((decay.on) && (protneut == NEUTRON))
	{	xend = - tau * log(mt_get_double(v_mt_state)); // time span till neutron decay
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
 * (2) Dicing the orientation of the decay proton according to the phase space and the neutron frame (and retriving the
 *     orientation of the virtual state).
 * (3) Analog 1a for electron.
 * (4) Analog 1b for electron (and neutrino) in respect to the orientation of the virtual state.
 * (5) Calculating the neutrino energy and momentum via center of mass energy.
 * (6) Lorentz boosting all 4-momentums in the moving neutron frame.
 * (7) Rotating all 3-momentums in respect to the orientaion of the neutrons 3-momentum.
 * 
 * Cross-check:
 * (8) Cross-check via centre of mass energy and determination of the decay error. The 'decay.error' shows the square root
 *     of the normed difference of the centre of mass energy, the sign indicates if energy has been (-) lost or (+) gained.   
 * 
 */
 	long double p[5], e[5], nue[5]; // 4-momentums [0..3] and norm [4]
 	long double beta1, beta2, delta1, delta2;
 	long double epx, epy, epz, pxy;
 	long double s, sign, ds; // center of mass energy, sign, difference in s
 	
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
		
		decay.PEnergie = Energie; // [eV]
/*
 * E   = E_kin + m*c^2
 * p^2 = (E/c)^2 - (m*c)^2
 *
 */
 		// energy of proton
		p[0] = Energie / c_0 + m_p * c_0;
		// momentum norm of proton
		p[4] = sqrtl(powl(p[0], 2) - powl(m_p * c_0, 2));

//-------- Step 2 ----------------------------------------------------------------------------------------------------------
		// isotropic emission characteristics (in the rest frame of the neutron with the z''-axis as its track)
		beta1 = 360 * (mt_get_double(v_mt_state));
		delta1 = acosl(1 - 2 * (mt_get_double(v_mt_state))) / conv;
		// spharical coordinates
		epx = cosl(delta1 * conv) * cosl(beta1 * conv);
		epy = cosl(delta1 * conv) * sinl(beta1 * conv);
		epz = sinl(delta1 * conv);
		// 3-momentum of the proton
		p[1] = p[4] * epx;
		p[2] = p[4] * epy;
		p[3] = p[4] * epz;
		
//-------- Step 3 ----------------------------------------------------------------------------------------------------------
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
/*
 * E   = E_kin + m*c^2
 * p^2 = (E/c)^2 - (m*c)^2
 *
 */
 		// energy of electron
		e[0] = Energie / c_0 + m_e * c_0;
		// momentum norm of electron
		e[4] = sqrtl(powl(e[0], 2) - powl(m_e * c_0, 2));

//-------- Step 4 ----------------------------------------------------------------------------------------------------------
		// isotropic emission characteristics (in the rest frame of the neutron with the z''-axis as its track)
		beta2 = 360 * (mt_get_double(v_mt_state));
		delta2 = acosl(1 - 2 * (mt_get_double(v_mt_state))) / conv;
		// spharical coordinates
		epx = cosl(delta2 * conv) * cosl(beta2 * conv);
		epy = cosl(delta2 * conv) * sinl(beta2 * conv);
		epz = sinl(delta2 * conv);
		// rotation in respect to the orientation of the virtual state
		zyROTROT((180 + beta1), (180 - delta1), &epx, &epy, &epz);		
		// 3-momentum of the electron
		e[1] = e[4] * epx;
		e[2] = e[4] * epy;
		e[3] = e[4] * epz;

//-------- Step 5 ----------------------------------------------------------------------------------------------------------
		// momentum norm of neutrino
/*
 * s = (m_n*c^2)^2 = (SUM_i (E_i, p_ix..z*c))^2   ,   i = {p, e, nue}
 * 
 * C = E*E_nue - P*p_nue 
 * 
 */
		long double M2 = powl(m_n, 2) - powl(m_p, 2) - powl(m_e, 2) - powl(m_nue, 2);
		long double EP = 2 * p[0] * e[0] - 2 * (p[1] * e[1] + p[2] * e[2] + p[3] * e[3]);
		long double C = M2 * powl(c_0, 2) - EP; 
		long double E = 2 * (p[0] + e[0]);
		long double P = 2 * ((p[1] * epx + p[2] * epy + p[3] * epz) + e[4]);
		nue[4] = (C * P + E * sqrtl(powl(C, 2) + (powl(P, 2) - powl(E, 2)) * powl(m_nue * c_0, 2))) / (powl(E, 2) - powl(P, 2));
		// 4-momentum of neutrino
		nue[0] = sqrtl(powl(nue[4], 2) + powl(m_nue * c_0, 2));
		nue[1] = nue[4] * epx;
		nue[2] = nue[4] * epy;
		nue[3] = nue[4] * epz;

//-------- Step 6 ----------------------------------------------------------------------------------------------------------
		// Lorentz boost in z''-direction with neutron velocity
		
		// for proton
		zBOOST(decay.Nv / c_0, &p[0], &p[3]);
		p[4] = sqrtl(powl(p[1], 2) + powl(p[2], 2) + powl(p[3], 2));
		
		// for electron
		zBOOST(decay.Nv / c_0, &e[0], &e[3]);
		e[4] = sqrtl(powl(e[1], 2) + powl(e[2], 2) + powl(e[3], 2));
		
		// for neutrino
		zBOOST(decay.Nv / c_0, &nue[0], &nue[3]);
		nue[4] = sqrtl(powl(nue[1], 2) + powl(nue[2], 2) + powl(nue[3], 2));

//-------- Step 7 ----------------------------------------------------------------------------------------------------------
		// rotation around z''-axis with (alpha_n) and rotation around y'-axis with (gamma_n)  
		
		// for proton
		epx = p[1] / p[4];
		epy = p[2] / p[4];
		epz = p[3] / p[4];
		zyROTROT(decay.Nalpha, decay.Ngamma, &epx, &epy, &epz);
		// determination of the protons starting angles !!
		pxy = sqrt(powl(epx, 2) + powl(epy, 2)); 
		decay.Pgamma = GetAngle(pxy, epz);
		decay.Palpha = GetAngle(epx / pxy, epy / pxy);
		// 3-momentum of proton
		p[1] = p[4] * epx;
		p[2] = p[4] * epy;
		p[3] = p[4] * epz;
		
		// for electron
		epx = e[1] / e[4];
		epy = e[2] / e[4];
		epz = e[3] / e[4];
		zyROTROT(decay.Nalpha, decay.Ngamma, &epx, &epy, &epz);
		// determination of the protons starting angles !!
		pxy = sqrt(powl(epx, 2) + powl(epy, 2)); 
		decay.Egamma = GetAngle(pxy, epz);
		decay.Ealpha = GetAngle(epx / pxy, epy / pxy);
		// 3-momentum of electron
		e[1] = e[4] * epx;
		e[2] = e[4] * epy;
		e[3] = e[4] * epz;
		
		// for neutrino
		epx = nue[1] / nue[4];
		epy = nue[2] / nue[4];
		epz = nue[3] / nue[4];
		zyROTROT(decay.Nalpha, decay.Ngamma, &epx, &epy, &epz);
		// 3-momentum of neutrino
		nue[1] = nue[4] * epx;
		nue[2] = nue[4] * epy;
		nue[3] = nue[4] * epz;

//-------- Step 8 ----------------------------------------------------------------------------------------------------------
		// cross-check via s
/*
 * s = (m_n*c^2)^2 = (SUM_i (E_i, p_ix..z*c))^2   ,   i = {p, e, nue}
 * 
 */
		s = 0;
		sign = 1;
		for(int i = 0; i < 4; i++)
		{	s = s + sign * powl((p[i] + e[i] + nue[i]) * c_0, 2);
			sign = -1;
		}
		ds = s - powl(m_n * powl(c_0, 2),2);
		if(ds != 0)
		{	sign = ds / fabs(ds);
			decay.error = sign * sqrtl(fabs(ds));
		}
		else
		{	decay.error = 0;
		}				
		
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

//======== Getting the angle out of the cosine and sine ====================================================================
extern long double GetAngle(long double cos, long double sin)
{	long double angle = 0;	
	if((cos >= 0) && (sin >=0))				// quadrant I
	{	angle = acosl(cos) / conv;
	}
	else if((cos < 0) && (sin >= 0))		// quadrant II
	{	angle = acosl(cos) / conv;
	}
	else if((cos < 0) && (sin < 0))			// quadrant III
	{	angle = 360 - (acosl(cos) / conv);
	}
	else if((cos >= 0) && (sin < 0))		// quadrant IV
	{	angle = (-1) * (acosl(cos) / conv);		
	}
	return angle;
}
//======== end of GetAngle =================================================================================================

//======== Lorentz boost in z''-direction with a given beta ================================================================ 
extern void zBOOST(long double BOOSTbeta, long double* Dp0, long double* Dp3)
{	long double BOOSTgamma = powl(1 - powl(BOOSTbeta, 2), -0.5);
	long double D2p0 = *Dp0;
	long double D2p3 = *Dp3;
	
	*Dp0 = BOOSTgamma * D2p0 - BOOSTbeta * BOOSTgamma * D2p3;
	*Dp3 = BOOSTgamma * D2p3 - BOOSTbeta * BOOSTgamma * D2p0;
	return;
}
//======== end of zBOOST ===================================================================================================

//======== Rotaion around z''-axis with (+beta) and rotation around y'-axis with (+delta) ==================================
extern void zyROTROT(long double zROTbeta, long double yROTdelta, long double* Dp1, long double* Dp2, long double* Dp3)
{	// Note: The minus sign for the rotation angle is already implemented in the (sign of the sines in the) calculation
	long double cb = cosl(zROTbeta * conv);
	long double sb = sinl(zROTbeta * conv);
	long double cd = cosl(yROTdelta * conv);
	long double sd = sinl(yROTdelta * conv);
	long double D0p1 = *Dp1;
	long double D0p2 = *Dp2;
	long double D0p3 = *Dp3;
	
	*Dp1 = 	      cd * cb * D0p1 - cd * sb * D0p2 + sd * D0p3;
	*Dp2 =             sb * D0p1 +      cb * D0p2;
	*Dp3 = (-1) * sd * cb * D0p1 + sd * sb * D0p2 + cd * D0p3;
	return;
}
//======== end of zyROTROT =================================================================================================
