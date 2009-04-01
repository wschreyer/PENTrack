
#include "main.h"
//#include "vars.h"

/*
Aus den Grenzen f�r die Startbedingungen wird ein zuf�lliger Startwert und eine zuf�llige
Zerfallszeit f�r das Teilchen berechnet
*/

void MCStartwerte(long double delx)
{	long double nEnerg; // dummy for energie of neutron
	long double crit;
	long double WktTMP;
	// examine correlation of fields with proton collection 
	 // BFeldSkalGlobal = mt_get_double(v_mt_state);
	 // EFeldSkal = mt_get_double(v_mt_state);
	 // EFeldSkal = 0;
	
	// dicing of spin orientation towards magnetic field
	if (polarisationsave==4)
	{
		long double poldice = mt_get_double(v_mt_state);
		if(poldice<=0.5)
		{
			polarisation=1;
			hfs = -1;
			mu_n= hfs * mu_nSI / ele_e;
			mumB= mu_n/M;  // [c^2/T]
		}
		else if(poldice>0.5) 
		{
			polarisation=2;
			hfs = 1;
			mu_n=hfs * mu_nSI / ele_e;		
			mumB= mu_n/M;  // [c^2/T]
		}
	}
	
	if(BFeldSkalGlobalSave==-1)
	{
		BFeldSkalGlobal=mt_get_double(v_mt_state);
		cout << "BFeldSkalGlobal:" << BFeldSkalGlobal << endl;
	}
	
	if(EFeldSkalSave==-1)
	{
		EFeldSkal=mt_get_double(v_mt_state);
		cout << "EFeldSkal:" << EFeldSkal << endl;
	}
	
	
	/*
	// neutron energy distribution for AbEx @ ILL
		
		long double p1=-77.6138,p2=-1.34704,p3=0.00739579,p4=0.00012494,p5=-1.88103e-6,p6=8.52798e-9;
		do
		{
			cout << "above distribution... dicing on..." << endl;
			x = EnergieS*1e9 + mt_get_double(v_mt_state)*(EnergieE*1e9-EnergieS*1e9);			
			WktTMP = mt_get_double(v_mt_state)*140;
					cout << "AbEx energy dist: E=" << x << " Wkt=" << WktTMP << endl; 
		}while(WktTMP>(-1*(p1+2*p2*x+3*p3*pow(x,2)+4*p4*pow(x,3)+5*p5*pow(x,4)+6*p6*pow(x,5))));
		cout << "below! Take this value E=" << x <<  endl;
		Energie = NeutEnergie = x/1e9;
		// END AbEX@ILL*/
		
		// check if lower neutron energy boundary is possible, if not set to possible value 
		if(z_ns>=0 && Emin_n>EnergieE)
			EnergieS=Emin_n+1e-9;
		
		// square root neutron energy distribution
		if ((protneut != PROTON)&&(protneut != ELECTRONS))
		{
			//Energie = powl(mt_get_double(v_mt_state) * (EnergieE*EnergieE - EnergieS*EnergieS) + EnergieS*EnergieS,0.5);    // quadratisch gewichtete Energie
			Energie = powl( (powl(EnergieE,1.5)-powl(EnergieS,1.5)) *mt_get_double(v_mt_state) + powl(EnergieS,1.5) ,2.0/3.0);
			NeutEnergie = Energie;		
			//cout << "Die gw�rfelte Energie ist :" << Energie << endl;
			//sleep(5);
		}		
		
		// proton energy distribution
		else if (protneut == PROTON)
		{            // energy distribution of protons 0 < Energie < 750eV			
			// proton recoil spectrum aus Diplomarbeit M. Simson
			long double Wahrschl;
			// dice as long as point is above curve
			do
			{
				//cout << "above distribution... dicing on..." << endl;
				Energie = EnergieS + mt_get_double(v_mt_state)*(EnergieE-EnergieS);	 // constant distribution betwwen 0 and 800
				WktTMP = mt_get_double(v_mt_state)*2;
				
				long double DeltaM=m_n-m_p;
				long double Xi=m_e/DeltaM;
				long double Sigma=1-2*Energie*m_n/powl(DeltaM,2)/powl(c_0,2);
				long double g1=powl(((Sigma-powl(Xi,2))/Sigma),2)*sqrtl(1-Sigma)*(4*(1+powl(Xi,2)/Sigma)-4/3*((Sigma-powl(Xi,2))/Sigma*(1-Sigma)));
				long double g2=powl(((Sigma-powl(Xi,2))/Sigma),2)*sqrtl(1-Sigma)*(4*(1+powl(Xi,2)/Sigma-2*Sigma)-4/3*((Sigma-powl(Xi,2))/Sigma*(1-Sigma)));
				long double littlea=-0.1017;
				Wahrschl=g1+littlea*g2;
				//cout << "Proton energy dist: E= " << Energie << " decay rate: " << Wahrschl << " diced value = " << WktTMP << endl; 
			
			}while(WktTMP>Wahrschl);			
			
			// temporary for energy dependence test of proton extraction efficiency 
			//NeutEnergie =  15e-9+mt_get_double(v_mt_state)*(110.0e-9-15e-9);    // gleich gewichtete Neutronen-Energie
			// temporary for energy dependence test of proton extraction efficiency END
			
			nEnerg = max(Emin_n, (nini.EnergieS * 1e-9));	
			NeutEnergie  = powl( (powl(nini.EnergieE*1e-9,1.5)-powl(nEnerg+1e-9,1.5)) *mt_get_double(v_mt_state) + powl(nEnerg+1e-9,1.5) ,2.0/3.0);
						
			
		}
		// electron
		else if (protneut == ELECTRONS)
		{            // energy distribution of electrons 0 < Energie < 782 keV			
			// from http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html
			long double Elvert;
			long double Qvalue = 782e3;
			// dice as long as point is above curve
			do
			{				
				Energie = EnergieS*1000 + mt_get_double(v_mt_state)*(EnergieE-EnergieS)*1000;	 // constant distribution between 0 and 800 keV
				WktTMP = mt_get_double(v_mt_state)*0.12;					
				
				Elvert = sqrtl(Energie/1e6*Energie/1e6+2*Energie/1e6*m_e*c_0*c_0/1e6) * powl(Qvalue/1e6-Energie/1e6,2) * (Energie/1e6+m_e*c_0*c_0/1e6);			
				cout << "Proton energy dist: E= " << Energie << " decay rate: " << Elvert << " diced value = " << WktTMP << endl; 
			}while(WktTMP>Elvert);			
			
			nEnerg = max(Emin_n, (nini.EnergieS * 1e-9));
			NeutEnergie  = powl( (powl(nini.EnergieE*1e-9,1.5)-powl(nEnerg+1e-9,1.5)) *mt_get_double(v_mt_state) + powl(nEnerg+1e-9,1.5) ,2.0/3.0);
			
		}
			
	
	// repeat dicing for starting point until particle is valid (energy is possible at this specific point)
	// because neutrons are filled in according to the energy spectrum above, the points in the trap that are reachable are determined by that
	while(1)
	{		
		
		// starting point
		r_n = powl(mt_get_double(v_mt_state)*(r_ne*r_ne - r_ns*r_ns) + r_ns*r_ns,0.5); // Wichtung da Volumenelement und damit Whk nach aussen mit r^2 steigt
		phi_n = phis + mt_get_double(v_mt_state) * (phie-phis);   // dice starting coordinate phi
		/* // z distribution according to rho0*sqrt(1-mgz/e)		
		z_n = NeutEnergie*1e9/mg* ( powl(mt_get_double(v_mt_state)*(-1)*(powl(1-mg*z_ne/(NeutEnergie*1e9),1.5) - powl(1-mg*z_ns/(NeutEnergie*1e9),1.5)) - powl(1-mg*z_ns/(NeutEnergie*1e9),1.5),2.0/3.0) +1 ); */
		z_n = z_ns + (mt_get_double(v_mt_state)) * (z_ne-z_ns);
		//cout << "Die gew�rfelte H�he ist :" << z_n << endl;
			//sleep(5);
		
		if(protneut == NEUTRON)
		{			
			// Test ob Neutron �berhaupt da hingekommen w�re (Energy)
			BFeld(r_n,phi_n*conv,z_n, 0.0);
			crit = NeutEnergie - m_n*gravconst*z_n + (mu_n)*Bws;
			//printf("Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie*1.0e9,m_n*gravconst*z_n*1.0e9,(mu_nSI/ele_e)*Bws*1.0e9,crit*1.0e9);
			//fprintf(LOGSCR,"Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie*1.0e9,m_n*gravconst*z_n*1.0e9,(mu_nSI/ele_e)*Bws*1.0e9,crit*1.0e9);
			if((crit<0.0))
			{   //  Anfangsen. - Pot Energie + Energie im B-Feld
				printf("Dice anew... \n");
				fprintf(LOGSCR,"Dice anew... \n");
			}
			else if (crit>=0) 
					break;		
		}
		// proton
		else if ((protneut == PROTON) ||(protneut == ELECTRONS))
		{      
		// Test ob Neutron �berhaupt da hingekommen w�re (NeutEnergie)
			BFeld(r_n,phi_n*conv,z_n, 0.0);
			crit = NeutEnergie - m_n*gravconst*z_n + (mu_n)*Bws;
			//printf("Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie*1.0e9,m_n*gravconst*z_n*1.0e9,(mu_nSI/ele_e)*Bws*1.0e9,crit*1.0e9);
			//fprintf(LOGSCR,"Energy of supposed neutron of decay:\n (E_ges)%LG - (E_grav)%LG - (E_B)%LG = (dE)%LG\n",NeutEnergie*1.0e9,m_n*gravconst*z_n*1.0e9,(mu_nSI/ele_e)*Bws*1.0e9,crit*1.0e9);
			if((crit<0.0))
			{   //  Anfangsen. - Pot Energie + Energie im B-Feld
				printf("Dice anew... \n");
				fprintf(LOGSCR,"Dice anew... \n");
			}
			else if (crit>=0) 
					break;		
		}
		
	}

	alpha = alphas + (mt_get_double(v_mt_state)) * (alphae-alphas);      // constant angular distribution
	gammaa = acosl( cosl(gammas*conv) - mt_get_double(v_mt_state) * (cosl(gammas*conv) - cosl(gammae*conv)) ) / conv;   // isotropic emission characteristics
	
	
	// for checking the dependence of spin flip on inner rod current
	if((protneut==NEUTRON)&&(DiceRodField==1))
		RodFieldMultiplicator = (mt_get_double(v_mt_state));
	
	
	
	
	// random start times from 0 to 2e-5 for protons
	if(protneut==(PROTON||ELECTRONS))
		xstart =  (mt_get_double(v_mt_state))*2e-5;
	
	// start randomly in time during filling period
	if((FillingTime>0)&&(r_ne<0.12)&&(protneut==NEUTRON))
	{
		xstart = FillingTime * (mt_get_double(v_mt_state));
	}
		
	if((decay)&&(protneut==NEUTRON))
	{
		xend = - tau * log(mt_get_double(v_mt_state)) +xstart;            // Zeit nach der das Neutron zerfallen soll
	}	
	
		
	return ;
}
