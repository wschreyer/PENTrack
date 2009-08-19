#include "main.h"
//#include "vars.h"


// write the current configuration of the programm to screen and to the logfile
void PrintConfig(void)
{
	// ONSCREEN
	printf("The following parameters were set for program execution: \n");
	printf("protneut: %i => (1) neutron, (2) proton tracking, (3) magnetic field evaluation \n",protneut);
	if (protneut == NEUTRON)
	{
		if (polarisation == 1)
		{
			printf("Low field seekers are tracked! \n");		
		}
		else if (polarisation == 2)
		{
			printf("High field seekers are tracked! \n");		
		}
		else if (polarisation == 3)
		{
			printf("No polarisation! \n");		
		}
		else if (polarisation == 4)
		{
			printf("Polarisation is diced! \n");		
		}
		
		if (fireflekt == 1)
			printf("Filling Time %LG \n Reflection is on, ", FillingTime );
		else if (fireflekt == 0)
			printf("Filling Time %LG \n Reflection is off, ",FillingTime);
		if (fiBruteForce == 1)
			printf("Brute Force on, ");
		else if (fiBruteForce == 0)
			printf("Brute Force off, ");		
		
		if (clreflekt == 1)
			printf("Cleaning Time %LG \n Reflection is on, ", CleaningTime );
		else if (clreflekt == 0)
			printf("Cleaning  Time %LG \n Reflection is off, ",CleaningTime);
		if (clBruteForce == 1)
			printf("Brute Force on, ");
		else if (clBruteForce == 0)
			printf("Brute Force off, ");
		
		if (rureflekt == 1)
			printf("Ramp Up Time %LG \n Reflection is on, ", RampUpTime );
		else if (rureflekt == 0)
			printf("Ramp Up  Time %LG \n Reflection is off, ",RampUpTime);
		if (ruBruteForce == 1)
			printf("Brute Force on, ");
		else if (ruBruteForce == 0)
			printf("Brute Force off, ");
		
		if (ffreflekt == 1)
			printf("Full Field Time %LG \n Reflection is on, ", FullFieldTime );
		else if (ffreflekt == 0)
			printf("Full Field Time %LG \n Reflection is off, ",FullFieldTime);
		if (ffBruteForce == 1)
			printf("Brute Force on, ");
		else if (ffBruteForce == 0)
			printf("Brute Force off, ");
		
		if (rdreflekt == 1)
			printf("Ramp Down Time  %LG \n Reflection is on, ", RampDownTime );
		else if (rdreflekt == 0)
			printf("Ramp Down Time  %LG \n Reflection is off, ",RampDownTime);
		if (rdBruteForce == 1)
			printf("Brute Force on, ");
		else if (rdBruteForce == 0)
			printf("Brute Force off, ");
		
		if (coreflekt == 1)
			printf("Counting Time  \n Reflection is on, ");
		else if (coreflekt == 0)
			printf("Coutning Time  \n Reflection is off, ");
		if (coBruteForce == 1)
			printf("Brute Force on, ");
		else if (coBruteForce == 0)
			printf("Brute Force off, ");
		
			
	}
	printf("The choice of Bfield is: %i  (0) interpolated field, (1) no field, (2) check interpolation routine \n",bfeldwahl);
	printf("The choice of reflection is: %i  (1) specular, (2) diffuse, (3) statistically specular or diffuse \n",diffuse);
	printf("Choice of output: %i  (1) Endpoints and track (2) only endpoints (3) Endpoints, track and spin \n (4) Endpoints and spin (5) nothing \n ", ausgabewunsch);
	
	// LOGSCREEN
	fprintf(LOGSCR,"The following parameters were set for program execution: \n");
	fprintf(LOGSCR,"protneut: %i => (1) neutron, (2) proton tracking, (3) magnetic field evaluation \n",protneut);
	if (protneut == NEUTRON)
	{
		if (polarisation == 1)
		{
			fprintf(LOGSCR,"Low field seekers are tracked! \n");		
		}
		else if (polarisation == 2)
		{
			fprintf(LOGSCR,"High field seekers are tracked! \n");		
		}
		else if (polarisation == 3)
		{
			fprintf(LOGSCR,"No polarisation! \n");		
		}
		else if (polarisation == 4)
		{
			fprintf(LOGSCR,"Polarisation is diced! \n");		
		}
		
		if (fireflekt == 1)
			fprintf(LOGSCR,"Filling Time %LG \n Reflection is on, ", FillingTime );
		else if (fireflekt == 0)
			fprintf(LOGSCR,"Filling Time %LG \n Reflection is off, ",FillingTime);
		if (fiBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (fiBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		if (clreflekt == 1)
			fprintf(LOGSCR,"Cleaning Time %LG \n Reflection is on, ", CleaningTime );
		else if (clreflekt == 0)
			fprintf(LOGSCR,"Cleaning  Time %LG \n Reflection is off, ",CleaningTime);
		if (clBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (clBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		if (rureflekt == 1)
			fprintf(LOGSCR,"Ramp Up Time %LG \n Reflection is on, ", RampUpTime );
		else if (rureflekt == 0)
			fprintf(LOGSCR,"Ramp Up  Time %LG \n Reflection is off, ",RampUpTime);
		if (ruBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (ruBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		if (ffreflekt == 1)
			fprintf(LOGSCR,"Full Field Time %LG \n Reflection is on, ", FullFieldTime );
		else if (ffreflekt == 0)
			fprintf(LOGSCR,"Full Field Time %LG \n Reflection is off, ",FullFieldTime);
		if (ffBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (ffBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		if (rdreflekt == 1)
			fprintf(LOGSCR,"Ramp Down Time %LG \n Reflection is on, ", RampDownTime );
		else if (rdreflekt == 0)
			fprintf(LOGSCR,"Ramp Down Time %LG \n Reflection is off, ",RampDownTime);
		if (rdBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (rdBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		if (coreflekt == 1)
			fprintf(LOGSCR,"Counting Time  \n Reflection is on, ");
		else if (coreflekt == 0)
			fprintf(LOGSCR,"Coutning Time  \n Reflection is off, ");
		if (coBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (coBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		
		
		
	}
	fprintf(LOGSCR,"The choice of Bfield is: %i  (0) interpolated field, (1) no field, (2) check interpolation routine \n",bfeldwahl);
	fprintf(LOGSCR,"The choice of reflection is: %i  (1) speculat, (2) diffuse, (3) statistically specular or diffuse \n",bfeldwahl);
	fprintf(LOGSCR,"Choice of output: %i  (1) Endpoints and track (2) only endpoints (3) Endpoints, track and spin \n (4) Endpoints and spin (5) nothing \n ", ausgabewunsch);
	
	
}


long double lastr = INFINITY, lastphi = INFINITY, lastz = INFINITY, lastt = INFINITY; // remember last call to BFeld to avoid multiple calculations of same value

void BFeld (long double rloc, long double philoc, long double zloc, long double t){      //B-Feld am Ort des Teilchens berechnen
	if (lastr == rloc && lastphi == philoc && lastz == zloc && lastt == t)
		return;
	else{
		lastr = rloc; lastphi = philoc; lastz = zloc; lastt = t;
	}
	//clock_t mytime1, mytime2;
	long double Brtemp,Bztemp,dBdrtemp,dBdztemp;
	
	if (protneut == NEUTRON) // switching of field only for neutrons valid
		SwitchField(t);
	else if(protneut == PROTON || protneut == ELECTRONS || protneut == BF_CUT)
		BFeldSkal=BFeldSkalGlobal;
		
	
	Bnull();
	switch (bfeldwahl)
	{
		
		case 0:	//mytime1 = clock(); 
						BInterpol(rloc, philoc, zloc);
						//mytime2 = clock();
						//timer1 =( ((long double)mytime2 - (long double)mytime1) / CLOCKS_PER_SEC ) + timer1;
						
						break;
		case 2:	//mytime1 = clock(); 
						BInterpol(rloc, philoc, zloc);
						//mytime2 = clock();
						//timer1 =( ((long double)mytime2 - (long double)mytime1) / CLOCKS_PER_SEC ) + timer1;
						Brtemp=Br; Bztemp=Bz;dBdrtemp=dBdr;dBdztemp=dBdz;
						Bwsomaple(rloc, philoc, zloc);
						printf("BrI BrM deltaBr %17LG %17LG %17LG \n", Brtemp, Br, (Brtemp-Br)/Br);
						printf("BzI BzM deltaBz %17LG %17LG %17LG \n", Bztemp, Bz, (Bztemp-Bz)/Bz);						
						Br=(Brtemp-Br)/Br; Bz=(Bztemp-Bz)/Bz; dBdr=(dBdrtemp-dBdr)/dBdr; dBdz=(dBdztemp-dBdz)/dBdz;
						
		case 3: Banalytic(rloc, philoc, zloc, t);
						break;
		case 4: BForbes(rloc, philoc, zloc, t);
						break;
	}   
	
	Bws = sqrtl(Br*Br+Bz*Bz+Bphi*Bphi);			
	if (Bws>1e-31)
	{			
		dBdr   = (Br*dBrdr + Bphi*dBphidr + Bz*dBzdr)  /Bws;
		dBdz   = (Br*dBrdz + Bphi*dBphidz + Bz*dBzdz)  /Bws;
		dBdphi = (Br*dBrdphi + Bphi*dBphidphi + Bz*dBzdphi)/Bws;
	}
	
	Feldcount++;
	return;
}


void EFeld(long double rloc, long double philoc, long double zloc){
	if(((protneut == PROTON) || (protneut == ELECTRONS))&(EFeldSkal!=0)){
		EInterpol(rloc, philoc, zloc);
		/*switch (Efeldwahl){
			case 0: EInterpol(rloc, philoc, zloc); break;
			case 2: ETrichter(rloc,philoc,zloc); break;
			case 3: EKegelpot(rloc,philoc,zloc); break;
            case 4: Sack(rloc,philoc,zloc); break;
            case 5: Flaechdet(rloc,philoc,zloc);
        }     */
		
    }
	 else if(EFeldSkal==0)
	 {		 
		 Er = 0;
		dErdr =0;
		dErdz = 0;
		Ephi = 0;
		dEphidr = 0;
		dEphidz = 0;
		Ez = 0;
		dEzdr = 0;
		dEzdz = 0;
  
 }
 Feldcount++;	
 return;
}

//-----------------------------------------------------------------------------
void SwitchField(long double t)
{        //B-Feld nach Wunsch skalieren, BFeldSkal wird an alle B Komp und deren Ableitungen multipliziert
	
		
		if ((t < FillingTime)&&(FillingTime>0))
		{      // filling in neutrons
			BFeldSkal = 0;
			BruteForce = fiBruteForce; // no spin tracking
			reflekt = fireflekt;
			spinflipcheck = fispinflipcheck;		
		}
		
		else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
		{      // spectrum cleaning
			BFeldSkal = 0;
			BruteForce = clBruteForce; // no spin tracking
			reflekt = clreflekt;
			spinflipcheck = clspinflipcheck;
		}
		else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
		{   // ramping up field
			//linear: BFeldSkal = (t-CleaningTime+1e-15)/RampUpTime; // ramp field
			//now smoothly with a cosine 
			BFeldSkal=(0.5-0.5*cosl(pi*(t-CleaningTime-FillingTime)/RampUpTime)) * BFeldSkalGlobal ;
			BruteForce = ruBruteForce; 
			reflekt = rureflekt;
			spinflipcheck = ruspinflipcheck;
		}
		else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
		{  // storage time
			BFeldSkal = BFeldSkalGlobal;
			if(protneut == NEUTRON)
			{
				BruteForce = ffBruteForce;  // start spin tracking
				reflekt = ffreflekt;
				spinflipcheck = ffspinflipcheck;
			}
		}
		else if ((t >= (RampUpTime+CleaningTime+FillingTime+FullFieldTime)) && (t < (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime)) && (RampDownTime != 0))
		{   // ramping down field
			//linear: BFeldSkal = 1- (t-(RampUpTime+CleaningTime+FullFieldTime)) / RampDownTime + 1e-15;  // ramp down field
			//now smoothly with a cosine 
			BFeldSkal=(0.5+0.5*cosl(pi*(t-(RampUpTime+CleaningTime+FillingTime+FullFieldTime)) / RampDownTime)) * BFeldSkalGlobal;
			BruteForce = rdBruteForce; // no spin tracking, neutrons shall stay in trap through reflection
			reflekt = rdreflekt;
			spinflipcheck=rdspinflipcheck;
		}
		else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
		{      // emptying of neutrons into detector
			BFeldSkal = 0;
			BruteForce = coBruteForce;
			spinflipcheck= cospinflipcheck;
		}
		else
			BFeldSkal = BFeldSkalGlobal;
	
	
	if (FieldOscillation==1){
		BFeldSkal = BFeldSkal* (1+(BFeldSkal*OscillationFraction*sinl(OscillationFrequency*2*pi*t)));		
	}
	return;
}

//-----------------------------------------------------------------------------
// Vektor W (kein Ortsvektor, sondern B-Feld o.�hnliches) von zyl in kart koord xyz umrechnen ## Wr, Wphi, Wz: Eingabewerte ## Wx0, Wy0, Wz0 werden ver�ndert von CylKartCoord
// phi ist dabei die winkelkoordinate des aktuellen ortsvektors
void CylKartCoord(long double Wr, long double Wphi, long double Wz, long double phi, long double *Wx0, long double *Wy0, long double *Wz0){
	(*Wx0) = cosl(phi) * Wr - sinl(phi) * Wphi;
	(*Wy0) = sinl(phi) * Wr + cosl(phi) * Wphi;
	(*Wz0) = Wz;
	return;
}


//-----------------------------------------------------------------------------
// Vektor W (kein Ortsvektor, sondern B-Feld o.�hnliches) von kart in zyl koord xyz umrechnen ## Wx, Wy, Wz: Eingabewerte ## Wr0, Wphi0, Wz0 werden ver�ndert von KartCylCoord
// phi ist dabei die winkelkoordinate des aktuellen ortsvektors
void KartCylCoord(long double Wx, long double Wy, long double Wz, long double phi, long double *Wr0, long double *Wphi0, long double *Wz0){
	(*Wr0) = cosl(phi) * Wx + sinl(phi) * Wy;
	(*Wphi0) = (-1) * sinl(phi) * Wx + cosl(phi) * Wy;
	(*Wz0) = Wz;
	return;
}


// return absolute value of a 3 vector in cartesian coordinates
long double AbsValueCart(long double x, long double y, long double z)
{
	return sqrt(x*x + y*y + z*z);
}

