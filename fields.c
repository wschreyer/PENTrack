#include "main.h"

list<TabField*> fields;


// read table files or conductor file
void PrepareFields(){
	if ((bfeldwahl == 0 || bfeldwahl == 2) && (BFeldSkalGlobal != 0 || EFeldSkal != 0))
	{
		TabField *tf;
		for (list<string>::iterator i = fieldvaltab.begin(); i != fieldvaltab.end(); i++){
			tf = new TabField((inpath + '/' + *i).c_str());
			fields.push_back(tf);
		}
	}
	else if(bfeldwahl==4)
	{
		ReadMagnets();
	
		printf("\n \n Test of integration\n");
		//	long double TestInt;
		BFeldSkal=1.0;
		for (int a = 0;a<1;a++)
		{
		
			BFeld(0.3,0,0.1, 500.0);
			cout << "T" << endl;
		}
		printf("Br = %.17LG \n",Br);
		printf("dBrdr = %.17LG \n",dBrdr);
		printf("dBrdz = %.17LG \n",dBrdz);
		printf("Bz = %.17LG \n",Bz);
		printf("dBzdr = %.17LG \n",dBzdr);
		printf("dBzdz = %.17LG \n",dBzdz);	
	}	
	if (Emin_n >= 1e20) Emin_n = 0; // if Emin_n was not set, set it to zero
	
}


// free memory
void FreeFields(){
	for (list<TabField*>::iterator i = fields.begin(); i != fields.end(); i++)
		delete (*i);	
}


// fill global B field variables with values
void BFeld (long double rloc, long double philoc, long double zloc, long double t){      //B-Feld am Ort des Teilchens berechnen
	timeval fieldstart, fieldend;
	gettimeofday(&fieldstart, NULL);	
	
	SwitchField(t); // set BFeldSkal for different experiment phases	
	
	Bnull();	// first, set all to zero
	long double Brtemp, Bztemp, dBdrtemp, dBdztemp;
	switch (bfeldwahl)
	{
		
		case 0:	if (BFeldSkal != 0){
					list<TabField*>::iterator i = fields.begin(); 
					while (i != fields.end()){						
						if ((*i)->BInterpol(rloc, zloc))
							break;
						else if (++i == fields.end()){
							Log("\nThe particle has left fieldval boundaries: r=%LG, z=%LG! Stopping particle...\n", rloc, zloc);
							kennz = KENNZAHL_LEFT_FIELD;
							stopall = 1;
						}		
					}				
				}
				RacetrackField(rloc,philoc,zloc);
				break;
				
		case 2:	if (BFeldSkal != 0){
					list<TabField*>::iterator i = fields.begin(); 
					while (i != fields.end()){						
						if ((*i)->BInterpol(rloc, zloc))
							break;
						else if (++i == fields.end()){
							Log("\nThe particle has left fieldval boundaries: r=%LG, z=%LG! Stopping particle...\n", rloc, zloc);
							kennz = KENNZAHL_LEFT_FIELD;
							stopall = 1;
						}		
					}				
				}
				RacetrackField(rloc,philoc,zloc);
				Brtemp = Br; Bztemp = Bz; dBdrtemp = dBdr; dBdztemp = dBdz;
				Bwsomaple(rloc, philoc, zloc);
				printf("BrI BrM deltaBr %17LG %17LG %17LG \n", Brtemp, Br, (Brtemp-Br)/Br);
				printf("BzI BzM deltaBz %17LG %17LG %17LG \n", Bztemp, Bz, (Bztemp-Bz)/Bz);						
				Br=(Brtemp-Br)/Br; Bz=(Bztemp-Bz)/Bz; dBdr=(dBdrtemp-dBdr)/dBdr; dBdz=(dBdztemp-dBdz)/dBdz;
				break;	
					
		case 3: Banalytic(rloc, philoc, zloc, t);
				break;
				
		case 4: BForbes(rloc, philoc, zloc, t);
				RacetrackField(rloc,philoc,zloc);
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
	
	gettimeofday(&fieldend, NULL);
	FieldTime += fieldend.tv_sec - fieldstart.tv_sec + (fieldend.tv_usec - fieldstart.tv_usec)/1e6;	
	
	return;
}


void SwitchField(long double t){        //B-Feld nach Wunsch skalieren, BFeldSkal wird an alle B Komp und deren Ableitungen multipliziert		
		if ((t < FillingTime)&&(FillingTime>0))
		{      // filling in neutrons
			BFeldSkal = 0;
			BruteForce = fiBruteForce; // no spin tracking
			reflekt = fireflekt;
			spinflipcheck = fispinflipcheck;		
			ExpPhase=1;
		}
		
		else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
		{      // spectrum cleaning
			BFeldSkal = 0;
			BruteForce = clBruteForce; // no spin tracking
			reflekt = clreflekt;
			spinflipcheck = clspinflipcheck;
			ExpPhase=2;
		}
		else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
		{   // ramping up field
			//linear: BFeldSkal = (t-CleaningTime+1e-15)/RampUpTime; // ramp field
			//now smoothly with a cosine 
			BFeldSkal=(0.5-0.5*cosl(pi*(t-CleaningTime-FillingTime)/RampUpTime)) * BFeldSkalGlobal ;
			BruteForce = ruBruteForce; 
			reflekt = rureflekt;
			spinflipcheck = ruspinflipcheck;
			ExpPhase=3;
		}
		else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
		{  // storage time
			BFeldSkal = BFeldSkalGlobal;
			if(protneut == NEUTRON)
			{
				BruteForce = ffBruteForce;  // start spin tracking
				reflekt = ffreflekt;
				spinflipcheck = ffspinflipcheck;
				ExpPhase=4;
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
			ExpPhase=5;
		}
		else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
		{      // emptying of neutrons into detector
			BFeldSkal = 0;
			BruteForce = coBruteForce;
			spinflipcheck= cospinflipcheck;
			ExpPhase=6;
		}
		else
			BFeldSkal = BFeldSkalGlobal;
	
	
	if (FieldOscillation==1){
		BFeldSkal = BFeldSkal* (1+(BFeldSkal*OscillationFraction*sinl(OscillationFrequency*2*pi*t)));		
	}
	return;
}


// produce zero field
long double Bnull(){
	Br = 0.0;
	Bphi = 0.0;
	Bz = 0.0;
	dBrdr=0.0;
	dBrdphi=0.0;
	dBrdz=0.0;
	dBphidr=0.0;
	dBphidphi=0.0;
	dBphidz=0.0;
	dBzdr=0.0;
	dBzdphi=0.0;
	dBzdz=0.0;

	Bws = 0.0;
	dBdr=0;
	dBdphi=0;
	dBdz=0;
	
	return 0;
}

// produce linearly rising B-field in z-direction and linearly changing in time
long double Banalytic(long double r,long double phi,long double z, long double t){
	long double dBx = 5;
	long double x0 = 1.1;
	long double alpha = 0.5;
	long double beta = 0.025;
	
	
	// transform Bfield to higher z-values
	long double offset = -0.0;
	z = z+offset;
	
	BFeldSkal = 1;
	
	Br = (-dBx*(x0*z-0.5*z*z)+6)*(alpha+beta*t);
	Bphi = 0;
	Bz = 0;
	dBrdr=0.0;
    dBrdphi=0.0;
    dBrdz=-dBx*(alpha+beta*t)*(x0-z);
	dBphidr=0.0;
    dBphidphi=0.0;
	dBphidz=0.0;
    dBzdr=0.0;
    dBzdphi=0.0;
    dBzdz=0.0;
        	
	Br *= BFeldSkal;
	dBrdr *= BFeldSkal;
	dBrdz *= BFeldSkal;
	Bphi *= RodFieldMultiplicator * BFeldSkal;
	dBphidr *= RodFieldMultiplicator * BFeldSkal;
	dBphidz *= RodFieldMultiplicator * BFeldSkal;
	dBphidphi = 0.0;
	Bz *= BFeldSkal;
	dBzdr *= BFeldSkal;
	dBzdz *= BFeldSkal;
	
	return 0;
}



// fill global E field variables with values
void EFeld(long double rloc, long double philoc, long double zloc){
	timeval fieldstart, fieldend;
	gettimeofday(&fieldstart, NULL);	
	
	Er = 0.0;
	Ephi = 0.0;
	Ez = 0.0;
	dErdr=0.0;
	dErdz=0.0;
	dEphidr=0.0;
	dEphidz=0.0;
	dEzdr=0.0;
	dEzdz=0.0;

	if (EFeldSkal != 0 && (bfeldwahl == 0 || bfeldwahl == 2)){
		list<TabField*>::iterator i = fields.begin(); 
		while (i != fields.end()){						
			if ((*i)->EInterpol(rloc, zloc))
				break;
			else if (++i == fields.end()){
				Log("\nThe particle has left fieldval boundaries: r=%LG, z=%LG! Stopping particle...\n", rloc, zloc);
				kennz = KENNZAHL_LEFT_FIELD;
				stopall = 1;
			}		
		}				
		/*switch (Efeldwahl){
			case 0: EInterpol(rloc, philoc, zloc); break;
			case 2: ETrichter(rloc,philoc,zloc); break;
			case 3: EKegelpot(rloc,philoc,zloc); break;
            case 4: Sack(rloc,philoc,zloc); break;
            case 5: Flaechdet(rloc,philoc,zloc);
        }     */
		
    }
	Feldcount++;	

	gettimeofday(&fieldend, NULL);
	FieldTime += fieldend.tv_sec - fieldstart.tv_sec + (fieldend.tv_usec - fieldstart.tv_usec)/1e6;	
	
	return;
}
