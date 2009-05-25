#include "main.h"
//#include "vars.h"


// read the conductor data from a .cond file written by Vectorfield Opera
int ReadMagnets(void)
{
	long double CoilAngle1=-1, CoilAngle2=-1, CoilAngle3=-1;
	long double dummy, CoilZOffset, Symmetry;
	long double CoilPoint1R, CoilPoint1Z,  CoilPoint2R,  CoilPoint2Z,  CoilPoint3R,  CoilPoint3Z,  CoilPoint4R,  CoilPoint4Z;
	
	
	FILE *cfg = NULL;
	
	char line[1024];
	string path(inpath + "/coils.cond");
	cfg = fopen(path.c_str(),mode_r);
	if(cfg == NULL){  
		printf("File coils.cond could not be opened!!!");
		fprintf(LOGSCR,"File coils.cond could not be opened!!!");
		return 1;
	}
	
	do // read the file
	{
		n=0;
		
		//for (i=1;i<=1024;i++)			
		//	line[i]='0'; // reset string
		
		fgets(line,1024,cfg); // get a line
		
		//if((line[0]=='C')||(line[0]=='Q')) continue; // if its a comment goto next line
			
		if((line[0]=='D')&&(line[7]=='G')&&(line[8]=='R')&&(line[9]=='A'))
		{cout << "RACETRACk COIL DETECTED" << endl; continue; }// discard racetrack coils
			
		if((line[0]=='D')&&(line[7]=='G')&&(line[8]=='S')&&(line[9]=='O'))
		{
			//printf(line);
			fgets(line,1024,cfg); // 1st line of coil definition: local coordinate system 1, only used for racetrack so far
		
			fgets(line,1024,cfg); // 2nd line of coil definition: local coordinat system 2 offsets
			sscanf(line,"%25Lf %25Lf %25Lf", &dummy, &dummy, &CoilZOffset);
		
			fgets(line,1024,cfg); // 3rd line of coil definition: local coordinat system 2 angles
			sscanf(line,"%25Lf %25Lf %25Lf", &CoilAngle1, &CoilAngle2, &CoilAngle3);
			if((CoilAngle1!=90)||(CoilAngle2!=0)||(CoilAngle3!=90))  // check if coils have right orientation: z as rotational axis
			{
				cout << "CoilAngle1 "	<< CoilAngle1 << "CoilAngle2 "<< CoilAngle2 << "CoilAngle3 " << CoilAngle3 << endl;
				cout << "Reading of coil number "	<< CoilNr << " failed!!! Wrong orientation!!! Aborting!!" << endl;
				return 1;
			}
		
			fgets(line,1024,cfg); // 4th line of coil definition: cross section points 1 and 2 
			sscanf(line,"%25Lf %25Lf %25Lf %25Lf", &CoilPoint1R, &CoilPoint1Z, &CoilPoint2R, &CoilPoint2Z);
			
			fgets(line,1024,cfg); // 5th line of coil definition: cross section points 3 and 4 
			sscanf(line,"%25Lf %25Lf %25Lf %25Lf", &CoilPoint3R, &CoilPoint3Z, &CoilPoint4R, &CoilPoint4Z);
			if (CoilPoint3Z!=CoilPoint1Z)
				aF[CoilNr]=fabsl((CoilPoint3Z-CoilPoint1Z)/2.0);
			else aF[CoilNr]=fabsl((CoilPoint2Z-CoilPoint1Z)/2.0);
				
			if (CoilPoint2R!=CoilPoint1R)
				bF[CoilNr]=fabsl((CoilPoint2R-CoilPoint1R)/2.0);
			else bF[CoilNr]=fabsl((CoilPoint3R-CoilPoint1R)/2.0);
				
			if (CoilPoint1R<CoilPoint2R)
				R_0[CoilNr]= CoilPoint1R+(CoilPoint2R-CoilPoint1R)/2.0;
			
			zoffset[CoilNr]=CoilZOffset+CoilPoint1Z+(CoilPoint3Z-CoilPoint1Z)/2.0;
			
			fgets(line,1024,cfg); // 6th line of coil definition: Curvatures -> not used
			
			fgets(line,1024,cfg); // 7th line of coil definition: current density, symmetry, label
			sscanf(line,"%25Lf %25Lf ", &J_0[CoilNr], &Symmetry);
			if ((Symmetry!=0)&&(Symmetry!=1))
			{
				cout << "Reading of coil number "	<< CoilNr << " failed!!! Wrong symmetry!!! Aborting!!" << endl;
				return 1;
			}
			
			fgets(line,1024,cfg); // 8th line of coil definition: reflection codes -> not used
			
			fgets(line,1024,cfg); // 9th line of coil definition: accuracy of field calculation in Opera -> not used		
			
			cout << "Coil " << CoilNr << " a " << aF[CoilNr] << " b " << bF[CoilNr] << " R_0 " << R_0[CoilNr] << " zoffset " << zoffset[CoilNr] << " J_0 " << J_0[CoilNr] << endl ;
			CoilNr++;
		}
			
			// printf(line);
		
	}while(!feof(cfg));
	
	printf("%i solenoids were read in!!! \n",CoilNr);
	fprintf(LOGSCR,"%i solenoids were read in!!! \n",CoilNr);
	
	return 0;
}


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
		if(fislit == 1)
			printf("slit on.\n");
		else if(fislit == 0)
			printf("slit off.\n");
		if(fiDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(fiDetOpen == 0)
			printf("The detector is hidden from UCN.\n");
		
		
		if (clreflekt == 1)
			printf("Cleaning Time %LG \n Reflection is on, ", CleaningTime );
		else if (clreflekt == 0)
			printf("Cleaning  Time %LG \n Reflection is off, ",CleaningTime);
		if (clBruteForce == 1)
			printf("Brute Force on, ");
		else if (clBruteForce == 0)
			printf("Brute Force off, ");
		if(clslit == 1)
			printf("slit on.\n");
		else if(clslit == 0)
			printf("slit off.\n");
		if(clDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(clDetOpen == 0)
			printf("The detector is hidden from UCN.\n");		
		
		if (rureflekt == 1)
			printf("Ramp Up Time %LG \n Reflection is on, ", RampUpTime );
		else if (rureflekt == 0)
			printf("Ramp Up  Time %LG \n Reflection is off, ",RampUpTime);
		if (ruBruteForce == 1)
			printf("Brute Force on, ");
		else if (ruBruteForce == 0)
			printf("Brute Force off, ");
		if(ruslit == 1)
			printf("slit on.\n");
		else if(ruslit == 0)
			printf("slit off.\n");
		if(clDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(fiDetOpen == 0)
			printf("The detector is hidden from UCN.\n");
		
		if (ffreflekt == 1)
			printf("Full Field Time %LG \n Reflection is on, ", FullFieldTime );
		else if (ffreflekt == 0)
			printf("Full Field Time %LG \n Reflection is off, ",FullFieldTime);
		if (ffBruteForce == 1)
			printf("Brute Force on, ");
		else if (ffBruteForce == 0)
			printf("Brute Force off, ");
		if(ffslit == 1)
			printf("slit on.\n");
		else if(ffslit == 0)
			printf("slit off.\n");
		if(ffDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(ffDetOpen == 0)
			printf("The detector is hidden from UCN.\n");		
		
		if (rdreflekt == 1)
			printf("Ramp Down Time  %LG \n Reflection is on, ", RampDownTime );
		else if (rdreflekt == 0)
			printf("Ramp Down Time  %LG \n Reflection is off, ",RampDownTime);
		if (rdBruteForce == 1)
			printf("Brute Force on, ");
		else if (rdBruteForce == 0)
			printf("Brute Force off, ");
		if(rdslit == 1)
			printf("slit on.\n");
		else if(rdslit == 0)
			printf("slit off.\n");
		if(rdDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(rdDetOpen == 0)
			printf("The detector is hidden from UCN.\n");		
		
		if (coreflekt == 1)
			printf("Counting Time  \n Reflection is on, ");
		else if (coreflekt == 0)
			printf("Coutning Time  \n Reflection is off, ");
		if (coBruteForce == 1)
			printf("Brute Force on, ");
		else if (coBruteForce == 0)
			printf("Brute Force off, ");
		if(coslit == 1)
			printf("slit on.\n");
		else if(coslit == 0)
			printf("slit off.\n");
		if(coDetOpen == 1)
			printf("The detector is visible to UCN.\n");
		else if(coDetOpen == 0)
			printf("The detector is hidden from UCN.\n");		
		
			
	}
	printf("The choice of Bfield is: %i  (0) interpolated field, (1) no field, (2) check interpolation routine \n",bfeldwahl);
	printf("The choice of reflection is: %i  (1) specular, (2) diffuse, (3) statistically specular or diffuse \n",diffuse);
	printf("The probability for diffuse reflection is: %LG\n", DiffProb);
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
		if(fislit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(fislit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(fiDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(fiDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");
		
		
		if (clreflekt == 1)
			fprintf(LOGSCR,"Cleaning Time %LG \n Reflection is on, ", CleaningTime );
		else if (clreflekt == 0)
			fprintf(LOGSCR,"Cleaning  Time %LG \n Reflection is off, ",CleaningTime);
		if (clBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (clBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		if(clslit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(clslit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(clDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(clDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");		
		
		if (rureflekt == 1)
			fprintf(LOGSCR,"Ramp Up Time %LG \n Reflection is on, ", RampUpTime );
		else if (rureflekt == 0)
			fprintf(LOGSCR,"Ramp Up  Time %LG \n Reflection is off, ",RampUpTime);
		if (ruBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (ruBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		if(ruslit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(ruslit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(clDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(fiDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");
		
		if (ffreflekt == 1)
			fprintf(LOGSCR,"Full Field Time %LG \n Reflection is on, ", FullFieldTime );
		else if (ffreflekt == 0)
			fprintf(LOGSCR,"Full Field Time %LG \n Reflection is off, ",FullFieldTime);
		if (ffBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (ffBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		if(ffslit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(ffslit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(ffDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(ffDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");		
		
		if (rdreflekt == 1)
			fprintf(LOGSCR,"Ramp Down Time %LG \n Reflection is on, ", RampDownTime );
		else if (rdreflekt == 0)
			fprintf(LOGSCR,"Ramp Down Time %LG \n Reflection is off, ",RampDownTime);
		if (rdBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (rdBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		if(rdslit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(rdslit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(rdDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(rdDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");		
		
		if (coreflekt == 1)
			fprintf(LOGSCR,"Counting Time  \n Reflection is on, ");
		else if (coreflekt == 0)
			fprintf(LOGSCR,"Coutning Time  \n Reflection is off, ");
		if (coBruteForce == 1)
			fprintf(LOGSCR,"Brute Force on, ");
		else if (coBruteForce == 0)
			fprintf(LOGSCR,"Brute Force off, ");
		if(coslit == 1)
			fprintf(LOGSCR,"slit on.\n");
		else if(coslit == 0)
			fprintf(LOGSCR,"slit off.\n");
		if(coDetOpen == 1)
			fprintf(LOGSCR,"The detector is visible to UCN.\n");
		else if(coDetOpen == 0)
			fprintf(LOGSCR,"The detector is hidden from UCN.\n");		
		
		
		
	}
	fprintf(LOGSCR,"The choice of Bfield is: %i  (0) interpolated field, (1) no field, (2) check interpolation routine \n",bfeldwahl);
	fprintf(LOGSCR,"The choice of reflection is: %i  (1) speculat, (2) diffuse, (3) statistically specular or diffuse \n",bfeldwahl);
	fprintf(LOGSCR,"The probability for diffuse reflection is: %LG\n", DiffProb);
	fprintf(LOGSCR,"Choice of output: %i  (1) Endpoints and track (2) only endpoints (3) Endpoints, track and spin \n (4) Endpoints and spin (5) nothing \n ", ausgabewunsch);
	
	
}



void Entkommen(long double *ystart, long double t, long double H)        //Klassifizierung des Endes der Neutronenbahn und Abbruch des aktuellen Teilchens
{	
	if((ystart[1]<rmin)&&(ystart[1]>rmax)&&(ystart[3]<zmin)&&(ystart[3]>zmax))
	{
		kennz=2;  
		stopall=1;
		printf("Particle has hit outer boundaries: Stopping it! t=%LG r=%LG z=%LG\n",t,ystart[1],ystart[3]);
	   fprintf(LOGSCR,"Particle has hit outer boundaries: Stopping it! t=%LG r=%LG z=%LG\n",t,ystart[1],ystart[3]);
		return;
	}
	   
	
	return;
}

void BFeld (long double rloc, long double philoc, long double zloc, long double t){      //B-Feld am Ort des Teilchens berechnen
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
//	else
//		Bnull();
	
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
			DetOpen=fiDetOpen; // detector hidden
			slit = fislit;     // emptying slit closed
			BruteForce = fiBruteForce; // no spin tracking
			reflekt = fireflekt;
			spinflipcheck = fispinflipcheck;		
		}
		
		else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
		{      // spectrum cleaning
			BFeldSkal = 0;
			DetOpen= clDetOpen; // detector hidden
			slit = clslit;     // emptying slit closed
			BruteForce = clBruteForce; // no spin tracking
			reflekt = clreflekt;
			spinflipcheck = clspinflipcheck;
		}
		else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
		{   // ramping up field
			//linear: BFeldSkal = (t-CleaningTime+1e-15)/RampUpTime; // ramp field
			//now smoothly with a cosine 
			BFeldSkal=(0.5-0.5*cosl(pi*(t-CleaningTime-FillingTime)/RampUpTime)) * BFeldSkalGlobal ;
			DetOpen=ruDetOpen; // detector hidden
			slit = ruslit;
			BruteForce = ruBruteForce; 
			reflekt = rureflekt;
			spinflipcheck = ruspinflipcheck;
		}
		else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
		{  // storage time
			BFeldSkal = BFeldSkalGlobal;
			DetOpen=ffDetOpen;  // detect high field seeker and marginal UCN
			slit = ffslit; // emptying slit open for depolarised neutrons
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
			DetOpen=rdDetOpen;  // detect UCN
			slit = rdslit; // slit open to empty neutrons already
			BruteForce = rdBruteForce; // no spin tracking, neutrons shall stay in trap through reflection
			reflekt = rdreflekt;
			spinflipcheck=rdspinflipcheck;
		}
		else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
		{      // emptying of neutrons into detector
			BFeldSkal = 0;
			DetOpen= coDetOpen;  // detect UCN
			slit = coslit;
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


//======== The endcodes of trajectories are added up here regarding the particle typ =======================================
void IncrementCodes(int kennz)
{	// protneut % 3 = {0, 1, 2} <<<	ELECTRONS (=6) % 3 = 0
	// 								NEUTRON   (=1) % 3 = 1	
	// 								PROTON    (=2) % 3 = 2
	switch (kennz)
	{	case  0:	kennz0[protneut % 3]++;
					break;
		case  1:	kennz1[protneut % 3]++;
					break;
		case  2:	kennz2[protneut % 3]++;
					break;
		case  3:	kennz3[protneut % 3]++;
					break;
		case  4:	kennz4[protneut % 3]++;
					break;
		case  5:	kennz5[protneut % 3]++;
					break;
		case  6:	kennz6[protneut % 3]++;
					break;
		case  7:	kennz7[protneut % 3]++;
					break;
		case  8:	kennz8[protneut % 3]++;
					break;
		case  9:	kennz9[protneut % 3]++;
					break;
		case 10:	kennz10[protneut % 3]++;
					break;
		case 11:	kennz11[protneut % 3]++;
					break;
		case 12:	kennz12[protneut % 3]++;
					break;
		case 99:	kennz99[protneut % 3]++;
					break;
	}
}
//======== end of IncrementCodes ===========================================================================================

//======== Output of the endcodes ==========================================================================================
void OutputCodes(int iMC)
{	if(decay.on == 2)
	{	printf("\nThe calculations of %li particle(s) yielded:\n"
		       "endcode    out of %i neutron(s)                            \tout of %li proton(s)                              \tout of %li electron(s)\n"
		       "   0       %*li (not categorized)                          \t%*li (not categorized)                          \t%*li (not categorized)\n"
		       "   1       %*li (did not finish)                           \t%*li (did not finish in %2.1LE s)              \t%*li (did not finish in %2.1LE s)\n"
		       "   2       %*li (hit outer boundaries)                     \t%*li (hit outer boundaries)                     \t%*li (hit outer boundaries)\n"
		       "   3       %*li (hit the walls)                            \t%*li (hit the walls)                            \t%*li (hit the walls)\n"
		       "   4       %*li (escaped to the top)                       \t%*li (escaped to the top)                       \t%*li (escaped to the top)\n"
		       "   5       %*li (escaped through the exit slit)            \t%*li (escaped through the exit slit)            \t%*li (escaped through the exit slit)\n"
		       "   6       %*li (hit the detector from the bottom)         \t%*li (hit the detector from the bottom)         \t%*li (hit the detector from the bottom)\n"
		       "   7       %*li (statistically absorbed)                   \t%*li (hit the detector from the top)            \t%*li (hit the detector from the top)\n"
		       "   8       %*li (decayed)                                  \t   -                                            	   -\n"
		       "   9       %*li (were eaten by the absorber)               \t%*li (were eaten by the absorber)               \t%*li (were eaten by the absorber)\n"
		       "  10       %*li (absorbed at the detector top)             \t%*li (absorbed at the detector top)             \t%*li (absorbed at the detector top)\n"
		       "  11       %*li (lost in a slit at the closure of the trap)\t%*li (lost in a slit at the closure of the trap)\t%*li (lost in a slit at the closure of the trap)\n"
		       "  12       %*li (reached the UCN detector)                 \t%*li (reached the UCN detector)                 \t%*li (reached the UCN detector)\n"
		       "  99       %*li (entered forbidden space)                  \t%*li (entered forbidden space)                  \t%*li (entered forbidden space)\n\n",		       
		       (iMC - 1 + 2 * decay.counter),
		       (iMC- 1), decay.counter, decay.counter,
		       4, kennz0[1], 4, kennz0[2], 4, kennz0[0],
		       4, kennz1[1], 4, kennz1[2], pini.xend, 4, kennz1[0], eini.xend,
		       4, kennz2[1], 4, kennz2[2], 4, kennz2[0],
		       4, kennz3[1], 4, kennz3[2], 4, kennz3[0],
		       4, kennz4[1], 4, kennz4[2], 4, kennz4[0],
		       4, kennz5[1], 4, kennz5[2], 4, kennz5[0],
		       4, kennz6[1], 4, kennz6[2], 4, kennz6[0],
		       4, kennz7[1], 4, kennz7[2], 4, kennz7[0],
		       4, kennz8[1],
		       4, kennz9[1], 4, kennz9[2], 4, kennz9[0],
		       4, kennz10[1], 4, kennz10[2], 4, kennz10[0],
		       4, kennz11[1], 4, kennz11[2], 4, kennz11[0],
		       4, kennz12[1], 4, kennz12[2], 4, kennz12[0],
		       4, kennz99[1], 4, kennz99[2], 4, kennz99[0]);
		fprintf(LOGSCR,"\nThe calculations of %li particle(s) yielded:\n"
		               "endcode    out of %i neutron(s)                            \tout of %li proton(s)                              \tout of %li electron(s)\n"
		               "   0       %*li (not categorized)                          \t%*li (not categorized)                          \t%*li (not categorized)\n"
		               "   1       %*li (did not finish)                           \t%*li (did not finish in %2.1LE s)              \t%*li (did not finish in %2.1LE s)\n"
		               "   2       %*li (hit outer boundaries)                     \t%*li (hit outer boundaries)                     \t%*li (hit outer boundaries)\n"
		               "   3       %*li (hit the walls)                            \t%*li (hit the walls)                            \t%*li (hit the walls)\n"
		               "   4       %*li (escaped to the top)                       \t%*li (escaped to the top)                       \t%*li (escaped to the top)\n"
		               "   5       %*li (escaped through the exit slit)            \t%*li (escaped through the exit slit)            \t%*li (escaped through the exit slit)\n"
		               "   6       %*li (hit the detector from the bottom)         \t%*li (hit the detector from the bottom)         \t%*li (hit the detector from the bottom)\n"
		               "   7       %*li (statistically absorbed)                   \t%*li (hit the detector from the top)            \t%*li (hit the detector from the top)\n"
		               "   8       %*li (decayed)                                  \t   -                                            	   -\n"
		               "   9       %*li (were eaten by the absorber)               \t%*li (were eaten by the absorber)               \t%*li (were eaten by the absorber)\n"
		               "  10       %*li (absorbed at the detector top)             \t%*li (absorbed at the detector top)             \t%*li (absorbed at the detector top)\n"
		               "  11       %*li (lost in a slit at the closure of the trap)\t%*li (lost in a slit at the closure of the trap)\t%*li (lost in a slit at the closure of the trap)\n"
		               "  12       %*li (reached the UCN detector)                 \t%*li (reached the UCN detector)                 \t%*li (reached the UCN detector)\n"
		               "  99       %*li (entered forbidden space)                  \t%*li (entered forbidden space)                  \t%*li (entered forbidden space)\n\n",		       
		               (iMC - 1 + 2 * decay.counter),
		               (iMC- 1), decay.counter, decay.counter,
		               4, kennz0[1], 4, kennz0[2], 4, kennz0[0],
		               4, kennz1[1], 4, kennz1[2], pini.xend, 4, kennz1[0], eini.xend,
		               4, kennz2[1], 4, kennz2[2], 4, kennz2[0],
		               4, kennz3[1], 4, kennz3[2], 4, kennz3[0],
		               4, kennz4[1], 4, kennz4[2], 4, kennz4[0],
		               4, kennz5[1], 4, kennz5[2], 4, kennz5[0],
		               4, kennz6[1], 4, kennz6[2], 4, kennz6[0],
		               4, kennz7[1], 4, kennz7[2], 4, kennz7[0],
		               4, kennz8[1],
		               4, kennz9[1], 4, kennz9[2], 4, kennz9[0],
		               4, kennz10[1], 4, kennz10[2], 4, kennz10[0],
		               4, kennz11[1], 4, kennz11[2], 4, kennz11[0],
		               4, kennz12[1], 4, kennz12[2], 4, kennz12[0],
		               4, kennz99[1], 4, kennz99[2], 4, kennz99[0]);
	}
	else // (decay.on != 2)
	{	int p3 = protneut % 3;
		// protneut % 3 = {0, 1, 2} <<<	ELECTRONS (=6) % 3 = 0
		// 								NEUTRON   (=1) % 3 = 1	
		// 								PROTON    (=2) % 3 = 2
		printf("The calculations yielded:\n"
		       "Out of %i particles:\n"
		       "  %li were not categorized (kennz=0).\n"
		       "  %li did not finish in %LG s (kennz=1).\n"
		       "  %li hit outer boundaries (kennz=2).\n"
		       "  %li hit the walls (kennz=3).\n"
		       "  %li escaped to the top (kennz=4).\n"
		       "  %li escaped through the exit slit (kennz=5).\n"
		       "  %li hit the detector from the bottom (kennz=6).\n",
		       (iMC - 1),
		       kennz0[p3],
		       kennz1[p3], xend,
		       kennz2[p3],
		       kennz3[p3],
		       kennz4[p3],
		       kennz5[p3],
		       kennz6[p3]);

		if(protneut==NEUTRON)
		{	printf("  %li were statistically absorbed (kennz=7). \n",kennz7[p3]);
		}
		else
		{	printf("  %li hit the detector from the top (kennz=7).\n",kennz7[p3]);
		}
	
		printf("  %li decayed (kennz=8).\n"
		       "  %li were eaten by the absorber (kennz=9).\n"
		       "  %li were absorbed at the detector top (kennz=10).\n"
		       "  %li were lost in a slit at the closure of the trap (kennz=11).\n"
		       "  %li reached the UCN detector (kennz=12). \n"
		       "  %li entered forbidden space.\n",
		       kennz8[p3],
		       kennz9[p3],
		       kennz10[p3],
		       kennz11[p3],
		       kennz12[p3],
		       kennz99[p3]);
	
		fprintf(LOGSCR,"The calculations yielded:\n"
		               "Out of %i particles:\n"
		               "  %li were not categorized (kennz=0).\n"
		               "  %li did not finish in %LG s (kennz=1).\n"
		               "  %li hit outer boundaries (kennz=2).\n"
		               "  %li hit the walls (kennz=3).\n"
		               "  %li escaped to the top (kennz=4).\n"
		               "  %li escaped through the exit slit (kennz=5).\n"
		               "  %li hit the detector from the bottom (kennz=6).\n",
		               (iMC - 1),
		               kennz0[p3],
		               kennz1[p3], xend,
		               kennz2[p3],
		               kennz3[p3],
		               kennz4[p3],
		               kennz5[p3],
		               kennz6[p3]);

		if(protneut==NEUTRON)
		{	fprintf(LOGSCR,"  %li were statistically absorbed (kennz=7). \n",kennz7[p3]);
		}
		else
		{	fprintf(LOGSCR,"  %li hit the detector from the top (kennz=7).\n",kennz7[p3]);
		}
	
		fprintf(LOGSCR,"  %li decayed (kennz=8).\n"
		               "  %li were eaten by the absorber (kennz=9).\n"
		               "  %li were absorbed at the detector top (kennz=10).\n"
		               "  %li were lost in a slit at the closure of the trap (kennz=11).\n"
		               "  %li reached the UCN detector (kennz=12). \n"
		               "  %li entered forbidden space.\n",
		               kennz8[p3],
		               kennz9[p3],
		               kennz10[p3],
		               kennz11[p3],
		               kennz12[p3],
		               kennz99[p3]);
	}
}
//======== end of OutputCodes ==============================================================================================


// This function takes a point, follows along a flux line in positive z-direction (up and in the r-z-plane),
// ends the line at the walls of after 5m and returns one for a detector hit and zero for no hit
int CalcFluxLine(long double r, long double phi, long double z, long double FluxStep){
	ystart[5]=0;
	BFeld(r,phi,z,500.0);
	EFeld(r,phi,z);
	Be0 = Bws;
	Bemax = Bws;
	Vflux=  0.0;                             // potential difference to the starting point
	long double AngleZ = tanl(Br/Bz), Vorzeichen = 1;        // angle of the magnetic field vector with respect to z-axis
	
	int i;
	for (i=1;i<=5000;i++){
		r = r + Vorzeichen * sinl(AngleZ) * FluxStep;
		z = z + Vorzeichen * cosl(AngleZ) * FluxStep;
		Vflux = Vflux + (sinl(AngleZ)*Er + cosl(AngleZ)*Ez)*FluxStep;    // linearised electric field strength integral along the step
		ystart[1]=r;
		ystart[3]=z;
		ystart[4]=cosl(AngleZ);
		Entkommen(ystart,  x2, H);
		if (stopall == 1)
		{
			if (kennz == KENNZAHL_DETECTOR_BOTTOM)    // 6
				return 1;
			else if ((kennz != KENNZAHL_DETECTOR_BOTTOM) && (Vorzeichen == 1))   // 6 
				Vorzeichen = -1;
			else if ((kennz != 6) && (Vorzeichen == -1))
				return 0;
		}
		BFeld(r,phi,z, 500.0);
		EFeld(r,phi,z);
		if (Bws>Bemax)
			Bemax = Bws;
		AngleZ = tanl(Br/Bz);
	}
	return 0;
}

// This function computes the critical angle for an electron with respect to the Bfield to hit the detector
long double CalcCritAngle(long double r_n, long double phi_n, long double z_n, long double Energie){
	long double frac1 =sqrtl((Bemax/Be0)-1), vpar = 0;
	gammarel =  Energie/m_e+1;
	v_n = c_0 * sqrtl(1-(1/(gammarel*gammarel)));                   // relatistic velocity from E_kin

	if (frac1 != 0 )
	vpar =v_n*sqrtl(1/(1+1/(frac1*frac1)));             // velocity parallel to flux line

	long double Epar =m_e*c_0*c_0*((1/sqrtl(1-vpar*vpar/(c_0*c_0)))-1) - Vflux; // Energy parallel to flux needed to hit det + electrostatic potential to overcome
	long double vparnew =c_0*sqrtl(1-((m_e*c_0*c_0)/(Epar+m_e*c_0*c_0))*((m_e*c_0*c_0)/(Epar+m_e*c_0*c_0)));  // corresponding velocity
	long double frac2 =fabsl((vparnew)/(v_n/sqrtl(1+frac1*frac1)));
	long double angle =atan2l(1,frac2) * 180 / pi;

	return angle;
}

// a function to calculate angle of incidence on the detector of an electron starting at r,z with vr, vz, vphi
long double CalcIncidentAngle(long double r, long double phi, long double z, long double vr, long double vphi, long double vz, long double Br0, long double Bphi0, long double Bz0, long double Bemax){
	long double B0 = sqrtl(Br0*Br0+Bphi0*Bphi0+Bz0*Bz0);
	long double vpar = (vr*Br0+vphi*Bphi0+vz*Bz0)/B0;     //projection onto flux vector
	if (vpar > 0){
        ElecAngleB = acosl((vr*Br0+vphi*Bphi0+vz*Bz0)/(B0*v_n))/conv;
        if (ElecAngleB < CritAngle) {
			gammarel =  Energie/m_e+1;
			long double vnew = c_0 * sqrtl(1-(1/(gammarel*gammarel)));
			long double vperp = sqrtl(Bemax/B0*(vnew*vnew-vpar*vpar));
			long double IncidentAngle = atan2l(vperp,sqrtl(vnew*vnew-vperp*vperp));
			return IncidentAngle/conv;
		}
	}
	ElecAngleB = 0.0;
	return 0.0;
}

// return absolute value of a 3 vector in cartesian coordinates
long double AbsValueCart(long double x, long double y, long double z)
{
	return sqrtl(powl(x,2)+powl(y,2)+powl(z,2));
}

