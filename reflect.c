// include code for one of the random number generators:
//#include "mersenne.cpp"         // members of class TRandomMersenne
#include "main.h"

// kill the particle if statistically absorbed
// v doesn't have to be positive, as it is cubed in the function
void StatAbsorbtion(long double v, long double FermiReal, long double FermiImag)
{
	long double trans = Transmission(0.5*m_n*powl(v,2)*1e9,FermiReal,FermiImag);
	long double dice = mt_get_double(v_mt_state);
	if(reflekt==0)
	{
		stopall = 1;
		kennz = 3;   // walls
		printf("Particle hit the walls (no reflection) at r=%LG z=%LG\n",ystart[1],ystart[3]);
		fprintf(LOGSCR,"Particle hit the walls (no reflection) at r=%LG z=%LG\n",ystart[1],ystart[3]);		
		return;
	}			
	else if (reflekt!=0)
	{
		if(dice<trans)
		{	
			stopall = 1;
			kennz = KENNZAHL_ABSORBED;   // 7
			printf("Statistical absorption!!! Transmission: %LG diced value %LG\n",trans,dice);
			fprintf(LOGSCR,"Statistical absorption!!! Transmission: %LG diced value %LG\n",trans,dice);
			return;
		}
	}
}

// transmission of a wall (loss of UCN) with Fermi potential  Mf (real) and PF (im) for a UCN with energy Er perpendicular to the wall (all in neV)
long double Transmission(long double Er, long double Mf, long double Pf){
return 0.2e1 * sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) / (Er + sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) + sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf));
}

// absorber with location absrmin, absrmax,abszmin,abszmax
int Absorber(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t)
{
   int AlreadyReflected;
	//printf("a");
   


   if((r>=absrmin) && (r<=absrmax))
	{
		if((z>=abszmin) && (z<=abszmax))
		{
			long double testphi = fmodl(phi/conv, 360.), prob;   // bring phi back to [0,360]
			if (testphi<0)
				testphi=360.0 + testphi;
			if((testphi>=absphimin) && (testphi<=absphimax))
			{
				long double Er = 1.0e9*0.5*m_n*(powl(vr,2));	// in  neV			
				prob = mt_get_double(v_mt_state);
				AbsProb = Transmission(Er,Mf,Pf);
				cout << "r: " << r << " z: " << z << endl; 
				cout << "vr: " << vr << " vz: " << vz << endl; 
				if((prob < AbsProb)&&(NoAbsorption==0))
				{
					AbsorberHits++;
					kennz=9;     // absorber hit
					stopall=1;
					printf("Absorber hit and absorption!\nvr = %LG m/s\nEr = %LG neV\nAbsProb = %LG\n", vr,Er, AbsProb);
					fprintf(LOGSCR,"Absorber hit and absorption!\nvr = %LG m/s\nEr = %LG neV\nAbsProb = %LG\n", vr,Er, AbsProb);
					return 1;
				}
				else if((prob >= AbsProb)&&(NoAbsorption==0))
				{
					printf("Absorber hit but no absorption!\nvr = %LG m/s\nEr = %LG neV\nAbsProb = %LG\n", vr,Er, AbsProb);
					fprintf(LOGSCR,"Absorber hit but no absorption!\nvr = %LG m/s\nEr = %LG neV\nAbsProb = %LG\n", vr,Er, AbsProb);
					//  reflexion at absorber
					if((z>=abszmin)&&(vr>0))  // outer reflection
						AlreadyReflected =ReflexOut(y,r,vr,z,vz,phi,vphi,t,DiffProb,FPrealNocado,FPimNocado);
					if((r>=absrmin+0.002)&&(vz>0)&&(!AlreadyReflected)) // reflection at bottom side
						AlreadyReflected =ReflexTop(y,r,vr,z,vz,phi,vphi,t,DiffProb,FPrealNocado,FPimNocado);
					if((r>=absrmin+0.002)&&(vz<0)&&(!AlreadyReflected)) // reflection at bottom side
						AlreadyReflected =ReflexBottom(y,r,vr,z,vz,phi,vphi,t,DiffProb,FPrealNocado,FPimNocado);				
					NoAbsorption = 1;
					AbsorberHits++;		
					return 1;
				}
			}
		}
	}
  
	NoAbsorption = 0;
	return 0;
     	
}
// Ende Absorber



// changes the vector ystart according to reflexion at a horizontal top boundary with FermiReal+i*FermiImag and probability for diffuse reflexion ProbDiffuse
int ReflexTop(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,long double ProbDiffuse, long double FermiReal, long double FermiImag)
{
	if((vz>0)) // only reflect particles moving up
	{					
	long double winkeben;     // Winkel in Spiegelebene
   long double winksenkr;  // Winkel senkrecht dazu, 0 in Ebene, sin verteilt, Lambert law around normal
   long double vabs; 
	
	diffuprob = mt_get_double(v_mt_state); // dice the diffuse reflection probability
	//printf("\nDIFFUPROB %.12LG\n",diffuprob);
	//fprintf(LOGSCR,"\nDIFFUPROB %.12LG\n",diffuprob);
	
	// specular reflexion
	if ((diffuse == 1)||(diffuse==3)&&(diffuprob>=ProbDiffuse))
	{
		StatAbsorbtion(vz,FermiReal,FermiImag);
		if (stopall==1) return 2;
	   y[4]= - vz;
		//cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl; 
		//cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		//cout << "r: " << r << "z: " << z << endl; 
		printf("pol %d t= %LG top Erefl= %LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z);
		fprintf(LOGSCR,"pol %d t= %LG top Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z);    
		nrefl++;
		if(reflektlog == 1)
			fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG  %LG 0 0 %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}
			
	// diffuse reflexion
	if ((diffuse == 2)||(diffuse == 3)&&(diffuprob<ProbDiffuse))
	{
		winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi);
		winksenkr = acosl( cosl(0.0) - mt_get_double(v_mt_state) * (cosl(0.0) - cosl(0.5 * pi)));
		vabs = sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2));
		StatAbsorbtion(vz,FPrealCu,FPimCu);
		if (stopall==1) return 2;
	
		y[2] = sinl(winkeben) * cosl(winksenkr) * vabs;
		y[4] = -sinl(winksenkr) * vabs;
		y[6] = cosl(winkeben) * cosl(winksenkr) * vabs / r;
		//cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;  
		//cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		// cout << "r: " << r << "z: " << z << endl; 
		printf("pol %d t= %LG top: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
		fprintf(LOGSCR,"pol %d t= %LG top: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
		nrefl++;
		if(reflektlog == 1)
			fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG  %LG %LG %LG %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,winkeben,winksenkr,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}		
}
	return 1;
}

// changes the vector ystart according to reflexion at a horizontal bottom boundary with FermiReal+i*FermiImag and probability for diffuse reflexion ProbDiffuse
int ReflexBottom(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,long double ProbDiffuse, long double FermiReal, long double FermiImag)
{
if((vz<0)) // only reflect particles moving down
	{	
	long double winkeben;     // Winkel in Spiegelebene
   long double winksenkr;  // Winkel senkrecht dazu, 0 in Ebene, sin verteilt, Lambert law around normal
   long double vabs; 
	
	diffuprob = mt_get_double(v_mt_state); // dice the diffuse reflection probability
	//printf("\nDIFFUPROB %.12LG\n",diffuprob);
	//fprintf(LOGSCR,"\nDIFFUPROB %.12LG\n",diffuprob);
	
	// specular reflexion
	if ((diffuse == 1)||(diffuse==3)&&(diffuprob>=ProbDiffuse))
	{
		StatAbsorbtion(vz,FermiReal,FermiImag);
		if (stopall==1) return 2;
		  
		y[4]= - vz;
		//cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
      // cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		// cout << "r: " << r << "z: " << z << endl; 
		printf("pol %d t= %LG bottom Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z);
		fprintf(LOGSCR,"pol %d t= %LG bottom Erefl=%LG neV r= %LG z=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z);    
		nrefl++;
       if(reflektlog == 1)
			fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG  %LG 0 0 %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}
			
	// diffuse reflexion
	if ((diffuse == 2)||(diffuse == 3)&&(diffuprob<ProbDiffuse))
	{
		winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi);
		winksenkr = acosl( cosl(0.0) - mt_get_double(v_mt_state) * (cosl(0.0) - cosl(0.5 * pi)));
		vabs = sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2));
		StatAbsorbtion(vz,FermiReal,FermiImag);
			if (stopall==1) return 2;
				
        y[2] = sinl(winkeben) * cosl(winksenkr) * vabs;
        y[4] = sinl(winksenkr) * vabs;
        y[6] = cosl(winkeben) * cosl(winksenkr) * vabs / r;
       //cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
       //cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		 //cout << "r: " << r << "z: " << z << endl; 
	   printf("pol %d t= %LG botttom: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       fprintf(LOGSCR,"pol %d t= %LG botttom: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vz,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       nrefl++;
       if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG  %LG %LG %LG %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,winkeben,winksenkr,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}	
}
	return 1;
}

// changes the vector ystart according to reflexion at a vertical outer boundary with FermiReal+i*FermiImag and probability for diffuse reflexion ProbDiffuse
int ReflexOut(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,long double ProbDiffuse, long double FermiReal, long double FermiImag)
{
	if((vr>0)) // only reflect particles moving outward
	{
	long double winkeben;     // Winkel in Spiegelebene
   long double winksenkr;  // Winkel senkrecht dazu, 0 in Ebene, sin verteilt, Lambert law around normal
   long double vabs; 
	
	diffuprob = mt_get_double(v_mt_state); // dice the diffuse reflection probability
	//printf("\nDIFFUPROB %.12LG\n",diffuprob);
	//fprintf(LOGSCR,"\nDIFFUPROB %.12LG\n",diffuprob);
	
	// specular reflexion
	if ((diffuse == 1)||(diffuse==3)&&(diffuprob>=ProbDiffuse))
	{
		StatAbsorbtion(vr,FermiReal,FermiImag);
		if (stopall==1) return 2;
		y[2]= - vr;
		// cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
      // cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		// cout << "r: " << r << "z: " << z << endl; 
		printf("pol %d t= %LG outer Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z);
	   fprintf(LOGSCR,"pol %d t= %LG outer Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z);
	   nrefl++;
	   if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG  %LG 0 0 %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}
			
	// diffuse reflexion
	if ((diffuse == 2)||(diffuse == 3)&&(diffuprob<ProbDiffuse))
	{
		winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi);
		winksenkr = acosl( cosl(0.0) - mt_get_double(v_mt_state) * (cosl(0.0) - cosl(0.5 * pi)));
		vabs = sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2));
		StatAbsorbtion(vr,FermiReal,FermiImag);
			if (stopall==1) return 2;
				
		 y[2] = - sinl(winksenkr) * vabs;
       y[4] = cosl(winkeben) * cosl(winksenkr) * vabs;
       y[6] = sinl(winkeben) * cosl(winksenkr) * vabs / r;
       // cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
       //cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		 //cout << "r: " << r << "z: " << z << endl; 
	   printf("pol %d t= %LG outer: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       fprintf(LOGSCR,"pol %d t= %LG outer: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       nrefl++;
	if(reflektlog == 1)
	  fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG  %LG %LG %LG %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,winkeben,winksenkr,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}	
}
	return 1;
}
	
// changes the vector ystart according to reflexion at a vertical inner boundary with FermiReal+i*FermiImag and probability for diffuse reflexion ProbDiffuse
int ReflexIn(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,long double ProbDiffuse, long double FermiReal, long double FermiImag)
{
	if((vr<0)) // only reflect particles moving innward
	{
	long double winkeben;     // Winkel in Spiegelebene
   long double winksenkr;  // Winkel senkrecht dazu, 0 in Ebene, sin verteilt, Lambert law around normal
	long double vabs; 
	
	diffuprob = mt_get_double(v_mt_state); // dice the diffuse reflection probability
	//printf("\nDIFFUPROB %.12LG\n",diffuprob);
	//fprintf(LOGSCR,"\nDIFFUPROB %.12LG\n",diffuprob);
	
	// specular reflexion
	if ((diffuse == 1)||(diffuse==3)&&(diffuprob>=ProbDiffuse))
	{
		StatAbsorbtion(vr,FermiReal,FermiImag);
		if (stopall==1) return 2;
       y[2]= - vr;
		//cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
       //cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		 //cout << "r: " << r << "z: " << z << endl; 
       printf("pol %d t= %LG outer Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z);
	   fprintf(LOGSCR,"pol %d t= %LG outer Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z);
	   nrefl++;
	   if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG  %LG 0 0 %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}
			
	// diffuse reflexion
	if ((diffuse == 2)||(diffuse == 3)&&(diffuprob<ProbDiffuse))
	{
		winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi);
		winksenkr = acosl( cosl(0.0) - mt_get_double(v_mt_state) * (cosl(0.0) - cosl(0.5 * pi)));
		vabs = sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2));
		StatAbsorbtion(vr,FermiReal,FermiImag);
			if (stopall==1) return 2;
		
		y[2] = sinl(winksenkr) * vabs;
       y[4] = cosl(winkeben) * cosl(winksenkr) * vabs;
       y[6] = sinl(winkeben) * cosl(winksenkr) * vabs /r;
		//cout << "vabsvorher" << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
      // cout << "vabsnachher" << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		// cout << "r: " << r << "z: " << z << endl; 
       printf("pol %d t= %LG inner: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       fprintf(LOGSCR,"pol %d t= %LG inner: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*powl(vr,2)*1e9,r,z,winkeben/conv,winksenkr/conv);
       nrefl++;
       if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG  %LG %LG %LG %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,winkeben,winksenkr,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}	
}
	return 1;
}

// changes the vector ystart according to reflexion at an inclined boundary (truncated cone) with FermiReal+i*FermiImag and probability for diffuse reflexion ProbDiffuse
int ReflexGeneral(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t,long double ProbDiffuse, long double FermiReal, long double FermiImag,long double rNorm,long double zNorm,long double phiNorm)
{
	long double rNormNorm=rNorm/sqrtl(powl(rNorm,2)+powl(zNorm,2)); // normalising
	long double zNormNorm=zNorm/sqrtl(powl(rNorm,2)+powl(zNorm,2)); // normalising
	long double vsenkr=(vr*rNormNorm+vz*zNormNorm);	
	//if(vsenkr>=0)
		//cin >> blankint;
	
	if((vsenkr<0)) // only reflect particles moving towards boundary
	{
	//cout << "r: " << r << " z: " << z << endl; 
	//cout << "rNorm: " << rNorm << " zNorm: " << zNorm << endl; 
	//cout << "vsenkr: " << vsenkr << " vr: " << vr << " vz: " << vz <<endl;  
	long double winkeben;     // Winkel in Spiegelebene
   long double winksenkr;  // Winkel senkrecht dazu, 0 in Ebene, sin verteilt, Lambert law around normal
	long double vabs;
	long double alphainc=atan2l(zNorm,rNorm)+pi;		// atan2l: return angle between -pi and pi, then plus pi, because first write-down of this function used wrong angles and +pi corrects this
		 cout << "alphainc: "<< alphainc << endl; 
		 //cin >> blankint;
		
	// specular reflexion
	if ((diffuse == 1)||(diffuse==3)&&(diffuprob>=ProbDiffuse))
	{
		StatAbsorbtion(vsenkr,FermiReal,FermiImag);
		if (stopall==1) return 2;
		y[2] = vr - 2 * vsenkr * rNormNorm;
		y[4] = vz - 2 * vsenkr * zNormNorm;
		//cout << "vabsvorher " << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;	
       //cout << "vabsnachher " << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		 //cout << "r: " << r << "z: " << z << endl; 
		//cin >> blankint;
       printf("pol %d t= %LG inclined Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*vsenkr*vsenkr*1e9,r,z);
	   fprintf(LOGSCR,"pol %d t= %LG inclined Erefl=%LG neV r=%LG z=%LG\n",polarisation,t,0.5*m_n*vsenkr*vsenkr*1e9,r,z);
	   nrefl++;
	   if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG  %LG 0 0 %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}
	
	// diffuse reflection
	if ((diffuse == 2)||(diffuse == 3)&&(diffuprob<ProbDiffuse))
	{
		winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi);
		winksenkr = acosl( cosl(0.0) - mt_get_double(v_mt_state) * (cosl(0.0) - cosl(0.5 * pi)));
		vabs = sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2));
		StatAbsorbtion(vsenkr,FermiReal,FermiImag);
			if (stopall==1) return 2;
		long double vrvirt = -1.0 * sinl(winksenkr) * vabs;     // values without rotating alphainc
       long double vzvirt = cosl(winkeben) * cosl(winksenkr) * vabs;  // values without rotating alphainc
		 y[2] = vrvirt *  cosl(alphainc) - vzvirt * sinl(alphainc);
       y[4] = vrvirt *  sinl(alphainc) + vzvirt * cosl(alphainc);
       y[6] = sinl(winkeben) * cosl(winksenkr) * vabs / r;
		//cout << "vabsvorher " << sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)) << endl;
       //cout << "vabsnachher " << sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2)) << endl;
		 //cout << "r: " << r << "z: " << z << endl; 
		 //cin >> blankint; 
       printf("pol %d t= %LG inclined: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*vsenkr*vsenkr*1e9,r,z,winkeben/conv,winksenkr/conv);
       fprintf(LOGSCR,"pol %d t= %LG inclined: Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,t,0.5*m_n*vsenkr*vsenkr*1e9,r,z,winkeben/conv,winksenkr/conv);
       nrefl++;
       if(reflektlog == 1)
		fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG  %LG %LG %LG %LG %LG %LG %LG %LG\n",t,r,z,phi,r*cosl(phi),r*sinl(phi),sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)),H,0.5*m_n*powl(vr,2)*1e9,winkeben,winksenkr,vr,vz,r*vphi,vphi,sqrtl(powl(y[2],2)+powl(y[1],2)*powl(y[6],2)+powl(y[4],2))-sqrtl(powl(vr,2)+powl(r,2)*powl(vphi,2)+powl(vz,2)));
	}	
	
							

	return 1;
	}
	else return 0;
}








// ALTE ROUTINE!!!!!!!!
//-------------------------------------------------------------------------------
//Reflektionsroutine für Besselfeld mit 90° Ecke
void ReflektEckig(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t)
{
   long double rmi, rma;                       //Abstand von der Ecke skaliert
   bool schlitz = false;   
   //cout << vr << vz << endl;
   rmi=rmin+wandinnen;
   rma=rmax-wanddicke;		
  //Absorber(y,r,vr,z,vz,phi,vphi,t);
  if (kennz == 9) // absorber hit
        return;

   if (slit){
	if((z < (r-rma+LueckeR)*(LueckeZ/LueckeR)) && (r <= rmax) && (z <= LueckeZ) && (r >= rma-LueckeR)){    // Spalt im Boden
		schlitz = true;		// do not reflect at walls
			if((vz<0)&&(vr>0)){    // only consider lost in schlitz when moving towards walls
				kennz=5;        // Teilchen im Füllschlitz                                                                              // Nachweis der
				stopall=1;
				return;
			}                                                     // dep. Neutronen
	}
    }
    //Innenwand
   if( (r<= rmi+epsi) && (vr<0.0) && (!schlitz)){		
		blankint = ReflexIn(y,r,vr,z,vz,phi,vphi,t,0,FPrealNocado,FPimNocado);
		cout << "NEWREFLEXION!";
   }
   //Außenwand
   if( (r>= rma-epsi) && (vr>0.0)  && (!schlitz)){		
		blankint = ReflexOut(y,r,vr,z,vz,phi,vphi,t,0,FPrealNocado,FPimNocado);
		cout << "NEWREFLEXION!";
    }
   //Boden
   if((z <= wanddicke+epsi) && (vz < 0)  && (!schlitz))
	{        
	   // simulates a slit at the closure of the trap at the bottom (AbsorberExp)
		if((r<0.035)&&(r>0.0345))
		{
			kennz=11;    // escape through bottom slit in AbEx
			stopall=1;
			printf("UCN was lost in a slit at the bottom !!\n");
			fprintf(LOGSCR,"UCN was lost in a slit at the bottom !!\n");
			return;	
		}		
		 blankint =  ReflexBottom(y,r,vr,z,vz,phi,vphi,t,0,FPrealNocado,FPimNocado);
		 cout << "NEWREFLEXION!";	
   }   
   // reflection at the (horizontal) lid
   if((z >= hlid-epsi) && (vz > 0)  && (!schlitz))
	{
		 blankint =  ReflexTop(y,r,vr,z,vz,phi,vphi,t,0,FPrealCu,FPimCu);
		 cout << "NEWREFLEXION!";	
	}
}
//-------------------------------------------------------------------------------
//diffuse reflection at the walls
void ReflektDiffuse(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t)
{
   long double rmi, rma;                       //Abstand von der Ecke skaliert
   bool schlitz = false;
   //cout << vr << vz << endl;
   rmi=rmin+wandinnen;
   rma=rmax-wanddicke;
   //Absorber(y,r,vr,z,vz,phi,vphi,t);
    if (kennz == 9)   // absorber hit
        return;   
  if (slit){
	if((z < (r-rma+LueckeR)*(LueckeZ/LueckeR)) && (r <= rmax) && (z <= LueckeZ) && (r >= rma-LueckeR)){    // Spalt im Boden
		schlitz = true;		// do not reflect at walls
		if((vz<0)&&(vr>0))    // only consider lost in schlitz when moving towards walls
		{
	                                                                                             // aussen für
			//if((r>rma)||(z<wanddicke))
			//{
			kennz=5;        // Teilchen im Füllschlitz                                                                              // Nachweis der
			stopall=1;
			//}                                                                       							 // dep. Neutronen
		}
	}
    }    
   if( (r<= rmi+epsi) && (vr<0.0) && (!schlitz)){           //Innenwand		 
		 blankint =  ReflexIn(y,r,vr,z,vz,phi,vphi,t,1,FPrealNocado,FPimNocado);
		 cout << "NEWREFLEXION!";	
   }	
	 //Außenwand
   if( (r>= rma-epsi) && (vr>0.0)  && (!schlitz)){          //Außenwand		
		blankint =  ReflexOut(y,r,vr,z,vz,phi,vphi,t,1,FPrealNocado,FPimNocado);
		 cout << "NEWREFLEXION!";		
		}
	//Boden
   if((z <= wanddicke+epsi) && (vz < 0) && (!schlitz))
	{    		
		// simulates a slit at the closure of the trap at the bottom (AbsorberExp)
		if((r<0.035)&&(r>0.0345))
		{
			kennz=11;    // escape through bottom slit in AbEx
			stopall=1;
			printf("UCN was lost in a slit at the bottom !!\n");
			fprintf(LOGSCR,"UCN was lost in a slit at the bottom !!\n");
			return;	
		}		
		 blankint = ReflexBottom(y,r,vr,z,vz,phi,vphi,t,1,FPrealNocado,FPimNocado);
		 cout << "NEWREFLEXION!";
   }	
	// reflection at the (horizontal) lid
   if((z >= hlid-epsi) && (vz > 0) && (!schlitz))
	{
		 blankint = ReflexTop(y,r,vr,z,vz,phi,vphi,t,1,FPrealNocado,FPimNocado);
		 cout << "NEWREFLEXION!";	
	}	 	
	return;
	}
