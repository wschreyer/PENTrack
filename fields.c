#include <list>
#include <string>
#include <cmath>

using namespace std;

#include "globals.h"
#include "fields.h"
#include "racetrack.h"
#include "2dinterpolfeld.h"
//#include "BForbes.h"
#include "maplefield.h"

void SwitchField(long double t);
void Banalytic(long double r,long double phi,long double z, long double t, long double Bi[3], long double dBidrj[3][3]);


int bfeldwahl = 0, Feldcount = 0;
long double BFeldSkalGlobal = 1.0, BFeldSkal = 1.0, EFeldSkal = 1.0;

// incorporate B-fieldoszillations into the code
int FieldOscillation = 0;        // turn field oscillation on if 1
long double OscillationFraction = 1e-4, OscillationFrequency = 1;    // Frequency in Hz

list<TabField*> fields;
long double dBrdr, dBrdz, dBrdphi=0.0, dBphidr, dBphidz, dBphidphi=0.0, dBzdr, dBzdz, dBzdphi=0.0;
long double Bws,dBdr,dBdz,dBdphi,Br,Bz,Bphi; //B-field: derivations in r, z, phi, Komponenten r, z, Phi
long double Ez,Er, Ephi, dErdr, dErdz, dEzdr, dEzdz, dEphidr, dEphidz;    // electric field

long double BCutPlanePoint[3], BCutPlaneNormalAlpha, BCutPlaneNormalGamma, BCutPlaneSampleDist;	// plane for B field cut
int BCutPlaneSampleCount;


// read table files or conductor file
void PrepareFields(list<string> &fieldvaltab){
	if ((bfeldwahl == 0 || bfeldwahl == 2) && (BFeldSkalGlobal != 0 || EFeldSkal != 0))
	{
		TabField *tf;
		for (list<string>::iterator i = fieldvaltab.begin(); i != fieldvaltab.end(); i++){
			tf = new TabField(i->c_str());
			fields.push_back(tf);
		}
	}
	else if(bfeldwahl==4)
	{
//		ReadMagnets((inpath+"/coils.cond").c_str());
	
		printf("\n \n Test of integration\n");
		//	long double TestInt;
		BFeldSkal=1.0;
		for (int a = 0;a<1;a++)
		{
		
			BFeld(0.3,0,0.1, 500.0);
			printf("T\n");
		}
		printf("Br = %.17LG \n",Br);
		printf("dBrdr = %.17LG \n",dBrdr);
		printf("dBrdz = %.17LG \n",dBrdz);
		printf("Bz = %.17LG \n",Bz);
		printf("dBzdr = %.17LG \n",dBzdr);
		printf("dBzdz = %.17LG \n",dBzdz);	
	}	
}


// free memory
void FreeFields(){
	for (list<TabField*>::iterator i = fields.begin(); i != fields.end(); i++)
		delete (*i);	
}


// fill global B field variables with values
void BFeld (long double rloc, long double philoc, long double zloc, long double t){      //B-Feld am Ort des Teilchens berechnen
	long double Bi[3] = {0.0, 0.0, 0.0};
	long double dBidrj[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
	
	SwitchField(t); // set BFeldSkal for different experiment phases	
		
	if (bfeldwahl != 1){
		long double Brtemp, Bztemp, dBdrtemp, dBdztemp;
		switch (bfeldwahl)
		{
			
			case 0:	if (BFeldSkal != 0){
						for (list<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){
							if ((*i)->BInterpol(rloc, zloc, Bi, dBidrj))
								break;
						}
						for (int i = 0; i < 3; i++){
							Bi[i] *= BFeldSkal;
							for (int j = 0; j < 3; j++){
								dBidrj[i][j] *= BFeldSkal;
							}
						}												
					}
					RacetrackField(rloc,philoc,zloc,Bi,dBidrj);
					break;
					
			case 2:	if (BFeldSkal != 0){
						for (list<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){ 
							if ((*i)->BInterpol(rloc, zloc, Bi, dBidrj))
								break;
						}				
						for (int i = 0; i < 3; i++){
							Bi[i] *= BFeldSkal;
							for (int j = 0; j < 3; j++){
								dBidrj[i][j] *= BFeldSkal;
							}
						}												
					}
					RacetrackField(rloc,philoc,zloc, Bi, dBidrj);
					Brtemp = Br; Bztemp = Bz; dBdrtemp = dBdr; dBdztemp = dBdz;
					Bwsomaple(rloc, philoc, zloc, Ibar, Bi, dBidrj);
					printf("BrI BrM deltaBr %17LG %17LG %17LG \n", Brtemp, Br, (Brtemp-Br)/Br);
					printf("BzI BzM deltaBz %17LG %17LG %17LG \n", Bztemp, Bz, (Bztemp-Bz)/Bz);						
					Br=(Brtemp-Br)/Br; Bz=(Bztemp-Bz)/Bz; dBdr=(dBdrtemp-dBdr)/dBdr; dBdz=(dBdztemp-dBdz)/dBdz;
					break;	
						
			case 3: Banalytic(rloc, philoc, zloc, t, Bi, dBidrj);
					break;
					
			case 4: //BForbes(rloc, philoc, zloc, t, Bi, dBidrj);
					RacetrackField(rloc,philoc,zloc, Bi, dBidrj);
					break;
		}   
	}
	
	Br = Bi[0];
	Bphi = Bi[1];
	Bz = Bi[2];
	dBrdr = dBidrj[0][0];
	dBrdphi = dBidrj[0][1];
	dBrdz = dBidrj[0][2];
	dBphidr = dBidrj[1][0];
	dBphidphi = dBidrj[1][1];
	dBphidz = dBidrj[1][2];
	dBzdr = dBidrj[2][0];
	dBzdphi = dBidrj[2][1];
	dBzdz = dBidrj[2][2];
	Bws = sqrt(Br*Br+Bz*Bz+Bphi*Bphi);			
	if (Bws>1e-31)
	{			
		dBdr   = (Br*dBrdr + Bphi*dBphidr + Bz*dBzdr)  /Bws;
		dBdz   = (Br*dBrdz + Bphi*dBphidz + Bz*dBzdz)  /Bws;
		dBdphi = (Br*dBrdphi + Bphi*dBphidphi + Bz*dBzdphi)/Bws;
	}
	else
	{
		dBdr = 0;
		dBdz = 0;
		dBdphi = 0;
	}	
	
	Feldcount++;
	
	return;
}


void SwitchField(long double t){        //B-Feld nach Wunsch skalieren, BFeldSkal wird an alle B Komp und deren Ableitungen multipliziert		
		if ((t < FillingTime)&&(FillingTime>0))
		{      // filling in neutrons
			BFeldSkal = 0;
		}
		
		else if ((t>=FillingTime)&&(t < (CleaningTime+FillingTime)) && (CleaningTime>0))
		{      // spectrum cleaning
			BFeldSkal = 0;
		}
		else if ((RampUpTime>0)&&(t >= CleaningTime+FillingTime) && (t < (RampUpTime+CleaningTime+FillingTime)) && (RampUpTime > 0))
		{   // ramping up field
			//linear: BFeldSkal = (t-CleaningTime+1e-15)/RampUpTime; // ramp field
			//now smoothly with a cosine 
			BFeldSkal=(0.5-0.5*cosl(pi*(t-CleaningTime-FillingTime)/RampUpTime)) * BFeldSkalGlobal ;
		}
		else if ((FullFieldTime>0)&&(t >= (RampUpTime+CleaningTime+FillingTime)) && (t < (RampUpTime+CleaningTime+FullFieldTime+FillingTime)))
		{  // storage time
			BFeldSkal = BFeldSkalGlobal;
		}
		else if ((t >= (RampUpTime+CleaningTime+FillingTime+FullFieldTime)) && (t < (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime)) && (RampDownTime != 0))
		{   // ramping down field
			//linear: BFeldSkal = 1- (t-(RampUpTime+CleaningTime+FullFieldTime)) / RampDownTime + 1e-15;  // ramp down field
			//now smoothly with a cosine 
			BFeldSkal=(0.5+0.5*cosl(pi*(t-(RampUpTime+CleaningTime+FillingTime+FullFieldTime)) / RampDownTime)) * BFeldSkalGlobal;
		}
		else if (t >=  (RampUpTime+CleaningTime+FillingTime+FullFieldTime+RampDownTime))
		{      // emptying of neutrons into detector
			BFeldSkal = 0;
		}
		else
			BFeldSkal = BFeldSkalGlobal;
	
	
	if (FieldOscillation==1){
		BFeldSkal = BFeldSkal* (1+(BFeldSkal*OscillationFraction*sinl(OscillationFrequency*2*pi*t)));		
	}
	return;
}

// produce linearly rising B-field in z-direction and linearly changing in time
void Banalytic(long double r,long double phi,long double z, long double t, long double Bi[3], long double dBidrj[3][3]){
	long double dBx = 5;
	long double x0 = 1.1;
	long double alpha = 0.5;
	long double beta = 0.025;
	
	
	// transform Bfield to higher z-values
	long double offset = -0.0;
	z = z+offset;
	
	Bi[0] += (-dBx*(x0*z-0.5*z*z)+6)*(alpha+beta*t);
    dBidrj[0][2] += -dBx*(alpha+beta*t)*(x0-z);
}



// fill global E field variables with values
void EFeld(long double rloc, long double philoc, long double zloc){
	long double Ei[3] = {0.0, 0.0, 0.0};
	long double dEidrj[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };

	if (EFeldSkal != 0 && (bfeldwahl == 0 || bfeldwahl == 2)){
		for (list<TabField*>::iterator i = fields.begin(); i != fields.end(); i++){ 
			if ((*i)->EInterpol(rloc, zloc, Ei, dEidrj))
				break;
/*			else if (++i == fields.end()){
				Log("\nThe particle has left fieldval boundaries: r=%LG, z=%LG! Stopping particle...\n", rloc, zloc);
				kennz = KENNZAHL_LEFT_FIELD;
				stopall = 1;
			}	*/	
			for (int i = 0; i < 3; i++){
				Ei[i] *= EFeldSkal;
				for (int j = 0; j < 3; j++){
					dEidrj[i][j] *= EFeldSkal;
				}
			}												
		}				
		
    }
	Er = Ei[0];
	Ephi = Ei[1];
	Ez = Ei[2];
	dErdr = dEidrj[0][0];
//	dErdphi = dEidrj[0][1];
	dErdz = dEidrj[0][2];
	dEphidr = dEidrj[1][0];
//	dEphidphi = dEidrj[1][1];
	dEphidz = dEidrj[1][2];
	dEzdr = dEidrj[2][0];
//	dEzdphi = dEidrj[2][1];
	dEzdz = dEidrj[2][2];

	Feldcount++;	
	return;
}


// print cut through BField to file
void PrintBFieldCut(const char *outfile){
	// transform plane parameters to cartesian coords
	long double P[3] = {BCutPlanePoint[0]*cosl(BCutPlanePoint[1]), BCutPlanePoint[0]*sinl(BCutPlanePoint[1]), BCutPlanePoint[2]};
	long double n[3] = {cosl(BCutPlanePoint[1]+BCutPlaneNormalAlpha)*sinl(BCutPlaneNormalGamma), 
						sinl(BCutPlanePoint[1]+BCutPlaneNormalAlpha)*sinl(BCutPlaneNormalGamma), 
						cosl(BCutPlaneNormalGamma)};

	// project x/y/z-axes on plane for direction vectors u,v
	long double u[3], v[3];
	if (n[0] < 0.9){		// n not parallel to x-axis
		u[0] = 1-n[0]*n[0];
		u[1] = -n[0]*n[1];
		u[2] = -n[0]*n[2];
	}
	else{					// if n parallel to x-axis use z-axis for u
		u[0] = -n[2]*n[0];
		u[1] = -n[2]*n[1];
		u[2] = 1-n[2]*n[1];
	}
	if (n[1] < 0.9){		// n not parallel to y-axis
		v[0] = -n[1]*n[0];
		v[1] = 1-n[1]*n[1];
		v[2] = -n[1]*n[2];
	}
	else{					// if n parallel to y-axis use z-axis for v
		v[0] = -n[2]*n[0];
		v[1] = -n[2]*n[1];
		v[2] = 1-n[2]*n[1];
	}
	
	int i,j,k;
	// normalize u,v
	long double uabs = sqrtl(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
	long double vabs = sqrtl(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	for (i = 0; i < 3; i++){
		u[i] /= uabs;
		v[i] /= vabs;
	}
	
	// print BField to file
	FILE *cutfile = fopen(outfile, "w");
	if (!cutfile){
		Log("Could not open %s!",outfile);
		exit(-1);
	}
	fprintf(cutfile, "r phi z x y Br Bphi Bz Bx By dBrdr dBrdphi dBrdz dBphidr dBphidphi dBphidz dBzdr dBzdphi dBzdz Babs dBdr dBdphi dBdz Er Ephi Ez Ex Ey dErdr dErdz dEphidr dEphidz dEzdr dEzdz\n");
	
	long double Pp[3],r,phi,Bx,By,Ex,Ey;
	float start = clock();
	for (i = -BCutPlaneSampleCount/2; i <= BCutPlaneSampleCount/2; i++) {
		for (j = -BCutPlaneSampleCount/2; j <= BCutPlaneSampleCount/2; j++){
			for (k = 0; k < 3; k++)
				Pp[k] = P[k] + i*BCutPlaneSampleDist*u[k] + j*BCutPlaneSampleDist*v[k];
			r = sqrtl(Pp[0]*Pp[0]+Pp[1]*Pp[1]);
			phi = atan2(Pp[1],Pp[0]);
			BFeld(r, phi, Pp[2], 0);
			Bx = Br*cosl(phi) - Bphi*sinl(phi);
			By = Br*sinl(phi) + Bphi*cosl(phi);
			fprintf(cutfile, "%LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG ",
							  r,phi,Pp[2],Pp[0],Pp[1],Br,Bphi,Bz,Bx,By,dBrdr,dBrdphi,dBrdz,dBphidr,dBphidphi,dBphidz,dBzdr,dBzdphi,dBzdz,Bws,dBdr,dBdphi,dBdz);
			EFeld(r,phi,Pp[2]);
			Ex = Er*cos(phi) - Ephi*sin(phi);
			Ey = Er*sin(phi) + Ephi*cos(phi);
			fprintf(cutfile, "%LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
							  Er,Ephi,Ez,Ex,Ey,dErdr,dErdz,dEphidr,dEphidz,dEzdr,dEzdz);
		}
	}
	start = (clock() - start)/CLOCKS_PER_SEC;
	fclose(cutfile);
	printf("Called BFeld and EFeld %u times in %fs (%fms per call)",BCutPlaneSampleCount*BCutPlaneSampleCount, start, start/BCutPlaneSampleCount/BCutPlaneSampleCount);
}


// calculate heating of the neutron gas by field rampup
void PrintBField(const char *outfile){
	// print BField to file
	FILE *bfile = fopen(outfile, "w");
	if (!bfile){
		Log("Could not open %s!",outfile);
		exit(-1);
	}

	fprintf(bfile,"r phi z Br Bphi Bz 0 0 Babs\n");
	BFeldSkal = 1.0;
	long double rmin = 0.12, rmax = 0.5, zmin = 0, zmax = 1.2;
	int E;
	const int Emax = 108;
	long double dr = 0.1, dz = 0.1;
	long double VolumeB[Emax];
	for (E = 0; E <= Emax; E++) VolumeB[E] = 0;
	
	long double EnTest;
	for (long double r = rmin; r <= rmax; r += dr){
		for (long double z = zmin; z <= zmax; z += dz){
			BFeld(r, 0, z, 500.0);
			fprintf(bfile,"%LG %G %LG %LG %LG %LG %G %G %LG \n",r,0.0,z,Br,Bphi,Bz,0.0,0.0,Bws);
			printf("r=%LG, z=%LG, Br=%LG T, Bz=%LG T\n",r,z,Br, Bz);
			
			// Ramp Heating Analysis
			for (E = 0; E <= Emax; E++){
				EnTest = E*1.0e-9 - m_n*gravconst*z - mu_nSI/ele_e * Bws;
				if (EnTest >= 0){
					// add the volume segment to the volume that is accessible to a neutron with energy Energie
					VolumeB[E] = VolumeB[E] + pi * dz * ((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
				}
			}
		}
	}

	// for investigating ramp heating of neutrons, volume accessible to neutrons with and
	// without B-field is calculated and the heating approximated by thermodynamical means
	Log("\nEnergie [neV], Volumen ohne B-Feld, mit B-Feld, 'Erwaermung'");
	long double Volume;
	for (E = 0; E <= Emax; E++) 
	{
		Volume = ((E * 1.0e-9 / (m_n * gravconst))) * pi * (rmax*rmax-rmin*rmin);
		// isentropische zustandsnderung, kappa=5/3
		Log("\n%i %.17LG %.17LG %.17LG",E,Volume,VolumeB[E],E * powl((Volume/VolumeB[E]),(2.0/3.0)) - E);
	}
}


