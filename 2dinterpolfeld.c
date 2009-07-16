// this file contains functions to read in E and B fields, interpolate
// them and calculate derivatives

#include "main.h"
FILE *FIN = NULL;
int fehler, indr, indz ;          // indr, indz: current indices for field interpolation
long double r_mi, r_ma, z_mi, z_ma; // minimum and maximum values, counter for field calls
long double *rind = NULL, *zind = NULL, **BrTab = NULL, **BzTab = NULL,**BphiTab = NULL,**BrTab1 = NULL,**BzTab1 = NULL,**BphiTab1 = NULL;  //B Arrays
long double **BrTab2 = NULL, **BzTab2 = NULL, **BphiTab2 = NULL, **BrTab12 = NULL, **BzTab12 = NULL, **BphiTab12 = NULL;          //B Arrays
long double *erind = NULL, *ezind = NULL, **ErTab = NULL, **EzTab = NULL, **EphiTab = NULL, **ErTab1 = NULL, **EzTab1 = NULL, **EphiTab1 = NULL;  //E Arrays
long double **ErTab2 = NULL, **EzTab2 = NULL, **EphiTab2 = NULL, **ErTab12 = NULL, **EzTab12 = NULL, **EphiTab12 = NULL;          //E Arrays
long double ****Brc  = NULL, ****Bphic  = NULL, ****Bzc  = NULL;
long double **ya=NULL, *rvec=NULL, *zvec=NULL;
long double rdist, zdist;
long double conv_rA, conv_rB, conv_zA, conv_zB; 

long double Emin_n = 1e30;  // minimum energy of neutrons in the B-field

// define racetrack current bars
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the straight wire 
long double Bars_1r[14],Bars_1phi[14],Bars_1z[14],Bars_2r[14],Bars_2phi[14],Bars_2z[14];  // lower horizontal 0 deg

void PrepareBField(){
	// current bar test
	// current from outside in
	Bars_1r[1]=0.60;Bars_1phi[1]=0.0;	Bars_1z[1]=-0.15;	Bars_2r[1]=0.0;	Bars_2phi[1]=0.0;	Bars_2z[1]=-0.15;   // lower horizontal 0 deg
	Bars_1r[2]=0.60;Bars_1phi[2]=pi/2.0;Bars_1z[2]=-0.15;	Bars_2r[2]=0.0;	Bars_2phi[2]=pi/2.0;Bars_2z[2]=-0.15; // lower horizontal 90 deg
	Bars_1r[3]=0.60;Bars_1phi[3]=pi;	Bars_1z[3]=-0.15;	Bars_2r[3]=0.0;	Bars_2phi[3]=pi;	Bars_2z[3]=-0.15; // lower horizontal 180 deg
	Bars_1r[4]=0.60;Bars_1phi[4]=pi*1.5;Bars_1z[4]=-0.15;	Bars_2r[4]=0.0;	Bars_2phi[4]=pi*1.5;Bars_2z[4]=-0.15; // lower horizontal 270 deg
	// current from inside out
	Bars_1r[5]=0.0;Bars_1phi[5]=0.0;	Bars_1z[5]=1.35;	Bars_2r[5]=0.6;	Bars_2phi[5]=0.0;	Bars_2z[5]=1.35;   // upper horizontal 0 deg
	Bars_1r[6]=0.0;Bars_1phi[6]=pi/2.0;	Bars_1z[6]=1.35;	Bars_2r[6]=0.6;	Bars_2phi[6]=pi/2.0;Bars_2z[6]=1.35;  // upper horizontal 90 deg
	Bars_1r[7]=0.0;Bars_1phi[7]=pi;		Bars_1z[7]=1.35;	Bars_2r[7]=0.6;	Bars_2phi[7]=pi;	Bars_2z[7]=1.35;  // upper horizontal 180 deg
	Bars_1r[8]=0.0;Bars_1phi[8]=pi*1.5;	Bars_1z[8]=1.35;	Bars_2r[8]=0.6;	Bars_2phi[8]=pi*1.5;Bars_2z[8]=1.35;  // upper horizontal 270 deg
	// current from high to low
	Bars_1r[9]=0.60;Bars_1phi[9]=0;		Bars_1z[9]=1.35;	Bars_2r[9]=0.6; Bars_2phi[9]=0;		Bars_2z[9]=-0.15;  //outer current 0 deg
	Bars_1r[10]=0.6;Bars_1phi[10]=pi/2;	Bars_1z[10]=1.35;	Bars_2r[10]=0.6;Bars_2phi[10]=pi/2;	Bars_2z[10]=-0.15; //outer current 90 deg
	Bars_1r[11]=0.6;Bars_1phi[11]=pi;	Bars_1z[11]=1.35;	Bars_2r[11]=0.6;Bars_2phi[11]=pi;	Bars_2z[11]=-0.15; //outer current 180 deg
	Bars_1r[12]=0.6;Bars_1phi[12]=pi*1.5;Bars_1z[12]=1.35;	Bars_2r[12]=0.6;Bars_2phi[12]=pi*1.5;Bars_2z[12]=-0.15; //outer current 270 deg
	// current from low to high
	Bars_1r[13]=0.0;Bars_1phi[13]=0.0;	Bars_1z[13]=-0.15;	Bars_2r[13]=0;	Bars_2phi[13]=0.0;	Bars_2z[13]=1.35;   // center current 4 TIMES THE CURRENT OF OTHERS!!!
  
	if (bfeldwahl == 0 || bfeldwahl == 2)
	{
        printf("\nPreparing the electromagnetic fields... \n");
        PrepIntpol(1);          // read input table file with E and B fields
		printf("allocating space for preinterpolation ... (about %.4LG MB)\n",(long double)n*m*12*16*3/1024/1024);

		// doing the preinterpolation
		Preinterpol(1);						
		
		/*
		// interpolation test
		wholetrackfile = (char*)malloc((outpathlength+20)*sizeof(char));
		sprintf(wholetrackfile, "%s/%06dinterpol.out", outpath, jobnumber);
		OUTFILE1 = fopen(wholetrackfile,mode_w);       // open outfile neut001.out
		Zeilencount=0;
		fprintf(OUTFILE1,"r z Br dBrdr dBrdphi dBrdz Bz dBzdr dBzdphi dBzdz Babs\n");
		//fprintf(OUTFILE1,"indr indz c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34 c41 c42 c43 c44 \n");
		BFeldSkal=1;
		long double rtst = 0.45;
		for(long double ztst = 0.137; ztst<=0.15; ztst=ztst+0.0001)
		{
			//BInterpol(rtst,0,ztst);
			BInterpol(rtst,0,ztst);
			fprintf(OUTFILE1,"%.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG %.17LG \n",
													rtst, ztst,Br,dBrdr,dBrdphi,dBrdz,Bz,dBzdr,dBzdphi,dBzdz,Bws);
			
		}
		return 0;
		// ENDE interpolation test*/
		
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
	if (Emin_n == 1e30) Emin_n = 0; // if Emin_n was not set, set it to zero
	
}

// get the size of the array ......................
void GetDim(int *m, int *n)
{
	int intmuell;
	long double muell, rtemp = 0, ztemp = 0, Brtemp = 0, Bphitemp=0, Bztemp=0, Ertemp=0, Ephitemp=0, Eztemp=0;	
	char str[1024];
	long double rtemptmp=-1e10, ztemptmp=-1e10;  // memory for r and z values in table
	int ri=1, zi=1;   // indizes
	
	if(FIN != NULL) 
		fclose(FIN);
	string path = inpath + "/fieldval.tab";
	FIN = fopen(path.c_str(),mode_r);
	
	// in der ersten zeile stehen die dimensionen auch drin
	sscanf(fgets(str,1024,FIN),"%d %d %d %LG",m,&intmuell,n,&muell);
	printf("The header says: \nThe arrays are %i by %i.\n",*m,*n);
	fprintf(LOGSCR,"The header says: \nThe arrays are %i by %i.\n",*m,*n);
	
	// discard next 10 lines because of header stuff
	int i;
	for(i=1; i<=10; i++){
		fgets(str,1024,FIN);
	}		
		
	// read all lines until the end of the file
	while (!feof(FIN) && (sscanf(fgets(str,1024,FIN),"%LG %LG %LG %LG %LG %LG %LG %LG %LG",&rtemp,&muell,&ztemp,&Brtemp,&Bphitemp,&Bztemp,&Ertemp,&Ephitemp,&Eztemp) == 9))
	{	
		rtemp = rtemp * lengthconv; ztemp = ztemp * lengthconv; Brtemp = Brtemp * Bconv; Bphitemp = Bphitemp * Bconv;    // Einheiten
		Bztemp = Bztemp * Bconv; Ertemp = Ertemp * Econv; Ephitemp = Ephitemp * Econv; Eztemp = Eztemp * Econv;          // ausgleichen
	
		/*if((zi<20)&&(ri==1))
			printf("%LG %LG %LG %LG %LG %LG %LG %LG %LG\n",rtemp,muell,ztemp,Brtemp,Bphitemp,Bztemp,Ertemp,Ephitemp,Eztemp);*/
		
		
		if((rtemptmp==-1e10)&&(ztemptmp==-1e10)) //first line
		{
			rtemptmp=r_mi=r_ma=rtemp;  // set maximum values to current
			ztemptmp=z_mi=z_ma=ztemp;		// set maximum values to current		
		}
	
			
		// Index ri nur erh�hen  wenn sich der r-Wert �ndert
		if (rtemp > rtemptmp){
			ri++;		
			if(rtemp>r_ma)
				r_ma=rtemp;
		}
		// Index zi nur erh�hen wenn sich der z-Wert �ndert
		if ( ztemp > ztemptmp){
			zi++;		
			if(ztemp>z_ma)
				z_ma=ztemp;
		}
		// set index back to one of it jumped down
		if (rtemp < rtemptmp){
			ri=1;	
			printf(".");
			fflush(stdout);
			if(rtemp<r_mi)
				r_mi=rtemp;
		}
		// set index back to one of it jumped down
		if (ztemp < ztemptmp){
			zi=1;	
			if(ztemp<z_mi)
				z_mi=ztemp;
		}
		
		rtemptmp = rtemp;  ztemptmp = ztemp;	 // memorize current r and z for later comparison
	}
	
	*m=ri; *n=zi; // set to maximum index values found in file
	
	printf("In reality they are %d by %d.\n",*m,*n);
	fprintf(LOGSCR,"In reality they are %d by %d.\n",*m,*n);
	printf("The r values go from %LG to %LG\n",r_mi,r_ma);
	fprintf(LOGSCR,"The r values go from %LG to %LG.\n",r_mi,r_ma);
	printf("The z values go from %LG to %LG.\n",z_mi,z_ma);
	fprintf(LOGSCR,"The z values go from %LG to %LG.\n",z_mi,z_ma);	
	return ;
}


// Zu den Arrays BrTab, BphiTab, BphiTab werden die entsprechenden Ableitungen berechnet:
// z.B: BrTab1 = dBrTab / dr
void CalcDeriv4th(long double rind[], long double zind[], long double **BrTab, long double **BzTab, long double **BphiTab, long double **BrTab1, long double **BrTab2,long double **BrTab12, long double **BzTab1, long double **BzTab2, long double **BzTab12, long double **BphiTab1, long double **BphiTab2,long double **BphiTab12,int m, int n)
{

	int j, k;

	for (j=3; j <= m-2; j++)
	{
			for (k=3; k <= n-2; k++)
			{
			BrTab1[j][k] = 1 / (12*(rind[j]-rind[j-1] )) *  ( (-1)*BrTab[j+2][k] + 8*BrTab[j+1][k] - 8*BrTab[j-1][k] + BrTab[j-2][k] ) ;   // derivative in r direction
			BrTab2[j][k] = 1 / (12*(zind[k]-zind[k-1] )) *  ( (-1)*BrTab[j][k+2] + 8*BrTab[j][k+1] - 8*BrTab[j][k-1] + BrTab[j][k-2] ) ;     // derivative in z direction
			BrTab12[j][k] = 1 / (12*(rind[j]-rind[j-1] )*12*(zind[k]-zind[k-1] )) 
								*(	-1*  ( (-1)*BrTab[j+2][k+2] + 8*BrTab[j+2][k+1] - 8*BrTab[j+2][k-1] + BrTab[j+2][k-2] ) 
										+8*  ( (-1)*BrTab[j+1][k+2] + 8*BrTab[j+1][k+1] - 8*BrTab[j+1][k-1] + BrTab[j+1][k-2] ) 
										-8*  ( (-1)*BrTab[j-1][k+2] + 8*BrTab[j-1][k+1] - 8*BrTab[j-1][k-1] + BrTab[j-1][k-2] )  
										+    ( (-1)*BrTab[j-2][k+2] + 8*BrTab[j-2][k+1] - 8*BrTab[j-2][k-1] + BrTab[j-2][k-2] )     );// cross derivative
			  
			BzTab1[j][k] = 1 / (12*(rind[j]-rind[j-1] )) *  ( (-1)*BzTab[j+2][k] + 8*BzTab[j+1][k] - 8*BzTab[j-1][k] + BzTab[j-2][k] ) ;   // derivative in r direction
			BzTab2[j][k] = 1 / (12*(zind[k]-zind[k-1] )) *  ( (-1)*BzTab[j][k+2] + 8*BzTab[j][k+1] - 8*BzTab[j][k-1] + BzTab[j][k-2] ) ;     // derivative in z direction
			BzTab12[j][k] = 1 / (12*(rind[j]-rind[j-1] )*12*(zind[k]-zind[k-1] )) 
									*(	-1*  ( (-1)*BzTab[j+2][k+2] + 8*BzTab[j+2][k+1] - 8*BzTab[j+2][k-1] + BzTab[j+2][k-2] ) 
										+8*  ( (-1)*BzTab[j+1][k+2] + 8*BzTab[j+1][k+1] - 8*BzTab[j+1][k-1] + BzTab[j+1][k-2] )
										-8*  ( (-1)*BzTab[j-1][k+2] + 8*BzTab[j-1][k+1] - 8*BzTab[j-1][k-1] + BzTab[j-1][k-2] )
										+    ( (-1)*BzTab[j-2][k+2] + 8*BzTab[j-2][k+1] - 8*BzTab[j-2][k-1] + BzTab[j-2][k-2] )        );// cross derivative
			  
			BphiTab1[j][k] = 1 / (12*(rind[j]-rind[j-1] )) *  ( (-1)*BphiTab[j+2][k] + 8*BphiTab[j+1][k] - 8*BphiTab[j-1][k] + BphiTab[j-2][k] ) ;   // derivative in r direction
			BphiTab2[j][k] = 1 / (12*(zind[k]-zind[k-1] )) *  ( (-1)*BphiTab[j][k+2] + 8*BphiTab[j][k+1] - 8*BphiTab[j][k-1] + BphiTab[j][k-2] ) ;     // derivative in z direction
			BphiTab12[j][k] = 1 / (12*(rind[j]-rind[j-1] )*12*(zind[k]-zind[k-1] )) 
									*(	-1*  ( (-1)*BphiTab[j+2][k+2] + 8*BphiTab[j+2][k+1] - 8*BphiTab[j+2][k-1] + BphiTab[j+2][k-2] ) 
										+8*  ( (-1)*BphiTab[j+1][k+2] + 8*BphiTab[j+1][k+1] - 8*BphiTab[j+1][k-1] + BphiTab[j+1][k-2] )
										-8*  ( (-1)*BphiTab[j-1][k+2] + 8*BphiTab[j-1][k+1] - 8*BphiTab[j-1][k-1] + BphiTab[j-1][k-2] )
										+    ( (-1)*BphiTab[j-2][k+2] + 8*BphiTab[j-2][k+1] - 8*BphiTab[j-2][k-1] + BphiTab[j-2][k-2] )        );// cross derivative
        }
	}
	
	printf("\nDerivatives calculated.  \n");
	fprintf(LOGSCR,"\nDerivatives calculated.  \n");
}


// Einlesen von B_r in die 2D Array BrTab[][] aus einer Textdatei folgenden Formats:
// x=r   y=egal z=z Bx=Br By=Bphi Bz Ex=Er Ey=Ephi Ez
// Es wird angenommen, dass die x und z Werte jeweils aufsteigend sortiert sind!!!!
int readWert(long double rind[], long double zind[], long double **BrTab,long double **BphiTab, long double **BzTab,long double **ErTab, long double **EphiTab,long double **EzTab, int *m, int *n){
	long double muell, rtemp = 1e-10, ztemp = 1e-10, Brtemp = 0, Bphitemp=0, Bztemp=0, Ertemp=0, Ephitemp=0, Eztemp=0;	
	int ri=1, zi=1; // indizes
	int ritemp, zitemp, zindex=2, perc=0, NrValueLines=0;     // temporary r,z indezes,   lineindex, percentage of file read
	char str[1024]; //retstr[1024];
	int mtmp=*m,ntmp=*n;
	
	printf("In readWert they are %d by %d.\n",*m,*n);
	fprintf(LOGSCR,"In readWert they are %d by %d.\n",*m,*n);

	
	// close file if it already open
	if(FIN != NULL)  
		fclose(FIN);
	string path(inpath + "/fieldval.tab");
	FIN = fopen(path.c_str(),mode_r);
	
	// discard first eleven lines with header
	int i;
	for(i=1; i<=11; i++)
	{
		fgets(str,1024,FIN);
	}

	
	printf("\nreading fieldval.tab \n");
	
   // Schleife �ber alle Zeilen der Input Datei
	do{
		
		// status if read is displayed
		if(100*zindex/(*n**m)>perc)
		{
			perc++;
			if(perc%10==0)
			{
				printf("%i%%",perc);
				fflush(stdout);
			}
			if((perc%10!=0)&&(perc%2==0)) 
			{
				printf(".");
				fflush(stdout);
			}
		}
		
		zindex++;		
		sscanf(fgets(str,1024,FIN),"%LG %LG %LG %LG %LG %LG %LG %LG %LG",&rtemp,&muell,&ztemp,&Brtemp,&Bphitemp,&Bztemp,&Ertemp,&Ephitemp,&Eztemp);
		NrValueLines++;
		rtemp = rtemp * lengthconv; ztemp = ztemp * lengthconv; Brtemp = Brtemp * Bconv; Bphitemp = Bphitemp * Bconv;    // Einheiten
		Bztemp = Bztemp * Bconv; Ertemp = Ertemp * Econv; Ephitemp = Ephitemp * Econv; Eztemp = Eztemp * Econv;          // ausgleichen
		ritemp = ri; zitemp = zi;  // memorize current r and z indexes for later comparison
		
		/*// set rind[1] and zind[1] with first line
		if(NrValueLines==1) 
		{
			rind[ri]=rtemp; zind[zi]=ztemp;
		}*/
		
		// Index ri nur erh�hen  wenn sich der r-Wert �ndert
		if ((rtemp > rind[ri])&&(NrValueLines!=1)&&(ri<mtmp)){
			ri++;
			zi=1;             // immer wenn sich der r-Wert �ndert springt z auch auf Eins
		}

		// Index zi nur erh�hen wenn sich der z-Wert �ndert
		if ((ztemp > zind[zi])&&(NrValueLines!=1)&&(zi<ntmp)){
			zi++;
		}
		
		// den Wert nur in die Array schreiben wenn sich der Index ge�ndert hat oder wenn nur eine zeile gelesen wurde bisher
		if ((ri != ritemp)||(zi!=zitemp)||(NrValueLines==1)) {
			rind[ri] = rtemp;
			zind[zi] = ztemp;
			BrTab[ri][zi] = Brtemp;
			BphiTab[ri][zi] = Bphitemp;
			BzTab[ri][zi] = Bztemp;
			if ((protneut==PROTON)||(protneut==ELECTRONS)){
				ErTab[ri][zi] = Ertemp;
				EphiTab[ri][zi] = Ephitemp;
				EzTab[ri][zi] = Eztemp;
			}
		}
		
		//if((NrValueLines<650)&&((zi%5==0)||(zi==1)||(zi==2))||((ri==mtmp)&&(zi>ntmp-3)))
		//	printf("%d %d %LG %LG %LG %LG\n",ri,zi,rind[ri],zind[zi],BrTab[ri][zi],BzTab[ri][zi]);	
		
	}while (!feof(FIN)); // maybe unsave

	printf("\n\nFinished reading field components from table file containing %d * %d values in %d lines",*m,*n,NrValueLines);
	fprintf(LOGSCR,"Finished reading field components from table file containing %d * %d values in %d lines",*m,*n,NrValueLines);
	
  return 0;
}





// Routine um in einer 1D Array xx den Index zu finden, so dass der Wert x in
// der Array zwischen den Indizes n und n+1 liegt  (jlo ist ein first guess des indexes und der r�ckgabewert)
//Given an array
// xx[1..n], and given a value x,returns a value jlo such that x is between
//xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or decreasing.
// jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as the
//initialguess for jlo on output.
 void hunt(long double xx[], int n, long double x, int *jlo){
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

// Koeffizienten f�r die sp�tere bikubische Interpolation berechnen
void bcucof(long double y[], long double y1[], long double y2[], long double y12[], long double d1, long double d2,long double **c){
	static int wt[16][16]=
		{{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
		{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
		{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
		{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
		{0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
		{0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
		{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
		{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
		{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
		{2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
		{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
		{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
	int l,k,j,i;
	long double xx,d1d2,cl[16],x[16];

	d1d2=d1*d2;
	for (i=1;i<=4;i++) {
		x[i-1]=y[i];
		x[i+3]=y1[i]*d1;
		x[i+7]=y2[i]*d2;
		x[i+11]=y12[i]*d1d2;
	}
	for (i=0;i<=15;i++) {
		xx=0.0;
		for (k=0;k<=15;k++) 
			xx += wt[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=1;i<=4;i++)
		for (j=1;j<=4;j++) 
			c[i][j]=cl[l++];
}


// Bikubische 2D Interpolation

void bcuint(long double y[],long double y1[] ,long double y2[],long double y12[],long double x1l, long double x1u, long double x2l, long double x2u, long double x1, long double x2, long double *ansy, long double *ansy1, long double *ansy2){
	int i;
	long double t,u,d1,d2,**c; //**dmatrix();
	//void bcucof(),nrerror(),free_matrix();

	c=dmatrix(1,4,1,4);
	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2,c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine BCUINT");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
		*ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
	}
	*ansy1 /= d1;
	*ansy2 /= d2;
	//fprintf(OUTFILE1,"%d %d  %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n"
		//																,0,0,	c[1][1],c[1][2],c[1][3],c[1][4],
			//																				c[2][1],c[2][2],c[2][3],c[2][4],
				//																			c[3][1],c[3][2],c[3][3],c[3][4],
					//																		c[4][1],c[4][2],c[4][3],c[4][4]);
	free_dmatrix(c,1,4,1,4);
}

// calculation of interpolated values from predetermined coefficients
//bcuint_new       (     indr,     indz,                     &Brc,                 rdist,               zdist,           rind[indr],       rind[indr +1],          zind[indz],      zind[indz +1],                r_n,                  z_n,                     &Br,                   &dBrdr,                   &dBrdz);
void bcuint_new(int indr, int indz,long double ****c, long double d1, long double d2, long double x1l, long double x1u, long double x2l, long double x2u, long double x1, long double x2, long double *ansy, long double *ansy1, long double *ansy2){
	int i;
	long double t,u; //**dmatrix();
	
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine BCUINT");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
		*ansy=t*(*ansy)+((c[indr][indz][i][4]*u+c[indr][indz][i][3])*u+c[indr][indz][i][2])*u+c[indr][indz][i][1];
		*ansy2=t*(*ansy2)+(3.0*c[indr][indz][i][4]*u+2.0*c[indr][indz][i][3])*u+c[indr][indz][i][2];
		*ansy1=u*(*ansy1)+(3.0*c[indr][indz][4][i]*t+2.0*c[indr][indz][3][i])*t+c[indr][indz][2][i];
	}
	*ansy1 /= d1;
	*ansy2 /= d2;
	//fprintf(OUTFILE1,"%d %d  %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n"
		//																,indr,indz,	c[indr][indz][1][1],c[indr][indz][1][2],c[indr][indz][1][3],c[indr][indz][1][4],
			//																				c[indr][indz][2][1],c[indr][indz][2][2],c[indr][indz][2][3],c[indr][indz][2][4],
				//																			c[indr][indz][3][1],c[indr][indz][3][2],c[indr][indz][3][3],c[indr][indz][3][4],
					//																		c[indr][indz][4][1],c[indr][indz][4][2],c[indr][indz][4][3],c[indr][indz][4][4]);
}




// Feld Interpolation vorbereiten
void PrepIntpol(int k){
  //long double muell;
// Gr��e der Arrays f�r das Einlesen der Tabelle bestimmen
   long double EnTemp;
	string path(inpath + "/fieldval.tab");
	FIN = fopen(path.c_str(),mode_r);
	if (FIN == NULL) 
		exit(-1);        // Fehlerbehandlung
	GetDim(&m, &n);


	// Vektoren der r und z Komponenten, Matrizen f�r B-Feld Komponenten und deren Ableitungen kreieren
	rind=dvector(1,m);
	zind=dvector(1,n);
	
	BrTab=dmatrix(1,m,1,n);
	BzTab=dmatrix(1,m,1,n);
	BphiTab=dmatrix(1,m,1,n);
	BrTab1=dmatrix(2,m-1,2,n-1);
	BzTab1=dmatrix(2,m-1,2,n-1);
	BphiTab1=dmatrix(2,m-1,2,n-1);
	BrTab2=dmatrix(2,m-1,2,n-1);
	BzTab2=dmatrix(2,m-1,2,n-1);
	BphiTab2=dmatrix(2,m-1,2,n-1);
	BrTab12=dmatrix(2,m-1,2,n-1);
	BzTab12=dmatrix(2,m-1,2,n-1);
	BphiTab12=dmatrix(2,m-1,2,n-1);

	// Vektoren der r und z Komponenten, Matrizen f�r B-Feld Komponenten und deren Ableitungen kreieren
	if ((protneut==PROTON)||(protneut==ELECTRONS)||(decay.on==2)){
		ErTab=dmatrix(1,m,1,n);
		EzTab=dmatrix(1,m,1,n);
		EphiTab=dmatrix(1,m,1,n);
		ErTab1=dmatrix(2,m-1,2,n-1);
		EzTab1=dmatrix(2,m-1,2,n-1);
		EphiTab1=dmatrix(2,m-1,2,n-1);
		ErTab2=dmatrix(2,m-1,2,n-1);
		EzTab2=dmatrix(2,m-1,2,n-1);
		EphiTab2=dmatrix(2,m-1,2,n-1);
		ErTab12=dmatrix(2,m-1,2,n-1);
		EzTab12=dmatrix(2,m-1,2,n-1);
		EphiTab12=dmatrix(2,m-1,2,n-1);
	}
	
	fehler = readWert(rind, zind, BrTab, BphiTab, BzTab, ErTab, EphiTab, EzTab, &m, &n);

	// Die Ableitungen der B-Werte in den Tabellen berechnen
	CalcDeriv4th(rind, zind, BrTab, BzTab, BphiTab, BrTab1, BrTab2, BrTab12, BzTab1, BzTab2, BzTab12, BphiTab1, BphiTab2, BphiTab12,m, n);
	// Die Ableitungen der E-Werte in den Tabellen berechnen

	if ((protneut==PROTON)||(protneut==ELECTRONS)||(decay.on==2)){
		CalcDeriv4th(rind, zind, ErTab, EzTab, EphiTab, ErTab1, ErTab2, ErTab12, EzTab1, EzTab2, EzTab12, EphiTab1, EphiTab2, EphiTab12,m, n);
	}
	
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB	
	conv_rA = r_mi; 													// y Abschnitt von r
	conv_rB = (r_ma-r_mi)/(m-1);						// Steigung von r
	conv_zA = z_mi; 													// y Abschnitt von z
	conv_zB = (z_ma-z_mi)/(n-1);						// Steigung von z
	printf("conv_rA =  %LG, conv_rB = %LG\n",conv_rA, conv_rB);
	printf("conv_zA =  %LG, conv_zB = %LG\n",conv_zA, conv_zB);	
	printf("r_mi = %LG, r_ma = %LG, z_mi = %LG, z_ma = %LG\n",r_mi,r_ma,z_mi,z_ma);
	fprintf(LOGSCR,"r_mi = %LG, r_ma = %LG, z_mi = %LG, z_ma = %LG\n",r_mi,r_ma,z_mi,z_ma);	
		
	long double Babsmax, Babsmin, Babsmaxtmp, Eabsmax, Eabsmin, Eabsmaxtmp;
	long double rBabsmin, zBabsmin;
	for (int j=1; j <= m; j++)
	{
		for (int k=1; k <= n; k++)
		{
			// Calculate maximum values
			Babsmaxtmp = sqrtl(BrTab[j][k]*BrTab[j][k]+BphiTab[j][k]*BphiTab[j][k]+
			BzTab[j][k]*BzTab[j][k]);
			EnTemp =  m_n*gravconst*zind[k] + (mu_nSI/ele_e)*Babsmaxtmp;
			if(Babsmax < Babsmaxtmp)
				Babsmax = Babsmaxtmp;
			if((Babsmin > Babsmaxtmp))
			{
				Babsmin = Babsmaxtmp;				
			}
			if(EnTemp < Emin_n && InSourceVolume(rind[j],0,zind[k]))
			{
				Emin_n = EnTemp;
				rBabsmin = rind[j];
				zBabsmin = zind[k];
			}
			if ((protneut==PROTON)||(protneut==ELECTRONS)||(decay.on==2))
			{
				Eabsmaxtmp = sqrtl(ErTab[j][k]*ErTab[j][k]+EphiTab[j][k]*EphiTab[j][k]+
				EzTab[j][k]*EzTab[j][k]);
				if(Eabsmax < Eabsmaxtmp)
					Eabsmax = Eabsmaxtmp;
				if(Eabsmin > Eabsmaxtmp)
					Eabsmin = Eabsmaxtmp;
			}
							
		}
	}
	
	//Emin_n =  m_n*gravconst*zBabsmin - (mu_n/ele_e)*Babsmin;
	
	printf("The interpolation can done from r = %.17LG to %.17LG and from n z = %.17LG to %.17LG\n",rind[3],rind[m-2],zind[3],zind[n-2]);
	printf("The input table file has values of |B| from %.17LG T to %.17LG T\n and values of |E| from %.17LG V/m to %.17LG V/m\n",Babsmin,Babsmax,Eabsmin,Eabsmax);
	printf("The minimum energy a low field seeking neutron has to have is %.3LG eV (at r:%.3LG, z: %.3LG).\n",Emin_n,rBabsmin,zBabsmin);
	fprintf(LOGSCR,"The interpolation can done from r = %.17LG to %.17LG and from n z = %.17LG to %.17LG\n",rind[3],rind[m-2],zind[3],zind[n-2]);
	fprintf(LOGSCR,"The input table file has values of |B| from %.17LG T to %.17LG T\n and values of |E| from %.17LG V/m to %.17LG V/m\n",Babsmin,Babsmax,Eabsmin,Eabsmax);
	fprintf(LOGSCR,"The minimum energy a low field seeking neutron has to have is %.3LG eV (at r:%.3LG, z: %.3LG).\n",Emin_n,rBabsmin,zBabsmin);
	
	return;
}



void BInterpol(long double r_n, long double phi, long double z_n){
	if (BFeldSkal != 0){
		indr = 1+(int) ((r_n- conv_rA)/(conv_rB));    // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!		
		indz = 1+(int) ((z_n-conv_zA) / (conv_zB));  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3)  ){			
			// bicubic interpolation
			bcuint_new(indr, indz, Brc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Br, &dBrdr, &dBrdz);
			//bcuint_new(indr, indz, Bphic, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bphi, &dBphidr, &dBphidz);
			bcuint_new(indr, indz, Bzc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bz, &dBzdr, &dBzdz);

			Br = BFeldSkal * Br;
			dBrdr = BFeldSkal * dBrdr;
			dBrdz = BFeldSkal * dBrdz;
			Bz = BFeldSkal * Bz;
			dBzdr = BFeldSkal * dBzdr;
			dBzdz = BFeldSkal * dBzdz;
		}
		else
		{
			printf("\n The array index has left fieldval boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n ,indz);
			fprintf(LOGSCR,"\n The array index has left fieldval boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n, indz);
			OutputState(ystart,1);
		}
	} 

	// add fields and field derivations of racetrack coils
	
	if((Ibar!=0)&&(RodFieldMultiplicator>0))
	{
		if(Racetracks==1)
		{
			// old implementation of racetracks as 5 infinitely long vertical current bars, fast
			Racetrack(r_n, phi, z_n, -0.424264, -0.424264, RodFieldMultiplicator*Ibar);
			Racetrack(r_n, phi, z_n, 0.424264, -0.424264, RodFieldMultiplicator*Ibar); 
			Racetrack(r_n, phi, z_n, -0.424264, 0.424264, RodFieldMultiplicator*Ibar);
			Racetrack(r_n, phi, z_n, 0.424264, 0.424264, RodFieldMultiplicator*Ibar); 
			CenterCurrent(r_n, phi, z_n, -4.0*RodFieldMultiplicator*Ibar);
		}
		else if(Racetracks==2)
		{			
			// new but slow implementation of racetracks including horizontal parts
			BarRaceTrack(r_n, phi, z_n, RodFieldMultiplicator*Ibar);   // here B field is only calculated for charged particles and derivatives only for neutrons
		}
		
	}  		
		
}



void EInterpol(long double r_n, long double phi, long double z_n){
	if (EFeldSkal != 0){
		yyy=dvector(1,4);
		yyy1=dvector(1,4);
		yyy2=dvector(1,4);
		yyy12=dvector(1,4);
	
		if (r_n >= r_mi && r_n <= r_ma && z_n >= z_mi && z_n <= z_ma ) 
		{
			
			
			//hunt(rind, m, r_n, &indr);
			//hunt(zind, n, z_n, &indz);
			
			// new calculation of indexes
			indr = 1 + (int) ((r_n - conv_rA) / (conv_rB));   // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
			indz = 1 + (int) ((z_n - conv_zA) / (conv_zB));  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
	
			// Abst�nde der Koordinaten voneinander
			//dr = rind[indr + 1] - rind [indr];
			//dz = zind[indz + 1] - zind [indz];
	
			// cout << "Indizes: " << indr << " " << indz << "\n";
	
			if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3)){
				// Rechteck mit den Werten und deren Ableitungen bef�llen und Interpolieren
				yyy[1] = ErTab[indr][indz];
				yyy[2] = ErTab[indr+1][indz];
				yyy[3] = ErTab[indr+1][indz+1];
				yyy[4] = ErTab[indr][indz+1];
	        
				yyy1[1] = ErTab1[indr][indz];
				yyy1[2] = ErTab1[indr+1][indz];
				yyy1[3] = ErTab1[indr+1][indz+1];
				yyy1[4] = ErTab1[indr][indz+1];
				
				yyy2[1] = ErTab2[indr][indz];
				yyy2[2] = ErTab2[indr+1][indz];
				yyy2[3] = ErTab2[indr+1][indz+1];
				yyy2[4] = ErTab2[indr][indz+1];
	        
				yyy12[1] = ErTab12[indr][indz];
				yyy12[2] = ErTab12[indr+1][indz];
				yyy12[3] = ErTab12[indr+1][indz+1];
				yyy12[4] = ErTab12[indr][indz+1];
				
				bcuint(yyy, yyy1, yyy2, yyy12, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Er, &dErdr, &dErdz);
	
				// Rechteck mit den Werten und deren Ableitungen bef�llen und Interpolieren
				yyy[1] = EzTab[indr][indz];
				yyy[2] = EzTab[indr+1][indz];
				yyy[3] = EzTab[indr+1][indz+1];
				yyy[4] = EzTab[indr][indz+1];
	        
				yyy1[1] = EzTab1[indr][indz];
				yyy1[2] = EzTab1[indr+1][indz];
				yyy1[3] = EzTab1[indr+1][indz+1];
				yyy1[4] = EzTab1[indr][indz+1];
				
				yyy2[1] = EzTab2[indr][indz];
				yyy2[2] = EzTab2[indr+1][indz];
				yyy2[3] = EzTab2[indr+1][indz+1];
				yyy2[4] = EzTab2[indr][indz+1];
				
				yyy12[1] = EzTab12[indr][indz];
				yyy12[2] = EzTab12[indr+1][indz];
				yyy12[3] = EzTab12[indr+1][indz+1];
				yyy12[4] = EzTab12[indr][indz+1];
				
				bcuint(yyy, yyy1, yyy2, yyy12, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Ez, &dEzdr, &dEzdz);
	
				// Rechteck mit den Werten und deren Ableitungen bef�llen und Interpolieren
				yyy[1] = EphiTab[indr][indz];
				yyy[2] = EphiTab[indr+1][indz];
				yyy[3] = EphiTab[indr+1][indz+1];
				yyy[4] = EphiTab[indr][indz+1];
				
				yyy1[1] = EphiTab1[indr][indz];
				yyy1[2] = EphiTab1[indr+1][indz];
				yyy1[3] = EphiTab1[indr+1][indz+1];
				yyy1[4] = EphiTab1[indr][indz+1];
				
				yyy2[1] = EphiTab2[indr][indz];
				yyy2[2] = EphiTab2[indr+1][indz];
				yyy2[3] = EphiTab2[indr+1][indz+1];
				yyy2[4] = EphiTab2[indr][indz+1];
				
				yyy12[1] = EphiTab12[indr][indz];
				yyy12[2] = EphiTab12[indr+1][indz];
				yyy12[3] = EphiTab12[indr+1][indz+1];
				yyy12[4] = EphiTab12[indr][indz+1];
				
				bcuint(yyy, yyy1, yyy2, yyy12, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Ephi, &dEphidr, &dEphidz);
	        }else{
				printf("\n The array index has left boundaries!!! Exiting!!! \n");
				fprintf(LOGSCR,"\n The array index has left boundaries!!! Exiting!!! \n");
				OutputState(ystart,1);
	        }
		}
		else	
		{
	        //char bloed;
	        printf("The particle has entered forbidden space!!! I will terminate it now! \n");
	        fprintf(LOGSCR,"The particle has entered forbidden space!!! I will terminate it now! \n");
				Er=0.0; Ephi=0.0; Ez=0.0;
				kennz=99;    // something's wrong
	        stopall=1;
		}
	
		free_dvector(yyy,1,4);
		free_dvector(yyy1,1,4);
		free_dvector(yyy2,1,4);
		free_dvector(yyy12,1,4);
		
		// scale electric field if desired
		Er = EFeldSkal * Er;
		dErdr = EFeldSkal * dErdr;
		dErdz = EFeldSkal * dErdz;
		Ephi = EFeldSkal * Ephi;
		dEphidr = EFeldSkal * dEphidr;
		dEphidz = EFeldSkal * dEphidz;
		Ez = EFeldSkal * Ez;
		dEzdr = EFeldSkal * dEzdr;
		dEzdz = EFeldSkal * dEzdz;
	}			
	return;
}


// produce zero field
long double Bnull()
{
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
	
	Er = 0;
	dErdr = 0;
	dErdz = 0;
	Ephi = 0;
	dEphidr = 0;
	dEphidz = 0;
	Ez = 0;
	dEzdr = 0;
	dEzdz = 0;	
	
	return 0;
}



void Preinterpol(int p)
{
	int indr, indz, perc=0; 
	long double d1, d2;
	long double **ctemp;
	
	// allocating space for the preinterpolation, we need to cast it to long double ****
	// The B*c are 4D arrays with m x n x 4 x 4 fields
	Brc = (long double ****) viertensor(1,m,1,n,1,4,1,4);
	Bphic = (long double ****) viertensor(1,m,1,n,1,4,1,4);
	Bzc = (long double ****) viertensor(1,m,1,n,1,4,1,4);	
	
	printf("\nStart Preinterpolation!\n");
	
	ctemp=dmatrix(1,4,1,4);	
	yyy=dvector(1,4);
	yyy1=dvector(1,4);
	yyy2=dvector(1,4);
	yyy12=dvector(1,4);
	
	
	// calculate distances of coordinates, normally 0.002 m
	rdist = d1 = rind[4] - rind[3];
	zdist = d2 = zind[4] - zind[3];
	
	cout << "rdist: " << rdist << " zdist: " << rdist << endl;

	for (indr=3; indr < m-2; indr++)
	{
		for (indz=3; indz < n-2; indz++)
		{
			
			
		// status if read is displayed
		if(100*indr/(m)>perc)
		{
			perc++;
			if(perc%10==0)
			{
				printf("%i%%",perc);
				fflush(stdout);
			}
			if((perc%10!=0)&&(perc%2==0)) 
			{
				printf(".");
				fflush(stdout);
			}
		}
			
			
			
			// fill rectancle with values and derivatives
			yyy[1] = BrTab[indr][indz];
			yyy[2] = BrTab[indr+1][indz];
			yyy[3] = BrTab[indr+1][indz+1];
			yyy[4] = BrTab[indr][indz+1];
			
			yyy1[1] = BrTab1[indr][indz];
			yyy1[2] = BrTab1[indr+1][indz];
			yyy1[3] = BrTab1[indr+1][indz+1];
			yyy1[4] = BrTab1[indr][indz+1];
			
			yyy2[1] = BrTab2[indr][indz];
			yyy2[2] = BrTab2[indr+1][indz];
			yyy2[3] = BrTab2[indr+1][indz+1];
			yyy2[4] = BrTab2[indr][indz+1];
			
			yyy12[1] = BrTab12[indr][indz];
			yyy12[2] = BrTab12[indr+1][indz];
			yyy12[3] = BrTab12[indr+1][indz+1];
			yyy12[4] = BrTab12[indr][indz+1];
			
															
			// determine koefficients of interpolation
			bcucof(yyy,yyy1,yyy2,yyy12,d1,d2,ctemp);
			
			// put them into major array
			for (int i=1;i<=4;i++)
			{
				for (int k=1;k<=4;k++)
				{
					Brc[indr][indz][i][k] = ctemp[i][k];
				}
			}
			
			
			yyy[1] = BphiTab[indr][indz];
			yyy[2] = BphiTab[indr+1][indz];
			yyy[3] = BphiTab[indr+1][indz+1];
			yyy[4] = BphiTab[indr][indz+1];
			
			yyy1[1] = BphiTab1[indr][indz];
			yyy1[2] = BphiTab1[indr+1][indz];
			yyy1[3] = BphiTab1[indr+1][indz+1];
			yyy1[4] = BphiTab1[indr][indz+1];
			
			yyy2[1] = BphiTab2[indr][indz];
			yyy2[2] = BphiTab2[indr+1][indz];
			yyy2[3] = BphiTab2[indr+1][indz+1];
			yyy2[4] = BphiTab2[indr][indz+1];
			
			yyy12[1] = BphiTab12[indr][indz];
			yyy12[2] = BphiTab12[indr+1][indz];
			yyy12[3] = BphiTab12[indr+1][indz+1];
			yyy12[4] = BphiTab12[indr][indz+1];
			
						
			bcucof(yyy,yyy1,yyy2,yyy12,d1,d2,ctemp);
			
			for (int i=1;i<=4;i++)
			{
				for (int k=1;k<=4;k++)
				{
					Bphic[indr][indz][i][k] = ctemp[i][k];
				}
			}
			
			
			yyy[1] = BzTab[indr][indz];
			yyy[2] = BzTab[indr+1][indz];
			yyy[3] = BzTab[indr+1][indz+1];
			yyy[4] = BzTab[indr][indz+1];
			
			yyy1[1] = BzTab1[indr][indz];
			yyy1[2] = BzTab1[indr+1][indz];
			yyy1[3] = BzTab1[indr+1][indz+1];
			yyy1[4] = BzTab1[indr][indz+1];
			
			yyy2[1] = BzTab2[indr][indz];
			yyy2[2] = BzTab2[indr+1][indz];
			yyy2[3] = BzTab2[indr+1][indz+1];
			yyy2[4] = BzTab2[indr][indz+1];
			
			yyy12[1] = BzTab12[indr][indz];
			yyy12[2] = BzTab12[indr+1][indz];
			yyy12[3] = BzTab12[indr+1][indz+1];
			yyy12[4] = BzTab12[indr][indz+1];
			
						
			bcucof(yyy,yyy1,yyy2,yyy12,d1,d2,ctemp);
			
			for (int i=1;i<=4;i++)
			{
				for (int k=1;k<=4;k++)
				{
					Bzc[indr][indz][i][k] = ctemp[i][k];
				}
			}

			//printf("Bis hier gehts: indr: %i indz: %i \n", indr, indz);
		}
	}
	
	free_dvector(yyy,1,4);
	free_dvector(yyy1,1,4);
	free_dvector(yyy2,1,4);
	free_dvector(yyy12,1,4);
	free_dmatrix(ctemp,1,4,1,4);
	
	printf("100%%\nDone with Preinterpolation!\n");
	printf("freeing the BField matrix ... (about %.4LG MB)\n",(long double)n*m*12*12/1024/1024);
	// now we don't need the BF matrix anymore
	free_dmatrix(BrTab,1,m,1,n);
	free_dmatrix(BzTab,1,m,1,n);
	free_dmatrix(BphiTab,1,m,1,n);
	
	free_dmatrix(BrTab1,2,m-1,2,n-1);
	free_dmatrix(BzTab1,2,m-1,2,n-1);
	free_dmatrix(BphiTab1,2,m-1,2,n-1);
	
	free_dmatrix(BrTab2,2,m-1,2,n-1);
	free_dmatrix(BzTab2,2,m-1,2,n-1);
	free_dmatrix(BphiTab2,2,m-1,2,n-1);
	
	free_dmatrix(BrTab12,2,m-1,2,n-1);
	free_dmatrix(BzTab12,2,m-1,2,n-1);
	free_dmatrix(BphiTab12,2,m-1,2,n-1);
	return;
}


// produce linearly rising B-field in z-direction and linearly changing in time
long double Banalytic(long double r,long double phi,long double z, long double t)
{
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

// produce   B-Field of one straight wire in z direction an positions lx, ly and current I_rt
void Racetrack(long double r,long double phi,long double z,  long double lx, long double ly, long double I_rt){
		long double Brv,dBrdrv,dBrdphiv,Bphiv,dBphidrv,dBphidphiv; 
	// cartesian coordinates of neutron
	//long double x = cosl(phi) * r, y = sinl(phi) *r;
	// cartesian coordinates of racetracks	
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	
	long double t2 = cosl(phi);    
  long double t3 = ly*r*t2;
  long double     t5 = sinl(phi);
  long double     t6 = lx*r*t5;
      long double t7 = t3-t6;
      long double t8 = vorfaktor*t7;
      long double t9 = r*r;
      long double t10 = t2*t2;
      long double t11 = t9*t10;
      long double t12 = r*t2;
      long double t13 = t12*lx;
      long double t15 = lx*lx;
      long double t16 = t5*t5;
      long double t17 = t9*t16;
      long double t18 = r*t5;
      long double t19 = t18*ly;
      long double t21 = ly*ly;
      long double t22 = t11-2.0*t13+t15+t17-2.0*t19+t21;
      long double t23 = 1/t22;
      long double t24 = 1/r;
      long double t25 = t23*t24;
      Brv = t8*t25;
      long double t31 = t22*t22;
      long double t33 = 1/t31*t24;
      long double t34 = r*t10;
      long double t35 = t2*lx;
      long double t36 = r*t16;
      long double t37 = t5*ly;
      long double t42 = 1/t9;
      dBrdrv = vorfaktor*(ly*t2-lx*t5)*t25-t8*t33*(2.0*t34-2.0*t35+2.0*t36-2.0*t37)-t8*t23*t42;
      dBrdphiv = vorfaktor*(-t19-t13)*t25+2.0*t8*t33*t7;
      long double t53 = vorfaktor*(-t17-t19-t11+t13);
      long double t54 = t12-lx;
      long double t55 = t54*t54;
      long double t56 = t18-ly;
      long double t57 = t56*t56;
      long double t58 = t55+t57;
      long double t59 = sqrt(t58);
      long double t60 = 1/t59;
      long double t61 = t60*t24;
      Bphiv = t53*t61;
      long double t69 = 1/t59/t58*t24;
      dBphidrv = vorfaktor*(-2.0*t36-t37-2.0*t34+t35)*t61-t53*t69*(2.0*t54*t2+2.0*t56*t5)/2.0-t53*t60*t42;
      dBphidphiv = vorfaktor*(-t3-t6)*t61-t53*t69*(-2.0*t54*r*t5+2.0*t56*r*t2)/2.0;


      // add values to current B-values (superposition principle)
      Br = Br + Brv;
      dBrdr = dBrdr + dBrdrv;
      dBrdphi = dBrdphi + dBrdphiv;
      Bphi = Bphi + Bphiv;
      dBphidr = dBphidr + dBphidrv;
      dBphidphi = dBphidphi + dBphidphiv;
      
      if((dBphidrv> 1) || (dBphidrv < -1))
         printf("Br = %17LG, %17LG, %17LG \n Bphi = %17LG, %17LG, %17LG \n r = %17LG, phi = %17LG, z = %17LG\n", Br,dBrdr,dBrdphi,Bphi ,dBphidr,dBphidphi, r, phi, z);
	
	
	return;
	
}

// produce center current rod field 
void CenterCurrent(long double r,long double phi,long double z, long double I_rt){
	
	long double vorfaktor = mu0 * I_rt / (2 * pi);	
		
	// only Bphi is present
	Bphi += vorfaktor/r;	
	
	// dBphidr
	dBphidr -= vorfaktor/(r*r);		
	
	// dBphidphi =0		
	
	return;	
}


// produce B-field at position r,phi,z of a straight wire with current I_rt 
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the ends of the straight wire 
// code generated by Maple ("FiniteWireDiffs.mw") with formula http://de.wikipedia.org/wiki/Biot-Savart#Gerader_Linienleiter
void StraightWireField(const long double r,const long double phi,const long double z, const long double I_rt, 
						const long double SW1r, const long double SW1phi, const long double SW1z, 
						const long double SW2r, const long double SW2phi, const long double SW2z)
{
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = cosl(SW1phi);
	long double t2 = t1 * SW1r;
	long double t3 = cosl(phi);
	long double t4 = t3 * r;
	long double t5 = t2 - t4;
	long double t6 = t5 * t5;
	long double t7 = sinl(SW1phi);
	long double t8 = t7 * SW1r;
	long double t9 = sin(phi);
	long double t10 = t9 * r;
	long double t11 = t8 - t10;
	long double t12 = t11 * t11;
	long double t13 = SW1z - z;
	long double t14 = t13 * t13;
	long double t15 = t6 + t12 + t14;
	long double t16 = sqrtl(t15);
	long double t17 = 0.1e1 / t16;
	long double t18 = t17 * vorfaktor;
	long double t19 = 0.1e1 / t15;
	long double t20 = SW1r * SW1r;
	long double t21 = t1 * t1;
	long double t23 = cosl(SW2phi);
	long double t24 = t23 * SW2r;
	long double t28 = t7 * t7;
	long double t30 = sinl(SW2phi);
	long double t31 = t30 * SW2r;
	long double t37 = SW2z - SW1z;
	long double t39 = t21 * t20 - t1 * (t24 + t4) * SW1r + t28 * t20 - t7 * (t10 + t31) * SW1r + t24 * t4 + t31 * t10 - t37 * t13;
	long double t40 = t39 * t39;
	long double t42 = t24 - t2;
	long double t43 = t42 * t42;
	long double t44 = t31 - t8;
	long double t45 = t44 * t44;
	long double t46 = t37 * t37;
	long double t47 = t43 + t45 + t46;
	long double t48 = 0.1e1 / t47;
	long double t50 = 0.1e1 - t48 * t40 * t19;
	long double t51 = sqrtl(t50);
	long double t52 = 0.1e1 / t51;
	long double t53 = t52 * t18;
	long double t54 = t24 - t4;
	long double t55 = t54 * t54;
	long double t56 = t31 - t10;
	long double t57 = t56 * t56;
	long double t58 = SW2z - z;
	long double t59 = t58 * t58;
	long double t60 = t55 + t57 + t59;
	long double t61 = sqrtl(t60);
	long double t62 = 0.1e1 / t61;
	long double t63 = SW2r * SW2r;
	long double t64 = t23 * t23;
	long double t69 = t30 * t30;
	long double t77 = -t64 * t63 + t23 * (t2 + t4) * SW2r - t69 * t63 + t30 * (t8 + t10) * SW2r - t4 * t2 - t10 * t8 - t37 * t58;
	long double t79 = sqrtl(t47);
	long double t80 = 0.1e1 / t79;
	long double t84 = -t80 * t77 * t62 + t80 * t39 * t17;
	long double t87 = -t37 * t11 + t44 * t13;
	long double t88 = t87 * t87;
	long double t91 = -t42 * t13 + t37 * t5;
	long double t92 = t91 * t91;
	long double t95 = -t44 * t5 + t42 * t11;
	long double t96 = t95 * t95;
	long double t97 = t88 + t92 + t96;
	long double t98 = sqrtl(t97);
	long double t99 = 0.1e1 / t98;
	long double t100 = t99 * t84;
	long double t101 = t3 * t87;
	long double t103 = t101 * t100 * t53;
	long double t104 = t9 * t91;
	long double t106 = t104 * t100 * t53;
	long double t109 = 0.1e1 / t16 / t15;
	long double t110 = t109 * vorfaktor;
	long double t111 = t84 * t52;
	long double t112 = t111 * t110;
	long double t113 = t87 * t99;
	long double t116 = -t3 * t5 - t9 * t11;
	long double t117 = 0.2e1 * t116 * t3;
	long double t122 = 0.1e1 / t51 / t50;
	long double t124 = t84 * t122 * t18;
	long double t125 = t15 * t15;
	long double t127 = t40 / t125;
	long double t130 = t39 * t19;
	long double t139 = -t1 * t3 * SW1r - t7 * t9 * SW1r + t23 * SW2r * t3 + t30 * SW2r * t9;
	long double t143 = 0.2e1 * t116 * t48 * t127 - 0.2e1 * t139 * t48 * t130;
	long double t144 = t143 * t3;
	long double t150 = t77 / t61 / t60;
	long double t159 = t39 * t109;
	long double t166 = t99 * ((-t3 * t54 - t9 * t56) * t80 * t150 - t80 * t139 * t62 - t116 * t80 * t159 + t80 * t139 * t17);
	long double t169 = t111 * t18;
	long double t171 = 0.1e1 / t98 / t97;
	long double t172 = t87 * t171;
	long double t173 = t9 * t87;
	long double t175 = t3 * t91;
	long double t178 = t42 * t9;
	long double t179 = t44 * t3 - t178;
	long double t181 = t37 * t173 - t37 * t175 + t179 * t95;
	long double t182 = 0.2e1 * t181 * t3;
	long double t186 = t91 * t99;
	long double t187 = 0.2e1 * t116 * t9;
	long double t191 = t143 * t9;
	long double t197 = t91 * t171;
	long double t198 = 0.2e1 * t181 * t9;
	long double t207 = t9 * r * t5 - t3 * r * t11;
	long double t208 = 0.2e1 * t207 * t3;
	long double t214 = SW1r * r;
	long double t221 = t1 * t9 * t214 - t7 * t3 * t214 - t24 * t10 + t31 * t4;
	long double t225 = 0.2e1 * t207 * t48 * t127 - 0.2e1 * t221 * t48 * t130;
	long double t226 = t225 * t3;
	long double t246 = t99 * ((t9 * r * t54 - t3 * r * t56) * t80 * t150 - t80 * t221 * t62 - t207 * t80 * t159 + t80 * t221 * t17);
	long double t257 = -t44 * t10 - t42 * t4;
	long double t259 = t37 * t3 * r * t87 + t37 * t9 * r * t91 + t257 * t95;
	long double t260 = 0.2e1 * t259 * t3;
	long double t264 = r * t99;
	long double t265 = t3 * t3;
	long double t266 = t37 * t265;
	long double t270 = t173 * t100 * t53;
	long double t271 = 0.2e1 * t207 * t9;
	long double t275 = t225 * t9;
	long double t281 = 0.2e1 * t259 * t9;
	long double t285 = t9 * t9;
	long double t286 = t37 * t285;
	long double t290 = t175 * t100 * t53;
	long double t291 = -t208 * t113 * t112 / 0.2e1 - t226 * t113 * t124 / 0.2e1 + t101 * t246 * t53 - t260 * t172 * t169 / 0.2e1 + t266 * t264 * t169 - t270 - t271 * t186 * t112 / 0.2e1 - t275 * t186 * t124 / 0.2e1 + t104 * t246 * t53 - t281 * t197 * t169 / 0.2e1 + t286 * t264 * t169 + t290;
	long double t292 = -0.2e1 * t13 * t3;
	long double t301 = -0.2e1 * t13 * t48 * t127 - 0.2e1 * t37 * t48 * t130;
	long double t302 = t301 * t3;
	long double t317 = t99 * (-t58 * t80 * t150 - t80 * t37 * t62 + t13 * t80 * t159 + t80 * t37 * t17);
	long double t322 = -t44 * t87 + t42 * t91;
	long double t323 = 0.2e1 * t322 * t3;
	long double t330 = -0.2e1 * t13 * t9;
	long double t334 = t301 * t9;
	long double t340 = 0.2e1 * t322 * t9;
	long double t429 = t52 * t110;
	long double t434 = t122 * t18;
	long double t441 = t171 * t84;
	Br += 		t103 + t106;
	dBrdr += 	-t117 * t113 * t112 / 0.2e1 - t144 * t113 * t124 / 0.2e1 + t101 * t166 * t53 - t182 * t172 * t169 / 0.2e1 - t187 * t186 * t112 / 0.2e1 - t191 * t186 * t124 / 0.2e1 + t104 * t166 * t53 - t198 * t197 * t169 / 0.2e1;
	dBrdphi += 	t291;
	dBrdz += 	-t292 * t113 * t112 / 0.2e1 - t302 * t113 * t124 / 0.2e1 + t101 * t317 * t53 - t323 * t172 * t169 / 0.2e1 - t44 * t3 * t100 * t53 - t330 * t186 * t112 / 0.2e1 - t334 * t186 * t124 / 0.2e1 + t104 * t317 * t53 - t340 * t197 * t169 / 0.2e1 + t178 * t100 * t53;
	Bphi += 	-t270 + t290;
	dBphidr += 	t187 * t113 * t112 / 0.2e1 + t191 * t113 * t124 / 0.2e1 - t173 * t166 * t53 + t198 * t172 * t169 / 0.2e1 - t286 * t100 * t53 - t117 * t186 * t112 / 0.2e1 - t144 * t186 * t124 / 0.2e1 + t175 * t166 * t53 - t182 * t197 * t169 / 0.2e1 - t266 * t100 * t53;
	dBphidphi += t271 * t113 * t112 / 0.2e1 + t275 * t113 * t124 / 0.2e1 - t173 * t246 * t53 + t281 * t172 * t169 / 0.2e1 - t103 - t208 * t186 * t112 / 0.2e1 - t226 * t186 * t124 / 0.2e1 + t175 * t246 * t53 - t260 * t197 * t169 / 0.2e1 - t106;
	dBphidz += 	t330 * t113 * t112 / 0.2e1 + t334 * t113 * t124 / 0.2e1 - t173 * t317 * t53 + t340 * t172 * t169 / 0.2e1 + t9 * t44 * t100 * t53 - t292 * t186 * t112 / 0.2e1 - t302 * t186 * t124 / 0.2e1 + t175 * t317 * t53 - t323 * t197 * t169 / 0.2e1 + t3 * t42 * t100 * t53;
	Bz += 		t95 * t100 * t53;
	dBzdr += 	-t116 * t95 * t100 * t429 - t143 * t95 * t100 * t434 / 0.2e1 + t95 * t166 * t53 - t181 * t95 * t441 * t53 + t179 * t100 * t53;
	dBzdphi += 	-t207 * t95 * t100 * t429 - t225 * t95 * t100 * t434 / 0.2e1 + t95 * t246 * t53 - t259 * t95 * t441 * t53 + t257 * t100 * t53;
	dBzdz += 	t13 * t95 * t100 * t429 - t301 * t95 * t100 * t434 / 0.2e1 + t95 * t317 * t53 - t322 * t95 * t441 * t53;

	return;	
}


// here the magnetic field for all current bars in the race track coils is generated and added to the current field
void BarRaceTrack(long double r_current, long double phi_current, long double z_current, long double I_bar)
{
	// calculate flux density at point r_current, phi_current, z_current
	// for all single current bars
	for(int cobar=1;cobar<=12;cobar++)
	{
		StraightWireField(r_current,phi_current,z_current, I_bar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar]);					
	}
	// then for the center bar (4x current)   
	StraightWireField(r_current,phi_current,z_current, 4.0*I_bar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13]);
}
