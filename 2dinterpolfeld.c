// this file contains functions to read in E and B fields, interpolate
// them and calculate derivatives

#include "main.h"


// get the size of the array ......................
void GetDim(int *m, int *n)
{
	int intmuell;
	long double muell, rtemp = 0, ztemp = 0, Brtemp = 0, Bphitemp=0, Bztemp=0, Ertemp=0, Ephitemp=0, Eztemp=0;	
	char str[1024];
	string string[1024];
	long double rtemptmp=-1e10, ztemptmp=-1e10;  // memory for r and z values in table
	int ri=1, zi=1;   // indizes
	
	if(FIN != NULL) 
		fclose(FIN);
	char *path;
	path=(char*)malloc((inpathlength+14)*sizeof(char));
	sprintf(path,"%s/fieldval.tab",inpath);
	FIN = fopen(path,mode_r);
	free(path);
	
	// in der ersten zeile stehen die dimensionen auch drin
	sscanf(fgets(str,1024,FIN),"%d %d %d %LG",m,&intmuell,n,&muell);
	printf("The header says: \nThe arrays are %i by %i.\n\n",*m,*n);
	fprintf(LOGSCR,"The header says: \nThe arrays are %i by %i.\n\n",*m,*n);
	
	// discard next 10 lines because of header stuff
	for(i=1; i<=10; i++){
		fgets(str,1024,FIN);
	}		
		
	// read all lines until the end of the file
	do
	{	
		sscanf(fgets(str,1024,FIN),"%LG %LG %LG %LG %LG %LG %LG %LG %LG",&rtemp,&muell,&ztemp,&Brtemp,&Bphitemp,&Bztemp,&Ertemp,&Ephitemp,&Eztemp);		
		rtemp = rtemp * lengthconv; ztemp = ztemp * lengthconv; Brtemp = Brtemp * Bconv; Bphitemp = Bphitemp * Bconv;    // Einheiten
		Bztemp = Bztemp * Bconv; Ertemp = Ertemp * Econv; Ephitemp = Ephitemp * Econv; Eztemp = Eztemp * Econv;          // ausgleichen
	
		/*if((zi<20)&&(ri==1))
			printf("%LG %LG %LG %LG %LG %LG %LG %LG %LG\n",rtemp,muell,ztemp,Brtemp,Bphitemp,Bztemp,Ertemp,Ephitemp,Eztemp);*/
		
		
		if((rtemptmp==-1e10)&&(ztemptmp==-1e10)) //first line
		{
			rtemptmp=r_mi=r_ma=rtemp;  // set maximum values to current
			ztemptmp=z_mi=z_ma=ztemp;		// set maximum values to current		
		}
	
			
		// Index ri nur erhöhen  wenn sich der r-Wert ändert
		if (rtemp > rtemptmp){
			ri++;		
			if(rtemp>r_ma)
				r_ma=rtemp;
		}
		// Index zi nur erhöhen wenn sich der z-Wert ändert
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
	}while (!feof(FIN));
	
	*m=ri; *n=zi; // set to maximum index values found in file
	
	printf("In reality they are %d by %d.\n\n",*m,*n);
	fprintf(LOGSCR,"In reality they are %d by %d.\n\n",*m,*n);
	printf("The r values go from %LG to %LG\n\n",r_mi,r_ma);
	fprintf(LOGSCR,"The r values go from %LG to %LG.\n\n",r_mi,r_ma);
	printf("The z values go from %LG to %LG.\n\n",z_mi,z_ma);
	fprintf(LOGSCR,"The z values go from %LG to %LG.\n\n",z_mi,z_ma);	
	rmin=r_mi;
	rmax=r_ma;
	zmin=z_mi;
	zmax=z_ma;
	return ;
}


// Zu den Arrays BrTab, BzTab, BphiTab werden die entsprechenden Ableitungen berechnet:
// z.B: BrTab1 = dBrTab / dr
void CalcDeriv(long double rind[], long double zind[], long double **BrTab, long double **BzTab, long double **BphiTab, long double **BrTab1, long double **BrTab2,long double **BrTab12, long double **BzTab1, long double **BzTab2, long double **BzTab12, long double **BphiTab1, long double **BphiTab2,long double **BphiTab12,int m, int n)
{

	int j, k;

	for (j=2; j < m; j++)
	{
			for (k=2; k < n; k++)
			{
			BrTab1[j][k] = ( BrTab[j+1][k] - BrTab[j-1][k] ) / ( rind[j+1] - rind[j-1] );   // derivative in r direction
			BrTab2[j][k] = ( BrTab[j][k+1] - BrTab[j][k-1] ) / ( zind[k+1] - zind[k-1] );  // derivative in z direction
			BrTab12[j][k] = ( BrTab[j+1][k+1] - BrTab[j+1][k-1] - BrTab[j-1][k+1] + BrTab[j-1][k-1]) / ((rind[j+1] - rind[j-1]) * (zind[k+1] - zind[k-1]));  // cross derivative
			  
			BzTab1[j][k] = ( BzTab[j+1][k] - BzTab[j-1][k] ) / ( rind[j+1] - rind[j-1] ); 
			BzTab2[j][k] = ( BzTab[j][k+1] - BzTab[j][k-1] ) / ( zind[k+1] - zind[k-1] );
			BzTab12[j][k] = ( BzTab[j+1][k+1] - BzTab[j+1][k-1] - BzTab[j-1][k+1] + BzTab[j-1][k-1] ) / ( (rind[j+1] - rind[j-1]) * (zind[k+1] - zind[k-1]) );
			  
			BphiTab1[j][k] = ( BphiTab[j+1][k] - BphiTab[j-1][k] ) / ( rind[j+1] - rind[j-1] );
			BphiTab2[j][k] = ( BphiTab[j][k+1] - BphiTab[j][k-1] ) / ( zind[k+1] - zind[k-1] );
			BphiTab12[j][k] = ( BphiTab[j+1][k+1] - BphiTab[j+1][k-1] - BphiTab[j-1][k+1] + BphiTab[j-1][k-1] ) / ( (rind[j+1] - rind[j-1]) * (zind[k+1] - zind[k-1]) );
        }
	}

	
	printf("\nDerivatives calculated.  \n");
	fprintf(LOGSCR,"\nDerivatives calculated.  \n");
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
	
	printf("In readWert they are %d by %d.\n\n",*m,*n);
	fprintf(LOGSCR,"In readWert they are %d by %d.\n\n",*m,*n);

	
	// close file if it already open
	if(FIN != NULL)  
		fclose(FIN);
	char *path;
	path=(char*)malloc((inpathlength+14)*sizeof(char));
	//snprintf(path,pathlen,"%s/fieldval.tab",inpath);
	sprintf(path,"%s/fieldval.tab",inpath);	
	FIN = fopen(path,mode_r);
	free(path);
	
	// discard first eleven lines with header
	for(i=1; i<=11; i++)
	{
		fgets(str,1024,FIN);
	}

	
	printf("\nreading fieldval.tab \n");
	
   // Schleife über alle Zeilen der Input Datei
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
		
		// Index ri nur erhöhen  wenn sich der r-Wert ändert
		if ((rtemp > rind[ri])&&(NrValueLines!=1)&&(ri<mtmp)){
			ri++;
			zi=1;             // immer wenn sich der r-Wert ändert springt z auch auf Eins
		}

		// Index zi nur erhöhen wenn sich der z-Wert ändert
		if ((ztemp > zind[zi])&&(NrValueLines!=1)&&(zi<ntmp)){
			zi++;
		}
		
		// den Wert nur in die Array schreiben wenn sich der Index geändert hat oder wenn nur eine zeile gelesen wurde bisher
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
// der Array zwischen den Indizes n und n+1 liegt  (jlo ist ein first guess des indexes und der rückgabewert)
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

// Koeffizienten für die spätere bikubische Interpolation berechnen
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
// Größe der Arrays für das Einlesen der Tabelle bestimmen
   long double EnTemp;
	char *path;
	path=(char*)malloc((inpathlength+14)*sizeof(char));
	sprintf(path,"%s/fieldval.tab",inpath);
	FIN = fopen(path,mode_r);
	free(path);
	if (FIN == NULL) 
		exit(-1);        // Fehlerbehandlung
	GetDim(&m, &n);


	// Vektoren der r und z Komponenten, Matrizen für B-Feld Komponenten und deren Ableitungen kreieren
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

	// Vektoren der r und z Komponenten, Matrizen für B-Feld Komponenten und deren Ableitungen kreieren
	if ((protneut==PROTON)||(protneut==ELECTRONS)){
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

	if ((protneut==PROTON)||(protneut==ELECTRONS)){
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
			//StorVolrmin=0.129,StorVolrmax=0.488,StorVolzmin=0.01, StorVolzmax=1.345
			if((Babsmin > Babsmaxtmp))
			{
				Babsmin = Babsmaxtmp;				
			}
			if((EnTemp < Emin_n)&&(rind[j]<StorVolrmax)&&(rind[j]>StorVolrmin)&&(zind[k]>0)&&(zind[k]<StorVolzmax))
			{
				Emin_n = EnTemp;
				rBabsmin = rind[j];
				zBabsmin = zind[k];
			}
			if ((protneut==PROTON)||(protneut==ELECTRONS))
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
	//long double dBrdphi = 0.0, dBzdphi = 0.0, dBphidphi = 0.0, x2, dy;
	//int ordnung = 4, ordleft = 1, ordright = 2;   // if ordnung even ordleft = (ordnung/2)-1 , if odd ordleft = (ordnung-1)/2 
	
	//fprintf(LOGSCR,"BInterpol got: r,phi,z (%.12LG / %.12LG / %.12LG)\n",r_n,phi,z_n);
	
	if ((r_n >= r_mi) && (r_n <= r_ma) && (z_n >= z_mi) && (z_n <= z_ma))
	{
		//Indizes indr und indz in der Array finden, der zu r_n und z_n passt
		//hunt(rind, m, r_n, &indr);
		//hunt(zind, n, z_n, &indz);		
		//printf("hunt: indr =  %i, indz = %i \n ",indr, indz);
		
		// new calculation of indexes
		indr = 1+(int) ((r_n- conv_rA)/(conv_rB));    // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!		
		indz = 1+(int) ((z_n-conv_zA) / (conv_zB));  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
		//printf("new: indr =  %i, indz = %i \n",indr, indz);
					
		// Abstände der Koordinaten voneinander
		//dr = rind[indr + 1] - rind [indr];
		//dz = zind[indz + 1] - zind [indz];
		//cout << "Indizes: " << indr << " " << indz << "\n";

		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3) ){			
									
			// 2D polynomial interpolation
			// for field interpolation tests 
			/*if(indr-ordleft<2) indr = 2+ordleft*2;
			if(indr+ordright>(m-1)) indr = m-(1+ordright)*2;
			if(indz-ordleft<2) indz = 2+ordleft*2;
			if(indz+ordright>(n-1)) indz = n-(1+ordright)*2;	*/
			
			if(BFeldSkal!=0)
			{
				// bicubic interpolation
				bcuint_new(indr, indz, Brc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Br, &dBrdr, &dBrdz);
				//bcuint_new(indr, indz, Bphic, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bphi, &dBphidr, &dBphidz);
				bcuint_new(indr, indz, Bzc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bz, &dBzdr, &dBzdz);
			}
			else if (BFeldSkal == 0)
			{
				Br = 0.0;
				dBrdr = 0.0;
				//dBrdphi = 0.0;
				dBrdz = 0;
				//Bphi=0.0;
				//dBphidr = 0.0;
				//dBphidphi = 0.0;
				//dBphidz = 0.0;
				Bz = 0.0;
				dBzdr = 0.0;
				//dBzdphi = 0.0;
				dBzdz = 0.0;
			}			
						
		} 
		else 
		{
			printf("\n The array index has left boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n ,indz);
			fprintf(LOGSCR,"\n The array index has left boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n, indz);
			OutputState(ystart,1);
		}
				
			Br = BFeldSkal * Br;
			dBrdr = BFeldSkal * dBrdr;
			dBrdphi = 0.0;
			dBrdz = BFeldSkal * dBrdz;
			//Bphi = RodFieldMultiplicator * BFeldSkal * Bphi;
			Bphi=0.0;
			//dBphidr = RodFieldMultiplicator * BFeldSkal * dBphidr;
			dBphidr = 0.0;
			dBphidphi = 0.0;
			//dBphidz = RodFieldMultiplicator * BFeldSkal * dBphidz;
			dBphidz = 0.0;
			Bz = BFeldSkal * Bz;
			dBzdr = BFeldSkal * dBzdr;
			dBzdphi = 0.0;
			dBzdz = BFeldSkal * dBzdz;
		
		
		// add fields and field derivations of racetrack coils
		
		if((Ibar>0)&&(RodFieldMultiplicator>0))
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
		if(Racetracks==2)
		{			
			// new but slow implementation of racetracks including horizontal parts
			BarRaceTrack(r_n, phi, z_n, RodFieldMultiplicator*Ibar);   // here B field is only calculated for charged particles and derivatives only for neutrons
		}
			
		}  		
		
		
		
		Bws = sqrtl(Br*Br+Bz*Bz+Bphi*Bphi);
		
		if (Bws>1e-31)
		{			
			dBdr   = (Br*dBrdr + Bphi*dBphidr + Bz*dBzdr)  /Bws;
			//dBdr = (Br*dBrdr +  Bz*dBzdr)  /Bws;
			dBdz   = (Br*dBrdz + Bphi*dBphidz + Bz*dBzdz)  /Bws;
			//dBdphi = 0.0; // intrinsic cylindrical symmetry => not any more!!!
			dBdphi = (Br*dBrdphi + Bphi*dBphidphi + Bz*dBzdphi)/Bws;
		}
		else
		{
			//if(BFeldSkal!=0&&RodFieldMultiplicator!=0)
			//	OutputState(ystart,1);
			Br = 0.;
			dBrdr = 0.;
			dBrdz = 0.;
			Bphi = 0.;
			dBphidr = 0.;
			dBphidz = 0.;
			Bz = 0.;
			dBzdr = 0.;
			dBzdz = 0.;
			dBdr   = 0.;
			dBdz   = 0.;
			dBdphi = 0.;
		}		
		
	}
	else
	{
		//char bloed;
		printf("The particle has entered forbidden space: r=%.17LG, z=%.17LG !!! I will terminate it now! \n", r_n, z_n);
		fprintf(LOGSCR,"The particle has entered forbidden space: r=%.17LG, z=%.17LG !!! I will terminate it now! \n", r_n, z_n);
		
		//OutputState(ystart,1);
		
		dBrdr=0.0;dBrdz=0.0;
		dBzdr=0.0;dBzdz=0.0;
		Br=0.0;Bphi=0.0;Bz=0.0;
		dBdr=0.; dBdphi=0.; dBdz=0.;
		kennz = 2;    // something's wrong
		stopall=1;		
	}
	//if(Bphi==0)
	//			cout << "hier";
			//	OutputState(ystart,1);
}


void BInterpolOld(long double r_n, long double phi, long double z_n){
	yyy=dvector(1,4);
	yyy1=dvector(1,4);
	yyy2=dvector(1,4);
	yyy12=dvector(1,4);
	
	if ((r_n >= r_mi) && (r_n <= r_ma) && (z_n >= z_mi) && (z_n <= z_ma)){
		//Indizes indr und indz in der Array finden, der zu r_n und z_n passt
		//hunt(rind, m, r_n, &indr);
		//hunt(zind, n, z_n, &indz);		
		//printf("hunt: indr =  %i, indz = %i \n ",indr, indz);
		
		// new calculation of indexes
		indr = 1+(int) ((r_n- conv_rA)/(conv_rB));    // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
		
		// for field interpolation tests if(indr % 2 != 0) indr-=1; 
		indz = 1+(int) ((z_n-conv_zA) / (conv_zB));  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
				// for field interpolation tests if(indz % 2 != 0) indz-=1;
		//printf("new: indr =  %i, indz = %i \n",indr, indz);
		//scanf("%LG",&muell);
		
		
				
		// Abstände der Koordinaten voneinander
		//dr = rind[indr + 1] - rind [indr];
		//dz = zind[indz + 1] - zind [indz];

		//cout << "Indizes: " << indr << " " << indz << "\n";

		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3) ){
			
					
						
			// 2D polynomial interpolation
			// for field interpolation tests 
			/*if(indr-ordleft<2) indr = 2+ordleft*2;
			if(indr+ordright>(m-1)) indr = m-(1+ordright)*2;
			if(indz-ordleft<2) indz = 2+ordleft*2;
			if(indz+ordright>(n-1)) indz = n-(1+ordright)*2;	*/
			
			if(BFeldSkal!=0)
			{
				
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
			// bicubic interpolation
			bcuint(yyy, yyy1, yyy2, yyy12, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Br, &dBrdr, &dBrdz);
			//bcuint_new(indr, indz, Brc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Br, &dBrdr, &dBrdz);
			//printf("\n Brold: %.17LG dBrdrold: %.17LG dBrdzold: %.17LG\n", Br, dBrdr, dBrdz);
			//bcuint_new(indr, indz, Bphic, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bphi, &dBphidr, &dBphidz);
			//printf("\n Brold: %.17LG dBdrold: %.17LG dBdzold: %.17LG\n", Bphi, dBphidr, dBphidz);
				
					
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
			bcuint(yyy, yyy1, yyy2, yyy12, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bz, &dBzdr, &dBzdz);
			//bcuint_new(indr, indz, Bzc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r_n, z_n, &Bz, &dBzdr, &dBzdz);
			//printf("\n Bzold: %.17LG dBdzold: %.17LG dBzdzold: %.17LG\n", Bz, dBzdr, dBzdz); 
			}
			else if (BFeldSkal == 0)
			{
				Br = 0.0;
				dBrdr = 0.0;
				//dBrdphi = 0.0;
				dBrdz = BFeldSkal * dBrdz;
				//Bphi=0.0;
				//dBphidr = 0.0;
				//dBphidphi = 0.0;
				//dBphidz = 0.0;
				Bz = 0.0;
				dBzdr = 0.0;
				//dBzdphi = 0.0;
				dBzdz = 0.0;
			}			
						
		} 
		else 
		{
			printf("\n The array index has left boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n ,indz);
			fprintf(LOGSCR,"\n The array index has left boundaries: r=%LG:%d, z=%LG:%i   Exiting!!! \n", r_n,indr,z_n, indz);
			OutputState(ystart,1);
		}
				
			Br = BFeldSkal * Br;
			dBrdr = BFeldSkal * dBrdr;
			dBrdphi = 0.0;
			dBrdz = BFeldSkal * dBrdz;
			//Bphi = RodFieldMultiplicator * BFeldSkal * Bphi;
			Bphi=0.0;
			//dBphidr = RodFieldMultiplicator * BFeldSkal * dBphidr;
			dBphidr = 0.0;
			dBphidphi = 0.0;
			//dBphidz = RodFieldMultiplicator * BFeldSkal * dBphidz;
			dBphidz = 0.0;
			Bz = BFeldSkal * Bz;
			dBzdr = BFeldSkal * dBzdr;
			dBzdphi = 0.0;
			dBzdz = BFeldSkal * dBzdz;
		
		
		// add fields and field derivations of racetrack coils
		
		if((Ibar>0)&&(RodFieldMultiplicator>0))
			{
			Racetrack(r_n, phi, z_n, -0.424264, -0.424264, Ibar);
			Racetrack(r_n, phi, z_n, 0.424264, -0.424264, Ibar); 
			Racetrack(r_n, phi, z_n, -0.424264, 0.424264, Ibar);
			Racetrack(r_n, phi, z_n, 0.424264, 0.424264, Ibar); 
			CenterCurrent(r_n, phi, z_n, -4.0l*Ibar);
			}  
		
		Bws = sqrtl(Br*Br+Bz*Bz+Bphi*Bphi);
		
		if (Bws>1e-31)
		{
			
			dBdr   = (Br*dBrdr + Bphi*dBphidr + Bz*dBzdr)  /Bws;
			//dBdr = (Br*dBrdr +  Bz*dBzdr)  /Bws;
			dBdz   = (Br*dBrdz + Bphi*dBphidz + Bz*dBzdz)  /Bws;
			//dBdphi = 0.0; // intrinsic cylindrical symmetry => not any more!!!
			dBdphi = (Br*dBrdphi + Bphi*dBphidphi + Bz*dBzdphi)/Bws;
		}
		else
		{
			Br = 0.;
			dBrdr = 0.;
			dBrdz = 0.;
			Bphi = 0.;
			dBphidr = 0.;
			dBphidz = 0.;
			Bz = 0.;
			dBzdr = 0.;
			dBzdz = 0.;
			dBdr   = 0.;
			dBdz   = 0.;
			dBdphi = 0.;
		}
		
		
	}
	else
	{
		//char bloed;
		printf("The particle has entered forbidden space: r=%.17LG, z=%.17LG !!! I will terminate it now! \n", r_n, z_n);
		fprintf(LOGSCR,"The particle has entered forbidden space: r=%.17LG, z=%.17LG !!! I will terminate it now! \n", r_n, z_n);
		
		//OutputState(ystart,1);
		
		dBrdr=0.0;dBrdz=0.0;
		dBzdr=0.0;dBzdz=0.0;
		Br=0.0;Bphi=0.0;Bz=0.0;
		dBdr=0.; dBdphi=0.; dBdz=0.;
		kennz = 2;    // something's wrong
		stopall=1;		
		
	}
	free_dvector(yyy,1,4);
	free_dvector(yyy1,1,4);
	free_dvector(yyy2,1,4);
	free_dvector(yyy12,1,4);
	
}


void EInterpol(long double r_n, long double phi, long double z_n){
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

		// Abstände der Koordinaten voneinander
		//dr = rind[indr + 1] - rind [indr];
		//dz = zind[indz + 1] - zind [indz];

		// cout << "Indizes: " << indr << " " << indz << "\n";

		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3)){
			// Rechteck mit den Werten und deren Ableitungen befüllen und Interpolieren
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

			// Rechteck mit den Werten und deren Ableitungen befüllen und Interpolieren
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

			// Rechteck mit den Werten und deren Ableitungen befüllen und Interpolieren
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
			csleep(10);
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
        if (z_n > zmax)
           kennz=4;    // escape to top
        if ((r_n < rmin+wandinnen) || (r_n > rmax+wanddicke) || (z_n < zmin+wanddicke) )
				if(!reflekt) kennz=3;   // walls
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
				
	return;
}


// TH: Dokumentation !!!!!!!
char *stripws(char *string,char *retstr){
	int pws_i=0,pws_n=0,pws_ws=0;
	
	for(pws_i=0;pws_i<strlen(string)-1;pws_i++){
		if((string[pws_i]==' ')||(string[pws_i]=='	')){       //????????
			if(pws_ws==0) retstr[pws_n++]=' ';
			pws_ws=1;
		}else{
			if(pws_ws==1){
				if(string[pws_i]!='-'){
					retstr[pws_n++]=' ';
				}
			}
			pws_ws=0;
			retstr[pws_n++]=string[pws_i];
		}
	}
	retstr[pws_n]='\0';
	fprintf(LOGSCR,"stringlaenge: %i\n",strlen(retstr));
	fprintf(LOGSCR,"%s\n",retstr);
	
	return retstr;
}


// produce zero Bfield
long double Bnull(long double r,long double phi,long double z)
{
	long double dBbdphi, dBphidr, dBphidz, dBphidphi;
  
	dBphidphi=0;
	R=0;


	
	dBrdr=0.0;
	dBrdphi=0.0;
	dBrdz=0.0;
	dBphidr=0.0;
	dBphidphi=0.0;
	dBphidz=0.0;
	dBzdr=0.0;
	dBzdphi=0.0;
	dBzdz=0.0;
	dBbdphi=0.0;

	dBphidr= 0.0;
	dBphidz=0.0;
	dBdphi= 0.0;
	Bws= 1e-31;

	if (Bws>1e-19)
	{
		dBdr   = (Bphi*dBphidr   + Br*dBrdr + Bz*dBzdr)  /Bws;
		dBdz   = (Bphi*dBphidz   + Br*dBrdz + Bz*dBzdz)  /Bws;
		dBdphi = (Bphi*dBphidphi )/Bws;
	}
	else
	{
		dBdr=0; dBdphi=0; dBdz=0;
	}
	
	return 0;
}



void Preinterpol(int p)
{
	int indr, indz, perc=0; 
	long double d1, d2;
	long double **ctemp;
	
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
	
	printf("\n Done with Preinterpolation!\n");
	return;
}

// 2 dimensional polynomial interpolation of order m-1, n-1
// **ya has to be the m*n array within which the point [x1,x2] lies
void polin2(long double x1a[],long double x2a[], long double **ya, int m, int n, long double x1, long double x2, long double *y, long double *dy)
{
	int j;
	long double *ymtmp,*vector(int nl,int nh);
	//void polint(),free_vector();

	ymtmp=vector(1,m);
	for (j=1;j<=m;j++) {
		BFpolint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	BFpolint(x1a,ymtmp,m,x1,y,dy);
	free_vector(ymtmp,1,m);
}

// produce linearly rising B-field in z-direction and linearly changing in time
long double Banalytic(long double r,long double phi,long double z, long double t)
{
	long double dBphidr, dBphidz, dBphidphi;
  
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
        	
	Bws= BFeldSkal * fabsl(Br);
	
		
	if (Bws>1e-31)
	{
		Br = BFeldSkal * Br;
		dBrdr = BFeldSkal * dBrdr;
		dBrdz = BFeldSkal * dBrdz;
		Bphi = RodFieldMultiplicator * BFeldSkal * Bphi;
		dBphidr = RodFieldMultiplicator * BFeldSkal * dBphidr;
		dBphidz = RodFieldMultiplicator * BFeldSkal * dBphidz;
		dBphidphi = 0.0;
		Bz = BFeldSkal * Bz;
		dBzdr = BFeldSkal * dBzdr;
		dBzdz = BFeldSkal * dBzdz;
		dBdr   = (Bphi*dBphidr   + Br*dBrdr   + Bz*dBzdr)  /Bws;
		dBdz   = (Bphi*dBphidz   + Br*dBrdz   + Bz*dBzdz)  /Bws;
		dBdphi = 0.0;
	}
	else
	{
		dBdr=0; dBdphi=0; dBdz=0;
	}
	
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
      dBrdrv = vorfaktor*(ly*t2-lx*t5)*t25-t8*t33*(2.0*t34-2.0*t35+2.0*t36-2.0*
t37)-t8*t23*t42;
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
      dBphidrv = vorfaktor*(-2.0*t36-t37-2.0*t34+t35)*t61-t53*t69*(2.0*t54*t2+
2.0*t56*t5)/2.0-t53*t60*t42;
      dBphidphiv = vorfaktor*(-t3-t6)*t61-t53*t69*(-2.0*t54*r*t5+2.0*t56*r*t2)/
2.0;


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
	Bphi = Bphi + vorfaktor/r;
	
	
	
	// dBphidr
	dBphidr = dBphidr - vorfaktor/(r*r);
		
	
	// dBphidphi =0
		
	
	return;
	
}

// produce B-field at position r,phi,z of a straight wire with current I_rt 
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the ends of the straight wire 
void StraightWireField(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z)
{
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	long double SW1x = cosl(SW1phi) * SW1r, SW1y = sinl(SW1phi) * SW1r;   // cart coord of first vector on line
	long double SW2x = cosl(SW2phi) * SW2r, SW2y = sinl(SW2phi) * SW2r;   // cart coord of second vector on line
	long double x = cosl(phi) * r, y = sinl(phi) * r;													// cart coord of position vector
	
	// cross product of vector (current position to vector on current 1) with (vector along current line)
	// determines direction of magnetic field vector
	long double BSWdirx,BSWdiry,BSWdirz;
	KartCrossProd(x-SW1x,y-SW1y,z-SW1z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z,&BSWdirx,&BSWdiry,&BSWdirz);
	
	long double absBSWdir = AbsValueCart(BSWdirx,BSWdiry,BSWdirz);
		if(absBSWdir==0) cout << "Problems with absBSWdir at " << r << " " << phi << " " << z << endl;
	long double BarLength = AbsValueCart(SW2x-SW1x,SW2y-SW1y,SW2z-SW1z);
	long double dabs =  absBSWdir / BarLength;   // absolute value of distance from position vector to current line  |(point-SW1)x(SW2-SW1)|/|(SW1-SW2)|
		if(dabs==0) cout << "Problems with dabs at " << r << " " << phi << " " << z << endl;
	long double Bxtmp =    vorfaktor * BSWdirx/absBSWdir / dabs;   
	long double Bytmp =  	vorfaktor * BSWdiry/absBSWdir / dabs;  
	long double BztmpC =  vorfaktor * BSWdirz/absBSWdir / dabs;   
	
	//  consider now the finiteness of the current bar through multiplying with (sin alpha2 - sin alpha 1), with alpha 1 and 2 the angles between the normal on the line
	// through the point r,phi,z and the ends of the current bar SW1, SW2
	long double finitecorr =  	1.0 /  BarLength * (  KartScalarProd(SW2x-x,SW2y-y,SW2z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW2x-x,SW2y-y,SW2z-z) 
																				- 	KartScalarProd(SW1x-x,SW1y-y,SW1z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW1x-x,SW1y-y,SW1z-z) );
		
		Bxtmp = Bxtmp * finitecorr;
	Bytmp = Bytmp * finitecorr;
	BztmpC = BztmpC * finitecorr;
	
	long double Brtmp, Bphitmp, Bztmp;
	KartCylCoord(Bxtmp, Bytmp,BztmpC,phi,&Brtmp,&Bphitmp,&Bztmp);
	
	Br = Br + Brtmp;
	Bphi = Bphi + Bphitmp;
	Bz = Bz + Bztmp;			
	
	return;	
}

// same as above only returning r component of B
long double StraightWireField_r(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z)
{
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	long double SW1x = cosl(SW1phi) * SW1r, SW1y = sinl(SW1phi) * SW1r;   // cart coord of first vector on line
	long double SW2x = cosl(SW2phi) * SW2r, SW2y = sinl(SW2phi) * SW2r;   // cart coord of second vector on line
	long double x = cosl(phi) * r, y = sinl(phi) * r;													// cart coord of position vector
	
	// cross product of vector (current position to vector on current 1) with (vector along current line)
	// determines direction of magnetic field vector
	long double BSWdirx,BSWdiry,BSWdirz;
	KartCrossProd(x-SW1x,y-SW1y,z-SW1z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z,&BSWdirx,&BSWdiry,&BSWdirz);
	long double absBSWdir = AbsValueCart(BSWdirx,BSWdiry,BSWdirz);
	long double BarLength = AbsValueCart(SW2x-SW1x,SW2y-SW1y,SW2z-SW1z);
	long double dabs =  absBSWdir / BarLength;   // absolute value of distance from position vector to current line  |(point-SW1)x(SW2-SW1)|/|(SW1-SW2)|
	long double Bxtmp = vorfaktor * BSWdirx/absBSWdir / dabs;   
	long double Bytmp =  	vorfaktor * BSWdiry/absBSWdir / dabs;  
	long double BztmpC = vorfaktor * BSWdirz/absBSWdir / dabs;   
	
	//  consider now the finiteness of the current bar through multiplying with (sin alpha2 - sin alpha 1), with alpha 1 and 2 the angles between the normal on the line
	// through the point r,phi,z and the ends of the current bar SW1, SW2
	long double finitecorr =  	1.0 /  BarLength * (  KartScalarProd(SW2x-x,SW2y-y,SW2z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW2x-x,SW2y-y,SW2z-z) 
																				- 	KartScalarProd(SW1x-x,SW1y-y,SW1z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW1x-x,SW1y-y,SW1z-z) );
		
	Bxtmp = Bxtmp * finitecorr;
	Bytmp = Bytmp * finitecorr;
	BztmpC = BztmpC * finitecorr;
	
	long double Brtmp, Bphitmp, Bztmp;
	KartCylCoord(Bxtmp, Bytmp,BztmpC,phi,&Brtmp,&Bphitmp,&Bztmp);
	
	return Brtmp;	
}

// same as above only returning phi component of B
long double StraightWireField_phi(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z)
{
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	long double SW1x = cosl(SW1phi) * SW1r, SW1y = sinl(SW1phi) * SW1r;   // cart coord of first vector on line
	long double SW2x = cosl(SW2phi) * SW2r, SW2y = sinl(SW2phi) * SW2r;   // cart coord of second vector on line
	long double x = cosl(phi) * r, y = sinl(phi) * r;													// cart coord of position vector
	
	// cross product of vector (current position to vector on current 1) with (vector along current line)
	// determines direction of magnetic field vector
	long double BSWdirx,BSWdiry,BSWdirz;
	KartCrossProd(x-SW1x,y-SW1y,z-SW1z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z,&BSWdirx,&BSWdiry,&BSWdirz);
	long double absBSWdir = AbsValueCart(BSWdirx,BSWdiry,BSWdirz);
	long double BarLength = AbsValueCart(SW2x-SW1x,SW2y-SW1y,SW2z-SW1z);
	long double dabs =  absBSWdir / BarLength;   // absolute value of distance from position vector to current line  |(point-SW1)x(SW2-SW1)|/|(SW1-SW2)|
	long double Bxtmp = vorfaktor * BSWdirx/absBSWdir / dabs;   
	long double Bytmp =  	vorfaktor * BSWdiry/absBSWdir / dabs;  
	
	//  consider now the finiteness of the current bar through multiplying with (sin alpha2 - sin alpha 1), with alpha 1 and 2 the angles between the normal on the line
	// through the point r,phi,z and the ends of the current bar SW1, SW2
	long double finitecorr =  	1.0 /  BarLength * (  KartScalarProd(SW2x-x,SW2y-y,SW2z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW2x-x,SW2y-y,SW2z-z) 
																				- 	KartScalarProd(SW1x-x,SW1y-y,SW1z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW1x-x,SW1y-y,SW1z-z) );
		
	Bxtmp = Bxtmp * finitecorr;
	Bytmp = Bytmp * finitecorr;	
	
	long double Bphitmp = (-1) * sinl(phi) * Bxtmp + cosl(phi) * Bytmp;
	
	return Bphitmp;	
}

// same as above only returning z component of B
long double StraightWireField_z(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z)
{
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	long double SW1x = cosl(SW1phi) * SW1r, SW1y = sinl(SW1phi) * SW1r;   // cart coord of first vector on line
	long double SW2x = cosl(SW2phi) * SW2r, SW2y = sinl(SW2phi) * SW2r;   // cart coord of second vector on line
	long double x = cosl(phi) * r, y = sinl(phi) * r;													// cart coord of position vector
	
	// cross product of vector (current position to vector on current 1) with (vector along current line)
	// determines direction of magnetic field vector
	long double BSWdirx,BSWdiry,BSWdirz;
	KartCrossProd(x-SW1x,y-SW1y,z-SW1z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z,&BSWdirx,&BSWdiry,&BSWdirz);
	long double absBSWdir = AbsValueCart(BSWdirx,BSWdiry,BSWdirz);
	long double BarLength = AbsValueCart(SW2x-SW1x,SW2y-SW1y,SW2z-SW1z);
	long double dabs =  absBSWdir / BarLength;   // absolute value of distance from position vector to current line  |(point-SW1)x(SW2-SW1)|/|(SW1-SW2)|
	
	long double BztmpC = vorfaktor * BSWdirz/absBSWdir / dabs;   
	
	//  consider now the finiteness of the current bar through multiplying with (sin alpha2 - sin alpha 1), with alpha 1 and 2 the angles between the normal on the line
	// through the point r,phi,z and the ends of the current bar SW1, SW2
	long double finitecorr =  	1.0 /  BarLength * (  KartScalarProd(SW2x-x,SW2y-y,SW2z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW2x-x,SW2y-y,SW2z-z) 
																				- 	KartScalarProd(SW1x-x,SW1y-y,SW1z-z,SW2x-SW1x,SW2y-SW1y,SW2z-SW1z) / AbsValueCart(SW1x-x,SW1y-y,SW1z-z) );
		
	
	BztmpC = BztmpC * finitecorr;	
	return BztmpC;	
}


// here the r component magnetic field for all current bars in the race track coils is returned
long double BarRaceTrack_Br(long double r_current, long double phi_current, long double z_current, long double I_bar)
{
	long double Brvorl = 0;
	
	// first calculate flux density at point r_current, phi_current, z_current
	// for all single current bars
	for(int cobar=1;cobar<=12;cobar++)
	{
		Brvorl= Brvorl + StraightWireField_r(r_current,phi_current,z_current, Ibar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar]);			
		
	}
	// then for the center bar (4x current)   
	return Brvorl+StraightWireField_r(r_current,phi_current,z_current, 4.0*Ibar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13]);
	
}

// the same just with a 4 vector
long double BarRaceTrack_Br(long double x[4])
{
		return BarRaceTrack_Br(x[0],x[1],x[2],x[3]);
}

// here the phi component magnetic field for all current bars in the race track coils is returned
long double BarRaceTrack_Bphi(long double r_current, long double phi_current, long double z_current, long double I_bar)
{
	long double Bphivorl = 0;
	
	// first calculate flux density at point r_current, phi_current, z_current
	// for all single current bars
	for(int cobar=1;cobar<=12;cobar++)
	{
		Bphivorl= Bphivorl + StraightWireField_phi(r_current,phi_current,z_current, Ibar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar]);			
		
	}
	// then for the center bar (4x current)   
	return Bphivorl+StraightWireField_phi(r_current,phi_current,z_current, 4.0*Ibar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13]);
	
}

// the same just with a 4 vector
long double BarRaceTrack_Bphi(long double x[4])
{
		return BarRaceTrack_Bphi(x[0],x[1],x[2],x[3]);
}

// here the phi component magnetic field for all current bars in the race track coils is returned
long double BarRaceTrack_Bz(long double r_current, long double phi_current, long double z_current, long double I_bar)
{
	long double Bzvorl = 0;
	
	// first calculate flux density at point r_current, phi_current, z_current
	// for all single current bars
	for(int cobar=1;cobar<=12;cobar++)
	{
		Bzvorl= Bzvorl + StraightWireField_z(r_current,phi_current,z_current, Ibar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar]);			
		
	}
	// then for the center bar (4x current)   
	return Bzvorl+StraightWireField_z(r_current,phi_current,z_current, 4.0*Ibar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13]);
	
}

// the same just with a 4 vector
long double BarRaceTrack_Bz(long double x[4])
{
		return BarRaceTrack_Bz(x[0],x[1],x[2],x[3]);
}


// here the magnetic field for all current bars in the race track coils is generated and added to the current field
void BarRaceTrack(long double r_current, long double phi_current, long double z_current, long double I_bar)
{
	long double dBrd[3]={0},dBphid[3]={0},dBzd[3]={0},position[4]={r_current,phi_current,z_current,I_bar},err;
	// first calculate flux density at point r_current, phi_current, z_current
	// for all single current bars
	
	// B field vector itself not necessary for neutron
	if(protneut!=NEUTRON)
	{
		for(int cobar=1;cobar<=12;cobar++)
		{
			//StraightWireField(r_current,phi_current,z_current, Ibar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar]);			
		}
		// then for the center bar (4x current)   
		StraightWireField(r_current,phi_current,z_current, 4.0*Ibar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13]);
	}
	
	// derivatives only necessary for neutron
	if(protneut==NEUTRON)
	{
		// differentiate all three field components with respect to r,phi,z
		for(int var=0;var<=2;var++)
		{
			dBrd[var] = dfridr3D(BarRaceTrack_Br, var,position, 0.01, &err);
			dBphid[var] = dfridr3D(BarRaceTrack_Bphi, var,position, 0.01, &err);
			dBzd[var] = dfridr3D(BarRaceTrack_Bz, var,position, 0.01, &err);
		}
	
		dBrdr = dBrdr + dBrd[0];
		dBrdphi = dBrdphi + dBrd[1];
		dBrdz = dBrdz + dBrd[2];
	
		dBphidr = dBphidr + dBphid[0];
		dBphidphi = dBphidphi + dBphid[1];
		dBphidz = dBphidz + dBphid[2];
	
		dBzdr = dBzdr + dBzd[0];
		dBzdphi = dBzdphi + dBzd[1];
		dBzdz = dBzdz + dBzd[2];
	
	}
		
}




// partial differentiation of a scalar function in 3 dimensions and a 4th variable used for straight bar current right now
// n: variable of differentialtion 0,1,2
#define CON 1.4 //Stepsize is decreased by CON at each iteration.
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10 //Sets maximum size of tableau.
#define SAFE 2.0 //Return when error is SAFE worse than the best so far.
long double dfridr3D(long double (*func)(long double x[4]), int n, long double x[4], long double h, long double *err)
//Returns the derivative of a function func at a point x by Ridders’ method of polynomial
//extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but
//rather should be an increment in x over which func changes substantially. An estimate of the
//error in the derivative is returned as err.
{
int i,j;
long double errt,fac,hh,**a,ans,xtmpp[4],xtmpm[4];
if (h == 0.0) nrerror("h must be nonzero in dfridr.");
a=matrix(1,NTAB,1,NTAB);
hh=h;
for(int co=0;co<=3;co++)
  xtmpp[co]=xtmpm[co]=x[co];

xtmpp[n]=x[n]+hh;
xtmpm[n]=x[n]-hh;
a[1][1]=((*func)(xtmpp)-(*func)(xtmpm))/(2.0*hh);
*err=BIG;
for (i=2;i<=NTAB;i++) {
//Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of
//extrapolation.
hh /= CON;
a[1][i]=((*func)(xtmpp)-(*func)(xtmpm))/(2.0*hh); //Try new, smaller stepsize
fac=CON2; 
for (j=2;j<=i;j++) { //Compute extrapolations of various orders, requiring
				//no new function evaluations.
a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
fac=CON2*fac;
errt=fmaxl(fabsl(a[j][i]-a[j-1][i]),fabsl(a[j][i]-a[j-1][i-1]));
	
//The error strategy is to compare each new extrapolation to one order lower, both
//at the present stepsize and the previous one.
if (errt <= *err) { //If error is decreased, save the improved answer.
*err=errt;
ans=a[j][i];
}
}
if (fabsl(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
//If higher order is worse by a significant factor SAFE, then quit early.
}
free_matrix(a,1,NTAB,1,NTAB);
return ans;
}
