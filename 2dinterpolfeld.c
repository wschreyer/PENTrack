// this file contains functions to read in E and B fields, interpolate
// them and calculate derivatives

#include "main.h"

long double lengthconv = 0.01 , Bconv = 1e-4, Econv = 1e2;    // Einheiten aus field Tabelle (cgs) und Programm (si) abgleichen 

// constructor
TabField::TabField(const char *tabfile){	
	rind = NULL;
	zind = NULL;
	Brc = NULL;
	Bzc = NULL;
	Erc = NULL;
	Ezc = NULL;
	long double **BrTab = NULL, **BzTab = NULL;	// Br/Bz values
	long double **ErTab = NULL, **EzTab = NULL;	// Er/Ez values
	ReadTabFile(tabfile,BrTab,BzTab,ErTab,EzTab);
	
	CheckTab(BrTab,BzTab,ErTab,EzTab);

	Log("\nStart Preinterpolation!\n");	
	if (BrTab){
		if (BFeldSkalGlobal != 0){
			PreInterpol(Brc,Bzc,BrTab,BzTab);
		}
		free_dmatrix(BrTab,1,m,1,n);
		free_dmatrix(BzTab,1,m,1,n);		
	}	
	if (ErTab){
		if (EFeldSkal != 0 && (protneut != NEUTRON || decay.on == 2)){ 
			PreInterpol(Erc,Ezc,ErTab,EzTab);
		}
		free_dmatrix(ErTab,1,m,1,n);
		free_dmatrix(EzTab,1,m,1,n);		
	}
	Log("Done with Preinterpolation!\n");
}


// read table file given in the constructor parameter
void TabField::ReadTabFile(const char *tabfile, long double **&BrTab, long double **&BzTab, long double **&ErTab, long double **&EzTab){
	ifstream FIN(tabfile, ifstream::in);
	Log("\nReading %s\n",tabfile);
	int intval;
	string line;
	FIN >> m >> intval >> n;
	rind = dvector(1,m);	// coordinates
	zind = dvector(1,n);
	
	getline(FIN,line);
	getline(FIN,line);
	getline(FIN,line);
	bool skipy = true;
	if (line.substr(0,12) == " 2 Y [LENGU]")  skipy = false;	
	getline(FIN,line);
	if (!skipy) getline(FIN,line);

	if (line.find("[FLUXU]") != string::npos){	// file contains B field?
		BrTab = dmatrix(1,m,1,n);
		BzTab = dmatrix(1,m,1,n);
		getline(FIN,line);
		if (!skipy) getline(FIN,line);
		getline(FIN,line);
	}

	if (line.find("[ELECU]") != string::npos){	// file contains E field?		
		ErTab = dmatrix(1,m,1,n);
		EzTab = dmatrix(1,m,1,n);
		getline(FIN,line);
		if (!skipy) getline(FIN,line);
		getline(FIN,line);		
	}

	if (!FIN || line.substr(0,2) != " 0"){
		Log("%s not found or corrupt! Exiting...\n",tabfile);	
		exit(-1);
	}
	
	int ri = 1,zi = 0, perc = 0;
	long double r, z, Br, Bz, Er, Ez, dmuell;
	while (FIN){
		FIN >> r;
		if (!skipy) FIN >> dmuell;
		FIN >> z;
		if (!FIN) break;
		r *= lengthconv;
		z *= lengthconv;
		if (zi > 0 && r != rind[ri]){
			ri++;
			zi = 1;
		}
		else zi++;
		
		// status if read is displayed
		if(100*zi*ri/(n*m)>perc)
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
		
		rind[ri] = r;
		zind[zi] = z;
		if (BrTab){
			FIN >> Br;
			if (!skipy) FIN >> dmuell;
			FIN >> Bz;
			BrTab[ri][zi] = Br*Bconv;
			BzTab[ri][zi] = Bz*Bconv;
		}
		if (ErTab){
			FIN >> Er;
			if (!skipy) FIN >> dmuell;
			FIN >> Ez;
			ErTab[ri][zi] = Er*Econv;
			EzTab[ri][zi] = Ez*Econv;
		}
	}
	printf("\n");
	if (ri != m || zi != n){
		Log("The header says the size is %i by %i, actually it is %i by %i! Exiting...\n", m, n, ri, zi);
		exit(-1);
	}
	FIN.close();
}


// print info about tablefile, calculate min/max fields and Emin_n
void TabField::CheckTab(long double **BrTab, long double **BzTab, long double **ErTab, long double **EzTab){
	//  calculate factors for conversion of coordinates to indexes  r = conv_rA + index * conv_rB	
	long double r_mi = rind[1];
	long double r_ma = rind[m];
	conv_rA = r_mi; 													// y Abschnitt von r
	conv_rB = (r_ma-r_mi)/(m-1);						// Steigung von r
	long double z_mi = zind[1];
	long double z_ma = zind[n];
	conv_zA = z_mi; 													// y Abschnitt von z
	conv_zB = (z_ma-z_mi)/(n-1);						// Steigung von z
	Log("The arrays are %d by %d.\n",m,n);
	Log("The r values go from %LG to %LG\n",r_mi,r_ma);
	Log("The z values go from %LG to %LG.\n",z_mi,z_ma);
	Log("conv_rA =  %LG, conv_rB = %LG\n",conv_rA, conv_rB);
	Log("conv_zA =  %LG, conv_zB = %LG\n",conv_zA, conv_zB);	
	// calculate distances of coordinates, normally 0.002 m
	rdist = rind[4] - rind[3];
	zdist = zind[4] - zind[3];
	Log("rdist = %LG zdist = %LG\n",rdist,zdist);
		
	long double Babsmax = -INFINITY, Babsmin = INFINITY, Babs;
	long double Eabsmax = -INFINITY, Eabsmin = INFINITY, Eabs;
	long double rBabsmin, zBabsmin, EnTemp;
	for (int j=1; j <= m; j++)
	{
		for (int k=1; k <= n; k++)
		{
			if (BrTab){
				// Calculate maximum values
				Babs = sqrt(BrTab[j][k]*BrTab[j][k] + BzTab[j][k]*BzTab[j][k]);
				EnTemp =  m_n*gravconst*zind[k] + (mu_nSI/ele_e)*Babs;
				Babsmax = max(Babs,Babsmax);
				Babsmin = min(Babs,Babsmin);
				if(EnTemp < Emin_n && InSourceVolume(rind[j],0,zind[k]))
				{
					Emin_n = EnTemp;
					rBabsmin = rind[j];
					zBabsmin = zind[k];
				}
			}
			if (ErTab)
			{
				Eabs = sqrt(ErTab[j][k]*ErTab[j][k] + EzTab[j][k]*EzTab[j][k]);
				Eabsmax = max(Eabs,Eabsmax);
				Eabsmin = min(Eabs,Eabsmin);
			}
							
		}
	}
	
	Log("The interpolation can be done from r = %LG to %LG and from z = %LG to %LG\n",rind[3],rind[m-2],zind[3],zind[n-2]);
	Log("The input table file has values of |B| from %LG T to %LG T\n and values of |E| from %LG V/m to %LG V/m\n",Babsmin,Babsmax,Eabsmin,Eabsmax);
	Log("The minimum energy a low field seeking neutron has to have is %.3LG eV (at r:%.3LG, z: %.3LG).\n",Emin_n,rBabsmin,zBabsmin);
}


// calculate derivatives of Tab and write them to Tab1/2/12
void TabField::CalcDeriv4th(long double **Tab, long double **Tab1, long double **Tab2,long double **Tab12)
{
	int j, k;
	for (j=3; j <= m-2; j++)
	{
		for (k=3; k <= n-2; k++)
		{
			Tab1[j][k] = 1 / (12*(rind[j]-rind[j-1] )) *  ( (-1)*Tab[j+2][k] + 8*Tab[j+1][k] - 8*Tab[j-1][k] + Tab[j-2][k] ) ;   // derivative in r direction
			Tab2[j][k] = 1 / (12*(zind[k]-zind[k-1] )) *  ( (-1)*Tab[j][k+2] + 8*Tab[j][k+1] - 8*Tab[j][k-1] + Tab[j][k-2] ) ;     // derivative in z direction
			Tab12[j][k] = 1 / (12*(rind[j]-rind[j-1] )*12*(zind[k]-zind[k-1] )) 
									*(	-1*  ( (-1)*Tab[j+2][k+2] + 8*Tab[j+2][k+1] - 8*Tab[j+2][k-1] + Tab[j+2][k-2] ) 
										+8*  ( (-1)*Tab[j+1][k+2] + 8*Tab[j+1][k+1] - 8*Tab[j+1][k-1] + Tab[j+1][k-2] ) 
										-8*  ( (-1)*Tab[j-1][k+2] + 8*Tab[j-1][k+1] - 8*Tab[j-1][k-1] + Tab[j-1][k-2] )  
										+    ( (-1)*Tab[j-2][k+2] + 8*Tab[j-2][k+1] - 8*Tab[j-2][k-1] + Tab[j-2][k-2] )     );// cross derivative			  
        }
	}
}


// calculate interpolation coefficients Brc, Bzc to speed up interpolation
void TabField::PreInterpol(long double ****&Brc, long double ****&Bzc, long double **BrTab, long double **BzTab){
	long double **BrTab1 = NULL, **BzTab1 = NULL;	// dBi/dr
	long double **BrTab2 = NULL, **BzTab2 = NULL;	// dBi/dz
	long double **BrTab12 = NULL, **BzTab12 = NULL;	// d²Bi/drdz
	Log("allocating memory for derivatives (%.4LG MB)\n",6*(long double)(m-2)*(n-2)*sizeof(long double)/1024/1024);	
	BrTab1=dmatrix(2,m-1,2,n-1);
	BzTab1=dmatrix(2,m-1,2,n-1);
	BrTab2=dmatrix(2,m-1,2,n-1);
	BzTab2=dmatrix(2,m-1,2,n-1);
	BrTab12=dmatrix(2,m-1,2,n-1);
	BzTab12=dmatrix(2,m-1,2,n-1);
	CalcDeriv4th(BrTab,BrTab1,BrTab2,BrTab12);
	CalcDeriv4th(BzTab,BzTab1,BzTab2,BzTab12);
	Log("Derivatives calculated.\n");
	
	Log("allocating memory for interpolation coefficients (%.4LG MB)\n",2*(long double)m*n*4*4*sizeof(long double)/1024/1024);	
	// allocating space for the preinterpolation
	// The B*c are 4D arrays with m x n x 4 x 4 fields
	Brc = viertensor(1,m,1,n,1,4,1,4);
	Bzc = viertensor(1,m,1,n,1,4,1,4);

	int indr, indz, perc=0; 
	long double *yyy = NULL, *yyy1 = NULL, *yyy2 = NULL, *yyy12 = NULL;        //rectangle with values at the 4 corners for interpolation
	yyy=dvector(1,4);
	yyy1=dvector(1,4);
	yyy2=dvector(1,4);
	yyy12=dvector(1,4);
	
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
			bcucof(yyy,yyy1,yyy2,yyy12,rdist,zdist,Brc[indr][indz]);
									
			
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
			
						
			bcucof(yyy,yyy1,yyy2,yyy12,rdist,zdist,Bzc[indr][indz]);
		}
	}
	printf("100%%\n");
	free_dvector(yyy,1,4);
	free_dvector(yyy1,1,4);
	free_dvector(yyy2,1,4);
	free_dvector(yyy12,1,4);
			
	free_dmatrix(BrTab1,2,m-1,2,n-1);
	free_dmatrix(BzTab1,2,m-1,2,n-1);		
	free_dmatrix(BrTab2,2,m-1,2,n-1);
	free_dmatrix(BzTab2,2,m-1,2,n-1);		
	free_dmatrix(BrTab12,2,m-1,2,n-1);
	free_dmatrix(BzTab12,2,m-1,2,n-1);
}


// interpolate B-field
bool TabField::BInterpol(long double r, long double z){
	if (Brc && Bzc){
		int indr = 1 + (int)((r - conv_rA)/conv_rB);    // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!		
		int indz = 1 + (int)((z - conv_zA)/conv_zB);  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3)  ){
			// bicubic interpolation
			long double Brtemp, dBrdrtemp, dBrdztemp, Bztemp, dBzdrtemp, dBzdztemp;
			bcuint_new(indr, indz, Brc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r, z, &Brtemp, &dBrdrtemp, &dBrdztemp);
			bcuint_new(indr, indz, Bzc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r, z, &Bztemp, &dBzdrtemp, &dBzdztemp);
			Br += BFeldSkal*Brtemp; Bz += BFeldSkal*Bztemp;
			dBrdr += BFeldSkal*dBrdrtemp; dBrdz += BFeldSkal*dBrdztemp;
			dBzdr += BFeldSkal*dBzdrtemp; dBzdz += BFeldSkal*dBzdztemp;
			return true;
		}
	}
	return false;
}


// interpolate E-field
bool TabField::EInterpol(long double r, long double z){
	if (Erc && Ezc){
		int indr = 1 + (int)((r - conv_rA)/conv_rB);    // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!		
		int indz = 1 + (int)((z - conv_zA)/conv_zB);  // 1 + ..., weil die werte nicht von 0, sondern von 1 beginnen!!!!!!
		if ((indr <= m-2) && (indr >= 3) && (indz <= n-2) && (indz >= 3)  ){			
			// bicubic interpolation
			long double Ertemp, dErdrtemp, dErdztemp, Eztemp, dEzdrtemp, dEzdztemp;
			bcuint_new(indr, indz, Erc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r, z, &Ertemp, &dErdrtemp, &dErdztemp);
			bcuint_new(indr, indz, Ezc, rdist, zdist, rind[indr], rind[indr +1], zind[indz], zind[indz +1], r, z, &Eztemp, &dEzdrtemp, &dEzdztemp);
			Er += EFeldSkal*Ertemp; Ez += EFeldSkal*Eztemp;
			dErdr += EFeldSkal*dErdrtemp; dErdz += EFeldSkal*dErdztemp;
			dEzdr += EFeldSkal*dEzdrtemp; dEzdz += EFeldSkal*dEzdztemp;
			return true;
		}
	}
	return false;
}


// destructor -> free memory
TabField::~TabField(){
	if (rind) free_dvector(rind,1,m);
	if (zind) free_dvector(zind,1,n);
	if (Brc) free_viertensor(Brc,1,m,1,n,1,4,1,4);
	if (Bzc) free_viertensor(Bzc,1,m,1,n,1,4,1,4);
	if (Erc) free_viertensor(Erc,1,m,1,n,1,4,1,4);
	if (Ezc) free_viertensor(Ezc,1,m,1,n,1,4,1,4);
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

