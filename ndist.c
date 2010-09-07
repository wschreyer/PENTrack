#include <cstdio>

#include "ndist.h"
#include "nrutil.h"
#include "globals.h"

long double *ndistr = NULL, *ndistz = NULL, **ndistW = NULL;    // matrix for probability of finding particle
int v=300,w=1200;                                       // dimension of matrix above, user choice for ndist

// Diese Funktionen werten die in yp[][] gespeichterten Werte aus und berechnen daraus die
// Neutronenverteilung in der r-z-Ebene. Dies wird �ber die Auswertung der Aufenthalts-
// wahrscheinlichkeit des Neutrons in Quadraten erreicht.
// Wie bei der Interpolation gibt es einen Vektor f�r die r-Werte der Ecken der Quadrate,
// eine f�r die z-Werte und eine Matrix in der die Verteilung abgespeichert wird

// Vorbereitung der r und z-Vektoren
void prepndist()
{
	long double dooo, increment;
	//v=300;       //  2mm steps in r-direction 0-60cm
	//w=1200;      //  z: -0.5 bis1.7
	increment = 0.6 / v;  // Schrittweite des Gitters (0.002m)
	
	ndistr=dvector(0,v);       // r-Werte
	ndistz=dvector(0,w);       // z-Werte
	ndistW=dmatrix(0,v,0,w);   // Verteilung

    // r-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
    dooo=0.0;
    for(int i=0;i<=v;i++)
    {
		ndistr[i] = dooo;
		dooo = dooo + increment;
    }

    // z-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
    dooo=-0.5;
    for(int i=0;i<=w;i++)
    {
		ndistz[i] = dooo;
		dooo = dooo + increment;
    }

    // Matrix der Verteilung auf 0 setzen
    for(int i=0;i<=v;i++)
    {
        for(int j=0;j<=w;j++)
        {
        ndistW[i][j]=0.0;
        }
    }
}

// Bef�llung der Verteilungsmatrix
void fillndist(long double x1, long double *y1, long double x2, long double *y2)
{	
	int ir1, iz1, ir2, iz2;

	hunt(ndistr, v, y1[1], &ir1);     // Index des ersten Punktes suchen
	hunt(ndistz, w, y2[3], &iz1);
	//ir1 = (yp[1][1] - conv_rA) / (conv_rB);
	//iz1 = (yp[3][1] - conv_zA) / (conv_zB);
	if((ir1>v)||(iz1>w)) printf("Ndist Error %d %d \n",ir1,iz1);

	hunt(ndistr, v, y2[1], &ir2);   // Index des n�chsten Punktes suchen
	hunt(ndistz, w, y2[3], &iz2);
	//ir2 = (yp[1][klauf] - conv_rA) / (conv_rB);
	//iz2 = (yp[3][klauf] - conv_zA) / (conv_zB);
	if((ir1>v)||(iz2>w)) printf("Ndist Error %d %d \n",ir2,iz2);
	
	if ((ir2==ir1) && (iz2==iz1))                     // sind die Pkte im gleichen Quadrat?
	{
		ndistW[ir1][iz1]=ndistW[ir1][iz1] + x2-x1;            // Zeit zum Quadrat dazuaddieren
	}
	else
	{
		ndistW[ir1][iz1]=ndistW[ir1][iz1] + (x2-x1)/2;        // H�lfte der Zeit zum ersten Quadrat
		ndistW[ir2][iz2]=ndistW[ir2][iz2] + (x2-x1)/2;        // H�lfte der Zeit zum anderen
	}
	
}

// Ausgabe in ndist.out
void outndist(const char *ndistfile)
{
	//printf("\nOutputting the particle spacial distribution... \n");
		
	FILE *NDIST=fopen(ndistfile,"w");
	int Treffer = 0;
	
	//FILE *NDIST=fopen("ndist.out",mode_w);
	
	fprintf(NDIST,"Rindex Zindex Rmtlpkt Zmtlpkt Whk Treffer\n");
	long double increment = ndistr[1] - ndistr[0];
	int i,j;
	for(i=0;i<=v;i++){
		for(j=0;j<=w;j++){
			if (ndistW[i][j] != 0) Treffer =1;
			if (ndistW[i][j] == 0) Treffer =0;	
			fprintf(NDIST,"%i %i %.5LG %.5LG %.17LG %i\n",i,j,(ndistr[i]+(increment/2.0)),(ndistz[j]+(increment/2.0)), ndistW[i][j], Treffer);
		}
	}
	fclose(NDIST);
	//printf("\nDone outputting the particle spacial distribution! \n");
	
}

