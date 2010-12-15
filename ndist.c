#include <cstdio>
#include <cmath>

#include "nr/nr3.h"

#include "ndist.h"

int neutdist = 0;
MatDoub ndist;	// matrix for probability of finding particle
int n_r = 310, n_z = 1000;	// number of bins
long double size = 0.002, rmin = 0, zmin = -0.8;	// size of bins, position of first bin

// Diese Funktionen werten die in yp[][] gespeichterten Werte aus und berechnen daraus die
// Neutronenverteilung in der r-z-Ebene. Dies wird �ber die Auswertung der Aufenthalts-
// wahrscheinlichkeit des Neutrons in Quadraten erreicht.
// Wie bei der Interpolation gibt es einen Vektor f�r die r-Werte der Ecken der Quadrate,
// eine f�r die z-Werte und eine Matrix in der die Verteilung abgespeichert wird

// Vorbereitung der r und z-Vektoren
void prepndist()
{
	ndist.assign(n_r,n_z,0);
}

// Bef�llung der Verteilungsmatrix
void fillndist(long double x1, VecDoub_I &y1, long double x2, VecDoub_I &y2)
{	
	long double r1 = sqrt(y1[0]*y1[0] + y1[1]*y1[1]);
	long double z1 = y1[2];
	long double r2 = sqrt(y2[0]*y2[0] + y2[1]*y2[1]);
	long double z2 = y2[2];
	int ir1 = int((r1 - rmin)/size);
	int iz1 = int((y1[2] - zmin)/size);
	int ir2 = int((r2 - rmin)/size);
	int iz2 = int((y2[2] - zmin)/size);
	while (ir1 != ir2 || iz1 != iz2){
		if((ir1>=n_r)||(iz1>=n_z)||(ir1>=n_r)||(iz2>=n_z)) printf("Ndist Error %d %d %d %d \n",ir1,iz1,ir2,iz2);
		long double s_left = (ir1*size + rmin - r1)/(r2 - r1);
		long double s_right = ((ir1+1)*size + rmin - r1)/(r2 - r1);
		long double s_bottom = (iz1*size + zmin - z1)/(z2 - z1);
		long double s_top = ((iz1+1)*size + zmin - z1)/(z2 - z1);
		long double s;
		int newir = ir1, newiz = iz1;
		if (s_left > 0){
			if (s_bottom > 0){
				if (s_left > s_bottom) newiz = iz1-1;
				else newir = ir1-1;
				s = min(s_left, s_bottom);
			}
			else{
				if (s_left > s_top) newiz = iz1+1;
				else newir = ir1-1;
				s = min(s_left, s_top);
			}
		}
		else{
			if (s_bottom > 0){
				if (s_right > s_bottom) newiz = iz1-1;
				else newir = ir1+1;
				s = min(s_right, s_bottom);
			}
			else{
				if (s_right > s_top) newiz = iz1+1;
				else newir = ir1+1;
				s = min(s_right, s_top);
			}
		}
		ndist[ir1][iz1] += s*(x2-x1);
		ir1 = newir;
		iz1 = newiz;
		r1 += (r2 - r1)*s;
		z1 += (z2 - z1)*s;
		x1 += (x2 - x1)*s;
	}
	ndist[ir2][iz2] += x2-x1;
}

// Ausgabe in ndist.out
void outndist(const char *ndistfile)
{
	//printf("\nOutputting the particle spacial distribution... \n");
		
	FILE *NDIST=fopen(ndistfile,"w");
	int Treffer = 0;
	
	//FILE *NDIST=fopen("ndist.out",mode_w);
	
	fprintf(NDIST,"Rindex Zindex Rmtlpkt Zmtlpkt Whk Treffer\n");
	for(int i=0; i<n_r; i++){
		for(int j=0; j<n_z; j++){
			if (ndist[i][j] != 0) Treffer =1;
			if (ndist[i][j] == 0) Treffer =0;	
			fprintf(NDIST,"%i %i %.5LG %.5LG %.17LG %i\n",i,j, i*size + rmin + size/2, j*size + zmin + size/2, ndist[i][j], Treffer);
		}
	}
	fclose(NDIST);
	ndist.resize(0,0);
}

