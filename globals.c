#include <cstdarg>
#include <cstdio>

#include "globals.h"

FILE *LOGSCR = NULL;

long double tau;

//set the duration of the experiment
long double FillingTime = 0;				// filling time, entrance open
long double CleaningTime = 0;				// cleaning without field
long double RampUpTime = 0;					// ramping up coils
long double FullFieldTime = 1000;			// storing in full field
long double RampDownTime = 5;				// ramping down coils
long double EmptyingTime = 0;				// emptying without field
long double StorageTime = 1500.0;			// time when ramping down shall start, if xend > storage time, let neutron decay
int ExpPhase = 0;							// current experiment phase

int jobnumber = 0;
string inpath = ".", outpath = ".";
int ausgabewunsch=5; // Ausgabewunsch

// works like printf, prints to logfile and, if jobnumber == 0, to stdout 
int Log(const char* format, ...){
	va_list args;
	va_start(args, format);
	int result;
	if (jobnumber == 0){
		result = vprintf(format, args);
		fflush(stdout);
	} 
	if (LOGSCR){
		result = vfprintf(LOGSCR, format, args);
		fflush(LOGSCR);
	}
	va_end(args);
	return result;
}

void percent(long double x, long double x1, long double len, int &perc){
	// write status to console
	// one point per 2 percent of endtime
	if((x - x1)*100/len > perc)
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



