#include "main.h"
#include "kdtree.h"

KDTree geometry, source;

void LoadGeometry(){
	geometry.ReadFile("geometry/heliumtankMay09.STL",KENNZAHL_HIT_WALL);
	geometry.ReadFile("geometry/absorber.STL",KENNZAHL_ABSORBED);
	//...	
	geometry.Init();
	
	/*
	source.ReadFile("..",0);
	long double ControlPoint = {0,0,0};
	source.Init(ControlPoint);
	*/
}

bool ReflectCheck(long double x1, long double *y1, long double x2, long double *y2){
	long double p1[3] = {y1[1]*cos(y1[5]), y1[1]*sin(y1[5]), y1[3]};
	long double p2[3] = {y2[1]*cos(y2[5]), y2[1]*sin(y2[5]), y2[3]};
	long double s, *n;
	int kennzahl;
	if (geometry.Collision(p1,p2,s,n,kennzahl)){
		
		return true;
	}
	return false;
}
