/**
 * \file
 * Manage all the fields.
 */

#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <cmath>

#include "fields.h"
#include "field.h"
#include "field_2d.h"
#include "field_3d.h"
#include "conductor.h"
#include "edmfields.h"

using namespace std;

TFieldManager::TFieldManager(TConfig &conf, int aFieldOscillation, double aOscillationFraction, double aOscillationFrequency):
			FieldOscillation(aFieldOscillation), OscillationFraction(aOscillationFraction), OscillationFrequency(aOscillationFrequency){
	for (map<string, string>::iterator i = conf["FIELDS"].begin(); i != conf["FIELDS"].end(); i++){
		string type;
		string ft;
		double Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime;
		double Ibar, p1, p2, p3, p4, p5, p6;
		TField *f = NULL;
		istringstream ss(i->second);

		if (i->first == "2Dtable"){
			ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime;
			if (ss)
				f = new TabField(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime);
		}
		else if (i->first == "3Dtable"){
			double BoundaryWidth;
			ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime >> BoundaryWidth;
			if (ss)
				f = new TabField3(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime, BoundaryWidth);
		}
		else if (i->first == "InfiniteWireZ"){
			ss >> Ibar >> p1 >> p2;
			if (ss)
				f = new TInfiniteWireZ(p1, p2, Ibar);
		}
		else if (i->first == "InfiniteWireZCenter"){
			ss >> Ibar;
			if (ss)
				f = new TInfiniteWireZCenter(Ibar);
		}
		else if (i->first == "FiniteWire"){
			ss >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6;
			if (ss)
				f = new TFiniteWire(p1, p2, p3, p4, p5, p6, Ibar);
		}
		else if (i->first == "FiniteWireX"){
			ss >> Ibar >> p1 >> p2 >> p3;
			if (ss)
				f = new TFiniteWireX(p1, p2, p3, Ibar);
		}
		else if (i->first == "FiniteWireY"){
			ss >> Ibar >> p1 >> p2 >> p3;
			if (ss)
				f = new TFiniteWireY(p1, p2, p3, Ibar);
		}
		else if (i->first == "FiniteWireZ"){
			ss >> Ibar >> p1 >> p2 >> p3 >> p4;
			if (ss)
				f = new TFiniteWireZ(p1, p2, p3, p4, Ibar);
		}
		else if (i->first == "FiniteWireZCenter"){
			ss >> Ibar >> p1 >> p2;
			if (ss)
				f = new TFiniteWireZCenter(p1, p2, Ibar);
		}
		else if (i->first == "FullRacetrack"){
			ss >> Ibar >> p1 >> p2 >> p3;
			if (ss)
				f = new TFullRacetrack(p1, p2, p3, Ibar);
		} 
		else if (i->first == "EDMStaticB0GradZField") {
			ss >> p1 >> p2; 
			if (ss) 
				f = new TEDMStaticB0GradZField(p1, p2);
		}
		else if (i->first == "EDMStaticEField") {
			ss >> p1;
			if (ss)
				f = new TEDMStaticEField (p1);
		}

		if (f)
			fields.push_back(f);
		else{
			cout << "\nCould not load field """ << i->first << """! Did you enter invalid parameters?\n";
			exit(-1);
		}
	}
}


TFieldManager::~TFieldManager(){
	for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++)
		delete (*i);
}


void TFieldManager::BField (double x, double y, double z, double t, double B[4][4]){
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			B[i][j] = 0;

	double BFeldSkal = BFieldScale(t);
	if (BFeldSkal != 0){
		for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++){
			(*i)->BField(x, y, z, t, B);
		}

		if (BFeldSkal != 1)
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 4; j++)
					B[i][j] *= BFeldSkal;

		B[3][0] = sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]); // absolute value of B-Vector
		if (B[3][0]>1e-31)
		{
			B[3][1] = (B[0][0]*B[0][1] + B[1][0]*B[1][1] + B[2][0]*B[2][1])/B[3][0]; // derivatives of absolute value
			B[3][2] = (B[0][0]*B[0][2] + B[1][0]*B[1][2] + B[2][0]*B[2][2])/B[3][0];
			B[3][3] = (B[0][0]*B[0][3] + B[1][0]*B[1][3] + B[2][0]*B[2][3])/B[3][0];
		}
	}
}


void TFieldManager::EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3]){
	V = 0;
	for (int i = 0; i < 3; i++){
		Ei[i] = 0;
		if (dEidxj != NULL){
			for (int j = 0; j < 3; j++)
				dEidxj[i][j] = 0;
		}
	}
	for (vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++){
		(*i)->EField(x, y, z, t, V, Ei, dEidxj);
	}
}


double TFieldManager::BFieldScale(double t){
	if (FieldOscillation==1)
		return 1 + OscillationFraction*sin(OscillationFrequency*2*pi*t);
	return 1;
}
