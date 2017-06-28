/**
 * \file
 * Manage all the fields.
 */

#include "fields.h"

#include <string>
#include <iostream>
#include <vector>

#include "field_2d.h"
#include "field_3d.h"
#include "conductor.h"
#include "edmfields.h"

TFieldManager::TFieldManager(TConfig &conf){
	for (std::map<std::string, std::string>::iterator i = conf["FIELDS"].begin(); i != conf["FIELDS"].end(); i++){
		std::string type;
		std::string ft;
		double Ibar, p1, p2, p3, p4, p5, p6, p7, freq, time1, time2, shift, bW, xma, xmi, yma, ymi, zma, zmi;
		std::string Bscale, Escale;
		TField *f = NULL;
		std::istringstream ss(i->second);
		ss >> type;

		if (type == "2Dtable"){
			ss >> ft >> Bscale >> Escale;
			if (ss)
				f = new TabField(ft.c_str(), Bscale, Escale);
		}
		else if (type == "3Dtable"){
			double BoundaryWidth;
			ss >> ft >> Bscale >> Escale >> BoundaryWidth;
			if (ss)
				f = new TabField3(ft.c_str(), Bscale, Escale, BoundaryWidth);
		}
		else if (type == "Conductor"){
			ss >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> Bscale;
			if (ss)
				f = new TConductorField(p1, p2, p3, p4, p5, p6, Ibar, Bscale);
		}
		else if (type == "EDMStaticB0GradZField") {
			ss >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> bW >> xma >> xmi >> yma >> ymi >> zma >> zmi >> Bscale;
			
			//conversion to radians
			p4*=pi/180;
			p5*=pi/180;
			
			
			if (ss)
				f = new TEDMStaticB0GradZField(p1, p2, p3, p4, p5, p6, p7, 0, 0, 0, 0, 0, bW, xma, xmi, yma, ymi, zma, zmi, Bscale);
		}
		else if (type == "EDM_AC_B1Field") {
			ss >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> freq >> time1 >> time2 >> shift >> bW >> xma >> xmi >> yma >> ymi >> zma >> zmi >> Bscale;
			
			//convert to radians
			p4*=pi/180;
			p5*=pi/180;
			
			
			if (ss)
				f = new TEDMStaticB0GradZField(p1, p2, p3, p4, p5, p6, p7, 1, freq, time1, time2, shift, bW, xma, xmi, yma, ymi, zma, zmi, Bscale);
		}
		else if (type == "EDMStaticEField") {
			ss >> p1 >> p2 >> p3 >> Bscale;
			if (ss)
				f = new TEDMStaticEField (p1, p2, p3, Bscale);
		}

		if (f)
			fields.push_back(f);
		else{
			std::cout << "\nCould not load field """ << type << """! Did you enter invalid parameters?\n";
			exit(-1);
		}
	}
	std::cout << "\n";
}


TFieldManager::~TFieldManager(){
	for (std::vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++)
		delete (*i);
}


void TFieldManager::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	for (int i = 0; i < 3; i++){
		B[i] = 0;
		if (dBidxj != NULL){
			for (int j = 0; j < 3; j++)
				dBidxj[i][j] = 0;
		}
	}

	for (std::vector<TField*>::const_iterator it = fields.begin(); it != fields.end(); ++it){
		double Btmp[3] = {0,0,0};
		double dBtmp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
		if (dBidxj != NULL)
			(*it)->BField(x, y, z, t, Btmp, dBtmp);
		else
			(*it)->BField(x, y, z, t, Btmp);

		for (int i = 0; i < 3; i++){
			B[i] += Btmp[i];
			if (dBidxj != NULL){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] += dBtmp[i][j];
			}
		}
	}
}


void TFieldManager::EField(const double x, const double y, const double z, const double t,
		double &V, double Ei[3], double dEidxj[3][3]) const{
	V = 0;
	for (int i = 0; i < 3; i++){
		Ei[i] = 0;
		if (dEidxj != NULL){
			for (int j = 0; j < 3; j++)
				dEidxj[i][j] = 0;
		}
	}
	for (std::vector<TField*>::const_iterator it = fields.begin(); it != fields.end(); ++it){
		double Vtmp = 0, Etmp[3] = {0,0,0};

		if (dEidxj != NULL){
			double dEtmp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

			(*it)->EField(x, y, z, t, Vtmp, Etmp, dEtmp);

			V += Vtmp;
			for (int i = 0; i < 3; i++){
				Ei[i] += Etmp[i];
				for (int j = 0; j < 3; j++)
					dEidxj[i][j] += dEtmp[i][j];
			}
		}
		else{
			(*it)->EField(x, y, z, t, Vtmp, Etmp);

			V += Vtmp;
			for (int i = 0; i < 3; i++)
				Ei[i] += Etmp[i];
		}
	}
}

