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
#include "harmonicfields.h"
#include "analyticFields.h"


TFieldManager::TFieldManager(TConfig &conf){
	for (const auto &i: conf["FIELDS"]){
		std::string type;
		boost::filesystem::path ft;
		double Ibar, p1, p2, p3, p4, p5, p6, p7;
		double bW, xma, xmi, yma, ymi, zma, zmi;
		double axis_x, axis_y, axis_z, angle, G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14, G15, G16, G17, G18, G19, G20, G21, G22, G23;
		std::string Bscale, Escale;
        std::unique_ptr<TField> f;
		std::istringstream ss(i.second);
		ss >> type;

        if (type == "OPERA2D" or type == "2Dtable"){
            fields.push_back(ReadOperaField2(i.second));
		}
        else if (type == "OPERA3D" or type == "3Dtable"){
            fields.push_back(ReadOperaField3(i.second));
		}
        else if (type == "COMSOL"){
            fields.push_back(ReadComsolField(i.second));
		}
        else if ((type == "Conductor") && (ss >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> Bscale)){
            fields.push_back(std::unique_ptr<TField>(new TConductorField(p1, p2, p3, p4, p5, p6, Ibar, Bscale)));
		}
        else if ((type == "EDMStaticB0GradZField") && (ss >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> bW >> xma >> xmi >> yma >> ymi >> zma >> zmi >> Bscale)){
			//conversion to radians
			p4*=pi/180;
			p5*=pi/180;
            fields.push_back(std::unique_ptr<TField>(new TEDMStaticB0GradZField(p1, p2, p3, p4, p5, p6, p7, bW, xma, xmi, yma, ymi, zma, zmi, Bscale)));
		}
		else if (type == "HarmonicExpandedBField" and
			ss >> p1 >> p2 >> p3 >> bW >> xma >> xmi >> yma >> ymi >> zma >> zmi >> Bscale >> axis_x >> axis_y >> axis_z >> angle and
			ss >> G0 >> G1 >> G2 >> G3 >> G4 >> G5 >> G6 >> G7 >> G8 >> G9 >> G10 >> G11 >> G12 >> G13 >> G14 >> G15 >> G16 >> G17 >> G18 >> G19 >> G20 >> G21 >> G22 >> G23){
			
			//conversion to radians
			p4*=pi/180;
			p5*=pi/180;
			
			fields.push_back(std::unique_ptr<TField>(new HarmonicExpandedBField(p1, p2, p3, bW, xma, xmi, yma, ymi, zma, zmi, Bscale, axis_x, axis_y, axis_z, angle,
			        G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14, G15, G16, G17, G18, G19, G20, G21, G22, G23)));
		}
        else if ((type == "ExponentialFieldX") && (ss >> p1 >> p2 >> p3 >> p4 >> p5 >> xma >> xmi >> yma >> ymi >> zma >> zmi)){
            fields.push_back(std::unique_ptr<TField>(new TExponentialFieldX(p1, p2, p3, p4, p5, xma, xmi, yma, ymi, zma, zmi)));
		}
        else if ((type == "LinearFieldZ") && (ss >> p1 >> p2 >> xma >> xmi >> yma >> ymi >> zma >> zmi)){
            fields.push_back(std::unique_ptr<TField>(new TLinearFieldZ(p1, p2, xma, xmi, yma, ymi, zma, zmi)));
		}
        else if ((type == "EDMStaticEField") && (ss >> p1 >> p2 >> p3 >> Bscale)){
            fields.push_back(std::unique_ptr<TField>(new TEDMStaticEField (p1, p2, p3, Bscale)));
		}
		else{
            throw std::runtime_error("Could not load field """ + type + """! Check config file for invalid field type or parameters.");
		}
	}
	std::cout << "\n";
}


void TFieldManager::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
	for (int i = 0; i < 3; i++){
		B[i] = 0;
		if (dBidxj != nullptr){
			for (int j = 0; j < 3; j++)
				dBidxj[i][j] = 0;
		}
	}

    for (const auto &it: fields){
		double Btmp[3] = {0,0,0};
		double dBtmp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
		if (dBidxj != nullptr)
            it->BField(x, y, z, t, Btmp, dBtmp);
		else
            it->BField(x, y, z, t, Btmp, nullptr);

		for (int i = 0; i < 3; i++){
			B[i] += Btmp[i];
			if (dBidxj != nullptr){
				for (int j = 0; j < 3; j++)
					dBidxj[i][j] += dBtmp[i][j];
			}
		}
	}
}


void TFieldManager::EField(const double x, const double y, const double z, const double t,
		double &V, double Ei[3]) const{
	V = 0;
	for (int i = 0; i < 3; i++){
		Ei[i] = 0;
		}
    for (const auto &it: fields){
		double Vtmp = 0, Etmp[3] = {0,0,0};

        it->EField(x, y, z, t, Vtmp, Etmp);

			V += Vtmp;
			for (int i = 0; i < 3; i++)
				Ei[i] += Etmp[i];
		}
}
