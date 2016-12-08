#include <iostream>
#include <sstream>

#include "geometry.h"
#include "globals.h"

using namespace std;

TGeometry::TGeometry(TConfig &geometryin){
	vector<material> materials;

	for (map<string, string>::iterator i = geometryin["MATERIALS"].begin(); i != geometryin["MATERIALS"].end(); i++){
		material mat;
		mat.name = i->first;
		istringstream ss(i->second);
		ss >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb >> mat.SpinflipProb >> mat.RMSRoughness >> mat.CorrelLength >> mat.UseMRModel;
		if (ss)
			materials.push_back(mat);
		else
			cout << "Could not load material " << i->first << '\n';
	}

	string line;
	string STLfile;
	string matname;
	char name[80];
	for (map<string, string>::iterator i = geometryin["GEOMETRY"].begin(); i != geometryin["GEOMETRY"].end(); i++){	// parse STLfile list
		solid model;
		istringstream(i->first) >> model.ID;
		if (model.ID < 0){
			cout << "You defined a solid with ID " << model.ID << " < 0! IDs have to be larger than zero!\n";
			exit(-1);
		}
		for (vector<solid>::iterator j = solids.begin(); j != solids.end(); j++){
			if (j->ID == model.ID){
				cout << "You defined solids with identical ID " << model.ID << "! IDs have to be unique!\n";
				exit(-1);
			}
		}

		istringstream ss(i->second);
		ss >> STLfile >> matname;
		if (ss){
			for (unsigned i = 0; i < materials.size(); i++){
				if (matname == materials[i].name){
					model.mat = materials[i];
					if (model.ID > 1){
						model.name = mesh.ReadFile(STLfile, solids.size());
					}
					else
						model.name = "default solid";
					break;
				}
				else if (i+1 == materials.size()){
					printf("Material %s used for %s but not defined in geometry.in!",matname.c_str(),name);
					exit(-1);
				}
			}
		}
		else{
			cout << "Could not load solid with ID " << model.ID << "! Did you define invalid parameters?\n";
			exit(-1);
		}

		double ignorestart, ignoreend;
		while (ss){
			ss >> ignorestart;
			if (!ss) // no more ignore times found
				break;
			char dash;
			ss >> dash;
			ss >> ignoreend;
			if (ss && dash == '-'){
				model.ignoretimes.push_back(ignorestart);
				model.ignoretimes.push_back(ignoreend);
			}
			else{
				cout << "Invalid ignoretimes in geometry.in" << '\n';
				exit(-1);
			}
		}

		if (model.ID == 1)
			defaultsolid = model;
		else
			solids.push_back(model);
	}
	cout << '\n';
}


bool TGeometry::GetCollisions(const double x1, const double p1[3], const double x2, const double p2[3], map<TCollision, bool> &colls) const{
	vector<TCollision> c = mesh.Collision(std::vector<double>{p1[0], p1[1], p1[2]}, std::vector<double>{p2[0], p2[1], p2[2]});
	colls.clear();
	for (vector<TCollision>::iterator it = c.begin(); it != c.end(); it++){
		colls[*it] = false;
		std::vector<double> times = solids.at(it->sldindex).ignoretimes;
		if (!times.empty()){
			double x = x1 + (x2 - x1)*it->s;
			for (unsigned int i = 0; i < times.size(); i += 2){
				if (x >= times[i] && x < times[i+1]){
					colls[*it] = true; // set ignored flag if collision should be ignored according to geometry.in
					break;
				}
			}
		}
	}
	return !colls.empty();
}


void TGeometry::GetSolids(const double t, const double p[3], std::map<solid, bool> &currentsolids) const{
	double p2[3] = {p[0], p[1], mesh.GetBoundingBox().zmin() - REFLECT_TOLERANCE};
	map<TCollision, bool> c;
	currentsolids.clear();
	currentsolids[defaultsolid] = false;
	if (GetCollisions(t,p,0,p2,c)){	// check for collisions of a vertical segment from p to lower border of bounding box
		for (map<TCollision, bool>::iterator i = c.begin(); i != c.end(); i++){
			solid sld = solids.at(i->first.sldindex);
			if (currentsolids.count(sld) > 0) // if there is a collision with a solid already in the list, remove it from list
				currentsolids.erase(sld);
			else
				currentsolids[sld] = i->second; // else add solid to list
		}
	}
}

solid TGeometry::GetSolid(const double t, const double p[3]) const{
	std::map<solid, bool> currentsolids;
	GetSolids(t, p, currentsolids);
	return GetSolid(t, p, currentsolids);
}

solid TGeometry::GetSolid(const double t, const double p[3], const map<solid, bool> &currentsolids) const{
	for (std::map<solid, bool>::const_iterator i = currentsolids.begin(); i != currentsolids.end(); i++)
		if (!i->second)
			return i->first;
	return defaultsolid;
}
