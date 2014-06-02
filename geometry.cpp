#include <iostream>
#include <sstream>

#include "globals.h"
#include "geometry.h"

using namespace std;

TGeometry::TGeometry(TConfig &geometryin){
	vector<material> materials;

	for (map<string, string>::iterator i = geometryin["MATERIALS"].begin(); i != geometryin["MATERIALS"].end(); i++){
		material mat;
		mat.name = i->first;
		istringstream ss(i->second);
		ss >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb >> mat.SpinflipProb;
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
						mesh.ReadFile(STLfile.c_str(),solids.size(),name);
						model.name = name;
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

		long double ignorestart, ignoreend;
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
	mesh.Init();
	cout << '\n';
}


bool TGeometry::GetCollisions(const long double x1, const long double p1[3], const long double h, const long double p2[3], map<TCollision, bool> &colls){
	set<TCollision> c;
	mesh.Collision(p1,p2,c);
	colls.clear();
	for (set<TCollision>::iterator it = c.begin(); it != c.end(); it++){
		colls[*it] = false;
		std::vector<long double> *times = &solids[it->sldindex].ignoretimes;
		if (!times->empty()){
			long double x = x1 + h*it->s;
			for (unsigned int i = 0; i < times->size(); i += 2){
				if (x >= (*times)[i] && x < (*times)[i+1]){
					colls[*it] = true; // set ignored flag if collision should be ignored according to geometry.in
					break;
				}
			}
		}
	}
	return !colls.empty();
}


void TGeometry::GetSolids(const long double t, const long double p[3], std::map<solid, bool> &currentsolids){
	long double p2[3] = {p[0], p[1], mesh.tree.bbox().zmin() - REFLECT_TOLERANCE};
	map<TCollision, bool> c;
	currentsolids.clear();
	currentsolids[defaultsolid] = false;
	if (GetCollisions(t,p,0,p2,c)){	// check for collisions of a vertical segment from p to lower border of bounding box
		for (map<TCollision, bool>::iterator i = c.begin(); i != c.end(); i++){
			solid *sld = &solids[i->first.sldindex];
			if (currentsolids.count(*sld) > 0) // if there is a collision with a solid already in the list, remove it from list
				currentsolids.erase(*sld);
			else
				currentsolids[*sld] = i->second; // else add solid to list
		}
	}
}

