#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <fstream>
#include <list>
#include <vector>

#include "mc.h"
#include "kdtree.h"
#include "fields.h"

using namespace std;

#define MAX_DICE_ROLL 42000000 // number of tries to find start point
#define REFLECT_TOLERANCE 1e-8 	// max distance of reflection point to actual surface collision point

struct material{
	string name;
	long double FermiReal, FermiImag, DiffProb;
};

struct solid{
	string name;
	material mat;
	unsigned kennz;
	set<int> ignoretimes;
};

extern solid defaultsolid;

// transmission probability of neutron hitting surface mit Fermi Potential Mf + i*Pf (all in neV)
long double Transmission(const long double Er, const long double Mf, const long double Pf);
// absorption probability of neutron flying distance l (in m) through Fermi potential Mf + i*Pf (all in neV)
long double Absorption(const long double E, const long double Mf, const long double Pf, const long double l);


struct TGeometry{
	public:
		int diffuse; // diffuse reflection switch, reflection switch
		int reflekt[7];	// reflection on/off in different experiment stages
		
		KDTree *geometry;
		vector<solid> solids;	// dynamic array to associate solids with material/kennzahl/name
		
		TGeometry(const char *geometryin, int areflekt[7], int adiffuse, vector<int> kennz_counter[3]){
			diffuse = adiffuse;
			for (int i = 0; i < 7; i++) reflekt[i] = areflekt[i];
			ifstream infile(geometryin);
			string line;
			vector<material> materials;	// dynamic array to store material properties
			char c;
			printf("\n");
			while (infile.good()){
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if ((infile.peek() == '[') && getline(infile,line).good()){	// parse geometry.in for section header
					if (line.compare(0,11,"[MATERIALS]") == 0){
						LoadMaterialsSection(infile, materials);
					}
					else if (line.compare(0,10,"[GEOMETRY]") == 0){
						LoadGeometrySection(infile, materials, kennz_counter);
					}
					else getline(infile,line);
				}
				else getline(infile,line);
			}
		};

		~TGeometry(){
			if (geometry) delete geometry;
		};
		
		// check if point is inside geometry volume
		bool CheckPoint(const long double x, const long double y[6]){
			if (!geometry->PointInBox(y)){ // particle no longer contained inside geometry-volume?
				printf("\nParticle has hit outer boundaries: Stopping it! t=%LG x=%LG y=%LG z=%LG\n",x,y[0],y[1],y[2]);
				return false;
			}	
			return true;
		};
		
		// get all collisions of ray y1->y2 with surfaces
		bool GetCollisions(const long double x1, const long double y1[6], const long double h, const long double y2[6], list<TCollision> &colls){
			if (geometry->Collision(y1,y2,colls)){ // search in the kdtree for collisions
				for (list<TCollision>::iterator it = colls.begin(); it != colls.end(); it++){ // go through all collisions
					if (!solids[(*it).ID].ignoretimes.empty()){
						long double x = x1 + (*it).s*h;
						if (solids[(*it).ID].ignoretimes.count(ExpPhase(x)) > 0)
							it = --colls.erase(it); // delete collision if it should be ignored according to geometry.in
					}
				}
				if (!colls.empty())
					return true; // return true if there is a collision
			}
			return false;		
		};
			
	private:	
		void LoadMaterialsSection(ifstream &infile, vector<material> &materials){
			char c;
			string line;
			do{	// parse material list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue; // skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				material mat;
				infile >> mat.name >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb;
				materials.push_back(mat);
			}while(infile.good() && getline(infile,line).good());	
		};
	
		void LoadGeometrySection(ifstream &infile, vector<material> &materials, vector<int> kennz_counter[3]){
			char c;
			string line;
			string STLfile;
			string matname;
			geometry = new KDTree;
			char name[80];
			int ignoretime;
			do{	// parse STLfile list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue;	// skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				solid model;
				infile >> model.kennz;
				if (kennz_counter[0].size() <= model.kennz){
					kennz_counter[0].resize(model.kennz+1,0);
					kennz_counter[1].resize(model.kennz+1,0);
					kennz_counter[2].resize(model.kennz+1,0);
				}
				infile >> STLfile;
				infile >> matname;
				for (unsigned i = 0; i < materials.size(); i++){
					if (matname == materials[i].name){
						model.mat = materials[i];
						geometry->ReadFile(STLfile.c_str(),solids.size(),name);
						model.name = name;
						break;
					}
					else if (i+1 == materials.size()){
						printf("Material %s used for %s but not defined in geometry.in!",matname.c_str(),name);
						exit(-1);
					}
				}
				while ((c = infile.peek()) == '\t' || c == ' ')
					infile.ignore();
				while (infile && c != '#' && c != '\n'){
					infile >> ignoretime;
					model.ignoretimes.insert(ignoretime);
					while ((c = infile.peek()) == '\t' || c == ' ')
						infile.ignore();
				}
				solids.push_back(model);
			}while(infile.good() && getline(infile,line).good());
			geometry->Init();
		};
		
};

struct TSource{
	public:
		KDTree *sourcevolume;
		TGeometry *geometry;
		string sourcemode; // volume/surface/customvol/customsurf
		long double r_min, r_max, phi_min, phi_max, z_min, z_max; // for customvol/customsurf
		long double E_normal, Emin_n;
		
		list<Triangle*> sourcetris;
		long double sourcearea;
		
		TSource(const char *geometryin, TGeometry &geom, TField &field){
			geometry = &geom;
			sourcevolume = NULL;
			E_normal = 0;
			ifstream infile(geometryin);
			string line;
			char c;
			while (infile.good()){
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if ((infile.peek() == '[') && getline(infile,line).good()){	// parse geometry.in for section header
					if (line.compare(0,8,"[SOURCE]") == 0){
						LoadSourceSection(infile);
					}
					else getline(infile,line);
				}
				else getline(infile,line);
			}
			
			if (sourcemode == "customvol"){ // get minimal possible energy for neutrons to avoid long inital value dicing for too-low-energy neutron
				Emin_n = 9e90;
				long double rBabsmin = NAN, zBabsmin = NAN, EnTemp, B[4][4];
				for (long double r = r_min; r <= r_max; r += 0.003){
					for (long double z = z_min; z <= z_max; z += 0.003){
						field.BFeld(r,0,z,0,B);
						EnTemp =  m_n*gravconst*z + (mu_nSI/ele_e)*B[3][0];
						if(EnTemp < Emin_n)
						{
							Emin_n = EnTemp;
							rBabsmin = r;
							zBabsmin = z;
						}
					}
				}
				printf("The minimum energy a low field seeking neutron has to have is %.3LG neV (at r:%.3LG, z: %.3LG).\n",Emin_n*1e9,rBabsmin,zBabsmin);
			}
			else Emin_n = -9e90;
		};
	
		~TSource(){
			if (sourcevolume) delete sourcevolume;
		};
	
		void RandomPointInSourceVolume(TMCGenerator &mc, long double &x, long double &y, long double &z, long double &alpha, long double &gamma, long double n[3]){
			long double p1[3];
			int count = 0;
			do{
				if (sourcemode == "volume"){ // random point inside a STL-volume
					p1[0] = mc.UniformDist(sourcevolume->lo[0],sourcevolume->hi[0]); // random point 
					p1[1] = mc.UniformDist(sourcevolume->lo[1],sourcevolume->hi[1]); // random point 
					p1[2] = mc.UniformDist(sourcevolume->lo[2],sourcevolume->hi[2]); // random point 
					if (!InSourceVolume(p1))
						continue;
					mc.IsotropicDist(0, alpha, gamma);
				}
				else if (sourcemode == "surface" || sourcemode == "customsurf"){ // random point on a surface inside custom or STL volume
					long double RandA = mc.UniformDist(0,sourcearea);
					long double CurrA = 0, SumA = 0;
					list<Triangle*>::iterator i;
					for (i = sourcetris.begin(); i != sourcetris.end(); i++){
						(*i)->CalcNormal(n);
						CurrA = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
						SumA += CurrA;
						if (RandA <= SumA) break; // select random triangle, weighted by triangle area
					}
					long double a = mc.UniformDist(0,1); // generate random point on triangle (see Numerical Recipes 3rd ed., p. 1114)
					long double b = mc.UniformDist(0,1);
					if (a+b > 1){
						a = 1 - a;
						b = 1 - b;
					}
					for (int j = 0; j < 3; j++){
						n[j] /= CurrA;
						p1[j] = (*i)->vertex[0][j] + a*((*i)->vertex[1][j] - (*i)->vertex[0][j]) + b*((*i)->vertex[2][j] - (*i)->vertex[0][j]);
						p1[j] += REFLECT_TOLERANCE*n[j];
					}
					if (!InSourceVolume(p1))
						continue;				
					long double winkeben = mc.UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
					long double winksenkr = mc.SinCosDist(0, 0.5*pi); // Lambert's law!
					long double v[3] = {cos(winkeben)*cos(winksenkr), sin(winkeben)*cos(winksenkr), sin(winksenkr)};
					RotateVector(v,n);
					
					alpha = atan2(v[1],v[0]) - atan2(p1[1],p1[0]);
					gamma = acos(v[2]);
				}
				else if (sourcemode == "customvol"){ // random point inside a volume defined by a custom parameter range
					long double r = mc.SquareDist(r_min, r_max); // weighting because of the volume element and a r^2 probability outwards
					long double phi = mc.UniformDist(phi_min,phi_max);
					z = mc.UniformDist(z_min,z_max);
					p1[0] = r*cos(phi);
					p1[1] = r*sin(phi);
					p1[2] = z;
					mc.IsotropicDist(0, alpha, gamma);
				}
				else{
					printf("Invalid sourcemode: %s",sourcemode.c_str());
					exit(-1);
				}
				
				long double p2[3] = {p1[0], p1[1], geometry->geometry->lo[2]};
				list<TCollision> c;
				if (geometry->geometry->Collision(p1,p2,c)){	// test if random point is inside solid
					for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
						if (i->normal[2] > 0) count++; // count surfaces whose normals point to random point
						else count--; // count surfaces whose normals point away from random point
					}
				}
			}while (count != 0);
			x = p1[0];
			y = p1[1];
			z = p1[2];
		};

	private:
		void LoadSourceSection(ifstream &infile){
			char c;
			string line;
			do{	// parse source line
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue; // skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				infile >> sourcemode;
				if (sourcemode == "customvol" || sourcemode == "customsurf"){
					infile >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max;
					phi_min *= conv;
					phi_max *= conv;
				}
				else if (sourcemode == "volume" || sourcemode == "surface"){
					string sourcefile;
					infile >> sourcefile;
					sourcevolume = new KDTree;
					sourcevolume->ReadFile(sourcefile.c_str(),0);
					sourcevolume->Init();
				}
				if (sourcemode == "customsurf" || sourcemode == "surface")
					infile >> E_normal;
			}while(infile.good() && getline(infile,line).good());
			
			if (sourcemode == "customsurf" || sourcemode == "surface"){
				sourcearea = 0; // add triangles, whose vertices are all in the source volume, to sourcetris list
				for (list<Triangle>::iterator i = geometry->geometry->alltris.begin(); i != geometry->geometry->alltris.end(); i++){
					if (InSourceVolume(i->vertex[0]) && InSourceVolume(i->vertex[1]) && InSourceVolume(i->vertex[2])){
						sourcetris.push_back(&*i);
						long double n[3];
						i->CalcNormal(n);
						sourcearea += sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
					}
				}
				printf("Source Area: %LG m²\n",sourcearea);
			}	
		};
		
		// check if a point is inside the source volume
		template<typename coord> bool InSourceVolume(const coord p[3]){
			if (sourcemode == "customvol" || sourcemode == "customsurf"){
				long double r = sqrt(p[0]*p[0] + p[1]*p[1]);
				long double phi = atan2(p[1],p[0]); 
				return (r >= r_min && r <= r_max && 
						phi >= phi_min && phi <= phi_max && 
						p[2] >= z_min && p[2] <= z_max); // check if point is in custom paramter range
			}
			else{ // shoot a ray to the bottom of the sourcevol bounding box and check for intersection with the sourcevol surface
				long double p1[3] = {p[0], p[1], p[2]};
				long double p2[3] = {p[0], p[1], sourcevolume->lo[2]}; 
				list<TCollision> c;
				if (sourcevolume->Collision(p1,p2,c)){
					int count = 0;
					for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
						if (i->normal[2] > 0) count++;  // surface normal pointing towards point?
						else count--;					// or away from point?
					}
					return (count == -1); // when there's exactly one more normal pointing away from than to it, the point is inside the sourcevol
				}
				else return false;
			}
		}
	
};

#endif /*GEOMETRY_H_*/
