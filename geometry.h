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
	vector<long double> ignoretimes;
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
		
		KDTree *kdtree;
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
			if (kdtree) delete kdtree;
		};
		
		// check if point is inside geometry volume
		bool CheckPoint(const long double x, const long double y[6]){
			if (!kdtree->PointInBox(y)){ // particle no longer contained inside geometry-volume?
				printf("\nParticle has hit outer boundaries: Stopping it! t=%LG x=%LG y=%LG z=%LG\n",x,y[0],y[1],y[2]);
				return false;
			}	
			return true;
		};
		
		// get all collisions of ray y1->y2 with surfaces
		bool GetCollisions(const long double x1, const long double p1[3], const long double h, const long double p2[3], list<TCollision> &colls){
			if (kdtree->Collision(p1,p2,colls)){ // search in the kdtree for collisions
				for (list<TCollision>::iterator it = colls.begin(); it != colls.end(); it++){ // go through all collisions
					vector<long double> *times = &solids[(*it).ID].ignoretimes;
					if (!times->empty()){
						long double x = x1 + (*it).s*h;
						for (unsigned int i = 0; i < times->size(); i += 2){
							if (x >= (*times)[i] && x < (*times)[i+1]){
								it = --colls.erase(it); // delete collision if it should be ignored according to geometry.in
								break;
							}
						}
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
			kdtree = new KDTree;
			char name[80];
			long double ignorestart, ignoreend;
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
						kdtree->ReadFile(STLfile.c_str(),solids.size(),name);
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
					infile >> ignorestart;
					if (infile.peek() == '-') infile.ignore();
					infile >> ignoreend;
					if (!infile.good()){
						cout << "Invalid ignoretimes in geometry.in" << endl;
						exit(-1);
					}
					model.ignoretimes.push_back(ignorestart);
					model.ignoretimes.push_back(ignoreend);
					while ((c = infile.peek()) == '\t' || c == ' ')
						infile.ignore();
				}
				solids.push_back(model);
			}while(infile.good() && getline(infile,line).good());
			kdtree->Init();
		};
		
};

struct TSource{
	public:
		KDTree *kdtree;
		TGeometry *geometry;
		string sourcemode; // volume/surface/customvol/customsurf
		long double r_min, r_max, phi_min, phi_max, z_min, z_max; // for customvol/customsurf
		long double E_normal, Hmin_lfs, Hmin_hfs;
		
		vector<Triangle*> sourcetris;
		long double sourcearea;
		
		TSource(const char *geometryin, TGeometry &geom, TField &field){
			geometry = &geom;
			kdtree = NULL;
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
			
			const long double d = 0.003;
			if (sourcemode == "customvol"){ // get minimal possible energy for neutrons to avoid long inital value dicing for too-low-energy neutron
				Hmin_lfs = Hmin_hfs = 9e99;
				for (long double r = r_min; r <= r_max; r += d){
					for (long double z = z_min; z <= z_max; z += d){
						long double p[3] = {r,0.0,z};
						CalcHmin(field,p,d*d*2*r*pi);
					}
				}
			}
			else if (sourcemode == "volume"){
				Hmin_lfs = Hmin_hfs = 9e99;
				long double p[3];
				for (p[0] = kdtree->lo[0]; p[0] <= kdtree->hi[0]; p[0] += d){
					for (p[1] = kdtree->lo[1]; p[1] <= kdtree->lo[1]; p[1] += d){
						for (p[2] = kdtree->lo[2]; p[2] <= kdtree->lo[2]; p[2] += d){
							if (InSourceVolume(p)) CalcHmin(field,p,d*d*d);
						}
					}
				}
			}
			else Hmin_lfs = Hmin_hfs = 0;
			printf("min. total energy a lfs/hfs neutron must have: %LG / %LG neV\n\n", Hmin_lfs, Hmin_hfs);
		};
	
		~TSource(){
			if (kdtree) delete kdtree;
		};
	
		void RandomPointInSourceVolume(TMCGenerator &mc, long double t, long double &NeutEnergie, long double &x, long double &y, long double &z, long double &alpha, long double &gamma){
			long double p1[3];
			int count = 0;
			for(;;){
				if (sourcemode == "volume"){ // random point inside a STL-volume
					p1[0] = mc.UniformDist(kdtree->lo[0],kdtree->hi[0]); // random point
					p1[1] = mc.UniformDist(kdtree->lo[1],kdtree->hi[1]); // random point
					p1[2] = mc.UniformDist(kdtree->lo[2],kdtree->hi[2]); // random point
					if (!InSourceVolume(p1)) continue;
					mc.IsotropicDist(alpha, gamma);
				}
				else if (sourcemode == "customvol"){ // random point inside a volume defined by a custom parameter range
					long double r = mc.SquareDist(r_min, r_max); // weighting because of the volume element and a r^2 probability outwards
					long double phi = mc.UniformDist(phi_min,phi_max);
					z = mc.UniformDist(z_min,z_max);
					p1[0] = r*cos(phi);
					p1[1] = r*sin(phi);
					p1[2] = z;
					mc.IsotropicDist(alpha, gamma);
				}
				else if (sourcemode == "surface" || sourcemode == "customsurf"){ // random point on a surface inside custom or STL volume
					long double *n = NULL;
					long double RandA = mc.UniformDist(0,sourcearea);
					long double CurrA = 0, SumA = 0;
					vector<Triangle*>::iterator i;
					for (i = sourcetris.begin(); i != sourcetris.end(); i++){
						n = (*i)->normal;
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
						n[j] /= CurrA; // normalize normal vector
						p1[j] = (*i)->vertex[0][j] + a*((*i)->vertex[1][j] - (*i)->vertex[0][j]) + b*((*i)->vertex[2][j] - (*i)->vertex[0][j]);
						p1[j] += REFLECT_TOLERANCE*n[j]; // add some tolerance to avoid starting inside solid
					}
					if (!InSourceVolume(p1)) continue;
					long double winkeben = mc.UniformDist(0, 2*pi); // generate random velocity angles in upper hemisphere
					long double winksenkr = mc.SinCosDist(0, 0.5*pi); // Lambert's law!
					if (E_normal > 0){
						long double vsenkr = sqrt(NeutEnergie*cos(winksenkr)*cos(winksenkr) + E_normal); // add E_normal to component normal to surface
						long double veben = sqrt(NeutEnergie)*sin(winksenkr);
						winksenkr = atan2(veben,vsenkr); // update angle
						NeutEnergie = vsenkr*vsenkr + veben*veben; // update energy
					}
					long double v[3] = {cos(winkeben)*sin(winksenkr), sin(winkeben)*sin(winksenkr), cos(winksenkr)};
					RotateVector(v,n);
					
					alpha = atan2(v[1],v[0]) - atan2(p1[1],p1[0]);
					gamma = acos(v[2]);
				}
				else{
					printf("Invalid sourcemode: %s",sourcemode.c_str());
					exit(-1);
				}
				
				long double p2[3] = {p1[0], p1[1], geometry->kdtree->lo[2] - REFLECT_TOLERANCE};
				list<TCollision> c;
				count = 0;
				if (geometry->GetCollisions(t,p1,0,p2,c)){	// test if random point is inside solid
					for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
						if (i->normal[2] > 0) count++; // count surfaces whose normals point to random point
						else count--; // count surfaces whose normals point away from random point
					}
				}
				if (count == 0) break;
			}
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
					kdtree = new KDTree;
					kdtree->ReadFile(sourcefile.c_str(),0);
					kdtree->Init();
				}
				if (sourcemode == "customsurf" || sourcemode == "surface")
					infile >> E_normal;
			}while(infile.good() && getline(infile,line).good());
			
			if (sourcemode == "customsurf" || sourcemode == "surface"){
				sourcearea = 0; // add triangles, whose vertices are all in the source volume, to sourcetris list
				for (vector<Triangle>::iterator i = geometry->kdtree->alltris.begin(); i != geometry->kdtree->alltris.end(); i++){
					if (InSourceVolume(i->vertex[0]) && InSourceVolume(i->vertex[1]) && InSourceVolume(i->vertex[2])){
						sourcetris.push_back(&*i);
						long double *n = i->normal;
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
				long double p2[3] = {p[0], p[1], kdtree->lo[2] - REFLECT_TOLERANCE};
				list<TCollision> c;
				if (kdtree->Collision(p1,p2,c)){
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
	
		void CalcHmin(TField &field, long double p[3], long double vol){
			long double p2[3] = {p[0], p[1], geometry->kdtree->lo[2] - REFLECT_TOLERANCE};
			list<TCollision> c;
			int count = 0;
			if (geometry->GetCollisions(0,p,0,p2,c)){	// test if random point is inside solid
				for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
					if (i->normal[2] > 0) count++; // count surfaces whose normals point to random point
					else count--; // count surfaces whose normals point away from random point
				}
			}
			if (count == 0){
				long double B[4][4];
				field.BFeld(p[0],p[1],p[2],0,B);
				long double Epot = m_n*gravconst*p[2];
				long double Emag = (mu_nSI/ele_e)*B[3][0];
				long double Hlfs = Epot + Emag;
				long double Hhfs = Epot - Emag;
				Hmin_lfs = min(Hmin_lfs,Hlfs);
				Hmin_hfs = min(Hmin_hfs,Hhfs);
			}
		}
};

#endif /*GEOMETRY_H_*/
