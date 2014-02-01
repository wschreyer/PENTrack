/**
 * \file
 * Contains classes to include experiment geometry and particle sources.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <limits>

#include "mc.h"
#include "trianglemesh.h"
#include "fields.h"

static const int MAX_DICE_ROLL = 42000000; ///< number of tries to find start point
static const double REFLECT_TOLERANCE = 1e-8;  ///< max distance of reflection point to actual surface collision point


/// Struct to store material properties (read from geometry.in, right now only for neutrons)
struct material{
	string name; ///< Material name
	long double FermiReal; ///< Real part of Fermi potential
	long double FermiImag; ///< Imaginary part of Fermi potential
	long double DiffProb; ///< Diffuse reflection probability
	long double SpinflipProb; ///< Probability for spin flip on reflection
};


/// Struct to store solid information (read from geometry.in)
struct solid{
	string name; ///< name of solid
	material mat; ///< material of solid
	int ID; ///< ID of solid
	vector<long double> ignoretimes; ///< pairs of times, between which the solid should be ignored

	/**
	 * Comparison operator used to sort solids by priority
	 */
	bool operator< (const solid s) const { return ID < s.ID; };
};


/**
 * Class to include experiment geometry.
 *
 * Loads solids and materials from geometry.in, maintains solids list, checks for collisions.
 */
struct TGeometry{
	public:
		TTriangleMesh mesh; ///< kd-tree structure containing triangle meshes from STL-files
		vector<solid> solids; ///< solids list
		solid defaultsolid; ///< "vacuum", this solid's properties are used when the particle is not inside any other solid
		
		/**
		 * Constructor, reads geometry configuration file, loads triangle meshes.
		 *
		 * @param geometryin Path to geometry configuration file
		 */
		TGeometry(const char *geometryin){
			ifstream infile(geometryin);
			string line;
			vector<material> materials;	// dynamic array to store material properties
			while (infile.good()){
				infile >> ws; // ignore whitespaces
				if ((infile.peek() == '[') && getline(infile,line).good()){	// parse geometry.in for section header
					if (line.compare(0,11,"[MATERIALS]") == 0){
						LoadMaterialsSection(infile, materials);
					}
					else if (line.compare(0,10,"[GEOMETRY]") == 0){
						LoadGeometrySection(infile, materials);
					}
					else getline(infile,line);
				}
				else getline(infile,line);
			}
		};

		/**
		 * Check if point is inside geometry bounding box.
		 *
		 * @param y1 Position vector of segment start
		 * @param y2 Position vector of segment end
		 *
		 * @return Returns true if point is inside the bounding box
		 */
		bool CheckSegment(const long double y1[3], const long double y2[3]){
			return CGAL::do_intersect(mesh.tree.bbox(), CSegment(CPoint(y1[0], y1[1], y1[2]), CPoint(y2[0], y2[1], y2[2])));
		};
		
		/**
		 * Checks if line segment p1->p2 collides with a surface.
		 *
		 * Calls KDTree::Collision to check for collisions and removes all collisions
		 * which should be ignored (given by ignore times in geometry configuration file).
		 *
		 * @param x1 Start time of line segment
		 * @param p1 Start point of line segment
		 * @param h Time length of line segment
		 * @param p2 End point of line segment
		 * @param colls List of collisions
		 *
		 * @return Returns true if line segment collides with a surface
		 */
		bool GetCollisions(const long double x1, const long double p1[3], const long double h, const long double p2[3], set<TCollision> &colls){
			if (mesh.Collision(p1,p2,colls)){ // search in the kdtree for collisions
				set<TCollision>::iterator it = colls.begin();
				while (it != colls.end()){ // go through all collisions
					vector<long double> *times = &solids[(*it).ID].ignoretimes;
					if (!times->empty()){
						long double x = x1 + (*it).s*h;
						for (unsigned int i = 0; i < times->size(); i += 2){
							if (x >= (*times)[i] && x < (*times)[i+1]){
								set<TCollision>::iterator del = it;
								it++;
								colls.erase(del); // delete collision if it should be ignored according to geometry.in
								break;
							}
							else if (i == times->size() - 2)
								it++;
						}
					}
					else it++;
				}
				if (!colls.empty())
					return true; // return true if there is a collision
			}
			return false;		
		};
			
		/**
		 * Get solids in which the point p lies
		 *
		 * @param t Time
		 * @param p Point to test
		 * @param currentsolids Set of solids in which the point is inside
		 */
		void GetSolids(const long double t, const long double p[3], std::set<solid> &currentsolids){
			long double p2[3] = {p[0], p[1], mesh.tree.bbox().zmin() - REFLECT_TOLERANCE};
			set<TCollision> c;
			currentsolids.clear();
			currentsolids.insert(defaultsolid);
			if (GetCollisions(t,p,0,p2,c)){	// check for collisions of a vertical segment from p to lower border of bounding box
				for (set<TCollision>::iterator i = c.begin(); i != c.end(); i++){
					solid sld = solids[i->ID];
					if (currentsolids.count(sld) > 0) // if there is a collision with a solid already in the list, remove it from list
						currentsolids.erase(sld);
					else
						currentsolids.insert(sld); // else add solid to list
				}
			}
		}

	private:	
		/**
		 * Read [MATERIALS] section in geometry configuration file.
		 *
		 * @param infile File stream to geometry configuration file
		 * @param materials Vector of material structs, returns list of materials in the configuration file
		 */
		void LoadMaterialsSection(ifstream &infile, vector<material> &materials){
			int c;
			string line;
			do{	// parse material list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue; // skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				material mat;
				infile >> mat.name >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb >> mat.SpinflipProb;
				materials.push_back(mat);
			}while(infile.good() && getline(infile,line).good());	
		};
	
		/**
		 * Read [GEOMETRY] section in geometry configuration file.
		 *
		 * @param infile File stream to geometry configuration file
		 * @param materials Vector of material structs given by LoadMaterialsSection
		 */
		void LoadGeometrySection(ifstream &infile, vector<material> &materials){
			int c;
			string line;
			string STLfile;
			string matname;
			char name[80];
			long double ignorestart, ignoreend;
			do{	// parse STLfile list
				infile >> ws; // ignore whitespaces
				c = infile.peek();
				if (c == '#') continue;	// skip comments
				else if (!infile.good() || c == '[') break;	// next section found
				solid model;
				infile >> model.ID;
				if (model.ID < 0){
					cout << "You defined a solid with ID " << model.ID << " < 0! IDs have to be larger than zero!\n";
					exit(-1);
				}
				for (vector<solid>::iterator i = solids.begin(); i != solids.end(); i++){
					if (i->ID == model.ID){
						cout << "You defined solids with identical ID " << model.ID << "! IDs have to be unique!\n";
						exit(-1);
					}
				}
				infile >> STLfile;
				infile >> matname;
				for (unsigned i = 0; i < materials.size(); i++){
					if (matname == materials[i].name){
						model.mat = materials[i];
						if (model.ID > 1)
							mesh.ReadFile(STLfile.c_str(),solids.size(),name);
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
						cout << "Invalid ignoretimes in geometry.in" << '\n';
						exit(-1);
					}
					model.ignoretimes.push_back(ignorestart);
					model.ignoretimes.push_back(ignoreend);
					while ((c = infile.peek()) == '\t' || c == ' ')
						infile.ignore();
				}
				if (model.ID == 1)
					defaultsolid = model;
				else
					solids.push_back(model);
			}while(infile.good() && getline(infile,line).good());
			mesh.Init();
		};
		
};

#endif /*GEOMETRY_H_*/
