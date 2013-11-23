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

#include "mc.h"
#include "kdtree.h"
#include "fields.h"

#define MAX_DICE_ROLL 42000000 ///< number of tries to find start point
#define REFLECT_TOLERANCE 1e-8  ///< max distance of reflection point to actual surface collision point


/// Struct to store material properties (read from geometry.in, right now only for neutrons)
struct material{
	string name; ///< Material name
	long double FermiReal; ///< Real part of Fermi potential
	long double FermiImag; ///< Imaginary part of Fermi potential
	long double DiffProb; ///< Diffuse reflection probability
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
		KDTree *kdtree; ///< kd-tree structure containing triangle meshes from STL-files
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

		~TGeometry(){
			if (kdtree) delete kdtree;
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
			return kdtree->SegmentInBox(y1,y2) ;
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
			if (kdtree->Collision(p1,p2,colls)){ // search in the kdtree for collisions
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
			long double p2[3] = {p[0], p[1], kdtree->lo[2] - REFLECT_TOLERANCE};
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
	
		/**
		 * Read [GEOMETRY] section in geometry configuration file.
		 *
		 * @param infile File stream to geometry configuration file
		 * @param materials Vector of material structs given by LoadMaterialsSection
		 */
		void LoadGeometrySection(ifstream &infile, vector<material> &materials){
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
			kdtree->Init();
		};
		
};

/**
 * Class which can produce random particle starting points from different sources.
 *
 * There are four source modes which can be specified in the [SOURCE] section of the configuration file:
 * "volume" - random points inside a STL solid (TSource::kdtree) are created;
 * "surface" - random points on triangles (part of TSource::geometry) COMPLETELY surrounded by a STL solid (TSource::kdtree) are created;
 * "customvol" - random points inside a cylindrical coordinate range ([TSource::r_min..TSource::r_max]:[TSource::phi_min..TSource::phi_max]:[TSource::z_min..TSource::z_max]) are created;
 * "customsurf"	- random points on triangles (part of TSource::geometry) COMPLETELY inside a cylindrical coordinate range  ([TSource::r_min..TSource::r_max]:[TSource::phi_min..TSource::phi_max]:[TSource::z_min..TSource::z_max]) are created
 */
struct TSource{
	public:
		KDTree *kdtree; ///< KDTree structure containing triangle meshes surrounding the source volume/surface (if sourcemode is volume or surface)
		TGeometry *geometry; ///< TGeometry structure to check the random points against, also used as source surface if sourcemode is "surface" or "customsurf"
		string sourcemode; ///< volume/surface/customvol/customsurf
		long double r_min; ///< minimum radial coordinate for customvol/customsurf source
		long double r_max; ///< maximum radial coordinate for customvol/customsurf source
		long double phi_min; ///< minimum azimuth for customvol/customsurf source
		long double phi_max; ///< maximum azimuth coordinate for customvol/customsurf source
		long double z_min; ///< minimum z coordinate for customvol/customsurf source
		long double z_max;  ///< maximum z coordinate for customvol/customsurf source
		long double E_normal; ///< For surface sources an additional energy "boost" perpendicular to the surface can be added (e.g. neutrons "dropping" from a Fermi potential >0 into vacuum)
		long double Hmin_lfs; ///< Minimum total energy of low field seeking neutrons inside the source volume
		long double Hmin_hfs; ///< Minimum total energy of high field seeking neutrons inside the source volume
		long double ActiveTime; ///< the source produces particles in the time range [0..ActiveTime]
		
		vector<Triangle*> sourcetris; ///< list of triangles (part of TSource::geometry) inside the source volume (only relevant if source mode is "surface" or "customsurf")
		long double sourcearea; ///< area of all triangles in TSource::sourcetris (needed for correct weighting)
		

		/**
		 * Constructor, loads [SOURCE] section of configuration file, calculates TSource::Hmin_lfs and TSource::Hmin_hfs
		 *
		 * @param geometryin Configuration file containing [SOURCE] section
		 * @param geom TGeometry class against which start points are checked and which contains source surfaces
		 * @param field TFieldManager class to calculate TSource::Hmin_lfs and TSource::Hmin_hfs
		 */
		TSource(const char *geometryin, TGeometry &geom, TFieldManager &field){
			geometry = &geom;
			kdtree = NULL;
			E_normal = 0;
			ifstream infile(geometryin);
			string line;
			while (infile.good()){
				infile >> ws; // ignore whitespaces
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
						CalcHmin(field,p);
					}
				}
			}
			else if (sourcemode == "volume"){
				Hmin_lfs = Hmin_hfs = 9e99;
				long double p[3];
				for (p[0] = kdtree->lo[0]; p[0] <= kdtree->hi[0]; p[0] += d){
					for (p[1] = kdtree->lo[1]; p[1] <= kdtree->lo[1]; p[1] += d){
						for (p[2] = kdtree->lo[2]; p[2] <= kdtree->lo[2]; p[2] += d){
							if (InSourceVolume(p)) CalcHmin(field,p);
						}
					}
				}
			}
			else Hmin_lfs = Hmin_hfs = 0;
			printf("min. total energy a lfs/hfs neutron must have: %LG / %LG eV\n\n", Hmin_lfs, Hmin_hfs);
		};
	

		/**
		 * Destructor, deletes TSource::kdtree
		 */
		~TSource(){
			if (kdtree) delete kdtree;
		};
	

		/**
		 * Create random point in source volume.
		 *
		 * "volume": Random points with isotropic velocity distribution inside the bounding box of TSource::kdtree are produced until TSource::InSourceVolume is true for a point.
		 * "customvol": Random points with isotropic velocity deistribution in the coordinate range ([TSource::r_min..TSource::r_max]:[TSource::phi_min..TSource::phi_max]:[TSource::z_min..TSource::z_max]) are produced.
		 * "surface/customsurf": A random point on a random triangle from the TSource::sourcetris list is chosen (weighted by triangle area).
		 * 						Starting angles are cosine-distributed to the triangle normal and then angles and Ekin are modified according to TSource::E_normal.
		 * This is repeated for all source modes until the point does not lie inside a solid in TSource::geometry.
		 *
		 * @param mc TMCGenerator to produce random number distributions
		 * @param t Start time of particle
		 * @param Ekin Kinetic start energy of particle, might be modified by E_normal
		 * @param x Returns starting x coordinate of particle
		 * @param y Returns starting y coordinate of particle
		 * @param z Returns starting z coordinate of particle
		 * @param phi Returns azimuth of velocity vector
		 * @param theta Returns polar angle of velocity vector
		 */
		void RandomPointInSourceVolume(TMCGenerator &mc, long double t, long double &Ekin, long double &x, long double &y, long double &z, long double &phi, long double &theta){
			long double p1[3];
			for(;;){
				if (sourcemode == "volume"){ // random point inside a STL-volume
					p1[0] = mc.UniformDist(kdtree->lo[0],kdtree->hi[0]); // random point
					p1[1] = mc.UniformDist(kdtree->lo[1],kdtree->hi[1]); // random point
					p1[2] = mc.UniformDist(kdtree->lo[2],kdtree->hi[2]); // random point
					if (!InSourceVolume(p1)) continue;
					mc.IsotropicDist(phi, theta);
				}
				else if (sourcemode == "customvol"){ // random point inside a volume defined by a custom parameter range
					long double r = mc.LinearDist(r_min, r_max); // weighting because of the volume element and a r^2 probability outwards
					long double phi_r = mc.UniformDist(phi_min,phi_max);
					z = mc.UniformDist(z_min,z_max);
					p1[0] = r*cos(phi_r);
					p1[1] = r*sin(phi_r);
					p1[2] = z;
					mc.IsotropicDist(phi, theta);
				}
				else if (sourcemode == "surface" || sourcemode == "customsurf"){ // random point on a surface inside custom or STL volume
					long double n[3];
					long double RandA = mc.UniformDist(0,sourcearea);
					long double CurrA = 0, SumA = 0;
					vector<Triangle*>::iterator i;
					for (i = sourcetris.begin(); i != sourcetris.end(); i++){
						n[0] = (*i)->normal[0];
						n[1] = (*i)->normal[1];
						n[2] = (*i)->normal[2];
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
						long double vsenkr = sqrt(Ekin*cos(winksenkr)*cos(winksenkr) + E_normal); // add E_normal to component normal to surface
						long double veben = sqrt(Ekin)*sin(winksenkr);
						winksenkr = atan2(veben,vsenkr); // update angle
						Ekin = vsenkr*vsenkr + veben*veben; // update energy
					}
					long double v[3] = {cos(winkeben)*sin(winksenkr), sin(winkeben)*sin(winksenkr), cos(winksenkr)};
					RotateVector(v,n);

					phi = atan2(v[1],v[0]);
					theta = acos(v[2]);
				}
				else{
					printf("Invalid sourcemode: %s",sourcemode.c_str());
					exit(-1);
				}
				
				break;
			}
			x = p1[0];
			y = p1[1];
			z = p1[2];
		};

	private:
		/**
		 * Loads [SOURCE] section from configuration file.
		 *
		 * Reads ::sourcemode, ::kdtree from given STL file, ::r_min, ::r_max, ::phi_min, ::phi_max, ::z_min, ::z_max, ::ActiveTime, ::E_normal.
		 * Determines ::sourcetris and ::sourcearea.
		 *
		 * @param infile File stream to configuration file containing [SOURCE] section
		 */
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
					infile >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime;
					phi_min *= conv;
					phi_max *= conv;
				}
				else if (sourcemode == "volume" || sourcemode == "surface"){
					string sourcefile;
					infile >> sourcefile >> ActiveTime;
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
				printf("Source Area: %LG m^2\n",sourcearea);
			}	
		};
		

		/**
		 * Checks if a point is inside the source volume.
		 *
		 * "customvol/customsurf": checks if point is in the coordinate range.
		 * "volume/surface": Shoots a ray to the bottom of the ::kdtree bounding box and counts collision with surfaces.
		 * Implemented as template to allow long double precision points (simulation) as well as float precision points (triangle vertices).
		 *
		 * @param p Point to check.
		 *
		 * @return Returns true if point lies inside source volume.
		 */
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
				set<TCollision> c;
				if (kdtree->Collision(p1,p2,c)){
					int count = 0;
					for (set<TCollision>::iterator i = c.begin(); i != c.end(); i++){
						if (i->normal[2] > 0) count++;  // surface normal pointing towards point?
						else count--;					// or away from point?
					}
					return (count == -1); // when there's exactly one more normal pointing away from than to it, the point is inside the sourcevol
				}
				else return false;
			}
		}
	

		/**
		 * Calculates energies of low and high field seeking neutrons at point p and updates Hmin_hfs/lfs accordingly.
		 *
		 * @param field TFieldManager to calculate magnet fields
		 * @param p Point to calculate energy at
		 */
		void CalcHmin(TFieldManager &field, long double p[3]){
			set<solid> solids;
			geometry->GetSolids(0, p, solids);

			long double B[4][4];
			field.BField(p[0],p[1],p[2],0,B);
			long double Epot = m_n*gravconst*p[2] + solids.rbegin()->mat.FermiReal*1e-9;
			long double Emag = -(mu_nSI/ele_e)*B[3][0];
			long double Hlfs = Epot + Emag;
			long double Hhfs = Epot - Emag;
			Hmin_lfs = min(Hmin_lfs,Hlfs);
			Hmin_hfs = min(Hmin_hfs,Hhfs);
		}
};

#endif /*GEOMETRY_H_*/
