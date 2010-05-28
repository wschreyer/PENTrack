#include "main.h"
#include "kdtree.h"
#include <complex>

#define REFLECT_TOLERANCE 1e-10 // if the integration point is farther from the reflection point, the integration will be repeated


int reflektlog = 0;
FILE *REFLECTLOG = NULL;

KDTree geometry(Log), sourcevolume(Log);
list<Triangle*> sourcetris;
long double sourcearea;

vector<solid> solids;	// dynamic array to associate solids with material/kennzahl/name
solid defaultsolid = {"Restgas",{"Restgas",0,0,0},0};
solid *currentsolid = &defaultsolid;

vector<int> kennz_counter[3];

long double r_min = INFINITY, r_max = -INFINITY, 
			phi_min = INFINITY, phi_max = -INFINITY, 
			z_min = INFINITY, z_max = -INFINITY;
			
string sourcemode; // volume/surface/customvol/customsurf


// transmission through a wall with Fermi potential  Mf (real) and PF (im) for a UCN with energy Er perpendicular to the wall (all in neV)
long double Transmission(long double Er, long double Mf, long double Pf){
	complex<long double> V(Mf,Pf), sEV = sqrt(Er - V);
	long double sE = sqrt(Er);
	return 1 - norm((sE - sEV)/(sE + sEV));
}

// absorption probability of neutron flying distance l (in m) through Fermi potential Mf + i*Pf (all in neV)
long double Absorption(long double E, long double Mf, long double Pf, long double l){
	complex<long double> V(Mf,Pf);
	return 1 - exp( 2*imag(sqrt(2*m_n*1e-9*(E - V))) / (hquer/ele_e) * l ); // absorption length 2*Im(sqrt(2m(E-V))/hquer)
}


// load STL-files and create kd-trees
void LoadGeometry(){
	cout << endl;
	ifstream infile((inpath+"/geometry.in").c_str());
	string line;
	vector<material> materials;	// dynamic array to store material properties
	char c;
	while (infile.good()){
		infile >> ws; // ignore whitespaces
		c = infile.peek();
		if ((infile.peek() == '[') && getline(infile,line).good()){	// parse geometry.in for section header
			if (line.compare(0,11,"[MATERIALS]") == 0){
				do{	// parse material list
					infile >> ws; // ignore whitespaces
					c = infile.peek();
					if (c == '#') continue; // skip comments
					else if (!infile.good() || c == '[') break;	// next section found
					material mat;
					infile >> mat.name >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb;
					materials.push_back(mat);
				}while(infile.good() && getline(infile,line).good());
			}
			else if (line.compare(0,10,"[GEOMETRY]") == 0){
				string STLfile;
				string matname;
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
					STLfile = inpath + '/' + STLfile;
					for (unsigned i = 0; i < materials.size(); i++){
						if (matname == materials[i].name){
							model.mat = materials[i];
							geometry.ReadFile(STLfile.c_str(),solids.size(),name);
							model.name = name;
							break;
						}
						else if (i+1 == materials.size()){
							Log("Material %s used for %s but not defined in geometry.in!",matname.c_str(),name);
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
				geometry.Init();
			}
			else if (line.compare(0,8,"[SOURCE]") == 0){
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
						sourcevolume.ReadFile((inpath + '/' + sourcefile).c_str(),0);
						sourcevolume.Init();
					}
				}while(infile.good() && getline(infile,line).good());
			}
			else if (line.compare(0,8,"[FIELDS]") == 0){
				do{	// parse STLfile list
					infile >> ws; // ignore whitespaces
					c = infile.peek();
					if (c == '#') continue;	// skip comments
					else if (!infile.good() || c == '[') break;	// next section found
					string ft;
					infile >> ft;
					fieldvaltab.push_back(ft);
				}while(infile.good() && getline(infile,line).good());
			}
			else getline(infile,line);
		}
		else getline(infile,line);
	}
	
	if (sourcemode == "customsurf" || sourcemode == "surface"){
		sourcearea = 0; // add all triangles with at least one vertex in the source volume to sourcetris list
		for (list<Triangle*>::iterator i = geometry.alltris.begin(); i != geometry.alltris.end(); i++){
			if (InSourceVolume((*i)->vertex[0]) || InSourceVolume((*i)->vertex[1]) || InSourceVolume((*i)->vertex[2])){
				sourcetris.push_back(*i);
				long double n[3];
				(*i)->CalcNormal(n);
				sourcearea += sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			}
		}
		Log("Source Area: %LG m²\n",sourcearea);
	}
	
	if(reflektlog == 1){
		ostringstream reflectlogfile;
		reflectlogfile << outpath << "/" << setw(8) << setfill('0') << jobnumber << setw(0) << "reflect.out";
		REFLECTLOG = fopen(reflectlogfile.str().c_str(),mode_w);
		fprintf(REFLECTLOG,"t r z phi x y diffuse vabs Eges Erefl winkeben winksenkr vr vz vphi phidot dvabs Transmprop\n"); // Header for Reflection File
	}	
}

// return a random point in sourcevolume
void RandomPointInSourceVolume(long double &r, long double &phi, long double &z, long double &alpha, long double &gamma){
	long double p1[3];
	bool valid;
	do{	// repeat until valid is set true
		valid = false;
		list<TCollision> c;
		if (sourcemode == "volume"){ // random point inside a STL-volume
			p1[0] = mt_get_double(v_mt_state)*(sourcevolume.hi[0] - sourcevolume.lo[0]) + sourcevolume.lo[0]; // random point 
			p1[1] = mt_get_double(v_mt_state)*(sourcevolume.hi[1] - sourcevolume.lo[1]) + sourcevolume.lo[1];
			p1[2] = mt_get_double(v_mt_state)*(sourcevolume.hi[2] - sourcevolume.lo[2]) + sourcevolume.lo[2];
			r = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
			phi = atan2(p1[1],p1[0]);
			z = p1[2];
			if (InSourceVolume(r,phi,z)){
				alpha = alphas + (mt_get_double(v_mt_state)) * (alphae - alphas); // constant angular distribution
				gamma = acosl(cosl(gammas) - mt_get_double(v_mt_state) * (cosl(gammas) - cosl(gammae))); // isotropic emission characteristics
				valid = true;
			}
		}
		else if (sourcemode == "surface" || sourcemode == "customsurf"){ // random point on a surface inside custom or STL volume
			long double RandA = mt_get_double(v_mt_state)*sourcearea;
			long double n[3], CurrA, SumA = 0;
			list<Triangle*>::iterator i;
			for (i = sourcetris.begin(); i != sourcetris.end(); i++){
				(*i)->CalcNormal(n);
				CurrA = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
				SumA += CurrA;
				if (RandA <= SumA) break;
			}
			long double a = mt_get_double(v_mt_state);
			long double b = mt_get_double(v_mt_state)*(1-a);
			for (int j = 0; j < 3; j++){
				n[j] /= CurrA;
				p1[j] = (*i)->vertex[0][j] + a*((*i)->vertex[1][j] - (*i)->vertex[0][j]) + b*((*i)->vertex[2][j] - (*i)->vertex[0][j]);
				p1[j] += REFLECT_TOLERANCE*n[j];
			}
			r = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
			phi = atan2(p1[1],p1[0]);
			z = p1[2];
			if (InSourceVolume(r, phi, z)){				
				long double winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi); // generate random reflection angles
				long double winksenkr = acos( cos(0.0) - mt_get_double(v_mt_state) * (cos(0.0) - cos(0.5 * pi)));
				long double v[3] = {cos(winkeben)*cos(winksenkr), sin(winkeben)*cos(winksenkr), sin(winksenkr)}; // new velocity with respect to z-axis
				RotateVector(v,n);
				
				alpha = acos((p1[0]*v[0] + p1[1]*v[1])/r);
				gamma = acos(v[2]);
				valid = true;
			}
		}
		else if (sourcemode == "customvol"){ // random point inside a volume defined by a custom parameter range
			r = powl(mt_get_double(v_mt_state) * (powl(r_max,3) - powl(r_min,3)) + powl(r_min,3),1.0/3.0); // weighting because of the volume element and a r^2 probability outwards
			phi = mt_get_double(v_mt_state)*(phi_max - phi_min) + phi_min;
			z = mt_get_double(v_mt_state)*(z_max - z_min) + z_min;
			p1[0] = r*cos(phi);
			p1[1] = r*sin(phi);
			p1[2] = z;
			alpha = alphas + (mt_get_double(v_mt_state)) * (alphae - alphas); // constant angular distribution
			gamma = acos(cos(gammas) - mt_get_double(v_mt_state) * (cos(gammas) - cos(gammae))); // isotropic emission characteristics
			valid = true;
		}
		else{
			Log("Invalid sourcemode: %s",sourcemode.c_str());
			exit(-1);
		}
		
		long double p2[3] = {p1[0], p1[1], geometry.lo[2]};
		if (valid && geometry.Collision(p1,p2,c)){	// test if random point is inside solid
			int count = 0;
			for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
				if (i->normal[2] > 0) count++; // count surfaces whose normals point to random point
				else count--; // count surfaces whose normals point away from random point
			}
			valid = (count == 0); // count is zero, if all surfaces are closed -> point does not lie in a solid
		}
	}while(!valid);
	currentsolid = &defaultsolid;
}


// check if a point is inside the source volume
bool InSourceVolume(const float p[3]){
	if (sourcemode == "customvol" || sourcemode == "customsurf"){
		long double r = sqrt(p[0]*p[0] + p[1]*p[1]);
		long double phi = atan2(p[1],p[0]); 
		return (r >= r_min && r <= r_max && 
				phi >= phi_min && phi <= phi_max && 
				p[2] >= z_min && p[2] <= z_max); // check if point is in custom paramter range
	}
	else{ // shoot a ray to the bottom of the sourcevol bounding box and check for intersection with the sourcevol surface
		long double p1[3] = {p[1], p[2], p[3]};
		long double p2[3] = {p[0], p[1], sourcevolume.lo[2]}; 
		list<TCollision> c;
		if (sourcevolume.Collision(p1,p2,c)){
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
	
// same as above for cyl. coordinates
bool InSourceVolume(long double r, long double phi, long double z){
	if (sourcemode == "customvol" || sourcemode == "customsurf")
		return (r >= r_min && r <= r_max && phi >= phi_min && phi <= phi_max && z >= z_min && z <= z_max);
	else{
		long double p1[3] = {r*cos(phi), r*sin(phi), z};
		long double p2[3] = {p1[0], p1[1], sourcevolume.lo[2]};
		list<TCollision> c;
		if (sourcevolume.Collision(p1,p2,c)){
			int count = 0;
			for (list<TCollision>::iterator i = c.begin(); i != c.end(); i++){
				if (i->normal[2] > 0) count++;
				else count--;
			}
			return (count == -1);
		}
		else return false;
	}
}


// delete ignored collision from colls list
void DeleteIgnored(long double x1, long double h, list<TCollision> &colls){
	list<TCollision>::iterator it;
	for (it = colls.begin(); it != colls.end(); it++){
		if (!solids[(*it).ID].ignoretimes.empty()){
			long double x = x1 + (*it).s*h;
			set<int> itimes = solids[(*it).ID].ignoretimes;
			if ((x < FillingTime && itimes.count(1) > 0) ||
				(x >= FillingTime && x < FillingTime + CleaningTime && itimes.count(2) > 0) ||
				(x >= FillingTime + CleaningTime && x < FillingTime + CleaningTime + RampUpTime && itimes.count(3) > 0) ||
				(x >= FillingTime + CleaningTime + RampUpTime && x < FillingTime + CleaningTime + RampUpTime + FullFieldTime && itimes.count(4) > 0) ||
				(x >= FillingTime + CleaningTime + RampUpTime + FullFieldTime && x < FillingTime + CleaningTime + RampUpTime + FullFieldTime + RampDownTime && itimes.count(5) > 0) ||
				(x >= FillingTime + CleaningTime + RampUpTime + FullFieldTime + RampDownTime && itimes.count(6) > 0))
				it = --colls.erase(it);
		}
	}
}

// check if p1 closer to collision than REFLECT_TOLERANCE (return true), else reduce h accordingly and return false
bool CloseEnough(long double p1[3], long double p2[3], long double *y1, long double *y2, long double &h, TCollision &coll){
	long double u[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
	long double distnormal = abs(u[0]*coll.normal[0] + u[1]*coll.normal[1] + u[2]*coll.normal[2]);
	long double s = coll.s;
	if (s*distnormal < REFLECT_TOLERANCE) return true;
	long double v[3];
	CylKartCoord(y1[2] + s*(y2[2] - y1[2]), y1[1]*y1[6] + s*(y2[1]*y2[6] - y1[1]*y1[6]), y1[4] + s*(y2[4] - y1[4]),
					y1[5] + s*(y2[5] - y1[5]), &v[0],&v[1],&v[2]);
	long double vnormal = abs(v[0]*coll.normal[0] + v[1]*coll.normal[1] + v[2]*coll.normal[2]);
	h = h*coll.s - REFLECT_TOLERANCE/1e3/vnormal;
	return false;
}

// check step y1->y2 for reflection, return true when step has to be repeated (with smaller stepsize or reflected velocity)
//									 return false for no reflection
bool ReflectCheck(long double x1, long double &h, long double *y1, long double *y2, int &itercount){
	long double p1[3] = {y1[1]*cos(y1[5]), y1[1]*sin(y1[5]), y1[3]}; // cart. coords
	
	if (!geometry.PointInBox(p1)){
		kennz=KENNZAHL_HIT_BOUNDARIES;  
		stopall=1;
		Log("\nParticle has hit outer boundaries: Stopping it! t=%LG r=%LG phi=%LG z=%LG\n",x1,y1[1],y1[5]/conv,y1[3]);
		return false;
	}
	
	long double p2[3] = {y2[1]*cos(y2[5]), y2[1]*sin(y2[5]), y2[3]};
	long double normal[3];
	list<TCollision> colls;
	if (geometry.Collision(p1,p2,colls)){ // search in the kdtree for reflection
		DeleteIgnored(x1,h,colls);
		if (colls.empty())
			return false;
			
		if (!CloseEnough(p1,p2,y1,y2,h,colls.front()))// if p1 too far from surface
			return true; // repeat integration step
		else{
			if (colls.size() > 1){
				if (CloseEnough(p1,p2,y1,y2,h,*++colls.begin())){
					cerr << "Solids closer than REFLECT_TOLERANCE!!\n";
					exit(-1);
				}
				else
					return true;
			}
			Log("\n#%d Particle hit %s (material %s) at t=%LG r=%LG z=%LG tries=%i: ",
				iMC,solids[colls.front().ID].name.c_str(),solids[colls.front().ID].mat.name.c_str(),x1,y1[1],y1[3],itercount);
			long double *normal_cart = colls.front().normal;
			normal[0] = normal_cart[0]*cos(y1[5]) + normal_cart[1]*sin(y1[5]); // transform normal into local cyl. coord.
			normal[1] = -normal_cart[0]*sin(y1[5]) + normal_cart[1]*cos(y1[5]); 
			normal[2] = normal_cart[2];
		}
	}
	else
		return false;

	unsigned ID = colls.front().ID;
	material *mat = &solids[ID].mat; // get material
	
	//************ absorption of electrons/protons ***********
	if (!reflekt||protneut!=1){
		stopall = 1;
		kennz = solids[ID].kennz;
		Log("Absorption!");
		return false;
	}			

	
//************ handle different absorption characteristics of materials ****************

	long double v[3] = {y1[2], y1[1]*y1[6], y1[4]}; // reflection velocity in local cylindrical coord.
	long double vabs = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	long double vnormal = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]; // velocity normal to reflection plane
	long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
	
	//************** statistical absorption ***********
	long double prob = mt_get_double(v_mt_state);
	long double Trans=Transmission(Enormal*1e9,mat->FermiReal,mat->FermiImag);	
	if (prob < Trans)
	{
/*		stopall = 1;
		kennz = solids[ID].kennz;*/
		Log("Transmission! pol %d Erefl=%LG neV Transprob=%LG\n",polarisation,Enormal*1e9,Trans);
		if (vnormal > 0)
			currentsolid = &defaultsolid;
		else
			currentsolid = &solids[ID];
		return false;
	}		
	
	//*************** specular reflexion ************
	prob = mt_get_double(v_mt_state);
	if ((diffuse == 1) || ((diffuse==3)&&(prob >= mat->DiffProb)))
	{
    	Log("Specular reflection! pol %d Erefl=%LG neV Transprop=%LG",polarisation,Enormal*1e9,Trans);
		nrefl++;
		v[0] -= 2*vnormal*normal[0]; // reflect velocity
		v[1] -= 2*vnormal*normal[1];
		v[2] -= 2*vnormal*normal[2];
		if(reflektlog == 1)
			fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %.20LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
								x1,y1[1],y1[3],y1[5],p1[0],p1[1],vabs,H,Enormal*1e9,(long double)0.0,(long double)0.0,y1[2],y1[4],y1[1]*y1[6],y1[6],sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])-vabs,Trans);
	}
	
	//************** diffuse reflection ************
	else if ((diffuse == 2) || ((diffuse == 3)&&(prob < mat->DiffProb)))
	{
		long double winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi); // generate random reflection angles
		long double winksenkr = acos( cos(0.0) - mt_get_double(v_mt_state) * (cos(0.0) - cos(0.5 * pi)));
		if (vnormal > 0) winksenkr += pi; // if normal points out of volume rotate by 180 degrees
		v[0] = vabs*cos(winkeben)*cos(winksenkr);	// new velocity with respect to z-axis
		v[1] = vabs*sin(winkeben)*cos(winksenkr);
		v[2] = vabs*sin(winksenkr);
		RotateVector(v,normal);
       	Log("Diffuse reflection! pol %d Erefl=%LG neV w_e=%LG w_s=%LG Transprop=%LG",polarisation,Enormal*1e9,winkeben/conv,winksenkr/conv,Trans);
       	nrefl++;
       	if(reflektlog == 1)
			fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %.20LG %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
								x1,y1[1],y1[3],y1[5],p1[0],p1[1],vabs,H,Enormal*1e9,winkeben,winksenkr,y1[2],y1[4],y1[1]*y1[6],y1[6],sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])-vabs,Trans);
	}
	y1[2] = v[0]; // write new velocity into y1-vector
	y1[4] = v[2];
	y1[6] = v[1]/y1[1];
	
	return true;
}


void AbsorbCheck(long double *y1, long double *y2){
	if (protneut == NEUTRON){
		long double E = (m_n*gravconst*y1[3] + 0.5*m_n*(y1[2]*y1[2] + y1[4]*y1[4] + y1[1]*y1[6]*y1[1]*y1[6]) - mu_n*Bws)*1e9;
		long double l = sqrt(pow(y2[1] - y1[1],2) + pow(y2[3] - y1[3],2) + pow(y2[5]*y2[1] - y1[5]*y1[1],2));
		if (mt_get_double(v_mt_state) < Absorption(E, currentsolid->mat.FermiReal, currentsolid->mat.FermiImag, l)){
			stopall = 1;
			kennz = currentsolid->kennz;
			Log("Absorption in %s (material %s)! Depth: %LG m\n",currentsolid->name.c_str(), currentsolid->mat.name.c_str(), mt_get_double(v_mt_state)*l);
		}
	}
}

// rotate coordinate system so, that new z-axis lies on NORMALIZED vector n
void RotateVector(long double v[3], long double n[3])
{
	long double cosalpha = n[2], sinalpha = sqrt(1 - cosalpha*cosalpha);	// rotation angle (angle between z and n)
	if (sinalpha > 1e-30){ // when normal not parallel to z-axis rotate new velocity into the coordinate system where the normal is the z-axis
		long double a[2] = {-n[1]/sinalpha, n[0]/sinalpha};	// rotation axis (z cross n), a[2] = 0
		long double vtemp[3] = {v[0],v[1],v[2]};
		// rotate velocity vector
		v[0] = (cosalpha + a[0]*a[0]*(1 - cosalpha))*	vtemp[0] +  a[0]*a[1]*(1 - cosalpha)*				vtemp[1] + a[1]*sinalpha*	vtemp[2];
		v[1] =  a[1]*a[0]*(1 - cosalpha)*				vtemp[0] + (cosalpha + a[1]*a[1]*(1 - cosalpha))*	vtemp[1] - a[0]*sinalpha*	vtemp[2];
		v[2] = -a[1]*sinalpha*							vtemp[0] +  a[0]*sinalpha*							vtemp[1] + cosalpha*		vtemp[2];
	}
}

//======== The endcodes of trajectories are added up here regarding the particle typ =======================================
void IncrementCodes(int kennz)
{	// protneut % 3 = {0, 1, 2} <<<	ELECTRONS (=6) % 3 = 0
	// 								NEUTRON   (=1) % 3 = 1	
	// 								PROTON    (=2) % 3 = 2
	kennz_counter[protneut % 3][kennz]++;
}
//======== end of IncrementCodes ===========================================================================================

//======== Output of the endcodes ==========================================================================================
void OutputCodes(int iMC){
	int ncount = accumulate(kennz_counter[1].begin(),kennz_counter[1].end(),0);
	int pcount = accumulate(kennz_counter[2].begin(),kennz_counter[2].end(),0);
	int ecount = accumulate(kennz_counter[0].begin(),kennz_counter[0].end(),0);
	Log("\nThe calculations of %li particle(s) yielded:\n"
	       "endcode:  of %4i neutron(s) ; of %4i proton(s) ; of %4i electron(s)\n"
	       "   0 %12i %20i %19i 		(were not categorized)\n"
	       "   1 %12i %20i %19i 		(did not finish)\n"
	       "   2 %12i %20i %19i 		(hit outer boundaries)\n"
	       "   3 %12i %20i %19i 		(left field boundaries)\n"
	       "   4 %12i %20i %19i 		(decayed)\n"
	       "   5 %12i %20i %19i 		(found no initial position)\n",
	       (iMC - 1 + 2 * decay.counter),
	       ncount, pcount, ecount,
	       kennz_counter[1][0], kennz_counter[2][0], kennz_counter[0][0],
	       kennz_counter[1][1], kennz_counter[2][1], kennz_counter[0][1],
	       kennz_counter[1][2], kennz_counter[2][2], kennz_counter[0][2],
	       kennz_counter[1][3], kennz_counter[2][3], kennz_counter[0][3],
	       kennz_counter[1][4], kennz_counter[2][4], kennz_counter[0][4],
	       kennz_counter[1][5], kennz_counter[2][5], kennz_counter[0][5]);
	for (unsigned i = 6; i < kennz_counter[0].size(); i++){
		string solidnames;
		for (vector<solid>::iterator it = solids.begin(); it != solids.end(); it++)
			if (it->kennz == i)
				solidnames += '/' + it->name;
		Log("  %2i %12i %20i %19i		(were statistically absorbed by %s)\n",
				i,kennz_counter[1][i],kennz_counter[2][i],kennz_counter[0][i],solidnames.c_str()+1);
	}
}
//======== end of OutputCodes ==============================================================================================


void Snapshooter(long double x2, long double *ystart, long double H)
{
	if((snapshots.count((int) x2)>0) && ((int) x2 != snapshotsdone))
	{
		
		vend    = sqrtl(fabsl(ystart[2]*ystart[2]+ystart[1]*ystart[1]*ystart[6]*ystart[6]+ystart[4]*ystart[4]));
		phiend  = fmodl(ystart[5], 2*pi);
		if (phiend<0)                    // from 0 to 360
			phiend=2*pi + phiend;

		if(protneut == NEUTRON)                // n
			H = (M*gravconst*ystart[3]+0.5*M*vend*vend-mu_n*Bws)*1E9 ;       // Energie in neV
		else if(protneut == PROTON)           // p			
				H= (0.5*m_p*vend*vend);           // Energie in eV for p		
		else if(protneut == ELECTRONS)           // p,e
			H= c_0*c_0  * M * (1/sqrtl(1-v_n*v_n/(c_0*c_0))-1);                                        // rel Energie in eV

	
	Log("\n Snapshot at %LG s \n", x2);
	
	long double dt = x2 - xstart; // simulation time dt

	if(vend>0) gammaend= acosl(ystart[4]/vend);
	else gammaend=0;
	alphaend= atan2l(ystart[6]*ystart[1],ystart[2]);
	
	// calculate spin flip lifetime tauSF and influence on lifetime measurement 
	long double tauSF = -x2/logl(1-BFflipprob);
	long double dtau=tau-1/(1/tau+1/tauSF) ;
	 
	if(tauSF < -9e99) tauSF = INFINITY; // "-INF"?
	
	
	// output of end values
	fprintf(SNAP,"%i %i %i %i "
	               "%LG %LG %LG %LG %LG "
	               "%LG %LG %LG "
	               "%LG %LG %LG %LG %LG "
	               "%LG %LG %LG %LG %LG "
	               "%LG %i %i %LG %LG "
	               "%li %LG  "
	               "%LG %LG %LG %LG %LG %LG\n",
	               jobnumber, iMC, protneut, polarisation,   
	               xstart, r_n, phi_n, z_n, NeutEnergie*1.0e9,
	               v_n, alpha, gammaa, 
	               ystart[1], phiend, ystart[3], ystart[1]*cosl(ystart[5]),ystart[1]*sinl(ystart[5]),
	               vend, alphaend, gammaend, x2, dt,
	               H, kennz, NSF, RodFieldMultiplicator, BFflipprob,
	               nrefl, trajlengthsum,
	               Hstart, Hmax, BFeldSkal, EFeldSkal, tauSF, dtau);
	
	fflush(SNAP);
	snapshotsdone=(int) x2;	
	}	
}
