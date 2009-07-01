#include "main.h"
#include "kdtree.h"
#include <vector>

FILE *REFLECTLOG = NULL;

KDTree geometry, sourcevolume;

vector<string> mat_name;	// dynamic arrays to store material properties
vector<long double> mat_FermiReal;
vector<long double> mat_FermiImag;
vector<long double> mat_DiffProb;

vector<int> mat_map;	// dynamic arrays to associate model ID with material/name
vector<string>name_map;


// transmission through a wall (loss of UCN) with Fermi potential  Mf (real) and PF (im) for a UCN with energy Er perpendicular to the wall (all in neV)
/*long double Transmission(long double Er, long double Mf, long double Pf){
return 0.2e1 * sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) / (Er + sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) + sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf));
}
*/

#define REFLECT_TOLERANCE 1e-6 // if the integration point is farther from the reflection point, the integration will be repeated

// load STL-files and create kd-trees
void LoadGeometry(){
	ifstream infile((inpath+"/geometry.in").c_str());
	string line;
	char c;
	while (infile){
		while (infile && (c = infile.peek()) == ' ') infile.ignore(); // ignore whitespaces
		if ((infile.peek() == '[') && getline(infile,line)){	// parse geometry.in for section header
			if (line.compare(0,11,"[MATERIALS]") == 0){
				string name;
				long double val;
				do{	// parse material list
					while (infile && (c = infile.peek()) == ' ') infile.ignore(); // ignore whitespaces
					if (c == '#' || c == '\n') continue; // skip comments and empty lines
					else if (c == '[') break;	// next section found
					infile >> name;
					mat_name.push_back(name);
					infile >> val;
					mat_FermiReal.push_back(val);
					infile >> val;
					mat_FermiImag.push_back(val);
					infile >> val;
					mat_DiffProb.push_back(val);
				}while(infile && getline(infile,line));
			}
			else if (line.compare(0,10,"[GEOMETRY]") == 0){
				string STLfile;
				unsigned ID;
				string matname;
				char name[80];
				do{	// parse STLfile list
					while (infile && (c = infile.peek()) == ' ') infile.ignore(); // ignore whitespaces
					if (c == '#' || c == '\n') continue;	// skip comments and empty lines
					else if (c == '[') break;	// next section found
					infile >> ID;
					infile >> STLfile;
					infile >> matname;
					STLfile = inpath + '/' + STLfile;
					geometry.ReadFile(STLfile.c_str(),ID,name);
					if (ID >= mat_map.size()){
						mat_map.resize(ID+1);
						name_map.resize(ID+1,"");
					}
					for (unsigned i = 0; i < mat_name.size(); i++){
						if (matname == mat_name[i]){
							mat_map[ID] = i;
							name_map[ID] = name;
							break;
						}
						else if (i+1 == mat_name.size()){
							fprintf(stderr,"Material %s used but not defined in geometry.in!",matname.c_str());
							exit(-1);
						}
					}
				}while(infile && getline(infile,line));
				geometry.Init();
			}
			else if (line.compare(0,8,"[SOURCE]") == 0){
				string name;
				long double p[3];
				do{	// parse source line
					while (infile && (c = infile.peek()) == ' ') infile.ignore(); // ignore whitespaces
					if (c == '#' || c == '\n') continue; // skip comments and empty lines
					else if (c == '[') break;	// next section found
					infile >> name;
					sourcevolume.ReadFile((inpath + '/' + name).c_str(),0);
					infile >> p[0];
					infile >> p[1];
					infile >> p[2];
					sourcevolume.Init(p);
				}while(infile && getline(infile,line));
			}
			else getline(infile,line);
		}
		else getline(infile,line);
	}
	
/*	
	if(reflektlog == 1){
		ostringstream reflectlogfile;
		reflectlogfile << outpath << "/" << jobnumber << "reflect.out";
		REFLECTLOG = fopen(reflectlogfile.str().c_str(),mode_w);
		fprintf(REFLECTLOG,"t r z phi x y diffuse vabs Eges Erefl winkeben winksenkr vr vz vtang phidot dvabs\n"); // Header for Reflection File
	}	*/
}

// return a random point in sourcevolume
void RandomPointInSourceVolume(long double &r, long double &phi, long double &z){
	long double p[3];
	do{	
		p[0] = mt_get_double(v_mt_state)*(sourcevolume.hi[0] - sourcevolume.lo[0]) + sourcevolume.lo[0];
		p[1] = mt_get_double(v_mt_state)*(sourcevolume.hi[1] - sourcevolume.lo[1]) + sourcevolume.lo[1];
		p[2] = mt_get_double(v_mt_state)*(sourcevolume.hi[2] - sourcevolume.lo[2]) + sourcevolume.lo[2];
	}while(!sourcevolume.PointInVolume(p));
	r = sqrt(p[0]*p[0] + p[1]*p[1]);
	phi = atan2(p[1],p[0]);
	z = p[2];
}

// check step y1->y2 for reflection, return -1 for failed reflection (step has to be repeated with smaller stepsize)
//											0 for no reflection
//											1 for reflection successful
short ReflectCheck(long double x1, long double *y1, long double &x2, long double *y2){
	long double p1[3] = {y1[1]*cos(y1[5]), y1[1]*sin(y1[5]), y1[3]}; // cart. coords
	long double p2[3] = {y2[1]*cos(y2[5]), y2[1]*sin(y2[5]), y2[3]};
	long double s, normal_cart[3];
	unsigned ID;
	if (geometry.Collision(p1,p2,s,normal_cart,ID)){ // search in the kdtree for reflection
		
		long double distnormal = abs((p2[0] - p1[0])*normal_cart[0] + (p2[1] - p1[1])*normal_cart[1] + (p2[2] - p1[2])*normal_cart[2]); // shortest distance of p1 to surface
		if (s*distnormal > REFLECT_TOLERANCE){ // if p1 too far from surface
			s -= REFLECT_TOLERANCE/distnormal/1e10; // decrease s by a small amount
			x2 = x1 + s*(x2-x1); // write smaller integration time into x2
			return -1; // return fail to repeat integration step
		}
		
		long double v[3] = {y1[2], y1[1]*y1[6], y1[4]}; // reflection velocity in local cylindrical coord.
		long double vabs = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		long double normal[3] = {normal_cart[0]*cos(y1[5]) + normal_cart[1]*sin(y1[5]), // transform normal into local cyl. coord.
								-normal_cart[0]*sin(y1[5]) + normal_cart[1]*cos(y1[5]), 
								normal_cart[2]};
		long double vnormal = v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]; // velocity normal to reflection plane
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		unsigned mat = mat_map[ID]; // get material-index from ID->material map
		
		
		//************ handle different absorption characteristics of materials ****************
		long double prob = mt_get_double(v_mt_state);	
		if(!reflekt){
			stopall = 1;
			kennz = ID;
			printf("Particle hit %s (no reflection) at r=%LG z=%LG\n",name_map[ID].c_str(),y1[1],y1[3]);
			fprintf(LOGSCR,"Particle hit %s (no reflection) at r=%LG z=%LG\n",name_map[ID].c_str(),y1[1],y1[3]);
			return 1;
		}
		else if (prob < Transmission(Enormal*1e9,mat_FermiReal[mat],mat_FermiImag[mat])) // statistical absorption
		{
			stopall = 1;
			kennz = ID;
			printf("Statistical absorption at %s (%s)!\n",name_map[ID].c_str(),mat_name[mat].c_str());
			fprintf(LOGSCR,"Statistical absorption at %s (%s)!\n",name_map[ID].c_str(),mat_name[mat].c_str());
			return 1;
		}		
		
		//*************** specular reflexion ************
		prob = mt_get_double(v_mt_state);
		if ((diffuse == 1) || ((diffuse==3)&&(prob >= mat_DiffProb[mat])))
		{
	    	printf(" pol %d t=%LG Erefl=%LG neV r=%LG z=%LG tol=%LG m ",polarisation,x1,Enormal*1e9,y1[1],y1[3],s*distnormal);
			fprintf(LOGSCR,"pol %d t=%LG Erefl=%LG neV r=%LG z=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3]);
			nrefl++;
			if(reflektlog == 1)
				fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
									x1,y1[1],y1[3],y1[5],p1[0],p1[1],vabs,H,Enormal*1e9,(long double)0.0,(long double)0.0,y1[2],y1[4],y1[1]*y1[6],y1[6]);
			v[0] -= 2*vnormal*normal[0]; // reflect velocity
			v[1] -= 2*vnormal*normal[1];
			v[2] -= 2*vnormal*normal[2];
		}
		
		//************** diffuse reflection ************
		else if ((diffuse == 2) || ((diffuse == 3)&&(prob < mat_DiffProb[mat])))
		{
			long double winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi); // generate random reflection angles
			long double winksenkr = acos( cos(0.0) - mt_get_double(v_mt_state) * (cos(0.0) - cos(0.5 * pi)));
			if (vnormal > 0) winksenkr += pi; // if normal points out of volume rotate by 180 degrees
			v[0] = vabs*cos(winkeben)*sin(winksenkr);	// new velocity with respect to z-axis
			v[1] = vabs*sin(winkeben)*sin(winksenkr);
			v[2] = vabs*cos(winksenkr);

			long double cosalpha = normal[2], sinalpha = sqrt(1 - cosalpha*cosalpha);	// rotation angle (angle between z and n)
			if (sinalpha > 1e-30){ // when normal not parallel to z-axis rotate new velocity into the coordinate system where the normal is the z-axis
				long double a[2] = {-normal[1]/sinalpha, normal[0]/sinalpha};	// rotation axis (z cross n), a[2] = 0
				long double vtemp[3] = {v[0],v[1],v[2]};
				// rotate velocity vector
				v[0] = (cosalpha + a[0]*a[0]*(1 - cosalpha))*	vtemp[0] +  a[0]*a[1]*(1 - cosalpha)*				vtemp[1] + a[1]*sinalpha*	vtemp[2];
				v[1] =  a[1]*a[0]*(1 - cosalpha)*				vtemp[0] + (cosalpha + a[1]*a[1]*(1 - cosalpha))*	vtemp[1] - a[0]*sinalpha*	vtemp[2];
				v[2] = -a[1]*sinalpha*							vtemp[0] +  a[0]*sinalpha*							vtemp[1] + cosalpha*		vtemp[2];
			}
	       	printf(" pol %d t= %LG Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG tol=%LG m ",polarisation,x1,Enormal*1e9,y1[1],y1[3],winkeben/conv,winksenkr/conv,s*distnormal);
	       	fprintf(LOGSCR,"pol %d t= %LG Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3],winkeben/conv,winksenkr/conv);
	       	nrefl++;
	       	if(reflektlog == 1)
				fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 2 %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
									x1,y1[1],y1[3],y1[5],p1[0],p1[1],vabs,H,Enormal*1e9,winkeben,winksenkr,y1[2],y1[4],y1[1]*y1[6],y1[6]);
		}
		y1[2] = v[0]; // write new velocity into y1-vector
		y1[4] = v[2];
		y1[6] = v[1]/y1[1];
		
		return 1;
	}
	return 0;
}
