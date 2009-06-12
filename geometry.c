#include "main.h"
#include "kdtree.h"

KDTree geometry, sourcevolume;

// transmission through a wall (loss of UCN) with Fermi potential  Mf (real) and PF (im) for a UCN with energy Er perpendicular to the wall (all in neV)
/*long double Transmission(long double Er, long double Mf, long double Pf){
return 0.2e1 * sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) / (Er + sqrtl(Er) * sqrtl(0.2e1) * sqrtl(sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf) + Er - Mf) + sqrtl(Er * Er - 0.2e1 * Er * Mf + Mf * Mf + Pf * Pf));
}
*/

#define SURFACE_WALL 0
#define SURFACE_ABSORBER 1
#define SURFACE_UCNDET 2
#define SURFACE_PROTDET 3

#define REFLECT_TOLERANCE 0.00001

void LoadGeometry(){
	string STLfile = inpath + "/heliumtankMay09.STL";
	geometry.ReadFile(STLfile.c_str(),SURFACE_WALL);
	STLfile = inpath + "/absorber.STL";
	geometry.ReadFile(STLfile.c_str(),SURFACE_ABSORBER);
	STLfile = inpath + "/UCNdet.STL";
	geometry.ReadFile(STLfile.c_str(),SURFACE_UCNDET);
	STLfile = inpath + "/protdet.STL";
	geometry.ReadFile(STLfile.c_str(),SURFACE_PROTDET);
	//...	
	geometry.Init();
	
	STLfile = inpath + "/source.STL";
	sourcevolume.ReadFile(STLfile.c_str(),0);
	long double ControlPoint[3] = {0.3,0,0.5};
	sourcevolume.Init(ControlPoint);	
}

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
	int surface;
	if (geometry.Collision(p1,p2,s,normal_cart,surface)){ // search in the kdtree for reflection
		
		long double distnormal = abs((p2[0] - p1[0])*normal_cart[0] + (p2[1] - p1[1])*normal_cart[1] + (p2[2] - p1[2])*normal_cart[2]);
		if (abs(s*distnormal) > REFLECT_TOLERANCE){ // if first point to far from surface
			s -= REFLECT_TOLERANCE/distnormal/100; // decrease step size
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
		
		// handle different absorption characteristics of materials
		long double prob = mt_get_double(v_mt_state);	
		long double MaterialDiffProb = DiffProb; // use DiffProb as default	
		switch (surface){
			case SURFACE_WALL: // particle hit the heliumtank wall
				if(!reflekt) // no reflection (e.g. for protons/electrons)
				{
					stopall = 1;
					kennz = KENNZAHL_HIT_WALL;
					printf("Particle hit the walls (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);
					fprintf(LOGSCR,"Particle hit the walls (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);	
					return 1;	
				}				
				else if (prob < Transmission(Enormal*1e9,FPrealNocado,FPimNocado)) // statistical absorption
				{
					stopall = 1;
					kennz = KENNZAHL_ABSORBED;
					printf("Statistical absorption!!!");
					fprintf(LOGSCR,"Statistical absorption!!!");
					return 1;
				}
				break;
			case SURFACE_ABSORBER: // particle hit the absorber
				if(!reflekt) // no reflection (e.g. for protons/electrons)
				{
					stopall = 1;
					kennz = KENNZAHL_HIT_WALL;
					printf("Particle hit the walls (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);
					fprintf(LOGSCR,"Particle hit the walls (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);	
					return 1;	
				}				
				else if(prob < Transmission(Enormal*1e9,Mf,Pf)){ // statistical absorption
					AbsorberHits++;
					stopall=1;
					kennz=KENNZAHL_EATEN;     // absorber hit
					printf("Absorber hit and absorption!\nvr = %LG m/s\nEr = %LG neV\n", vnormal,Enormal*1e9);
					fprintf(LOGSCR,"Absorber hit and absorption!\nvr = %LG m/s\nEr = %LG neV\n", vnormal,Enormal*1e9);
					return 1;
				}
				else{
					AbsorberHits++;		
					printf("Absorber hit but no absorption!\nvr = %LG m/s\nEr = %LG neV\n", vnormal,Enormal*1e9);
					fprintf(LOGSCR,"Absorber hit but no absorption!\nvr = %LG m/s\nEr = %LG neV\n", vnormal,Enormal*1e9);
				}			
				break;
			case SURFACE_UCNDET: // particle hit the UCN detector, 100% efficiency is assumed
				stopall = 1;
				kennz = 12;
				return 1;
			case SURFACE_PROTDET: // particle hit the proton detector
				if(!reflekt) // no reflection (e.g. for protons/electrons)
				{
					stopall = 1;
					if (v[2] < 0)
						kennz = KENNZAHL_DETECTOR_TOP;
					else
						kennz = KENNZAHL_DETECTOR_BOTTOM;
					printf("Particle hit the proton detector (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);
					fprintf(LOGSCR,"Particle hit the proton detector (no reflection) at r=%LG z=%LG\n",y1[1],y1[3]);	
					return 1;	
				}				
				else if(prob < Transmission(Enormal*1e9,FPrealCsI,FPimCsI)){ // statistical absorption
					if (v[2] < 0)
						kennz = KENNZAHL_DETECTOR_TOP;
					else
						kennz = KENNZAHL_DETECTOR_BOTTOM;
					if (vabs > 0) gammaend = acos(v[2]/vabs)/conv;
					else gammaend=0;
					alphaend= atan2(v[1],v[0])/conv;			
					return 1;						
				}
				MaterialDiffProb = 0.5; // reflection on CsI, rough -> probability of diffuse reflection 0.5				
		}
		
		// specular reflexion
		prob = mt_get_double(v_mt_state);
		if ((diffuse == 1) || ((diffuse==3)&&(prob >= MaterialDiffProb)))
		{
	    	printf("\npol %d t=%LG Erefl=%LG neV r=%LG z=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3]);
			fprintf(LOGSCR,"\npol %d t=%LG Erefl=%LG neV r=%LG z=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3]);
			nrefl++;
			if(reflektlog == 1)
				fprintf(REFLECTLOG,"%LG %LG %LG %LG %LG %LG 1 %LG %LG %LG %LG %LG %LG %LG %LG %LG\n",
									x1,y1[1],y1[3],y1[5],p1[0],p1[1],vabs,H,Enormal*1e9,(long double)0.0,(long double)0.0,y1[2],y1[4],y1[1]*y1[6],y1[6]);
			v[0] -= 2*vnormal*normal[0]; // reflect velocity
			v[1] -= 2*vnormal*normal[1];
			v[2] -= 2*vnormal*normal[2];
		}
		
		// diffuse reflection
		else if ((diffuse == 2) || ((diffuse == 3)&&(prob < MaterialDiffProb)))
		{
			long double winkeben = 0.0 + (mt_get_double(v_mt_state)) * (2 * pi); // generate random reflection angles
			long double winksenkr = acos( cos(0.0) - mt_get_double(v_mt_state) * (cos(0.0) - cos(0.5 * pi)));
			if (vnormal > 0) winksenkr += pi; // if normal points out of volume rotate by 180 degrees
			v[0] = vabs*cos(winkeben)*sin(winksenkr);	// new velocity with respect to z-axis
			v[1] = vabs*sin(winkeben)*sin(winksenkr);
			v[2] = vabs*cos(winksenkr);

			long double cosalpha = normal[2], sinalpha = sqrt(1 - cosalpha*cosalpha);	// rotation angle (angle between z and n)
			if (sinalpha > 1e-30){ // when normal not parallel to z-axis rotate new velocity into the coordinate system where the normal is the z-axis
				long double a[3] = {-normal[1]/sinalpha, normal[0]/sinalpha, 0};	// rotation axis (z cross n)
				long double vtemp[3] = {v[0],v[1],v[2]};
				// rotate velocity vector
				v[0] = (cosalpha + a[0]*a[0]*(1 - cosalpha))*	vtemp[0] +  a[0]*a[1]*(1 - cosalpha)*				vtemp[1] + a[1]*sinalpha*	vtemp[2];
				v[1] =  a[1]*a[0]*(1 - cosalpha)*				vtemp[0] + (cosalpha + a[1]*a[1]*(1 - cosalpha))*	vtemp[1] - a[0]*sinalpha*	vtemp[2];
				v[2] = -a[1]*sinalpha*							vtemp[0] +  a[0]*sinalpha*							vtemp[1] + cosalpha*		vtemp[2];
			}
	       	printf("\npol %d t= %LG Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3],winkeben/conv,winksenkr/conv);
	       	fprintf(LOGSCR,"\npol %d t= %LG Erefl=%LG neV r=%LG z=%LG w_e=%LG w_s=%LG\n",polarisation,x1,Enormal*1e9,y1[1],y1[3],winkeben/conv,winksenkr/conv);
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
