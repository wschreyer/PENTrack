#include <cmath>

#include "racetrack.h"
#include "globals.h"

void Racetrack(long double r,long double phi,long double z,  long double lx, long double ly, long double I_rt, long double Bi[3], long double dBidrj[3][3]);
void CenterCurrent(long double r,long double phi,long double z,  long double I_rt, long double Bi[3], long double dBidrj[3][3]);
void StraightWireField(long double r,long double phi,long double z, long double I_rt, long double SW1r, long double SW1phi, long double SW1z, long double SW2r, long double SW2phi, long double SW2z, long double Bi[3], long double dBidrj[3][3]);

void StraightWireField_z(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SWx, const long double SWy, const long double SW1z, const long double SW2z,
						long double Bi[3], long double dBidrj[3][3]);
void StraightWireField_x(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SW1x, const long double SW2x, const long double SWz,
						long double Bi[3], long double dBidrj[3][3]);
void StraightWireField_y(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SW1y, const long double SW2y, const long double SWz,
						long double Bi[3], long double dBidrj[3][3]);
void FullRacetrack(long double r, long double phi, long double z, long double I_rt, long double SWz1, long double SWz2, long double SWr,
					long double Bi[3], long double dBidrj[3][3]);
void StraightWireField_Center(const long double r, const long double phi, const long double z, const long double I_rt,
							const long double SW1z, const long double SW2z, long double Bi[3], long double dBidrj[3][3]);

int Racetracks = 2;	// type of racetrack calculation (1: infinite wires, 2: finite wires, 3: central rod)
long double Ibar = 2250.; // current through rods

// define racetrack current bars
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the straight wire 
								//		current from outside in		from inside out			from high to low			from low to high
const long double Bars_1r[14]	= {0,	0.6, 0.6, 0.6, 0.6, 		0, 0, 0, 0, 			0.6, 0.6, 0.6, 0.6, 		0};
const long double Bars_1phi[14]	= {0,	0, pi/2, pi, 1.5*pi,		0, pi/2, pi, 1.5*pi,	0, pi/2, pi, 1.5*pi,		0};
const long double Bars_1z[14]	= {0,	-0.15, -0.15, -0.15, -0.15, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 1.35, 	-0.15};
const long double Bars_2r[14]	= {0,	0, 0, 0, 0,					0.6, 0.6, 0.6, 0.6,		0.6, 0.6, 0.6, 0.6,			0};
const long double Bars_2phi[14]	= {0,	0, pi/2, pi, 1.5*pi,		0, pi/2, pi, 1.5*pi,	0, pi/2, pi, 1.5*pi,		0};
const long double Bars_2z[14]	= {0,	-0.15, -0.15, -0.15, -0.15, 1.35, 1.35, 1.35, 1.35, -0.15, -0.15, -0.15, -0.15, 1.35};


void RacetrackField(long double rloc, long double philoc, long double zloc, long double Bi[3], long double dBidrj[3][3]){
	// add fields and field derivations of racetrack coils
	
	if(Ibar!=0)
	{
		if(Racetracks==1)
		{
			// old implementation of racetracks as 5 infinitely long vertical current bars, fast
			Racetrack(rloc, philoc, zloc, -0.424264, -0.424264, Ibar, Bi, dBidrj);
			Racetrack(rloc, philoc, zloc, 0.424264, -0.424264, Ibar, Bi, dBidrj); 
			Racetrack(rloc, philoc, zloc, -0.424264, 0.424264, Ibar, Bi, dBidrj);
			Racetrack(rloc, philoc, zloc, 0.424264, 0.424264, Ibar, Bi, dBidrj); 
			CenterCurrent(rloc, philoc, zloc, -4.0*Ibar, Bi, dBidrj);
		}
		else if(Racetracks==2)
		{			
			// new but slow implementation of racetracks including horizontal parts
			// calculate flux density at point r_current, phi_current, z_current
			// for all single current bars
/*
			for(int cobar=1;cobar<=12;cobar++)
			{
				StraightWireField(rloc,philoc,zloc, Ibar, Bars_1r[cobar],Bars_1phi[cobar], Bars_1z[cobar],Bars_2r[cobar], Bars_2phi[cobar], Bars_2z[cobar], Bi, dBidrj);					
			}
			// then for the center bar (4x current)   
			StraightWireField(rloc,philoc,zloc, 4*Ibar, Bars_1r[13],Bars_1phi[13], Bars_1z[13],Bars_2r[13], Bars_2phi[13], Bars_2z[13], Bi, dBidrj);
*/
/*
			StraightWireField_z(rloc,philoc,zloc,4*Ibar, 0, 0, -0.15, 1.35, Bi,dBidrj);

			StraightWireField_z(rloc,philoc,zloc,Ibar, 0.6, 0, 1.35, -0.15, Bi,dBidrj);
			StraightWireField_z(rloc,philoc,zloc,Ibar, 0, 0.6, 1.35, -0.15, Bi,dBidrj);
			StraightWireField_z(rloc,philoc,zloc,Ibar, -0.6, 0, 1.35, -0.15, Bi,dBidrj);
			StraightWireField_z(rloc,philoc,zloc,Ibar, 0, -0.6, 1.35, -0.15, Bi,dBidrj);
			
			StraightWireField_x(rloc,philoc,zloc,Ibar, 0.6, 0, -0.15, Bi,dBidrj);
			StraightWireField_x(rloc,philoc,zloc,Ibar, -0.6, 0, -0.15, Bi,dBidrj);
			StraightWireField_y(rloc,philoc,zloc,Ibar, 0.6, 0, -0.15, Bi,dBidrj);
			StraightWireField_y(rloc,philoc,zloc,Ibar, -0.6, 0, -0.15, Bi,dBidrj);

			StraightWireField_x(rloc,philoc,zloc,Ibar, 0, 0.6, 1.35, Bi,dBidrj);
			StraightWireField_x(rloc,philoc,zloc,Ibar, 0, -0.6, 1.35, Bi,dBidrj);
			StraightWireField_y(rloc,philoc,zloc,Ibar, 0, 0.6, 1.35, Bi,dBidrj);
			StraightWireField_y(rloc,philoc,zloc,Ibar, 0, -0.6, 1.35, Bi,dBidrj);
*/
			FullRacetrack(rloc, philoc, zloc, Ibar, -0.15, 1.35, 0.6, Bi, dBidrj);
		}
		else if (Racetracks == 3)
			StraightWireField_Center(rloc,philoc,zloc,4*Ibar, -1, 2, Bi,dBidrj);
		
	}  		
}	


// produce   B-Field of one straight wire in z direction an positions lx, ly and current I_rt
void Racetrack(long double r,long double phi,long double z,  long double lx, long double ly, long double I_rt, long double Bi[3], long double dBidrj[3][3]){
	// cartesian coordinates of neutron
	//long double x = cosl(phi) * r, y = sinl(phi) *r;
	// cartesian coordinates of racetracks	
	long double vorfaktor = mu0 * I_rt / (2 * pi);
	
	long double t2 = cos(phi);    
  long double t3 = ly*r*t2;
  long double     t5 = sin(phi);
  long double     t6 = lx*r*t5;
      long double t7 = t3-t6;
      long double t8 = vorfaktor*t7;
      long double t9 = r*r;
      long double t10 = t2*t2;
      long double t11 = t9*t10;
      long double t12 = r*t2;
      long double t13 = t12*lx;
      long double t15 = lx*lx;
      long double t16 = t5*t5;
      long double t17 = t9*t16;
      long double t18 = r*t5;
      long double t19 = t18*ly;
      long double t21 = ly*ly;
      long double t22 = t11-2.0*t13+t15+t17-2.0*t19+t21;
      long double t23 = 1/t22;
      long double t24 = 1/r;
      long double t25 = t23*t24;
      Bi[0] += t8*t25;
      long double t31 = t22*t22;
      long double t33 = 1/t31*t24;
      long double t34 = r*t10;
      long double t35 = t2*lx;
      long double t36 = r*t16;
      long double t37 = t5*ly;
      long double t42 = 1/t9;
      dBidrj[0][0] += vorfaktor*(ly*t2-lx*t5)*t25-t8*t33*(2.0*t34-2.0*t35+2.0*t36-2.0*t37)-t8*t23*t42;
      dBidrj[0][1] += vorfaktor*(-t19-t13)*t25+2.0*t8*t33*t7;
      long double t53 = vorfaktor*(-t17-t19-t11+t13);
      long double t54 = t12-lx;
      long double t55 = t54*t54;
      long double t56 = t18-ly;
      long double t57 = t56*t56;
      long double t58 = t55+t57;
      long double t59 = sqrt(t58);
      long double t60 = 1/t59;
      long double t61 = t60*t24;
      Bi[1] += t53*t61;
      long double t69 = 1/t59/t58*t24;
      dBidrj[1][0] += vorfaktor*(-2.0*t36-t37-2.0*t34+t35)*t61-t53*t69*(2.0*t54*t2+2.0*t56*t5)/2.0-t53*t60*t42;
      dBidrj[1][1] += vorfaktor*(-t3-t6)*t61-t53*t69*(-2.0*t54*r*t5+2.0*t56*r*t2)/2.0;

	return;
	
}

// produce center current rod field 
void CenterCurrent(long double r,long double phi,long double z, long double I_rt, long double Bi[3], long double dBidrj[3][3]){
	
	long double vorfaktor = mu0 * I_rt / (2 * pi);	
		
	// only Bphi is present
	Bi[1] += vorfaktor/r;	
	
	// dBphidr
	dBidrj[1][0] -= vorfaktor/(r*r);		
	
	return;	
}

// produce B-field at position r,phi,z of a straight wire with current I_rt 
// defined by two position vectors (SW1r,SW1phi,SW1z and SW2r,SW2phi,SW2z) lying on the ends of the straight wire 
// code generated by Maple ("FiniteWireDiffs.mw") with formula http://de.wikipedia.org/wiki/Biot-Savart#Gerader_Linienleiter
void StraightWireField(const long double r,const long double phi,const long double z, const long double I_rt, 
						const long double SW1r, const long double SW1phi, const long double SW1z, 
						const long double SW2r, const long double SW2phi, const long double SW2z,
						long double Bi[3], long double dBidrj[3][3])
{
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = cos(SW1phi);
	long double t2 = t1 * SW1r;
	long double t3 = cosl(phi);
	long double t4 = t3 * r;
	long double t5 = t2 - t4;
	long double t6 = t5 * t5;
	long double t7 = sin(SW1phi);
	long double t8 = t7 * SW1r;
	long double t9 = sin(phi);
	long double t10 = t9 * r;
	long double t11 = t8 - t10;
	long double t12 = t11 * t11;
	long double t13 = SW1z - z;
	long double t14 = t13 * t13;
	long double t15 = t6 + t12 + t14;
	long double t16 = sqrt(t15);
	long double t17 = 0.1e1 / t16;
	long double t18 = t17 * vorfaktor;
	long double t19 = 0.1e1 / t15;
	long double t20 = SW1r * SW1r;
	long double t21 = t1 * t1;
	long double t23 = cosl(SW2phi);
	long double t24 = t23 * SW2r;
	long double t28 = t7 * t7;
	long double t30 = sinl(SW2phi);
	long double t31 = t30 * SW2r;
	long double t37 = SW2z - SW1z;
	long double t39 = t21 * t20 - t1 * (t24 + t4) * SW1r + t28 * t20 - t7 * (t10 + t31) * SW1r + t24 * t4 + t31 * t10 - t37 * t13;
	long double t40 = t39 * t39;
	long double t42 = t24 - t2;
	long double t43 = t42 * t42;
	long double t44 = t31 - t8;
	long double t45 = t44 * t44;
	long double t46 = t37 * t37;
	long double t47 = t43 + t45 + t46;
	long double t48 = 0.1e1 / t47;
	long double t50 = 0.1e1 - t48 * t40 * t19;
	long double t51 = sqrtl(t50);
	long double t52 = 0.1e1 / t51;
	long double t53 = t52 * t18;
	long double t54 = t24 - t4;
	long double t55 = t54 * t54;
	long double t56 = t31 - t10;
	long double t57 = t56 * t56;
	long double t58 = SW2z - z;
	long double t59 = t58 * t58;
	long double t60 = t55 + t57 + t59;
	long double t61 = sqrtl(t60);
	long double t62 = 0.1e1 / t61;
	long double t63 = SW2r * SW2r;
	long double t64 = t23 * t23;
	long double t69 = t30 * t30;
	long double t77 = -t64 * t63 + t23 * (t2 + t4) * SW2r - t69 * t63 + t30 * (t8 + t10) * SW2r - t4 * t2 - t10 * t8 - t37 * t58;
	long double t79 = sqrtl(t47);
	long double t80 = 0.1e1 / t79;
	long double t84 = -t80 * t77 * t62 + t80 * t39 * t17;
	long double t87 = -t37 * t11 + t44 * t13;
	long double t88 = t87 * t87;
	long double t91 = -t42 * t13 + t37 * t5;
	long double t92 = t91 * t91;
	long double t95 = -t44 * t5 + t42 * t11;
	long double t96 = t95 * t95;
	long double t97 = t88 + t92 + t96;
	long double t98 = sqrtl(t97);
	long double t99 = 0.1e1 / t98;
	long double t100 = t99 * t84;
	long double t101 = t3 * t87;
	long double t103 = t101 * t100 * t53;
	long double t104 = t9 * t91;
	long double t106 = t104 * t100 * t53;
	long double t109 = 0.1e1 / t16 / t15;
	long double t110 = t109 * vorfaktor;
	long double t111 = t84 * t52;
	long double t112 = t111 * t110;
	long double t113 = t87 * t99;
	long double t116 = -t3 * t5 - t9 * t11;
	long double t117 = 0.2e1 * t116 * t3;
	long double t122 = 0.1e1 / t51 / t50;
	long double t124 = t84 * t122 * t18;
	long double t125 = t15 * t15;
	long double t127 = t40 / t125;
	long double t130 = t39 * t19;
	long double t139 = -t1 * t3 * SW1r - t7 * t9 * SW1r + t23 * SW2r * t3 + t30 * SW2r * t9;
	long double t143 = 0.2e1 * t116 * t48 * t127 - 0.2e1 * t139 * t48 * t130;
	long double t144 = t143 * t3;
	long double t150 = t77 / t61 / t60;
	long double t159 = t39 * t109;
	long double t166 = t99 * ((-t3 * t54 - t9 * t56) * t80 * t150 - t80 * t139 * t62 - t116 * t80 * t159 + t80 * t139 * t17);
	long double t169 = t111 * t18;
	long double t171 = 0.1e1 / t98 / t97;
	long double t172 = t87 * t171;
	long double t173 = t9 * t87;
	long double t175 = t3 * t91;
	long double t178 = t42 * t9;
	long double t179 = t44 * t3 - t178;
	long double t181 = t37 * t173 - t37 * t175 + t179 * t95;
	long double t182 = 0.2e1 * t181 * t3;
	long double t186 = t91 * t99;
	long double t187 = 0.2e1 * t116 * t9;
	long double t191 = t143 * t9;
	long double t197 = t91 * t171;
	long double t198 = 0.2e1 * t181 * t9;
	long double t207 = t9 * r * t5 - t3 * r * t11;
	long double t208 = 0.2e1 * t207 * t3;
	long double t214 = SW1r * r;
	long double t221 = t1 * t9 * t214 - t7 * t3 * t214 - t24 * t10 + t31 * t4;
	long double t225 = 0.2e1 * t207 * t48 * t127 - 0.2e1 * t221 * t48 * t130;
	long double t226 = t225 * t3;
	long double t246 = t99 * ((t9 * r * t54 - t3 * r * t56) * t80 * t150 - t80 * t221 * t62 - t207 * t80 * t159 + t80 * t221 * t17);
	long double t257 = -t44 * t10 - t42 * t4;
	long double t259 = t37 * t3 * r * t87 + t37 * t9 * r * t91 + t257 * t95;
	long double t260 = 0.2e1 * t259 * t3;
	long double t264 = r * t99;
	long double t265 = t3 * t3;
	long double t266 = t37 * t265;
	long double t270 = t173 * t100 * t53;
	long double t271 = 0.2e1 * t207 * t9;
	long double t275 = t225 * t9;
	long double t281 = 0.2e1 * t259 * t9;
	long double t285 = t9 * t9;
	long double t286 = t37 * t285;
	long double t290 = t175 * t100 * t53;
	long double t291 = -t208 * t113 * t112 / 0.2e1 - t226 * t113 * t124 / 0.2e1 + t101 * t246 * t53 - t260 * t172 * t169 / 0.2e1 + t266 * t264 * t169 - t270 - t271 * t186 * t112 / 0.2e1 - t275 * t186 * t124 / 0.2e1 + t104 * t246 * t53 - t281 * t197 * t169 / 0.2e1 + t286 * t264 * t169 + t290;
	long double t292 = -0.2e1 * t13 * t3;
	long double t301 = -0.2e1 * t13 * t48 * t127 - 0.2e1 * t37 * t48 * t130;
	long double t302 = t301 * t3;
	long double t317 = t99 * (-t58 * t80 * t150 - t80 * t37 * t62 + t13 * t80 * t159 + t80 * t37 * t17);
	long double t322 = -t44 * t87 + t42 * t91;
	long double t323 = 0.2e1 * t322 * t3;
	long double t330 = -0.2e1 * t13 * t9;
	long double t334 = t301 * t9;
	long double t340 = 0.2e1 * t322 * t9;
	long double t429 = t52 * t110;
	long double t434 = t122 * t18;
	long double t441 = t171 * t84;
	Bi[0] += 		t103 + t106;
	dBidrj[0][0] += -t117 * t113 * t112 / 0.2e1 - t144 * t113 * t124 / 0.2e1 + t101 * t166 * t53 - t182 * t172 * t169 / 0.2e1 - t187 * t186 * t112 / 0.2e1 - t191 * t186 * t124 / 0.2e1 + t104 * t166 * t53 - t198 * t197 * t169 / 0.2e1;
	dBidrj[0][1] += t291;
	dBidrj[0][2] += -t292 * t113 * t112 / 0.2e1 - t302 * t113 * t124 / 0.2e1 + t101 * t317 * t53 - t323 * t172 * t169 / 0.2e1 - t44 * t3 * t100 * t53 - t330 * t186 * t112 / 0.2e1 - t334 * t186 * t124 / 0.2e1 + t104 * t317 * t53 - t340 * t197 * t169 / 0.2e1 + t178 * t100 * t53;
	Bi[1] += 		-t270 + t290;
	dBidrj[1][0] += t187 * t113 * t112 / 0.2e1 + t191 * t113 * t124 / 0.2e1 - t173 * t166 * t53 + t198 * t172 * t169 / 0.2e1 - t286 * t100 * t53 - t117 * t186 * t112 / 0.2e1 - t144 * t186 * t124 / 0.2e1 + t175 * t166 * t53 - t182 * t197 * t169 / 0.2e1 - t266 * t100 * t53;
	dBidrj[1][1] += t271 * t113 * t112 / 0.2e1 + t275 * t113 * t124 / 0.2e1 - t173 * t246 * t53 + t281 * t172 * t169 / 0.2e1 - t103 - t208 * t186 * t112 / 0.2e1 - t226 * t186 * t124 / 0.2e1 + t175 * t246 * t53 - t260 * t197 * t169 / 0.2e1 - t106;
	dBidrj[1][2] += t330 * t113 * t112 / 0.2e1 + t334 * t113 * t124 / 0.2e1 - t173 * t317 * t53 + t340 * t172 * t169 / 0.2e1 + t9 * t44 * t100 * t53 - t292 * t186 * t112 / 0.2e1 - t302 * t186 * t124 / 0.2e1 + t175 * t317 * t53 - t323 * t197 * t169 / 0.2e1 + t3 * t42 * t100 * t53;
	Bi[2] += 		t95 * t100 * t53;
	dBidrj[2][0] += -t116 * t95 * t100 * t429 - t143 * t95 * t100 * t434 / 0.2e1 + t95 * t166 * t53 - t181 * t95 * t441 * t53 + t179 * t100 * t53;
	dBidrj[2][1] += -t207 * t95 * t100 * t429 - t225 * t95 * t100 * t434 / 0.2e1 + t95 * t246 * t53 - t259 * t95 * t441 * t53 + t257 * t100 * t53;
	dBidrj[2][2] += t13 * t95 * t100 * t429 - t301 * t95 * t100 * t434 / 0.2e1 + t95 * t317 * t53 - t322 * t95 * t441 * t53;

	return;	
}

// B field of a finite wire in z direction
void StraightWireField_z(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SWx, const long double SWy, const long double SW1z, const long double SW2z,
						long double Bi[3], long double dBidrj[3][3])
{
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = SWx * SWx;
	long double t2 = SWx * r;
	long double t3 = cos(phi);
	long double t5 = 2 * t3 * t2;
	long double t6 = r * r;
	long double t7 = t3 * t3;
	long double t8 = t7 * t6;
	long double t9 = SWy * SWy;
	long double t10 = r * SWy;
	long double t11 = sin(phi);
	long double t13 = 2 * t11 * t10;
	long double t14 = t11 * t11;
	long double t15 = t14 * t6;
	long double t16 = SW2z * SW2z;
	long double t19 = z * z;
	long double t20 = t1 - t5 + t8 + t9 - t13 + t15 + t16 - 2 * SW2z * z + t19;
	long double t21 = sqrt(t20);
	long double t22 = 0.1e1 / t21;
	long double t23 = -SW2z + z;
	long double t25 = SW1z * SW1z;
	long double t28 = t1 - t5 + t8 + t9 - t13 + t15 + t25 - 2 * SW1z * z + t19;
	long double t29 = sqrt(t28);
	long double t30 = 0.1e1 / t29;
	long double t31 = -SW1z + z;
	long double t33 = t23 * t22 - t31 * t30;
	long double t34 = fabs(t33);
	long double t35 = t34 * vorfaktor;
	long double t36 = -SW2z + SW1z;
	long double t37 = t36 * t35;
	long double t40 = -SWy * t3 + SWx * t11;
	long double t41 = t1 - t5 + t8 + t9 - t13 + t15;
	long double t42 = sqrt(t41);
	long double t43 = 0.1e1 / t42;
	long double t44 = t43 * t40;
	long double t45 = t36 * t36;
	long double t46 = t45 * t41;
	long double t47 = sqrt(t46);
	long double t48 = 0.1e1 / t47;
	long double t50 = t48 * t44 * t37;
	long double t51 = fabs(t33) / t33;
	long double t52 = t51 * vorfaktor;
	long double t55 = t23 / t21 / t20;
	long double t56 = SWy * t11;
	long double t58 = SWx * t3;
	long double t60 = -t56 + t14 * r - t58 + t7 * r;
	long double t64 = t31 / t29 / t28;
	long double t67 = (-2 * t60 * t55 + 2 * t60 * t64) * t52 / 2;
	long double t69 = t48 * t43;
	long double t70 = t69 * t40 * t36;
	long double t73 = 0.1e1 / t42 / t41;
	long double t74 = t73 * t40;
	long double t75 = 2 * t60 * t48;
	long double t80 = t45 * t36 * t35;
	long double t82 = 0.1e1 / t47 / t46;
	long double t83 = 2 * t60 * t82;
	long double t90 = t11 * t2 - t3 * t10;
	long double t94 = (-2 * t90 * t55 + 2 * t90 * t64) * t52 / 2;
	long double t100 = 2 * t90 * t48;
	long double t104 = 2 * t90 * t82;
	long double t114 = (-t23 * t55 + t22 + t31 * t64 - t30) * t52;
	long double t116 = t43 * t60;
	long double t120 = t69 * t60 * t36;
	long double t126 = t73 * t60;
	Bi[0] += t50;
	dBidrj[0][0] += t70 * t67 - t75 * t74 * t37 / 2 - t83 * t44 * t80 / 2;
	dBidrj[0][1] += t70 * t94 + t48 * t43 * (t56 + t58) * t37 - t100 * t74 * t37 / 2 - t104 * t44 * t80 / 2;
	dBidrj[0][2] += t70 * t114;
	Bi[1] += -t48 * t116 * t37;
	dBidrj[1][0] += -t120 * t67 - t48 * t43 * (t7 + t14) * t37 + t75 * t126 * t37 / 2 + t83 * t116 * t80 / 2;
	dBidrj[1][1] += -t120 * t94 - t50 + t100 * t126 * t37 / 2 + t104 * t116 * t80 / 2;
	dBidrj[1][2] += -t120 * t114;
}

// B field of a finite wire in x-z-plane and parallel to x axis
void StraightWireField_x(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SW1x, const long double SW2x, const long double SWz,
						long double Bi[3], long double dBidrj[3][3])
{
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = cos(phi);
	long double t2 = t1 * r;
	long double t3 = -SW2x + t2;
	long double t4 = r * r;
	long double t5 = t1 * t1;
	long double t6 = t5 * t4;
	long double t7 = SW2x * r;
	long double t10 = sin(phi);
	long double t11 = t10 * t10;
	long double t12 = t11 * t4;
	long double t13 = SW2x * SW2x;
	long double t14 = -SWz + z;
	long double t15 = t14 * t14;
	long double t16 = t6 - 2 * t1 * t7 + t12 + t13 + t15;
	long double t17 = sqrt(t16);
	long double t18 = 0.1e1 / t17;
	long double t20 = -SW1x + t2;
	long double t21 = SW1x * r;
	long double t24 = SW1x * SW1x;
	long double t25 = t6 - 2 * t1 * t21 + t12 + t24 + t15;
	long double t26 = sqrt(t25);
	long double t27 = 0.1e1 / t26;
	long double t29 = t18 * t3 - t27 * t20;
	long double t30 = fabs(t29);
	long double t31 = t30 * vorfaktor;
	long double t32 = t14 * t31;
	long double t33 = -SW2x + SW1x;
	long double t35 = t12 + t15;
	long double t36 = sqrt(t35);
	long double t37 = 0.1e1 / t36;
	long double t38 = t33 * t33;
	long double t39 = t38 * t35;
	long double t40 = sqrt(t39);
	long double t41 = 0.1e1 / t40;
	long double t42 = t41 * t37;
	long double t43 = t42 * t10 * t33;
	long double t44 = t43 * t32;
	long double t45 = fabs(t29) / t29;
	long double t46 = t45 * vorfaktor;
	long double t50 = 0.1e1 / t17 / t16 * t3;
	long double t51 = t5 * r;
	long double t53 = t11 * r;
	long double t60 = 0.1e1 / t26 / t25 * t20;
	long double t65 = t1 * t18 - (t51 - t1 * SW2x + t53) * t50 - t1 * t27 + (t51 - t1 * SW1x + t53) * t60;
	long double t67 = t14 * t65 * t46;
	long double t70 = t33 * t14 * t31;
	long double t71 = t11 * t10;
	long double t73 = 0.1e1 / t36 / t35;
	long double t75 = r * t41;
	long double t78 = t38 * t33;
	long double t80 = t78 * t14 * t31;
	long double t83 = 0.1e1 / t40 / t39;
	long double t84 = r * t83;
	long double t88 = t10 * r;
	long double t95 = -t18 * t88 - t10 * t7 * t50 + t27 * t88 + t10 * t21 * t60;
	long double t97 = t14 * t95 * t46;
	long double t100 = t42 * t1 * t33;
	long double t101 = t100 * t32;
	long double t103 = t4 * t41;
	long double t108 = t4 * t83;
	long double t115 = -2 * t14 * t50 + 2 * t14 * t60;
	long double t117 = t14 * t115 * t46 / 2;
	long double t119 = t33 * t31;
	long double t120 = t37 * t10;
	long double t122 = t41 * t120 * t119;
	long double t124 = 2 * t14 * t41;
	long double t128 = 2 * t14 * t83;
	long double t134 = t73 * t1;
	long double t138 = t37 * t1;
	long double t163 = r * t31;
	long double t168 = t4 * t31;
	long double t184 = t11 * t4 * r * t31;
	long double t185 = t73 * t33;
	long double t189 = t37 * t78;
	long double t197 = t88 * t31;
	Bi[0] += t44;
	dBidrj[0][0] += t43 * t67 - t75 * t73 * t71 * t70 - t84 * t37 * t71 * t80;
	dBidrj[0][1] += t43 * t97 + t101 - t1 * t103 * t73 * t11 * t70 - t1 * t108 * t37 * t11 * t80;
	dBidrj[0][2] += t43 * t117 + t122 - t124 * t73 * t10 * t70 / 2 - t128 * t120 * t80 / 2;
	Bi[1] += t101;
	dBidrj[1][0] += t100 * t67 - t11 * t75 * t134 * t70 - t11 * t84 * t138 * t80;
	dBidrj[1][1] += t100 * t97 - t44 - t10 * t103 * t73 * t5 * t70 - t10 * t108 * t37 * t5 * t80;
	dBidrj[1][2] += t100 * t117 + t41 * t138 * t119 - t124 * t134 * t70 / 2 - t128 * t138 * t80 / 2;
	Bi[2] += -t43 * t163;
	dBidrj[2][0] += -t43 * r * t65 * t46 - t122 + t41 * t73 * t33 * t71 * t168 + t83 * t37 * t78 * t71 * t168;
	dBidrj[2][1] += -t43 * r * t95 * t46 - t100 * t163 + t1 * t41 * t185 * t184 + t1 * t83 * t189 * t184;
	dBidrj[2][2] += -t43 * r * t115 * t46 / 2 + t124 * t185 * t197 / 2 + t128 * t189 * t197 / 2;
}

// B field of a finite wire in y-z-plane and parallel to y axis
void StraightWireField_y(const long double r, const long double phi, const long double z, const long double I_rt,
						const long double SW1y, const long double SW2y, const long double SWz,
						long double Bi[3], long double dBidrj[3][3])
{
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = r * r;
	long double t2 = cos(phi);
	long double t3 = t2 * t2;
	long double t4 = t3 * t1;
	long double t5 = -SWz + z;
	long double t6 = t5 * t5;
	long double t7 = t4 + t6;
	long double t8 = sqrt(t7);
	long double t10 = 0.1e1 / t8 * vorfaktor;
	long double t11 = sin(phi);
	long double t12 = t11 * r;
	long double t13 = -SW2y + t12;
	long double t14 = t11 * t11;
	long double t15 = t14 * t1;
	long double t16 = r * SW2y;
	long double t19 = SW2y * SW2y;
	long double t20 = t15 - 2 * t11 * t16 + t4 + t19 + t6;
	long double t21 = sqrt(t20);
	long double t22 = 0.1e1 / t21;
	long double t24 = -SW1y + t12;
	long double t25 = r * SW1y;
	long double t28 = SW1y * SW1y;
	long double t29 = t15 - 2 * t11 * t25 + t4 + t28 + t6;
	long double t30 = sqrt(t29);
	long double t31 = 0.1e1 / t30;
	long double t33 = t22 * t13 - t31 * t24;
	long double t34 = fabs(t33);
	long double t35 = t34 * t10;
	long double t36 = -SW2y + SW1y;
	long double t37 = -t36 * t5;
	long double t38 = t36 * t36;
	long double t39 = t7 * t38;
	long double t40 = sqrt(t39);
	long double t41 = 0.1e1 / t40;
	long double t42 = t2 * t41;
	long double t43 = t42 * t37;
	long double t47 = 0.1e1 / t8 / t7 * vorfaktor;
	long double t48 = -t5 * t34;
	long double t49 = t48 * t47;
	long double t50 = t41 * t36;
	long double t51 = t3 * t2;
	long double t52 = r * t51;
	long double t55 = fabs(t33) / t33;
	long double t59 = 0.1e1 / t21 / t20 * t13;
	long double t60 = t14 * r;
	long double t62 = t3 * r;
	long double t69 = 0.1e1 / t30 / t29 * t24;
	long double t76 = (t11 * t22 - (t60 - t11 * SW2y + t62) * t59 - t11 * t31 + (t60 - t11 * SW1y + t62) * t69) * t55 * t10;
	long double t78 = t48 * t10;
	long double t79 = t38 * t36;
	long double t81 = 0.1e1 / t40 / t39;
	long double t82 = t81 * t79;
	long double t86 = t11 * t4;
	long double t89 = t2 * r;
	long double t98 = (t22 * t89 + t2 * t16 * t59 - t31 * t89 - t2 * t25 * t69) * t55 * t10;
	long double t102 = t11 * t41;
	long double t106 = 2 * t5 * t2;
	long double t107 = t106 * t50;
	long double t114 = (-2 * t5 * t59 + 2 * t5 * t69) * t55 * t10 / 2;
	long double t117 = t2 * t50 * t35;
	long double t118 = t106 * t82;
	long double t122 = t36 * t5;
	long double t123 = t102 * t122;
	long double t125 = t5 * t34;
	long double t126 = t125 * t47;
	long double t127 = t3 * t12;
	long double t131 = t125 * t10;
	long double t135 = t2 * t15;
	long double t144 = 2 * t5 * t11;
	long double t155 = t50 * t89;
	long double t158 = t51 * t1;
	long double t166 = t1 * r * t34;
	long double t180 = r * t34;
	Bi[0] += t43 * t35;
	dBidrj[0][0] += -t52 * t50 * t49 + t43 * t76 - t52 * t82 * t78;
	dBidrj[0][1] += t86 * t50 * t49 + t43 * t98 + t86 * t82 * t78 - t102 * t37 * t35;
	dBidrj[0][2] += -t107 * t49 / 2 + t43 * t114 - t117 - t118 * t78 / 2;
	Bi[1] += t123 * t35;
	dBidrj[1][0] += -t127 * t50 * t126 + t123 * t76 - t127 * t82 * t131;
	dBidrj[1][1] += t135 * t50 * t126 + t123 * t98 + t135 * t82 * t131 + t42 * t122 * t35;
	dBidrj[1][2] += -t144 * t50 * t126 / 2 + t123 * t114 + t11 * t50 * t35 - t144 * t82 * t131 / 2;
	Bi[2] += t155 * t35;
	dBidrj[2][0] += -t50 * t158 * t34 * t47 + t155 * t76 + t117 - t82 * t158 * t35;
	dBidrj[2][1] += t102 * t36 * t3 * t166 * t47 + t155 * t98 - t50 * t12 * t35 + t11 * t81 * t79 * t3 * t166 * t10;
	dBidrj[2][2] += -t107 * t180 * t47 / 2 + t155 * t114 - t118 * t180 * t10 / 2;
}

// B field of four race track coils with bottom SWz1, top SWz2 in x,-x,y,-y directions with length SWr
void FullRacetrack(long double r, long double phi, long double z, long double I_rt, long double SWz1, long double SWz2, long double SWr,
					long double Bi[3], long double dBidrj[3][3]){
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = r * r;
	long double t2 = cos(phi);
	long double t3 = t2 * t2;
	long double t4 = t3 * t1;
	long double t5 = sin(phi);
	long double t6 = t5 * t5;
	long double t7 = t6 * t1;
	long double t8 = SWz2 * SWz2;
	long double t10 = 2 * SWz2 * z;
	long double t11 = z * z;
	long double t12 = t4 + t7 + t8 - t10 + t11;
	long double t13 = sqrt(t12);
	long double t14 = 0.1e1 / t13;
	long double t15 = SWz2 - z;
	long double t17 = -SWz1 + z;
	long double t18 = SWz1 * SWz1;
	long double t20 = 2 * SWz1 * z;
	long double t21 = t4 + t7 + t18 - t20 + t11;
	long double t22 = sqrt(t21);
	long double t23 = 0.1e1 / t22;
	long double t25 = -t15 * t14 - t23 * t17;
	long double t26 = fabs(t25);
	long double t27 = t26 * vorfaktor;
	long double t28 = r * t27;
	long double t29 = -SWz2 + SWz1;
	long double t30 = t29 * t5;
	long double t31 = t7 + t4;
	long double t32 = sqrt(t31);
	long double t33 = 0.1e1 / t32;
	long double t34 = t29 * t29;
	long double t35 = t34 * t31;
	long double t36 = sqrt(t35);
	long double t37 = 0.1e1 / t36;
	long double t38 = t37 * t33;
	long double t39 = t38 * t30;
	long double t42 = SWr * SWr;
	long double t43 = SWr * r;
	long double t44 = t2 * t43;
	long double t45 = 2 * t44;
	long double t46 = t7 + t42 - t45 + t4;
	long double t47 = sqrt(t46);
	long double t49 = vorfaktor / t47;
	long double t50 = t42 - t45 + t4 + t7 + t8 - t10 + t11;
	long double t51 = sqrt(t50);
	long double t52 = 0.1e1 / t51;
	long double t54 = t42 - t45 + t4 + t7 + t18 - t20 + t11;
	long double t55 = sqrt(t54);
	long double t56 = 0.1e1 / t55;
	long double t58 = -t15 * t52 - t17 * t56;
	long double t59 = fabs(t58);
	long double t60 = t59 * t49;
	long double t61 = t46 * t34;
	long double t62 = sqrt(t61);
	long double t63 = 0.1e1 / t62;
	long double t64 = r * t63;
	long double t65 = t30 * t64;
	long double t66 = t65 * t60;
	long double t67 = t7 + t42 + t45 + t4;
	long double t68 = sqrt(t67);
	long double t70 = vorfaktor / t68;
	long double t71 = t42 + t45 + t4 + t7 + t8 - t10 + t11;
	long double t72 = sqrt(t71);
	long double t73 = 0.1e1 / t72;
	long double t75 = t42 + t45 + t4 + t7 + t18 - t20 + t11;
	long double t76 = sqrt(t75);
	long double t77 = 0.1e1 / t76;
	long double t79 = -t15 * t73 - t17 * t77;
	long double t80 = fabs(t79);
	long double t81 = t80 * t70;
	long double t82 = t67 * t34;
	long double t83 = sqrt(t82);
	long double t84 = 0.1e1 / t83;
	long double t85 = r * t84;
	long double t86 = t30 * t85;
	long double t87 = t86 * t81;
	long double t88 = t8 - t10 + t11 + t4;
	long double t89 = sqrt(t88);
	long double t91 = vorfaktor / t89;
	long double t92 = t5 * t43;
	long double t93 = 2 * t92;
	long double t94 = t4 + t42 - t93 + t7 + t8 - t10 + t11;
	long double t95 = sqrt(t94);
	long double t96 = 0.1e1 / t95;
	long double t97 = t5 * r;
	long double t98 = SWr - t97;
	long double t100 = r * t14;
	long double t101 = t5 * t100;
	long double t102 = t98 * t96 + t101;
	long double t103 = fabs(t102);
	long double t104 = t103 * t91;
	long double t105 = t88 * t42;
	long double t106 = sqrt(t105);
	long double t107 = 0.1e1 / t106;
	long double t108 = -t15 * t107;
	long double t111 = t42 - t93 + t7 + t4;
	long double t112 = sqrt(t111);
	long double t114 = vorfaktor / t112;
	long double t116 = t4 + t42 - t93 + t7 + t18 - t20 + t11;
	long double t117 = sqrt(t116);
	long double t118 = 0.1e1 / t117;
	long double t120 = -t15 * t96 - t118 * t17;
	long double t121 = fabs(t120);
	long double t122 = t121 * t114;
	long double t123 = t111 * t34;
	long double t124 = sqrt(t123);
	long double t125 = 0.1e1 / t124;
	long double t126 = t98 * t125;
	long double t129 = t18 - t20 + t11 + t4;
	long double t130 = sqrt(t129);
	long double t132 = vorfaktor / t130;
	long double t134 = r * t23;
	long double t135 = t5 * t134;
	long double t136 = t98 * t118 + t135;
	long double t137 = fabs(t136);
	long double t138 = t137 * t132;
	long double t139 = t129 * t42;
	long double t140 = sqrt(t139);
	long double t141 = 0.1e1 / t140;
	long double t142 = -t17 * t141;
	long double t145 = t4 + t42 + t93 + t7 + t8 - t10 + t11;
	long double t146 = sqrt(t145);
	long double t147 = 0.1e1 / t146;
	long double t148 = -SWr - t97;
	long double t150 = 1;//csgn(SWr);
	long double t152 = t150 * t5;
	long double t153 = t152 * t100;
	long double t154 = -t150 * t148 * t147 - t153;
	long double t155 = fabs(t154);
	long double t156 = t155 * t91;
	long double t157 = t15 * t107;
	long double t160 = t42 + t93 + t7 + t4;
	long double t161 = sqrt(t160);
	long double t163 = vorfaktor / t161;
	long double t165 = t4 + t42 + t93 + t7 + t18 - t20 + t11;
	long double t166 = sqrt(t165);
	long double t167 = 0.1e1 / t166;
	long double t169 = -t15 * t147 - t167 * t17;
	long double t170 = fabs(t169);
	long double t171 = t170 * t163;
	long double t172 = t160 * t34;
	long double t173 = sqrt(t172);
	long double t174 = 0.1e1 / t173;
	long double t175 = t148 * t174;
	long double t178 = t152 * t134;
	long double t181 = -t178 - t150 * t148 * t167;
	long double t182 = fabs(t181);
	long double t183 = t182 * t132;
	long double t184 = t17 * t141;
	long double t187 = 4 * t39 * t28 - t66 - t87 + SWr * t108 * t104 + t29 * t126 * t122 + SWr * t142 * t138 + SWr * t157 * t156 + t29 * t175 * t171 + SWr * t184 * t183;
	long double t188 = t2 * t187;
	long double t189 = vorfaktor * t33;
	long double t190 = t26 * t189;
	long double t191 = r * t37;
	long double t193 = -t29 * t2 * t191;
	long double t196 = t8 - t10 + t11 + t7;
	long double t197 = sqrt(t196);
	long double t199 = vorfaktor / t197;
	long double t200 = t2 * r;
	long double t201 = SWr - t200;
	long double t203 = t2 * t100;
	long double t204 = t201 * t52 + t203;
	long double t205 = fabs(t204);
	long double t206 = t205 * t199;
	long double t207 = t196 * t42;
	long double t208 = sqrt(t207);
	long double t209 = 0.1e1 / t208;
	long double t210 = t15 * t209;
	long double t213 = -t201 * t63;
	long double t216 = t18 - t20 + t11 + t7;
	long double t217 = sqrt(t216);
	long double t219 = vorfaktor / t217;
	long double t221 = t2 * t134;
	long double t222 = t201 * t56 + t221;
	long double t223 = fabs(t222);
	long double t224 = t223 * t219;
	long double t225 = t216 * t42;
	long double t226 = sqrt(t225);
	long double t227 = 0.1e1 / t226;
	long double t228 = t17 * t227;
	long double t231 = -SWr - t200;
	long double t234 = t150 * t2;
	long double t235 = t234 * t100;
	long double t236 = -t150 * t231 * t73 - t235;
	long double t237 = fabs(t236);
	long double t238 = t237 * t199;
	long double t239 = -t15 * t209;
	long double t242 = -t29 * t231;
	long double t245 = t234 * t134;
	long double t248 = -t245 - t150 * t231 * t77;
	long double t249 = fabs(t248);
	long double t250 = t249 * t219;
	long double t251 = -t17 * t227;
	long double t254 = r * t125;
	long double t255 = t29 * t2;
	long double t256 = t255 * t254;
	long double t257 = t256 * t122;
	long double t258 = r * t174;
	long double t259 = t255 * t258;
	long double t260 = t259 * t171;
	long double t261 = 4 * t193 * t190 + SWr * t210 * t206 + t29 * t213 * t60 + SWr * t228 * t224 + SWr * t239 * t238 + t84 * t242 * t81 + SWr * t251 * t250 + t257 + t260;
	long double t262 = t5 * t261;
	long double t264 = t5 * t187;
	long double t265 = t2 * t261;
	long double t267 = SWr * t209;
	long double t268 = t97 * t267;
	long double t270 = SWr * t227;
	long double t271 = t97 * t270;
	long double t275 = SWr * t107;
	long double t276 = t200 * t275;
	long double t278 = SWr * t141;
	long double t279 = t200 * t278;
	long double t285 = 0.1e1 / t83 / t82;
	long double t287 = t285 * t80 * t70;
	long double t288 = t34 * t29;
	long double t289 = t6 * r;
	long double t290 = t2 * SWr;
	long double t291 = t3 * r;
	long double t292 = t289 + t290 + t291;
	long double t299 = vorfaktor / t130 / t129;
	long double t300 = t141 * t137;
	long double t301 = t300 * t299;
	long double t302 = -SWr * t17;
	long double t306 = 0.1e1 / t140 / t139;
	long double t308 = t306 * t137 * t132;
	long double t309 = t42 * SWr;
	long double t310 = -t309 * t17;
	long double t315 = vorfaktor / t89 / t88;
	long double t316 = t107 * t155;
	long double t317 = t316 * t315;
	long double t318 = SWr * t15;
	long double t322 = 0.1e1 / t106 / t105;
	long double t324 = t322 * t155 * t91;
	long double t325 = t309 * t15;
	long double t328 = t141 * t182;
	long double t329 = t328 * t299;
	long double t330 = SWr * t17;
	long double t334 = t306 * t182 * t132;
	long double t335 = t309 * t17;
	long double t338 = t107 * t103;
	long double t339 = t338 * t315;
	long double t340 = -SWr * t15;
	long double t344 = t322 * t103 * t91;
	long double t345 = -t309 * t15;
	long double t348 = fabs(t79) / t79;
	long double t350 = 0.1e1 / t72 / t71;
	long double t351 = t15 * t350;
	long double t354 = 0.1e1 / t76 / t75;
	long double t355 = t17 * t354;
	long double t357 = 2 * t292 * t351 + 2 * t292 * t355;
	long double t362 = 0.1e1 / t173 / t172;
	long double t364 = t5 * SWr;
	long double t365 = t291 + t364 + t289;
	long double t366 = 2 * t365 * t288;
	long double t370 = fabs(t181) / t181;
	long double t371 = t370 * t132;
	long double t373 = 0.1e1 / t22 / t21;
	long double t374 = r * t373;
	long double t375 = t291 + t289;
	long double t376 = 2 * t375 * t152;
	long double t379 = t5 * t23;
	long double t382 = 0.1e1 / t166 / t165;
	long double t383 = t148 * t382;
	long double t384 = 2 * t365 * t150;
	long double t389 = t376 * t374 / 2 - t150 * t379 + t384 * t383 / 2 + t150 * t5 * t167;
	long double t393 = fabs(t120) / t120;
	long double t394 = t393 * t114;
	long double t396 = 0.1e1 / t95 / t94;
	long double t397 = t15 * t396;
	long double t398 = t291 - t364 + t289;
	long double t401 = 0.1e1 / t117 / t116;
	long double t402 = t401 * t17;
	long double t404 = 2 * t398 * t397 + 2 * t398 * t402;
	long double t406 = t29 * t98;
	long double t409 = fabs(t136) / t136;
	long double t410 = t409 * t132;
	long double t411 = t98 * t401;
	long double t415 = 2 * t375 * t5;
	long double t418 = -t398 * t411 - t5 * t118 - t415 * t374 / 2 + t379;
	long double t422 = fabs(t102) / t102;
	long double t423 = t422 * t91;
	long double t424 = t98 * t396;
	long double t429 = 0.1e1 / t13 / t12;
	long double t430 = r * t429;
	long double t433 = t5 * t14;
	long double t434 = -t398 * t424 - t5 * t96 - t415 * t430 / 2 + t433;
	long double t440 = vorfaktor / t112 / t111;
	long double t442 = 2 * t398 * t29;
	long double t446 = t292 * t288 * t97 * t287 - t291 * t302 * t301 - t291 * t310 * t308 - t291 * t318 * t317 - t291 * t325 * t324 - t291 * t330 * t329 - t291 * t335 * t334 - t291 * t340 * t339 - t291 * t345 * t344 - t86 * t357 * t348 * t70 / 2 - t366 * t148 * t362 * t171 / 2 + t330 * t141 * t389 * t371 + t406 * t125 * t404 * t394 / 2 + t302 * t141 * t418 * t410 + t340 * t107 * t434 * t423 - t442 * t126 * t121 * t440 / 2;
	long double t447 = fabs(t154) / t154;
	long double t448 = t447 * t91;
	long double t450 = 0.1e1 / t146 / t145;
	long double t451 = t148 * t450;
	long double t459 = t384 * t451 / 2 + t150 * t5 * t147 + t376 * t430 / 2 - t150 * t433;
	long double t465 = vorfaktor / t161 / t160;
	long double t467 = 2 * t365 * t29;
	long double t471 = fabs(t169) / t169;
	long double t472 = t471 * t163;
	long double t473 = t15 * t450;
	long double t475 = t382 * t17;
	long double t477 = 2 * t365 * t473 + 2 * t365 * t475;
	long double t479 = t29 * t148;
	long double t483 = 0.1e1 / t124 / t123;
	long double t485 = 2 * t398 * t288;
	long double t506 = t97 * t27;
	long double t508 = 0.1e1 / t32 / t31;
	long double t514 = fabs(t25) / t25;
	long double t515 = t514 * vorfaktor;
	long double t516 = t15 * t429;
	long double t518 = t373 * t17;
	long double t520 = 2 * t375 * t516 + 2 * t375 * t518;
	long double t527 = 0.1e1 / t36 / t35;
	long double t534 = vorfaktor / t47 / t46;
	long double t536 = t63 * t59 * t534;
	long double t537 = t289 - t290 + t291;
	long double t538 = 2 * t537 * t29;
	long double t542 = fabs(t58) / t58;
	long double t544 = 0.1e1 / t51 / t50;
	long double t545 = t15 * t544;
	long double t548 = 0.1e1 / t55 / t54;
	long double t549 = t17 * t548;
	long double t551 = 2 * t537 * t545 + 2 * t537 * t549;
	long double t556 = 0.1e1 / t62 / t61;
	long double t558 = t556 * t59 * t49;
	long double t559 = 2 * t537 * t288;
	long double t565 = vorfaktor / t68 / t67;
	long double t567 = t84 * t80 * t565;
	long double t572 = t318 * t107 * t459 * t448 - t467 * t175 * t170 * t465 / 2 + t479 * t174 * t477 * t472 / 2 - t485 * t98 * t483 * t122 / 2 - t29 * t5 * t63 * t60 - t29 * t5 * t84 * t81 + 4 * t37 * t33 * t29 * t5 * t27 - t29 * t5 * t125 * t122 - t29 * t5 * t174 * t171 - 4 * t375 * t37 * t508 * t29 * t506 + 2 * t39 * r * t520 * t515 - 4 * t375 * t527 * t33 * t288 * t506 + t538 * t97 * t536 / 2 - t65 * t551 * t542 * t49 / 2 + t559 * t97 * t558 / 2 + t292 * t29 * t97 * t567;
	long double t573 = t446 + t572;
	long double t580 = fabs(t248) / t248;
	long double t581 = t580 * t219;
	long double t582 = 2 * t375 * t234;
	long double t585 = t2 * t23;
	long double t587 = t231 * t354;
	long double t588 = 2 * t292 * t150;
	long double t593 = t582 * t374 / 2 - t150 * t585 + t588 * t587 / 2 + t150 * t2 * t77;
	long double t597 = fabs(t204) / t204;
	long double t598 = t597 * t199;
	long double t599 = t201 * t544;
	long double t603 = 2 * t375 * t2;
	long double t606 = t2 * t14;
	long double t607 = -t537 * t599 - t2 * t52 - t603 * t430 / 2 + t606;
	long double t615 = t542 * t49;
	long double t617 = -t29 * t201;
	long double t620 = fabs(t222) / t222;
	long double t621 = t620 * t219;
	long double t622 = t201 * t548;
	long double t628 = -t537 * t622 - t2 * t56 - t603 * t374 / 2 + t585;
	long double t632 = fabs(t236) / t236;
	long double t633 = t632 * t199;
	long double t634 = t231 * t350;
	long double t642 = t588 * t634 / 2 + t150 * t2 * t73 + t582 * t430 / 2 - t150 * t606;
	long double t651 = t348 * t70;
	long double t653 = t84 * t29;
	long double t677 = vorfaktor / t217 / t216;
	long double t678 = t227 * t249;
	long double t679 = t678 * t677;
	long double t682 = t292 * t285 * t288 * t231 * t81 + t302 * t227 * t593 * t581 + t318 * t209 * t607 * t598 - t538 * t213 * t59 * t534 / 2 + t617 * t63 * t551 * t615 / 2 + t330 * t227 * t628 * t621 + t340 * t209 * t642 * t633 - t292 * t84 * t242 * t80 * t565 - t653 * t231 * t357 * t651 / 2 + t559 * t201 * t556 * t60 / 2 + t29 * t2 * t63 * t60 - 4 * t29 * t2 * t37 * t190 + t84 * t255 * t81 + t29 * t2 * t125 * t122 + t29 * t2 * t174 * t171 - t289 * t302 * t679;
	long double t684 = 0.1e1 / t226 / t225;
	long double t686 = t684 * t249 * t219;
	long double t690 = t125 * t121 * t440;
	long double t698 = t483 * t121 * t114;
	long double t703 = t174 * t170 * t465;
	long double t711 = t362 * t170 * t163;
	long double t717 = vorfaktor / t197 / t196;
	long double t718 = t209 * t237;
	long double t719 = t718 * t717;
	long double t723 = 0.1e1 / t208 / t207;
	long double t725 = t723 * t237 * t199;
	long double t731 = -2 * t375 * t29;
	long double t745 = t227 * t223;
	long double t746 = t745 * t677;
	long double t750 = t684 * t223 * t219;
	long double t753 = t209 * t205;
	long double t754 = t753 * t717;
	long double t758 = t723 * t205 * t199;
	long double t761 = -t289 * t310 * t686 - t442 * t200 * t690 / 2 + t256 * t404 * t393 * t114 / 2 - t485 * t200 * t698 / 2 - t467 * t200 * t703 / 2 + t259 * t477 * t471 * t163 / 2 - t366 * t200 * t711 / 2 - t289 * t340 * t719 - t289 * t345 * t725 - 2 * t731 * t200 * t37 * t26 * vorfaktor * t508 + 2 * t193 * t520 * t514 * t189 - 2 * t34 * t731 * t200 * t527 * t26 * t189 - t289 * t330 * t746 - t289 * t335 * t750 - t289 * t318 * t754 - t289 * t325 * t758;
	long double t762 = t682 + t761;
	long double t769 = t5 * t2 * t1;
	long double t770 = t769 * t302;
	long double t772 = t769 * t310;
	long double t774 = t769 * t318;
	long double t776 = t769 * t325;
	long double t783 = t769 * t330;
	long double t785 = t769 * t335;
	long double t794 = t44 * t288 * t98 * t698 + t770 * t301 + t772 * t308 + t774 * t317 + t776 * t324 - t44 * t479 * t703 - t44 * t288 * t148 * t711 + t783 * t329 + t785 * t334 - t257 - t260 + 4 * t38 * t255 * t28 - t255 * t64 * t60 - t255 * t85 * t81;
	long double t798 = t44 * t424 - t2 * r * t96 + t203;
	long double t804 = -t44 * t397 - t44 * t402;
	long double t811 = t44 * t411 - t2 * r * t118 + t221;
	long double t819 = t44 * t150 * t451 + t234 * r * t147 - t235;
	long double t825 = t44 * t473 + t44 * t475;
	long double t833 = -t245 + t44 * t150 * t383 + t234 * r * t167;
	long double t837 = SWr * t29;
	long double t838 = t837 * t7;
	long double t842 = t92 * t545 + t92 * t549;
	long double t846 = SWr * t288;
	long double t847 = t846 * t7;
	long double t852 = -t92 * t351 - t92 * t355;
	long double t857 = t769 * t340;
	long double t859 = t769 * t345;
	long double t863 = t340 * t107 * t798 * t423 + t406 * t125 * t804 * t394 + t302 * t141 * t811 * t410 + t318 * t107 * t819 * t448 + t479 * t174 * t825 * t472 + t330 * t141 * t833 * t371 + t838 * t536 - t65 * t842 * t542 * t49 + t847 * t558 - t838 * t567 - t86 * t852 * t348 * t70 - t847 * t287 + t857 * t339 + t859 * t344 + t44 * t406 * t690;
	long double t864 = t794 + t863;
	long double t868 = -t231 * t80;
	long double t878 = t837 * t4;
	long double t883 = t846 * t4;
	long double t892 = -t857 * t719 - t859 * t725 + t92 * t653 * t868 * t565 + t92 * t285 * t288 * t868 * t70 - t770 * t679 - t772 * t686 + t878 * t690 + t256 * t804 * t393 * t114 + t883 * t698 - t878 * t703 + t259 * t825 * t471 * t163 - t883 * t711 - t774 * t754 - t776 * t758;
	long double t907 = -t92 * t599 + t5 * r * t52 - t101;
	long double t917 = -t92 * t622 + t5 * r * t56 - t135;
	long double t925 = -t92 * t150 * t634 - t152 * r * t73 + t153;
	long double t936 = t178 - t92 * t150 * t587 - t152 * r * t77;
	long double t944 = -t92 * t617 * t536 + t92 * t288 * t201 * t558 - t783 * t746 - t785 * t750 - t66 - t87 + 4 * t29 * t5 * t191 * t190 + t318 * t209 * t907 * t598 + t617 * t63 * t842 * t615 + t330 * t227 * t917 * t621 + t340 * t209 * t925 * t633 - t653 * t231 * t852 * t651 + t302 * t227 * t936 * t581 - t30 * t254 * t122 - t30 * t258 * t171;
	long double t945 = t892 + t944;
	long double t952 = -t15 * t516 + t14 - t23 + t17 * t518;
	long double t961 = -t15 * t545 + t52 + t17 * t549 - t56;
	long double t969 = -t15 * t351 + t73 + t17 * t355 - t77;
	long double t973 = t103 * t315;
	long double t974 = -2 * SWr * t15;
	long double t979 = -2 * t15 * t5;
	long double t981 = 2 * t15 * t424 - t979 * t430;
	long double t986 = -2 * t309 * t15;
	long double t996 = -t15 * t397 + t96 - t118 + t17 * t402;
	long double t1000 = t137 * t299;
	long double t1001 = 2 * SWr * t17;
	long double t1006 = 2 * t17 * t5;
	long double t1008 = -2 * t17 * t411 - t1006 * t374;
	long double t1014 = 2 * t309 * t17;
	long double t1020 = t155 * t315;
	long double t1024 = -2 * t15 * t150;
	long double t1028 = t1024 * t451 - 2 * t15 * t152 * t430;
	long double t1042 = -t15 * t473 + t147 - t167 + t17 * t475;
	long double t1046 = t182 * t299;
	long double t1052 = 2 * t17 * t150;
	long double t1054 = 2 * t17 * t152 * t374 + t1052 * t383;
	long double t1064 = t1014 * t17 * t306 * t138 / 2 - SWr * t300 * t132 - t974 * t157 * t1020 / 2 + t318 * t107 * t1028 * t448 / 2 - t986 * t15 * t322 * t156 / 2 - SWr * t316 * t91 + t479 * t174 * t1042 * t472 - t1001 * t184 * t1046 / 2 + t330 * t141 * t1054 * t371 / 2 - t1014 * t17 * t306 * t183 / 2 + SWr * t328 * t132;
	long double t1065 = 4 * t39 * r * t952 * t515 - t65 * t961 * t542 * t49 - t86 * t969 * t348 * t70 - t974 * t108 * t973 / 2 + t340 * t107 * t981 * t423 / 2 + t986 * t15 * t322 * t104 / 2 + SWr * t338 * t91 + t406 * t125 * t996 * t394 - t1001 * t142 * t1000 / 2 + t302 * t141 * t1008 * t410 / 2 + t1064;
	long double t1071 = t205 * t717;
	long double t1076 = -2 * t15 * t2;
	long double t1078 = 2 * t15 * t599 - t1076 * t430;
	long double t1091 = t223 * t677;
	long double t1096 = 2 * t17 * t2;
	long double t1098 = -2 * t17 * t622 - t1096 * t374;
	long double t1109 = t237 * t717;
	long double t1116 = t1024 * t634 - 2 * t15 * t234 * t430;
	long double t1129 = t249 * t677;
	long double t1136 = 2 * t17 * t234 * t374 + t1052 * t587;
	long double t1152 = -t974 * t239 * t1109 / 2 + t340 * t209 * t1116 * t633 / 2 + t986 * t15 * t723 * t238 / 2 + SWr * t718 * t199 - t653 * t231 * t969 * t651 - t1001 * t251 * t1129 / 2 + t302 * t227 * t1136 * t581 / 2 + t1014 * t17 * t684 * t250 / 2 - SWr * t678 * t219 + t256 * t996 * t393 * t114 + t259 * t1042 * t471 * t163;
	long double t1153 = 4 * t193 * t952 * t514 * t189 - t974 * t210 * t1071 / 2 + t318 * t209 * t1078 * t598 / 2 - t986 * t15 * t723 * t206 / 2 - SWr * t753 * t199 + t617 * t63 * t961 * t615 - t1001 * t228 * t1091 / 2 + t330 * t227 * t1098 * t621 / 2 - t1014 * t17 * t684 * t224 / 2 + SWr * t745 * t219 + t1152;
	long double t1187 = t6 * t5 * t1;
	long double t1188 = t1187 * t267;
	long double t1191 = t1187 * t309 * t723;
	long double t1193 = t1187 * t270;
	long double t1196 = t1187 * t309 * t684;
	long double t1203 = t3 * t2 * t1;
	long double t1204 = t1203 * t275;
	long double t1206 = t268 * t607 * t597 * t199 - t271 * t628 * t620 * t219 - t268 * t642 * t632 * t199 + t271 * t593 * t580 * t219 - t276 * t434 * t422 * t91 + t276 * t459 * t447 * t91 - t279 * t389 * t370 * t132 - t1188 * t1071 - t1191 * t206 + t1193 * t1091 + t1196 * t224 + t1188 * t1109 + t1191 * t238 - t1193 * t1129 - t1196 * t250 + t1204 * t973;
	long double t1208 = t1203 * t309 * t322;
	long double t1210 = t1203 * t278;
	long double t1213 = t1203 * t309 * t306;
	long double t1222 = t5 * t267;
	long double t1224 = t5 * t270;
	long double t1228 = t2 * t275;
	long double t1230 = t2 * t278;
	long double t1234 = t1208 * t104 - t1210 * t1000 - t1213 * t138 - t1204 * t1020 - t1208 * t156 + t1210 * t1046 + t1213 * t183 + t279 * t418 * t409 * t132 + t1222 * t206 - t1224 * t224 - t1222 * t238 + t1224 * t250 - t1228 * t104 + t1230 * t138 + t1228 * t156 - t1230 * t183;
	long double t1239 = t1 * r;
	long double t1240 = t1239 * t309;
	long double t1241 = t5 * t3;
	long double t1242 = t1241 * t1240;
	long double t1244 = t1239 * SWr;
	long double t1245 = t2 * t6;
	long double t1246 = t1245 * t1244;
	long double t1248 = t200 * t267;
	long double t1250 = t200 * t270;
	long double t1255 = t1245 * t1240;
	long double t1272 = -t279 * t833 * t370 * t132 - t1242 * t334 - t1246 * t754 + t1248 * t206 - t1250 * t224 + t268 * t907 * t597 * t199 - t1255 * t758 + t1246 * t746 - t271 * t917 * t620 * t219 + t1255 * t750 + t1246 * t719 - t268 * t925 * t632 * t199 + t1255 * t725 - t1246 * t679 + t271 * t936 * t580 * t219 - t1255 * t686;
	long double t1273 = t1241 * t1244;
	long double t1292 = t97 * t275;
	long double t1294 = t97 * t278;
	long double t1298 = -t1273 * t339 - t276 * t798 * t422 * t91 - t1242 * t344 + t1273 * t301 + t279 * t811 * t409 * t132 + t1242 * t308 + t1273 * t317 + t276 * t819 * t447 * t91 + t1242 * t324 - t1273 * t329 - t1248 * t238 + t1250 * t250 + t1292 * t104 - t1294 * t138 - t1292 * t156 + t1294 * t183;
	long double t1300 = t979 * t43;
	long double t1306 = r * t309;
	long double t1307 = t979 * t1306;
	long double t1310 = t1006 * t43;
	long double t1316 = t1006 * t1306;
	long double t1333 = -t1300 * t754 / 2 + t268 * t1078 * t597 * t199 / 2 - t1307 * t758 / 2 + t1310 * t746 / 2 - t271 * t1098 * t620 * t219 / 2 + t1316 * t750 / 2 + t1300 * t719 / 2 - t268 * t1116 * t632 * t199 / 2 + t1307 * t725 / 2 - t1310 * t679 / 2 + t271 * t1136 * t580 * t219 / 2 - t1316 * t686 / 2;
	long double t1334 = t1076 * t43;
	long double t1340 = t1076 * t1306;
	long double t1343 = t1096 * t43;
	long double t1349 = t1096 * t1306;
	long double t1366 = t1334 * t339 / 2 - t276 * t981 * t422 * t91 / 2 + t1340 * t344 / 2 - t1343 * t301 / 2 + t279 * t1008 * t409 * t132 / 2 - t1349 * t308 / 2 - t1334 * t317 / 2 + t276 * t1028 * t447 * t91 / 2 - t1340 * t324 / 2 + t1343 * t329 / 2 - t279 * t1054 * t370 * t132 / 2 + t1349 * t334 / 2;
	Bi[0] += t188 + t262;
	Bi[1] += -t264 + t265;
	Bi[2] += t268 * t206 - t271 * t224 - t268 * t238 + t271 * t250 - t276 * t104 + t279 * t138 + t276 * t156 - t279 * t183;
	dBidrj[0][0] += t2 * t573 + t5 * t762;
	dBidrj[0][1] += t2 * t864 - t264 + t5 * t945 + t265;
	dBidrj[0][2] += t2 * t1065 + t5 * t1153;
	dBidrj[1][0] += -t5 * t573 + t2 * t762;
	dBidrj[1][1] += -t5 * t864 - t188 + t2 * t945 - t262;
	dBidrj[1][2] += -t5 * t1065 + t2 * t1153;
	dBidrj[2][0] += t1206 + t1234;
	dBidrj[2][1] += t1272 + t1298;
	dBidrj[2][2] += t1333 + t1366;				
}

// B field of a finite wire in the center
void StraightWireField_Center(const long double r, const long double phi, const long double z, const long double I_rt,
							const long double SW1z, const long double SW2z, long double Bi[3], long double dBidrj[3][3]){
	long double vorfaktor = mu0 * I_rt / (4 * pi);

	long double t1 = r * r;
	long double t2 = cos(phi);
	long double t3 = t2 * t2;
	long double t4 = t3 * t1;
	long double t5 = sin(phi);
	long double t6 = t5 * t5;
	long double t7 = t6 * t1;
	long double t8 = -SW2z + z;
	long double t9 = t8 * t8;
	long double t10 = t4 + t7 + t9;
	long double t11 = sqrt(t10);
	long double t12 = 0.1e1 / t11;
	long double t14 = -SW1z + z;
	long double t15 = t14 * t14;
	long double t16 = t4 + t7 + t15;
	long double t17 = sqrt(t16);
	long double t18 = 0.1e1 / t17;
	long double t20 = t8 * t12 - t18 * t14;
	long double t21 = fabs(t20);
	long double t22 = t21 * vorfaktor;
	long double t24 = -SW2z + SW1z;
	long double t25 = t3 + t6;
	long double t27 = t1 * t25;
	long double t28 = sqrt(t27);
	long double t29 = 0.1e1 / t28;
	long double t30 = t24 * t24;
	long double t31 = t30 * t27;
	long double t32 = sqrt(t31);
	long double t33 = 0.1e1 / t32;
	long double t35 = t33 * t29 * t25 * t24;
	long double t37 = fabs(t20) / t20;
	long double t38 = t37 * vorfaktor;
	long double t41 = t8 / t11 / t10;
	long double t44 = t3 * r + t6 * r;
	long double t48 = 0.1e1 / t17 / t16 * t14;
	long double t58 = t1 * t22;
	long double t59 = t25 * t25;
	Bi[1] += -t35 * r * t22;
	dBidrj[1][0] += -t35 * r * (-2 * t44 * t41 + 2 * t44 * t48) * t38 / 2 - t33 * t29 * t25 * t24 * t22 + t33 / t28 / t27 * t59 * t24 * t58 + 0.1e1 / t32 / t31 * t29 * t59 * t30 * t24 * t58;
	dBidrj[1][2] += -t35 * r * (-t8 * t41 + t12 - t18 + t14 * t48) * t38;								
}


