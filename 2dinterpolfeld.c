#include "nr/nr3.h"

// coefficients for bicubic interpolation
void bcucof(VecDoub_I &y, VecDoub_I &y1, VecDoub_I &y2, VecDoub_I &y12,
	const Doub d1, const Doub d2, Doub c[4][4]) {
	static Int wt_d[16*16]=
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
		2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
		0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
		-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
		9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
		-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
		2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
		-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
		4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
	Int l,k,j,i;
	Doub xx,d1d2=d1*d2;
	VecDoub cl(16),x(16);
	static MatInt wt(16,16,wt_d);
	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++) {
		xx=0.0;
		for (k=0;k<16;k++) xx += wt[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) c[i][j]=cl[l++];
}

// bicubic interpolation with predetermined coefficients
void bcuint_new(Doub c[4][4],
	const Doub x1l, const Doub x1u, const Doub x2l, const Doub x2u,
	const Doub x1, const Doub x2, Doub &ansy, Doub &ansy1, Doub &ansy2) {
	Int i;
	Doub t,u,d1=x1u-x1l,d2=x2u-x2l;
	if (x1u == x1l || x2u == x2l)
		throw("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}
	ansy1 /= d1;
	ansy2 /= d2;
}


void CylToCart(long double v_r, long double v_phi, long double phi, long double &v_x, long double &v_y){
	v_x = v_r*cos(phi) - v_phi*sin(phi);
	v_y = v_r*sin(phi) + v_phi*cos(phi);
}
