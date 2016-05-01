/**
 * \file
 * Several optimized routines to calculate magnetic fields produced by (in)finite straight wires.
 */

#include <cmath>

#include "conductor.h"
#include "globals.h"

TConductorField::TConductorField(double aI, std::string Bscale): TField(Bscale, "1"), I(aI){ };


TFiniteWire::TFiniteWire(double SW1xx, double SW1yy, double SW1zz, double SW2xx, double SW2yy, double SW2zz, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SW1x(SW1xx), SW1y(SW1yy), SW1z(SW1zz), SW2x(SW2xx), SW2y(SW2yy), SW2z(SW2zz){

};

void TFiniteWire::BField(const double x, const double y, const double z, double t, double B[4][4])
{
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);

	double t1 = SW2z * SW2z;
	double t2 = z * z;
	double t3 = SW2y * SW2y;
	double t4 = SW2z * z;
	double t5 = 2 * t4;
	double t6 = y * y;
	double t7 = SW2y * y;
	double t8 = 2 * t7;
	double t10 = SW1x * SW1x;
	double t12 = -SW2y + y;
	double t13 = -SW2x + x;
	double t17 = -SW2z + z;
	double t23 = t17 * SW2z;
	double t29 = x * x;
	double t30 = SW2x * x;
	double t31 = 2 * t30;
	double t32 = SW2x * SW2x;
	double t34 = SW1y * SW1y;
	double t40 = 2 * y * t32;
	double t41 = y + SW2y;
	double t47 = SW2z * y;
	double t53 = SW1z * SW1z;
	double t55 = z * t32;
	double t56 = z + SW2z;
	double t59 = z * t3;
	double t62 = t29 + t6;
	double t68 = t4 + t7;
	double t74 = SW2y * z;
	double t78 = t10 * (t1 + t2 + t3 - t5 + t6 - t8) + SW1x * (-2 * SW1y * t13 * t12 - 2 * SW1z * t13 * t17 + 2 * SW2x * (t4 - t6 + t7 - t2) + 2 * (-t3 + t7 + t23) * x) + t34 * (t29 + t1 + t2 - t5 - t31 + t32) + SW1y * (-2 * SW1z * t12 * t17 - t40 + 2 * SW2x * t41 * x + 2 * SW2y * (-t2 + t4 - t29) + 2 * t17 * t47) + t53 * (t29 - t8 + t32 + t6 - t31 + t3) + 2 * SW1z * (-t55 + SW2x * t56 * x - t59 + SW2y * t56 * y - t62 * SW2z) + t32 * (t2 + t6) - 2 * SW2x * t68 * x + t3 * (t2 + t29) - 2 * t74 * t47 + t62 * t1;
	double t79 = sqrt(t78);
	double t80 = 0.1e1 / t79;
	double t81 = t80 * vorfaktor;
	double t82 = SW1y - y;
	double t90 = sqrt(t32 - 2 * SW2x * SW1x + t10 + t3 - 2 * SW2y * SW1y + t34 + t1 - 2 * SW2z * SW1z + t53);
	double t91 = 0.1e1 / t90;
	double t92 = t91 * t82;
	double t93 = SW2z - SW1z;
	double t95 = SW1z - z;
	double t96 = t91 * t95;
	double t97 = SW2y - SW1y;
	double t99 = t93 * t92 - t97 * t96;
	double t100 = t99 * t99;
	double t101 = SW2x - SW1x;
	double t103 = SW1x - x;
	double t104 = t91 * t103;
	double t106 = t101 * t96 - t93 * t104;
	double t107 = t106 * t106;
	double t110 = t97 * t104 - t101 * t92;
	double t111 = t110 * t110;
	double t112 = t100 + t107 + t111;
	double t113 = sqrt(t112);
	double t114 = 0.1e1 / t113;
	double t115 = t32 - t31 + t29 + t3 - t8 + t6 + t1 - t5 + t2;
	double t116 = sqrt(t115);
	double t117 = 0.1e1 / t116;
	double t119 = t101 * t91;
	double t122 = t97 * t91;
	double t125 = t93 * t91;
	double t133 = t10 - 2 * SW1x * x + t29 + t34 - 2 * SW1y * y + t6 + t53 - 2 * SW1z * z + t2;
	double t134 = sqrt(t133);
	double t135 = 0.1e1 / t134;
	double t142 = -t119 * t13 * t117 - t122 * t12 * t117 - t125 * t17 * t117 - t119 * t103 * t135 - t122 * t82 * t135 - t125 * t95 * t135;
	double t143 = fabs(t142);
	double t144 = t143 * t114;
	double t146 = SW1z * t12;
	double t147 = -SW1y * t17 + t146 + t74 - t47;
	double t151 = 0.1e1 / t79 / t78;
	double t153 = t114 * t151 * vorfaktor;
	double t154 = t147 * t143;
	double t156 = SW1z * t17;
	double t162 = x * SW2y;
	double t169 = SW2z * x;
	double t179 = 2 * SW1x * (-SW1y * t12 - t156 - t3 + t7 + t23) + 2 * t34 * t13 + SW1y * (2 * SW2x * t41 - 4 * t162) + 2 * t53 * t13 + SW1z * (2 * SW2x * t56 - 4 * t169) - 2 * SW2x * t68 + 2 * x * t3 + 2 * x * t1;
	double t184 = 0.1e1 / t113 / t112;
	double t185 = t184 * t81;
	double t186 = t91 * t106;
	double t188 = t91 * t110;
	double t190 = t93 * t186 - t97 * t188;
	double t194 = t114 * t81;
	double t195 = fabs(t142) / t142;
	double t197 = 0.1e1 / t116 / t115;
	double t198 = -t13 * t197;
	double t202 = t91 * t117;
	double t204 = -t12 * t197;
	double t208 = -t17 * t197;
	double t213 = 0.1e1 / t134 / t133;
	double t214 = t103 * t213;
	double t218 = t91 * t135;
	double t220 = t82 * t213;
	double t224 = t95 * t213;
	double t228 = -t13 * t119 * t198 - t101 * t202 - t13 * t122 * t204 - t13 * t125 * t208 - t103 * t119 * t214 + t101 * t218 - t103 * t122 * t220 - t103 * t125 * t224;
	double t229 = t228 * t195;
	double t234 = SW1y * t13;
	double t257 = 2 * t10 * t12 + SW1x * (-2 * t234 + (-4 * y + 2 * SW2y) * SW2x + 2 * t162) + 2 * SW1y * (-t156 - t32 + t30 + t23) + 2 * t53 * t12 + SW1z * (2 * SW2y * t56 - 4 * t47) + t40 - 2 * SW2x * t162 - 2 * SW2y * t4 + 2 * y * t1;
	double t261 = t91 * t99;
	double t264 = -t93 * t261 + t101 * t188;
	double t288 = -t12 * t119 * t198 - t12 * t122 * t204 - t97 * t202 - t12 * t125 * t208 - t82 * t119 * t214 - t82 * t122 * t220 + t97 * t218 - t82 * t125 * t224;
	double t289 = t288 * t195;
	double t293 = -t93 * t144 * t81;
	double t296 = SW1z * t13;
	double t300 = 2 * SW2z - 4 * z;
	double t319 = 2 * t10 * t17 + SW1x * (-2 * t296 + SW2x * t300 + 2 * t169) + 2 * t34 * t17 + SW1y * (-2 * t146 + SW2y * t300 + 2 * t47) + 2 * SW1z * (-t32 + t30 - t3 + t7) + 2 * t55 - 2 * SW2z * t30 + 2 * t59 - 2 * SW2y * t47;
	double t325 = t97 * t261 - t101 * t186;
	double t349 = -t17 * t119 * t198 - t17 * t122 * t204 - t17 * t125 * t208 - t93 * t202 - t95 * t119 * t214 - t95 * t122 * t220 - t95 * t125 * t224 + t93 * t218;
	double t350 = t349 * t195;
	double t359 = (-SW1x * t17 + t296 + SW2x * z - t169) * vorfaktor;
	double t361 = t143 * t114 * t80;
	double t363 = t151 * t359;
	double t367 = t80 * t359;
	double t368 = t143 * t184;
	double t372 = t195 * t114;
	double t386 = t361 * t101 * vorfaktor;
	double t398 = -SW1x * t12 + t234 + y * SW2x - t162;
	double t401 = t398 * t143;
	B[0][0] = t147 * t144 * t81;
	B[0][1] = -t179 * t154 * t153 / 2 - t190 * t154 * t185 + t147 * t229 * t194;
	B[0][2] = -t257 * t154 * t153 / 2 - t264 * t154 * t185 + t147 * t289 * t194 + t293;
	B[0][3] = -t319 * t154 * t153 / 2 - t325 * t154 * t185 + t147 * t350 * t194 + t97 * t144 * t81;
	B[1][0] = -t361 * t359;
	B[1][1] = -t293 + t179 * t144 * t363 / 2 + t190 * t368 * t367 - t228 * t372 * t367;
	B[1][2] = t257 * t144 * t363 / 2 + t264 * t368 * t367 - t288 * t372 * t367;
	B[1][3] = -t386 + t319 * t144 * t363 / 2 + t325 * t368 * t367 - t349 * t372 * t367;
	B[2][0] = t398 * t144 * t81;
	B[2][1] = -t179 * t401 * t153 / 2 - t190 * t401 * t185 + t398 * t229 * t194 - t97 * t144 * t81;
	B[2][2] = -t257 * t401 * t153 / 2 - t264 * t401 * t185 + t398 * t289 * t194 + t386;
	B[2][3] = -t319 * t401 * t153 / 2 - t325 * t401 * t185 + t398 * t350 * t194;

}


TFiniteWireX::TFiniteWireX(double SW1xx, double SW2xx, double SWzz, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SW1x(SW1xx), SW2x(SW2xx), SWz(SWzz){

};


void TFiniteWireX::BField(const double x, const double y, const double z, double t, double B[4][4])
{
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);


	double t1 = SW2x * SW2x;
	double t4 = x * x;
	double t5 = y * y;
	double t6 = SWz * SWz;
	double t8 = 2 * SWz * z;
	double t9 = z * z;
	double t10 = t1 - 2 * SW2x * x + t4 + t5 + t6 - t8 + t9;
	double t11 = sqrt(t10);
	double t12 = 0.1e1 / t11;
	double t13 = SW2x - x;
	double t15 = -SW1x + x;
	double t16 = SW1x * SW1x;
	double t19 = t16 - 2 * SW1x * x + t4 + t5 + t6 - t8 + t9;
	double t20 = sqrt(t19);
	double t21 = 0.1e1 / t20;
	double t23 = -t13 * t12 - t21 * t15;
	double t24 = fabs(t23);
	double t25 = t24 * vorfaktor;
	double t26 = -SWz + z;
	double t27 = t26 * t25;
	double t28 = -SW2x + SW1x;
	double t29 = t6 - t8 + t9 + t5;
	double t30 = sqrt(t29);
	double t31 = 0.1e1 / t30;
	double t33 = t28 * t28;
	double t34 = t33 * t29;
	double t35 = sqrt(t34);
	double t36 = 0.1e1 / t35;
	double t37 = t36 * t31 * t28;
	double t39 = fabs(t23) / t23;
	double t40 = t39 * vorfaktor;
	double t43 = t13 / t11 / t10;
	double t48 = 0.1e1 / t20 / t19 * t15;
	double t51 = -t13 * t43 + t12 - t21 + t15 * t48;
	double t55 = t36 * t31 * t28 * t26;
	double t59 = y * t43 + y * t48;
	double t63 = 0.1e1 / t30 / t29;
	double t64 = t63 * t28;
	double t65 = y * t36;
	double t69 = t31 * t33 * t28;
	double t71 = 0.1e1 / t35 / t34;
	double t78 = 2 * t26 * t43 + 2 * t26 * t48;
	double t91 = t31 * vorfaktor;
	double t92 = t24 * t91;
	double t95 = t39 * t91;
	double t97 = -t28 * y;
	double t101 = t24 * t63 * vorfaktor;
	B[1][0] = t37 * t27;
	B[1][1] = t55 * t51 * t40;
	B[1][2] = t55 * t59 * t40 - t65 * t64 * t27 - y * t71 * t69 * t27;
	B[1][3] = t55 * t78 * t40 / 2 + t37 * t25 - t26 * t36 * t64 * t27 - t26 * t71 * t69 * t27;
	B[2][0] = -t28 * t65 * t92;
	B[2][1] = t97 * t36 * t51 * t95;
	B[2][2] = t28 * t5 * t36 * t101 + t97 * t36 * t59 * t95 + t33 * t28 * t5 * t71 * t92 - t28 * t36 * t24 * t91;
	B[2][3] = t28 * t26 * t65 * t101 + t97 * t36 * t78 * t95 / 2 - t33 * t26 * t97 * t71 * t24 * t91;
}


TFiniteWireY::TFiniteWireY(double SW1yy, double SW2yy, double SWzz, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SW1y(SW1yy), SW2y(SW2yy), SWz(SWzz){

};


void TFiniteWireY::BField(const double x, const double y, const double z, double t, double B[4][4])
{
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);

	double t1 = SWz * SWz;
	double t3 = 2 * SWz * z;
	double t4 = z * z;
	double t5 = x * x;
	double t6 = t1 - t3 + t4 + t5;
	double t7 = sqrt(t6);
	double t9 = 0.1e1 / t7 * vorfaktor;
	double t10 = SW2y * SW2y;
	double t13 = y * y;
	double t14 = t5 + t10 - 2 * SW2y * y + t13 + t1 - t3 + t4;
	double t15 = sqrt(t14);
	double t16 = 0.1e1 / t15;
	double t17 = SW2y - y;
	double t19 = -SW1y + y;
	double t20 = SW1y * SW1y;
	double t23 = t5 + t20 - 2 * SW1y * y + t13 + t1 - t3 + t4;
	double t24 = sqrt(t23);
	double t25 = 0.1e1 / t24;
	double t27 = -t17 * t16 - t25 * t19;
	double t28 = fabs(t27);
	double t29 = t28 * t9;
	double t30 = SWz - z;
	double t31 = -SW2y + SW1y;
	double t32 = t31 * t30;
	double t33 = t31 * t31;
	double t34 = t33 * t6;
	double t35 = sqrt(t34);
	double t36 = 0.1e1 / t35;
	double t42 = t28 / t7 / t6 * vorfaktor;
	double t46 = fabs(t27) / t27;
	double t47 = t46 * t9;
	double t50 = t17 / t15 / t14;
	double t54 = 0.1e1 / t24 / t23 * t19;
	double t56 = x * t50 + x * t54;
	double t58 = t36 * t31;
	double t61 = t33 * t31;
	double t62 = t61 * t30;
	double t64 = 0.1e1 / t35 / t34;
	double t73 = -t17 * t50 + t16 - t25 + t19 * t54;
	double t77 = -2 * t30 * t36;
	double t83 = -2 * t30 * t50 - 2 * t30 * t54;
	double t89 = t36 * t31 * t28 * t9;
	double t90 = -2 * t30 * t64;
	double t95 = t31 * x;
	B[0][0] = t36 * t32 * t29;
	B[0][1] = -x * t36 * t32 * t42 + t58 * t30 * t56 * t47 - x * t64 * t62 * t29;
	B[0][2] = t58 * t30 * t73 * t47;
	B[0][3] = -t77 * t32 * t42 / 2 + t58 * t30 * t83 * t47 / 2 - t89 - t90 * t62 * t29 / 2;
	B[2][0] = t36 * t95 * t29;
	B[2][1] = -t36 * t31 * t5 * t42 + t58 * x * t56 * t47 + t89 - t64 * t61 * t5 * t29;
	B[2][2] = t58 * x * t73 * t47;
	B[2][3] = -t77 * t95 * t42 / 2 + t58 * x * t83 * t47 / 2 - t90 * t61 * x * t29 / 2;
}


TFiniteWireZ::TFiniteWireZ(double SWxx, double SWyy, double SW1zz, double SW2zz, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SWx(SWxx), SWy(SWyy), SW1z(SW1zz), SW2z(SW2zz){

};


void TFiniteWireZ::BField(const double x, const double y, const double z, double t, double B[4][4])
{
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);

	double t1 = SWx * SWx;
	double t3 = 2 * SWx * x;
	double t4 = x * x;
	double t5 = SWy * SWy;
	double t7 = 2 * SWy * y;
	double t8 = y * y;
	double t9 = SW1z * SW1z;
	double t12 = z * z;
	double t13 = t1 - t3 + t4 + t5 - t7 + t8 + t9 - 2 * SW1z * z + t12;
	double t14 = sqrt(t13);
	double t15 = 0.1e1 / t14;
	double t16 = SW1z - z;
	double t18 = SW2z * SW2z;
	double t21 = t1 - t3 + t4 + t5 - t7 + t8 + t18 - 2 * SW2z * z + t12;
	double t22 = sqrt(t21);
	double t23 = 0.1e1 / t22;
	double t24 = -SW2z + z;
	double t26 = t16 * t15 + t24 * t23;
	double t27 = fabs(t26);
	double t28 = t27 * vorfaktor;
	double t29 = -SWy + y;
	double t30 = t29 * t28;
	double t31 = -SW2z + SW1z;
	double t32 = t1 - t3 + t4 + t5 - t7 + t8;
	double t33 = sqrt(t32);
	double t34 = 0.1e1 / t33;
	double t36 = t31 * t31;
	double t37 = t36 * t32;
	double t38 = sqrt(t37);
	double t39 = 0.1e1 / t38;
	double t40 = t39 * t34 * t31;
	double t42 = fabs(t26) / t26;
	double t43 = t42 * vorfaktor;
	double t46 = t16 / t14 / t13;
	double t47 = -SWx + x;
	double t51 = t24 / t22 / t21;
	double t54 = (-2 * t47 * t46 - 2 * t47 * t51) * t43 / 2;
	double t56 = t39 * t34;
	double t57 = t56 * t31 * t29;
	double t61 = 0.1e1 / t33 / t32 * t31;
	double t63 = 2 * t47 * t39 * t61;
	double t67 = t34 * t36 * t31;
	double t69 = 0.1e1 / t38 / t37;
	double t71 = 2 * t47 * t69 * t67;
	double t78 = (-2 * t29 * t46 - 2 * t29 * t51) * t43 / 2;
	double t80 = t40 * t28;
	double t82 = 2 * t29 * t39 * t61;
	double t86 = 2 * t29 * t69 * t67;
	double t95 = (t16 * t46 - t15 - t24 * t51 + t23) * t43;
	double t97 = -t47 * t28;
	double t100 = -t56 * t31 * t47;
	B[0][0] = t40 * t30;
	B[0][1] = t57 * t54 - t63 * t30 / 2 - t71 * t30 / 2;
	B[0][2] = t57 * t78 + t80 - t82 * t30 / 2 - t86 * t30 / 2;
	B[0][3] = t57 * t95;
	B[1][0] = t40 * t97;
	B[1][1] = t100 * t54 - t80 - t63 * t97 / 2 - t71 * t97 / 2;
	B[1][2] = t100 * t78 - t82 * t97 / 2 - t86 * t97 / 2;
	B[1][3] = t100 * t95;
}



TFiniteWireZCenter::TFiniteWireZCenter(double SW1zz, double SW2zz, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SW1z(SW1zz), SW2z(SW2zz){

};


void TFiniteWireZCenter::BField(const double x, const double y, const double z, double t, double B[4][4]){
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);

	double t1 = x * x;
	double t2 = y * y;
	double t3 = SW2z * SW2z;
	double t6 = z * z;
	double t7 = t1 + t2 + t3 - 2 * SW2z * z + t6;
	double t8 = sqrt(t7);
	double t9 = 0.1e1 / t8;
	double t10 = SW2z - z;
	double t12 = -SW1z + z;
	double t13 = SW1z * SW1z;
	double t16 = t1 + t2 + t13 - 2 * SW1z * z + t6;
	double t17 = sqrt(t16);
	double t18 = 0.1e1 / t17;
	double t20 = -t10 * t9 - t18 * t12;
	double t21 = fabs(t20);
	double t22 = t21 * vorfaktor;
	double t23 = y * t22;
	double t24 = -SW2z + SW1z;
	double t25 = t2 + t1;
	double t26 = sqrt(t25);
	double t27 = 0.1e1 / t26;
	double t29 = t24 * t24;
	double t30 = t29 * t25;
	double t31 = sqrt(t30);
	double t32 = 0.1e1 / t31;
	double t33 = t32 * t27 * t24;
	double t35 = fabs(t20) / t20;
	double t36 = t35 * vorfaktor;
	double t39 = t10 / t8 / t7;
	double t43 = 0.1e1 / t17 / t16 * t12;
	double t45 = x * t39 + x * t43;
	double t49 = t32 * t27 * t24 * y;
	double t52 = 0.1e1 / t26 / t25;
	double t53 = t52 * t24;
	double t54 = x * t32;
	double t58 = t27 * t29 * t24;
	double t60 = 0.1e1 / t31 / t30;
	double t67 = y * t39 + y * t43;
	double t71 = t2 * t22;
	double t81 = -t10 * t39 + t9 - t18 + t12 * t43;
	double t84 = t27 * vorfaktor;
	double t85 = t21 * t84;
	double t89 = t21 * t52 * vorfaktor;
	double t93 = t35 * t84;
	double t95 = -t24 * x;
	B[0][0] = t33 * t23;
	B[0][1] = t49 * t45 * t36 - t54 * t53 * t23 - x * t60 * t58 * t23;
	B[0][2] = t49 * t67 * t36 + t33 * t22 - t32 * t53 * t71 - t60 * t58 * t71;
	B[0][3] = t49 * t81 * t36;
	B[1][0] = -t24 * t54 * t85;
	B[1][1] = t24 * t1 * t32 * t89 + t95 * t32 * t45 * t93 + t29 * t24 * t1 * t60 * t85 - t24 * t32 * t21 * t84;
	B[1][2] = t24 * y * t54 * t89 + t95 * t32 * t67 * t93 - t29 * y * t95 * t60 * t21 * t84;
	B[1][3] = t95 * t32 * t81 * t93;
}



TFullRacetrack::TFullRacetrack(double SW1zz, double SW2zz, double SWrr, double aI, std::string Bscale):
		TConductorField(aI, Bscale), SW1z(SW1zz), SW2z(SW2zz), SWr(SWrr){

};


void TFullRacetrack::BField(double x, double y, double z, double t, double B[4][4]){
	double vorfaktor = mu0 * I / (4 * pi)*BScaling(t);

	double t1 = x * x;
	double t2 = y * y;
	double t3 = SW2z * SW2z;
	double t5 = 2 * SW2z * z;
	double t6 = z * z;
	double t7 = t1 + t2 + t3 - t5 + t6;
	double t8 = sqrt(t7);
	double t9 = 0.1e1 / t8;
	double t10 = SW2z - z;
	double t12 = -SW1z + z;
	double t13 = SW1z * SW1z;
	double t15 = 2 * SW1z * z;
	double t16 = t1 + t2 + t13 - t15 + t6;
	double t17 = sqrt(t16);
	double t18 = 0.1e1 / t17;
	double t20 = -t10 * t9 - t18 * t12;
	double t21 = fabs(t20);
	double t22 = t21 * vorfaktor;
	double t23 = y * t22;
	double t24 = -SW2z + SW1z;
	double t25 = t2 + t1;
	double t26 = sqrt(t25);
	double t27 = 0.1e1 / t26;
	double t29 = t24 * t24;
	double t30 = t29 * t25;
	double t31 = sqrt(t30);
	double t32 = 0.1e1 / t31;
	double t33 = t32 * t27 * t24;
	double t36 = SWr * SWr;
	double t37 = SWr * x;
	double t38 = 2 * t37;
	double t39 = t2 + t36 - t38 + t1;
	double t40 = sqrt(t39);
	double t42 = vorfaktor / t40;
	double t43 = t36 - t38 + t1 + t2 + t3 - t5 + t6;
	double t44 = sqrt(t43);
	double t45 = 0.1e1 / t44;
	double t47 = t36 - t38 + t1 + t2 + t13 - t15 + t6;
	double t48 = sqrt(t47);
	double t49 = 0.1e1 / t48;
	double t51 = -t10 * t45 - t12 * t49;
	double t52 = fabs(t51);
	double t53 = t52 * t42;
	double t54 = t39 * t29;
	double t55 = sqrt(t54);
	double t56 = 0.1e1 / t55;
	double t57 = y * t56;
	double t60 = t2 + t36 + t38 + t1;
	double t61 = sqrt(t60);
	double t63 = vorfaktor / t61;
	double t64 = t36 + t38 + t1 + t2 + t3 - t5 + t6;
	double t65 = sqrt(t64);
	double t66 = 0.1e1 / t65;
	double t68 = t36 + t38 + t1 + t2 + t13 - t15 + t6;
	double t69 = sqrt(t68);
	double t70 = 0.1e1 / t69;
	double t72 = -t10 * t66 - t12 * t70;
	double t73 = fabs(t72);
	double t74 = t73 * t63;
	double t75 = t60 * t29;
	double t76 = sqrt(t75);
	double t77 = 0.1e1 / t76;
	double t78 = y * t77;
	double t81 = t3 - t5 + t6 + t1;
	double t82 = sqrt(t81);
	double t84 = vorfaktor / t82;
	double t85 = SWr * y;
	double t86 = 2 * t85;
	double t87 = t1 + t36 - t86 + t2 + t3 - t5 + t6;
	double t88 = sqrt(t87);
	double t89 = 0.1e1 / t88;
	double t90 = SWr - y;
	double t92 = y * t9;
	double t93 = t90 * t89 + t92;
	double t94 = fabs(t93);
	double t95 = t94 * t84;
	double t96 = t81 * t36;
	double t97 = sqrt(t96);
	double t98 = 0.1e1 / t97;
	double t99 = -t10 * t98;
	double t102 = t36 - t86 + t2 + t1;
	double t103 = sqrt(t102);
	double t105 = vorfaktor / t103;
	double t107 = t1 + t36 - t86 + t2 + t13 - t15 + t6;
	double t108 = sqrt(t107);
	double t109 = 0.1e1 / t108;
	double t111 = -t10 * t89 - t12 * t109;
	double t112 = fabs(t111);
	double t113 = t112 * t105;
	double t114 = t102 * t29;
	double t115 = sqrt(t114);
	double t116 = 0.1e1 / t115;
	double t117 = t90 * t116;
	double t120 = t13 - t15 + t6 + t1;
	double t121 = sqrt(t120);
	double t123 = vorfaktor / t121;
	double t124 = y * t18;
	double t126 = t124 + t90 * t109;
	double t127 = fabs(t126);
	double t128 = t127 * t123;
	double t129 = t120 * t36;
	double t130 = sqrt(t129);
	double t131 = 0.1e1 / t130;
	double t132 = -t12 * t131;
	double t135 = t1 + t36 + t86 + t2 + t3 - t5 + t6;
	double t136 = sqrt(t135);
	double t137 = 0.1e1 / t136;
	double t138 = -SWr - y;
	double t140 = SWr >= 0 ? 1 : -1;
	double t143 = -t140 * t138 * t137 - t140 * t92;
	double t144 = fabs(t143);
	double t145 = t144 * t84;
	double t146 = t10 * t98;
	double t149 = t36 + t86 + t2 + t1;
	double t150 = sqrt(t149);
	double t152 = vorfaktor / t150;
	double t154 = t1 + t36 + t86 + t2 + t13 - t15 + t6;
	double t155 = sqrt(t154);
	double t156 = 0.1e1 / t155;
	double t158 = -t10 * t137 - t12 * t156;
	double t159 = fabs(t158);
	double t160 = t159 * t152;
	double t161 = t149 * t29;
	double t162 = sqrt(t161);
	double t163 = 0.1e1 / t162;
	double t164 = t138 * t163;
	double t170 = -t140 * t124 - t140 * t138 * t156;
	double t171 = fabs(t170);
	double t172 = t171 * t123;
	double t173 = t12 * t131;
	double t177 = vorfaktor * t27;
	double t178 = t21 * t177;
	double t179 = x * t32;
	double t183 = t3 - t5 + t6 + t2;
	double t184 = sqrt(t183);
	double t186 = vorfaktor / t184;
	double t187 = SWr - x;
	double t189 = x * t9;
	double t190 = t187 * t45 + t189;
	double t191 = fabs(t190);
	double t192 = t191 * t186;
	double t193 = t183 * t36;
	double t194 = sqrt(t193);
	double t195 = 0.1e1 / t194;
	double t196 = t10 * t195;
	double t199 = -t187 * t56;
	double t202 = t13 - t15 + t6 + t2;
	double t203 = sqrt(t202);
	double t205 = vorfaktor / t203;
	double t207 = x * t18;
	double t208 = t187 * t49 + t207;
	double t209 = fabs(t208);
	double t210 = t209 * t205;
	double t211 = t202 * t36;
	double t212 = sqrt(t211);
	double t213 = 0.1e1 / t212;
	double t214 = t12 * t213;
	double t217 = -SWr - x;
	double t221 = -t140 * t217 * t66 - t140 * t189;
	double t222 = fabs(t221);
	double t223 = t222 * t186;
	double t224 = -t10 * t195;
	double t227 = -t24 * t217;
	double t233 = -t140 * t207 - t140 * t217 * t70;
	double t234 = fabs(t233);
	double t235 = t234 * t205;
	double t236 = -t12 * t213;
	double t239 = x * t116;
	double t242 = x * t163;
	double t246 = SWr * t195;
	double t247 = y * t246;
	double t249 = SWr * t213;
	double t250 = y * t249;
	double t254 = SWr * t98;
	double t255 = x * t254;
	double t257 = SWr * t131;
	double t258 = x * t257;
	double t264 = 0.1e1 / t76 / t75;
	double t265 = y * t264;
	double t266 = t29 * t24;
	double t274 = t159 * vorfaktor / t150 / t149;
	double t275 = x * t24;
	double t278 = fabs(t158) / t158;
	double t279 = t278 * t152;
	double t281 = 0.1e1 / t136 / t135;
	double t282 = t10 * t281;
	double t285 = 0.1e1 / t155 / t154;
	double t286 = t12 * t285;
	double t289 = t163 * (x * t282 + x * t286);
	double t290 = t24 * t138;
	double t294 = 0.1e1 / t162 / t161;
	double t295 = t138 * t294;
	double t296 = x * t266;
	double t301 = vorfaktor / t121 / t120;
	double t302 = t171 * t301;
	double t305 = fabs(t170) / t170;
	double t306 = t305 * t123;
	double t308 = 0.1e1 / t17 / t16;
	double t309 = y * t308;
	double t310 = x * t140;
	double t311 = t310 * t309;
	double t312 = t138 * t285;
	double t315 = t131 * (t311 + t310 * t312);
	double t316 = SWr * t12;
	double t320 = 0.1e1 / t130 / t129;
	double t321 = t12 * t320;
	double t322 = t36 * SWr;
	double t323 = x * t322;
	double t326 = fabs(t51) / t51;
	double t327 = t326 * t42;
	double t329 = 0.1e1 / t44 / t43;
	double t330 = t10 * t329;
	double t333 = 0.1e1 / t48 / t47;
	double t334 = t12 * t333;
	double t337 = t56 * (-2 * t187 * t330 - 2 * t187 * t334) / 2;
	double t338 = t24 * y;
	double t344 = t73 * vorfaktor / t61 / t60;
	double t351 = vorfaktor / t82 / t81;
	double t352 = t144 * t351;
	double t355 = fabs(t143) / t143;
	double t356 = t355 * t84;
	double t357 = t138 * t281;
	double t360 = 0.1e1 / t8 / t7;
	double t361 = y * t360;
	double t362 = t310 * t361;
	double t364 = t98 * (t310 * t357 + t362);
	double t365 = SWr * t10;
	double t369 = 0.1e1 / t97 / t96;
	double t370 = t10 * t369;
	double t376 = t112 * vorfaktor / t103 / t102;
	double t379 = -t217 * t266 * t265 * t74 - t275 * t164 * t274 + t290 * t289 * t279 - t296 * t295 * t160 - t37 * t173 * t302 + t316 * t315 * t306 - t323 * t321 * t172 - t338 * t337 * t327 - t24 * t217 * t78 * t344 - t37 * t146 * t352 + t365 * t364 * t356 - t323 * t370 * t145 - t275 * t117 * t376;
	double t380 = fabs(t111) / t111;
	double t381 = t380 * t105;
	double t383 = 0.1e1 / t88 / t87;
	double t384 = t10 * t383;
	double t387 = 0.1e1 / t108 / t107;
	double t388 = t12 * t387;
	double t391 = t116 * (x * t384 + x * t388);
	double t392 = t24 * t90;
	double t396 = 0.1e1 / t115 / t114;
	double t397 = t90 * t396;
	double t400 = t127 * t301;
	double t403 = fabs(t126) / t126;
	double t404 = t403 * t123;
	double t405 = x * t309;
	double t406 = t90 * t387;
	double t409 = t131 * (-t405 - x * t406);
	double t410 = -SWr * t12;
	double t413 = -t12 * t320;
	double t416 = fabs(t20) / t20;
	double t417 = t416 * vorfaktor;
	double t418 = t10 * t360;
	double t420 = t308 * t12;
	double t422 = x * t418 + x * t420;
	double t425 = t32 * t27 * t338;
	double t429 = 0.1e1 / t26 / t25;
	double t430 = t429 * t24;
	double t434 = t27 * t266;
	double t436 = 0.1e1 / t31 / t30;
	double t444 = t52 * vorfaktor / t40 / t39;
	double t445 = -2 * t187 * t24;
	double t450 = 0.1e1 / t55 / t54;
	double t452 = -2 * t187 * t266;
	double t456 = t94 * t351;
	double t459 = fabs(t93) / t93;
	double t460 = t459 * t84;
	double t461 = t90 * t383;
	double t463 = x * t361;
	double t465 = t98 * (-x * t461 - t463);
	double t466 = -SWr * t10;
	double t469 = -t10 * t369;
	double t472 = fabs(t72) / t72;
	double t473 = t472 * t63;
	double t475 = 0.1e1 / t65 / t64;
	double t476 = t10 * t475;
	double t479 = 0.1e1 / t69 / t68;
	double t480 = t12 * t479;
	double t482 = -2 * t217 * t476 - 2 * t217 * t480;
	double t486 = t392 * t391 * t381 - t296 * t397 * t113 - t37 * t132 * t400 + t410 * t409 * t404 - t323 * t413 * t128 + 4 * t425 * t422 * t417 - 4 * t179 * t430 * t23 - 4 * x * t436 * t434 * t23 + t445 * t57 * t444 / 2 + t452 * y * t450 * t53 / 2 - t37 * t99 * t456 + t466 * t465 * t460 - t323 * t469 * t95 - t338 * t77 * t482 * t473 / 2;
	double t490 = y * t418 + y * t420;
	double t496 = t2 * t22;
	double t509 = t56 * (y * t330 + y * t334);
	double t517 = t24 * t56 * t52 * t42;
	double t523 = y * t476 + y * t480;
	double t532 = t24 * t77 * t73 * t63;
	double t533 = 4 * t425 * t490 * t417 + 4 * t33 * t22 - 4 * t32 * t430 * t496 - 4 * t436 * t434 * t496 + t24 * t2 * t56 * t444 - t338 * t509 * t327 + t266 * t2 * t450 * t53 - t517 + t24 * t2 * t77 * t344 - t338 * t77 * t523 * t473 + t266 * t2 * t264 * t74 - t532;
	double t536 = t2 * t360;
	double t538 = t98 * (t90 * t461 - t89 - t536 + t9);
	double t541 = -2 * t24 * t90;
	double t548 = t116 * (-2 * t90 * t384 - 2 * t90 * t388) / 2;
	double t551 = -2 * t90 * t266;
	double t557 = t24 * t116 * t112 * t105;
	double t558 = t2 * t308;
	double t562 = t131 * (-t558 + t18 + t90 * t406 - t109);
	double t565 = -2 * t138 * t140;
	double t570 = t140 * t9;
	double t572 = t98 * (t565 * t357 / 2 + t140 * t137 + t140 * t536 - t570);
	double t575 = -2 * t24 * t138;
	double t582 = t163 * (-2 * t138 * t282 - 2 * t138 * t286) / 2;
	double t585 = -2 * t138 * t266;
	double t591 = t24 * t163 * t159 * t152;
	double t593 = t140 * t18;
	double t598 = t131 * (t140 * t558 - t593 + t565 * t312 / 2 + t140 * t156);
	double t601 = t466 * t538 * t460 - t541 * t117 * t376 / 2 + t392 * t548 * t381 - t551 * t397 * t113 / 2 - t557 + t410 * t562 * t404 + t365 * t572 * t356 - t575 * t164 * t274 / 2 + t290 * t582 * t279 - t585 * t295 * t160 / 2 - t591 + t316 * t598 * t306;
	double t607 = -t10 * t418 + t9 - t18 + t12 * t420;
	double t616 = t56 * (-t10 * t330 + t45 + t12 * t334 - t49);
	double t623 = -t10 * t476 + t66 + t12 * t480 - t70;
	double t627 = -2 * SWr * t10;
	double t634 = t98 * (2 * t10 * t461 + 2 * t10 * t361) / 2;
	double t637 = -2 * t10 * t322;
	double t643 = SWr * t98 * t94 * t84;
	double t649 = t116 * (-t10 * t384 + t89 + t12 * t388 - t109);
	double t652 = 2 * SWr * t12;
	double t659 = t131 * (-2 * t12 * t309 - 2 * t12 * t406) / 2;
	double t663 = 2 * t12 * t322;
	double t669 = SWr * t131 * t127 * t123;
	double t673 = -2 * t10 * t140;
	double t677 = t98 * (t673 * t357 + t673 * t361) / 2;
	double t685 = SWr * t98 * t144 * t84;
	double t691 = t163 * (-t10 * t282 + t137 + t12 * t286 - t156);
	double t697 = 2 * t12 * t140;
	double t701 = t131 * (t697 * t309 + t697 * t312) / 2;
	double t709 = SWr * t131 * t171 * t123;
	double t710 = -t663 * t413 * t128 / 2 - t669 - t627 * t146 * t352 / 2 + t365 * t677 * t356 - t637 * t370 * t145 / 2 - t685 + t290 * t691 * t279 - t652 * t173 * t302 / 2 + t316 * t701 * t306 - t663 * t321 * t172 / 2 + t709;
	double t713 = t21 * vorfaktor * t429;
	double t718 = t416 * t177;
	double t720 = -x * t24;
	double t733 = fabs(t190) / t190;
	double t734 = t733 * t186;
	double t735 = t187 * t329;
	double t738 = t1 * t360;
	double t740 = t195 * (t187 * t735 - t45 - t738 + t9);
	double t746 = -t187 * t24;
	double t749 = -t187 * t450;
	double t753 = fabs(t208) / t208;
	double t754 = t753 * t205;
	double t755 = t187 * t333;
	double t758 = t1 * t308;
	double t760 = t213 * (t187 * t755 - t49 - t758 + t18);
	double t763 = fabs(t221) / t221;
	double t764 = t763 * t186;
	double t765 = t217 * t475;
	double t766 = -2 * t217 * t140;
	double t772 = t195 * (t766 * t765 / 2 + t140 * t66 + t140 * t738 - t570);
	double t779 = 4 * t24 * t1 * t32 * t713 + 4 * t720 * t32 * t422 * t718 + 4 * t29 * t24 * t1 * t436 * t178 - 4 * t24 * t32 * t21 * t177 + t365 * t740 * t734 - t445 * t199 * t444 / 2 + t746 * t337 * t327 - t452 * t749 * t53 / 2 + t517 + t316 * t760 * t754 + t466 * t772 * t764 + t217 * t77 * t227 * t344;
	double t781 = t77 * t24;
	double t784 = -t217 * t266;
	double t789 = fabs(t233) / t233;
	double t790 = t789 * t205;
	double t792 = t217 * t479;
	double t797 = t213 * (t140 * t758 - t593 + t766 * t792 / 2 + t140 * t70);
	double t816 = -t781 * t217 * t482 * t473 / 2 + t532 + t217 * t264 * t784 * t74 + t410 * t797 * t790 - t24 * t1 * t116 * t376 + t275 * t391 * t381 - t266 * t1 * t396 * t113 + t557 - t24 * t1 * t163 * t274 + t275 * t289 * t279 - t266 * t1 * t294 * t160 + t591;
	double t825 = 0.1e1 / t194 / t193;
	double t826 = -t10 * t825;
	double t827 = y * t322;
	double t840 = vorfaktor / t184 / t183;
	double t841 = t191 * t840;
	double t846 = t195 * (-y * t735 - t463);
	double t849 = t10 * t825;
	double t856 = y * t140;
	double t859 = t195 * (t856 * t765 + t362);
	double t862 = t222 * t840;
	double t870 = vorfaktor / t203 / t202;
	double t871 = t209 * t870;
	double t874 = -4 * t29 * y * t720 * t436 * t21 * t177 - t827 * t826 * t223 + 4 * t24 * y * t179 * t713 + 4 * t720 * t32 * t490 * t718 - t85 * t196 * t841 + t365 * t846 * t734 - t827 * t849 * t192 - t338 * t199 * t444 + t746 * t509 * t327 + t466 * t859 * t764 - t85 * t224 * t862 - y * t266 * t749 * t53 - t85 * t214 * t871;
	double t877 = t213 * (-y * t755 - t405);
	double t881 = 0.1e1 / t212 / t211;
	double t882 = t12 * t881;
	double t892 = t234 * t870;
	double t897 = t213 * (t311 + t856 * t792);
	double t900 = -t12 * t881;
	double t921 = t316 * t877 * t754 - t827 * t882 * t210 - t78 * t227 * t344 - t781 * t217 * t523 * t473 - t265 * t784 * t74 - t85 * t236 * t892 + t410 * t897 * t790 - t827 * t900 * t235 - t541 * t239 * t376 / 2 + t275 * t548 * t381 - t575 * t242 * t274 / 2 + t275 * t582 * t279 - t551 * x * t396 * t113 / 2 - t585 * x * t294 * t160 / 2;
	double t931 = x * t360;
	double t934 = t195 * (2 * t10 * t735 + 2 * t10 * t931) / 2;
	double t942 = SWr * t195 * t191 * t186;
	double t949 = x * t308;
	double t952 = t213 * (-2 * t12 * t755 - 2 * t12 * t949) / 2;
	double t960 = SWr * t213 * t209 * t205;
	double t968 = t195 * (t673 * t765 + t673 * t931) / 2;
	double t976 = SWr * t195 * t222 * t186;
	double t986 = t213 * (t697 * t949 + t697 * t792) / 2;
	double t994 = SWr * t213 * t234 * t205;
	double t999 = -t627 * t224 * t862 / 2 + t466 * t968 * t764 - t637 * t826 * t223 / 2 + t976 - t781 * t217 * t623 * t473 - t652 * t236 * t892 / 2 + t410 * t986 * t790 - t663 * t900 * t235 / 2 - t994 + t275 * t649 * t381 + t275 * t691 * t279;
	double t1009 = t1 * t254;
	double t1013 = t322 * t369;
	double t1014 = t1 * t1013;
	double t1016 = t1 * t257;
	double t1020 = t322 * t320;
	double t1021 = t1 * t1020;
	double t1031 = t85 * t740 * t734 - t85 * t760 * t754 - t85 * t772 * t764 + t85 * t797 * t790 + t1009 * t456 - t37 * t465 * t460 + t1014 * t95 - t643 - t1016 * t400 + t37 * t409 * t404 - t1021 * t128 + t669 - t1009 * t352 + t37 * t364 * t356 - t1014 * t145 + t685 + t1016 * t302 - t37 * t315 * t306 + t1021 * t172 - t709;
	double t1032 = t2 * t246;
	double t1036 = t322 * t825;
	double t1037 = t2 * t1036;
	double t1039 = t2 * t249;
	double t1043 = t322 * t881;
	double t1044 = t2 * t1043;
	double t1062 = -t1032 * t841 + t85 * t846 * t734 - t1037 * t192 + t942 + t1039 * t871 - t85 * t877 * t754 + t1044 * t210 - t960 + t1032 * t862 - t85 * t859 * t764 + t1037 * t223 - t976 - t1039 * t892 + t85 * t897 * t790 - t1044 * t235 + t994 - t37 * t538 * t460 + t37 * t562 * t404 + t37 * t572 * t356 - t37 * t598 * t306;
	double t1063 = -2 * t10 * y;
	double t1064 = t1063 * t246;
	double t1069 = t1063 * t1036;
	double t1072 = 2 * t12 * y;
	double t1073 = t1072 * t249;
	double t1078 = t1072 * t1043;
	double t1093 = -t1064 * t841 / 2 + t85 * t934 * t734 - t1069 * t192 / 2 + t1073 * t871 / 2 - t85 * t952 * t754 + t1078 * t210 / 2 + t1064 * t862 / 2 - t85 * t968 * t764 + t1069 * t223 / 2 - t1073 * t892 / 2 + t85 * t986 * t790 - t1078 * t235 / 2;
	double t1094 = -2 * t10 * x;
	double t1095 = t1094 * t254;
	double t1100 = t1094 * t1013;
	double t1103 = 2 * t12 * x;
	double t1104 = t1103 * t257;
	double t1109 = t1103 * t1020;
	double t1124 = t1095 * t456 / 2 - t37 * t634 * t460 + t1100 * t95 / 2 - t1104 * t400 / 2 + t37 * t659 * t404 - t1109 * t128 / 2 - t1095 * t352 / 2 + t37 * t677 * t356 - t1100 * t145 / 2 + t1104 * t302 / 2 - t37 * t701 * t306 + t1109 * t172 / 2;
	B[0][0] = 4 * t33 * t23 - t24 * t57 * t53 - t24 * t78 * t74 + SWr * t99 * t95 + t24 * t117 * t113 + SWr * t132 * t128 + SWr * t146 * t145 + t24 * t164 * t160 + SWr * t173 * t172;
	B[0][1] = -4 * t24 * t179 * t178 + SWr * t196 * t192 + t24 * t199 * t53 + SWr * t214 * t210 + SWr * t224 * t223 + t77 * t227 * t74 + SWr * t236 * t235 + t24 * t239 * t113 + t24 * t242 * t160;
	B[0][2] = t247 * t192 - t250 * t210 - t247 * t223 + t250 * t235 - t255 * t95 + t258 * t128 + t255 * t145 - t258 * t172;
	B[0][3] = t379 + t486;
	B[1][0] = t533 + t601;
	B[1][1] = 4 * t425 * t607 * t417 - t338 * t616 * t327 - t338 * t77 * t623 * t473 - t627 * t99 * t456 / 2 + t466 * t634 * t460 - t637 * t469 * t95 / 2 + t643 + t392 * t649 * t381 - t652 * t132 * t400 / 2 + t410 * t659 * t404 + t710;
	B[1][2] = t779 + t816;
	B[1][3] = t874 + t921;
	B[2][0] = 4 * t720 * t32 * t607 * t718 - t627 * t196 * t841 / 2 + t365 * t934 * t734 - t637 * t849 * t192 / 2 - t942 + t746 * t616 * t327 - t652 * t214 * t871 / 2 + t316 * t952 * t754 - t663 * t882 * t210 / 2 + t960 + t999;
	B[2][1] = t1031;
	B[2][2] = t1062;
	B[2][3] = t1093 + t1124;
}


TInfiniteWireZ::TInfiniteWireZ(double lxx, double lyy, double aI, std::string Bscale):
		TConductorField(aI, Bscale), lx(lxx), ly(lyy){

};


void TInfiniteWireZ::BField(double x,double y,double z, double t, double B[4][4]){
	// cartesian coordinates of neutron
	// cartesian coordinates of racetracks
	double vorfaktor = mu0 * I / (2 * pi)*BScaling(t);

	double t1 = ly - y;
	double t2 = vorfaktor * t1;
	double t3 = lx * lx;
	double t6 = x * x;
	double t7 = ly * ly;
	double t10 = y * y;
	double t11 = t3 - 2 * lx * x + t6 + t7 - 2 * ly * y + t10;
	double t12 = 1 / t11;
	double t14 = t11 * t11;
	double t15 = 1 / t14;
	double t16 = -lx + x;
	double t17 = 2 * t15 * t16;
	double t19 = vorfaktor * t12;
	double t20 = -2 * t15 * t1;
	double t23 = -vorfaktor * t16;
	B[0][0] = t2 * t12;
	B[0][1] = -t2 * t17;
	B[0][2] = -t19 - t2 * t20;
	B[1][0] = -t23 * t12;
	B[1][1] = t19 + t23 * t17;
	B[1][2] = t23 * t20;
}


TInfiniteWireZCenter::TInfiniteWireZCenter(double aI, std::string Bscale):
		TConductorField(aI, Bscale){

};

void TInfiniteWireZCenter::BField(const double x, const double y, const double z, double t, double B[4][4]){
	double vorfaktor = mu0 * I / (2 * pi)*BScaling(t);
	double r2 = x*x + y*y;
	B[0][0] = vorfaktor*y/r2;
	B[0][1] = -2*vorfaktor*x*y/r2/r2;
	B[0][2] = vorfaktor*(x*x - y*y)/r2/r2;

	B[1][0] = -vorfaktor*x/r2;
	B[1][1] = vorfaktor*(x*x - y*y)/r2/r2;
	B[1][2] = 2*vorfaktor*x*y/r2/r2;
}
