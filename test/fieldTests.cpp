/**
 * This file contains unit tests for (eventually all) field classes
 */

#include <random>

#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>

#include "analyticFields.h"
#include "conductor.h"
#include "globals.h"
#include "edmfields.h"
#include "fields.h"
#include "config.h"

#include <iostream>

std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> uni(-2., 2.);


/**
 * Compare magnetic fields calculated using two different field routines
 * The field template classes need to expose a BField(x, y, z, t, B[3], dBidxj[3][3]) routine
 */
template<class Field1, class Field2>
void compareMagneticFields(const Field1 &f1, const Field2 &f2,const double x, const double y, const double z, const double t = 0.){
    double field_tolerance = 1e-10; // difference in field components should only be due to rounding errors
    double deriv_tolerance = 1e-5; // use higher tolerance for derivatives, since they are approximated numerically
    double B1[3], dB1idxj[3][3], B2[3], dB2idxj[3][3];
    f1.BField(x, y, z, t, B1, nullptr); // check that function can be called without calculating derivatives
    f2.BField(x, y, z, t, B2, nullptr);

    for (int i = 0; i < 3; ++i){
        BOOST_TEST_INFO("i = " << i);
        BOOST_CHECK_SMALL(B1[i] - B2[i], field_tolerance); // check that difference between field components is smaller than field_tolerance
    }

    f1.BField(x, y, z, t, B1, dB1idxj);
    f2.BField(x, y, z, t, B2, dB2idxj);

    for (int i = 0; i < 3; ++i){
        BOOST_TEST_INFO("i = " << i);
        BOOST_CHECK_SMALL(B1[i] - B2[i], field_tolerance);
        for (int j = 0; j < 3; ++j){
        BOOST_TEST_INFO("i = " << i << ", j = " << j);
            BOOST_CHECK_SMALL(dB1idxj[i][j] - dB2idxj[i][j], deriv_tolerance); // check that difference between field derivatives is smaller than deriv_tolerance
        }
    }

}

/**
 * Check that the field routine returns a zero magnetic field
 * The field template classes need to expose a BField(x, y, z, t, B[3], dBidxj[3][3]) routine
 */
template<class Field>
void checkMagneticFieldZero(const Field &f, const double x, const double y, const double z){
    double B[3] = {0., 0., 0.}, dBidxj[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    f.BField(x, y, z, 0., B, dBidxj);
    for (int i = 0; i < 3; ++i){
        BOOST_TEST_INFO("i = " << i);
        BOOST_CHECK_EQUAL(B[i], 0.);
        for (int j = 0; j < 3; ++j){
            BOOST_TEST_INFO("i = " << i << ", j = " << j);
            BOOST_CHECK_EQUAL(dBidxj[i][j], 0.);
        }
    }
}

/**
 * Check that the field routine returns a zero electric field
 * The field template classes need to expose an EField(x, y, z, t, B[3], dBidxj[3][3]) routine
 */
template<class Field>
void checkElectricFieldZero(const Field &f, const double x, const double y, const double z){
    double V = 0, Ei[3] = {0, 0, 0};
    f.EField(x, y, z, 0., V, Ei);
    BOOST_CHECK_EQUAL(V, 0.);
    for (int i = 0; i < 3; ++i){
    BOOST_TEST_INFO("i = " << i);
        BOOST_CHECK_EQUAL(Ei[i], 0.);
    }
}


// compare field calculated from TExponentialFieldX to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TExponentialFieldXTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), a2 = uni(rng), a3 = uni(rng), c1 = uni(rng), c2 = uni(rng);
        TExponentialFieldX f1(a1, a2, a3, c1, c2);
        auto Bx = boost::format("%1$.20g * exp(- %2$.20g * x + %3$.20g) + %4$.20g") % a1 % a2 % a3 % c1;
        auto By = boost::format("y * %1$.20g * %2$.20g / 2 * exp(- %2$.20g*x + %3$.20g) + %4$.20g") % a1 % a2 % a3 % c2;
        auto Bz = boost::format("z * %1$.20g * %2$.20g / 2 * exp(- %2$.20g*x + %3$.20g) + %4$.20g") % a1 % a2 % a3 % c2;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << ", c1 = " << c1 << ", c2 = " << c2 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TLinearFieldZ to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TLinearFieldZTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), a2 = uni(rng);
        TLinearFieldZ f1(a1, a2);
        auto Bx = boost::format("0");
        auto By = boost::format("0");
        auto Bz = boost::format("%1$.20g * x + %2$.20g") % a1 % a2;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", a2 = " << a2 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TB0GradZ to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TB0GradZTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), a2 = uni(rng), z0 = uni(rng);
        TB0GradZ f1(a1, a2, z0);
        auto Bx = boost::format("(- %1$.20g / 2 * z - %2$.20g/2) * x") % a1 % a2;
        auto By = boost::format("(- %1$.20g / 2 * z - %2$.20g/2) * y") % a1 % a2;
        auto Bz = boost::format("%1$.20g / 2 * z*z + %2$.20g * z + %3$.20g") % a1 % a2 % z0;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", a2 = " << a2 << ", z0 = " << z0 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TB0GradZ to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TB0GradX2Test){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), a2 = uni(rng), a3 = uni(rng), z0 = uni(rng);
        TB0GradX2 f1(a1, a2, a3, z0);
        auto Bx = boost::format(" - %1$.20g / 6 * x^3 - %2$.20g/4 * x^2 - %3$.20g/2 * x") % a1 % a2 % a3;
        auto By = boost::format("(- %1$.20g / 2 * x^2 - %2$.20g/2 * x   - %3$.20g/2) * y") % a1 % a2 % a3;
        auto Bz = boost::format("(%1$.20g * x*x + %2$.20g * x + %3$.20g)*z + %4$.20g") % a1 % a2 % a3 % z0;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << ", z0 = " << z0 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TB0GradXY to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TB0GradXYTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), a2 = uni(rng), z0 = uni(rng);
        TB0GradXY f1(a1, a2, z0);
        auto Bx = boost::format("- %1$.20g / 4 * x^2*y - %2$.20g/2 * x") % a1 % a2;
        auto By = boost::format("- %1$.20g / 4 * x*y^2 - %2$.20g/2 * y") % a1 % a2;
        auto Bz = boost::format("%1$.20g * x*y*z + %2$.20g * z + %3$.20g") % a1 % a2 % z0;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", a2 = " << a2 << ", z0 = " << z0 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TB0_XY to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TB0_XYTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double a1 = uni(rng), z0 = uni(rng);
        TB0_XY f1(a1, z0);
        auto Bx = boost::format("%1$.20g *y*z") % a1;
        auto By = boost::format("%1$.20g *x*z") % a1;
        auto Bz = boost::format("%1$.20g *x*y + %2$.20g") % a1 % z0;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: a1 = " << a1 << ", z0 = " << z0 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// check that TCustomBField correctly identifies invalid formulas
BOOST_AUTO_TEST_CASE(TCustomBFieldTest){
    BOOST_CHECK_NO_THROW(TCustomBField("10.*x*t", "y*z", "x*y*z*t"));
    BOOST_CHECK_THROW(TCustomBField("a", "b", "c"), std::runtime_error);
    TCustomBField f("0", "0", "0");
    checkMagneticFieldZero(f, 1., 2., 3.);
}

// check that TFieldScaler correctly identifies invalid formulas and returns expected scaling factor
BOOST_AUTO_TEST_CASE(TFieldScalerTest){
    BOOST_CHECK_THROW(TFieldScaler("asgd"), std::runtime_error);
    BOOST_CHECK_NO_THROW(TFieldScaler());
    BOOST_CHECK_NO_THROW(TFieldScaler("10.*t+1."));
    TFieldScaler scaler;
    BOOST_CHECK_EQUAL(scaler.scalingFactor(0.), 0.);
    TFieldScaler scaler2("10*t + 1");
    BOOST_CHECK_EQUAL(scaler2.scalingFactor(0.), 1.);
    BOOST_CHECK_EQUAL(scaler2.scalingFactor(1.), 11.);
    double F[3] = {1., 1., 1.};
    double dFidxj[3][3] = {{1., 1., 1.}, {1., 1., 1.}, {1., 1., 1.}};
    scaler2.scaleVectorField(1., F, dFidxj);
    for (int i = 0; i < 3; ++i){
        BOOST_CHECK_EQUAL(F[i], scaler2.scalingFactor(1.));
        for (int j = 0; j < 3; ++j){
            BOOST_CHECK_EQUAL(dFidxj[i][j], scaler2.scalingFactor(1.));
        }
    }
}

// check that TFieldBoundaryBox correctly identifies invalid parameters and correctly scales within and outside of boundary (without smoothing)
BOOST_AUTO_TEST_CASE(TFieldBoundaryBoxTest){
    BOOST_CHECK_THROW(TFieldBoundaryBox(0., 1., 0., 0., 0., 0., 0.), std::runtime_error); // check that invalid parameter combinations throw exception
    BOOST_CHECK_THROW(TFieldBoundaryBox(0., 0., 0., 1., 0., 0., 0.), std::runtime_error);
    BOOST_CHECK_THROW(TFieldBoundaryBox(0., 0., 0., 0., 0., 1., 0.), std::runtime_error);
    BOOST_CHECK_THROW(TFieldBoundaryBox(0., 0., 0., 0., 0., 0., -1.), std::runtime_error);
    BOOST_CHECK_THROW(TFieldBoundaryBox(0., 0., 0., 0., 0., 0., 0.5), std::runtime_error);
    BOOST_CHECK_THROW(TFieldBoundaryBox(1., 0., 1., 0., 1., 0., 0.6), std::runtime_error);
    BOOST_CHECK_NO_THROW(TFieldBoundaryBox(1., 0., 1., 0., 1., 0., 0.5)); // valid parameter combination should not throw exception

    TFieldBoundaryBox box0;
    BOOST_CHECK(not box0.hasBounds()); // box constructed with default constructor should not have bounds
    BOOST_CHECK(box0.inBounds(0., 0., 0.));
    BOOST_CHECK(box0.inBounds(1., 1., 1.));

    TFieldBoundaryBox box(1., 0., 1., 0., 1., 0., 0.);
    BOOST_CHECK(box.hasBounds());
    BOOST_CHECK(box.inBounds(0.5, 0.5, 0.5)); // check that bounds are evaluated correctly
    BOOST_CHECK(not box.inBounds(1.01, 1.01, 1.01));

    double F[3] = {1., 1., 1.}, dFidxj[3][3] = {{1., 1., 1.}, {1., 1., 1.}, {1., 1., 1.}}; // set up test vector and derivates

    box.scaleVectorFieldAtBounds(0.5, 0.5, 0.5, F, dFidxj); // scale field within bounds without smoothing, field should not be scaled
    for (int i = 0; i < 3; ++i){
        BOOST_CHECK_EQUAL(F[i], 1.);
        for (int j = 0; j < 3; ++j){
            BOOST_CHECK_EQUAL(dFidxj[i][j], 1.);
        }
    }

    box.scaleVectorFieldAtBounds(1.01, 1.01, 1.01, F, dFidxj); // scale field outside bounds, field should be scaled to zero
    for (int i = 0; i < 3; ++i){
        BOOST_CHECK_EQUAL(F[i], 0.);
        for (int j = 0; j < 3; ++j){
            BOOST_CHECK_EQUAL(dFidxj[i][j], 0.);
        }
    }
}

// check that TFieldContainer properly applies TFieldBoundaryBox smoothing and TFieldScaler scaling
BOOST_AUTO_TEST_CASE(TFieldContainerTest){
    // set up field resembling the smoothing curve provided by a boundary box x/y/z = -1 .. 1 with boundary width 1, apply sin(t) oscillation
    std::string Bx = "abs(x) > 1 or abs(y) > 1 or abs(z) > 1 ? 0 : (6*(1 - abs(x))^5 - 15*(1 - abs(x))^4 + 10*(1 - abs(x))^3)*(6*(1 - abs(y))^5 - 15*(1 - abs(y))^4 + 10*(1 - abs(y))^3)*(6*(1 - abs(z))^5 - 15*(1 - abs(z))^4 + 10*(1 - abs(z))^3)*sin(t)";
    std::string By = Bx;
    std::string Bz = Bx;
    TCustomBField smoothed(Bx, By, Bz);
    std::unique_ptr<TField> f(new TCustomBField("1", "1", "1")); // set up homogeneous field
    TFieldContainer c(std::move(f), "sin(t)", "0", 1., -1., 1., -1., 1., -1., 1.); // combine homogeneous field with boundary box and sin(t) oscillation
 
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double x = uni(rng), y = uni(rng), z = uni(rng), t = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: x = " << x << ", y = " << y << ", z = " << z << ", t = " << t){
            compareMagneticFields(smoothed, c, x, y, z, t);
            checkElectricFieldZero(c, x, y, z);
            double B[3];
            c.BField(0., 0., 0., t, B, nullptr);
            BOOST_CHECK_EQUAL(sin(t), B[0]);
            BOOST_CHECK_EQUAL(sin(t), B[1]);
            BOOST_CHECK_EQUAL(sin(t), B[2]);
        }
    }

    for (int n = 0; n < nTests; ++n){
        double Ex = uni(rng), Ey = uni(rng), Ez = uni(rng);
        double x = uni(rng), y = uni(rng), z = uni(rng), t = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: Ex = " << Ex << ", Ey = " << Ey << ", Ez = " << Ez << ", x = " << x << ", y = " << y << ", z = " << z << ", t = " << t){
            TFieldContainer cE(std::unique_ptr<TField>(new TEDMStaticEField(Ex, Ey, Ez)), "0", "sin(t)", 0., 0., 0., 0., 0., 0., 0.);
            double V = 0, Ei[3] = {0., 0., 0.};
            cE.EField(x, y, z, t, V, Ei);
            BOOST_CHECK_EQUAL(Ex*sin(t), Ei[0]);
            BOOST_CHECK_EQUAL(Ey*sin(t), Ei[1]);
            BOOST_CHECK_EQUAL(Ez*sin(t), Ei[2]);
        }
    }
}

// compare field calculated from a conductor along the z axis to a TCustomBField with same field calculation formula, using randomly selected parameters and positions
BOOST_AUTO_TEST_CASE(TConductorFieldTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double z1 = uni(rng), z2 = uni(rng);
        TConductorField f1(0., 0., z1, 0., 0., z2, 4.*pi/mu0); // conductor along z axis
        auto sintheta = boost::format("(%1$.20g - z)/sqrt(x*x + y*y + (%1$.20g - z)^2)"); // angle between line connecting position and conductor end point and line perpendicular to conductor going through position
        std::string B = "1./sqrt(x*x + y*y)*(" + (sintheta % z2).str() + " - " + (sintheta % z1).str() + ")"; // absolute value of magnetic field B = mu0 I / (4 pi R)*(sin(theta2) - sin(theta1))
        TCustomBField f2(B + "*(-y)/sqrt(x*x + y*y)", B + "*x/sqrt(x*x + y*y)", "0."); // calculate x and y components from absolute field value (circular around conductor)
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: z1 = " << z1 << ", z2 = " << z2 << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// compare field calculated from TEDMStaticB0GradZField along the y axis to a TCustomBField with same field calculation formula, using randomly selected offsets, parameters and positions
BOOST_AUTO_TEST_CASE(TEDMStaticB0GradZFieldTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double xoff = uni(rng), yoff = uni(rng), zoff = uni(rng), abz = uni(rng), adB0zdz = uni(rng);
        TEDMStaticB0GradZField f1(xoff, yoff, zoff, pi/2, pi/2, abz, adB0zdz); // field with gradient along y axis
        auto By = boost::format("%1$.20g + %2$.20g*(y - %3$.20g)") % abz % adB0zdz % yoff;
        auto Bx = boost::format("-%1$.20g/2*(x - %2$.20g)") % adB0zdz % xoff;
        auto Bz = boost::format("-%1$.20g/2*(z - %2$.20g)") % adB0zdz % zoff;
        TCustomBField f2(Bx.str(), By.str(), Bz.str());
        double x = uni(rng), y = uni(rng), z = uni(rng);
        BOOST_TEST_CONTEXT("Parameters: xoff = " << xoff << ", yoff = " << yoff << ", zoff = " << zoff << ", abz = " << abz << ", adB0zdz = " << adB0zdz << ", x = " << x << ", y = " << y << ", z = " << z){
            compareMagneticFields(f1, f2, x, y, z);
            checkElectricFieldZero(f1, x, y, z);
        }
    }
}

// check that TEDMStaticEField returns a constant electric field at random positions
BOOST_AUTO_TEST_CASE(TEDMStaticEFieldTest){
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double Ex = uni(rng), Ey = uni(rng), Ez = uni(rng);
        TEDMStaticEField f1(Ex, Ey, Ez);
        double x = uni(rng), y = uni(rng), z = uni(rng), t = uni(rng);
        double V, Ei[3];
        f1.EField(x, y, z, t, V, Ei);
        BOOST_CHECK_EQUAL(Ei[0], Ex);
        BOOST_CHECK_EQUAL(Ei[1], Ey);
        BOOST_CHECK_EQUAL(Ei[2], Ez);
        checkMagneticFieldZero(f1, x, y, z);
    }
}

// integration test checking that TFieldManager properly handles scaling and boundaries
BOOST_AUTO_TEST_CASE(TFieldManagerTest){
    TConfig config({{"FIELDS", {{"0", "LinearFieldZ 0 1 1 -1 1 -1 1 -1 sin(t)"}}} , {"FORMULAS",{} }}); // homogeneous, oscillating field with hard boundaries at +/-1
    TFieldManager m(config);
    int nTests = 100;
    for (int n = 0; n < nTests; ++n){
        double x = uni(rng), y = uni(rng), z = uni(rng), t = uni(rng);
        double B[3];
        m.BField(x, y, z, t, B, nullptr);
        BOOST_TEST_CONTEXT("Parameters: x = " << x << ", y = " << y << ", z = " << z << ", t = " << t){
            BOOST_CHECK_EQUAL(B[0], 0.);
            BOOST_CHECK_EQUAL(B[1], 0.);
            if (abs(x) >= 1 or abs(y) >= 1 or abs(z) >= 1){
                BOOST_CHECK_EQUAL(B[2], 0.);
            }
            else{
                BOOST_CHECK_EQUAL(B[2], sin(t));
            }
        }
    }
}


/*****************************************************************************
 * MORE TO COME --- tests for TabField, TabField3, HarmonicExpandedBField, ...
 ****************************************************************************/
