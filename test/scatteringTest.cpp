#include <random>
#include <boost/test/unit_test.hpp>

#include "scattering.h"

template<typename T, class ScatterDistribution, class Random>
std::tuple<T, T, T, T> integrateScatteringDistribution(ScatterDistribution &dist, const T theta_i, Random &rng, int N = 400000){
    int diffuse_reflections = 0;
    int diffuse_transmissions = 0;
    int specular_reflections = 0;
    int specular_transmissions = 0;
    std::uniform_real_distribution uni(0., 1.);
    std::array<double, 3> n{uni(rng), uni(rng), uni(rng)};
    normalize(n);
    auto n_n = normalized(constructoOrthogonalVector(n));
    assert(dot(n, n_n) == 0.);
    std::array<double, 3> v = std::cos(theta_i)*n + std::sin(theta_i)*n_n;
    for (int i = 0; i < N; ++i){
        
        auto in_normal = dot(v, n);
        std::array<double, 3> out = scattered_vector(v, n, dist, rng);
        auto out_normal = dot(out, n);
        if (in_normal*out_normal < 0){
            auto r = v - 2*in_normal*n;
            if (mag(r - out) < 1e-10)
                ++specular_reflections;
            else
                ++diffuse_reflections;
        }
        else{
            auto v_plane = v - in_normal*n;
            auto t = std::copysign(std::sqrt(mag_sqr(v)/std::pow(std::real(dist.indexOfRefractionRatio()), 2) - mag_sqr(v_plane)), in_normal)*n + v_plane;
            if (mag(t - out) < 1e-10)
                ++specular_transmissions;
            else
                ++diffuse_transmissions;
        }
    }
    return {(double)specular_reflections/N, (double)specular_transmissions/N, (double)diffuse_reflections/N, (double)diffuse_transmissions/N};
}

BOOST_AUTO_TEST_CASE(lambertScatteringDistributionTest, *boost::unit_test::tolerance(1e-15)){
    std::random_device rd;
    std::mt19937 rng(rd());

    lambert_scattering_distribution lambert(0., {1., 0.}, 0.);
    auto [theta, phi] = lambert(M_PI/4, rng);
    BOOST_TEST(theta == M_PI - M_PI/4);
    BOOST_TEST(phi == 0);

    std::array<double, 3> v{-1., 0., 0.};
    RotateVector(v, {1., 0., 0.}, {0, 1, 0});
    BOOST_TEST(v[0] == 0);
    BOOST_TEST(v[1] == -1);
    BOOST_TEST(v[2] == 0);

    auto out = scattered_vector(std::array<double, 3>({1, 1, 1}), {-1, 0, 0}, lambert, rng);
    BOOST_TEST(out[0] == 1);
    BOOST_TEST(out[1] == 1);
    BOOST_TEST(out[2] == 1);

    lambert = lambert_scattering_distribution(0., {0., 0.}, 0.);
    auto [theta2, phi2] = lambert(M_PI/4, rng);
    BOOST_TEST(theta2 == M_PI/4);
    BOOST_TEST(phi2 == 0);

    out = scattered_vector(std::array<double, 3>({1, 1, 1}), {-1, 0, 0}, lambert, rng);
    BOOST_TEST(out[0] == -1);
    BOOST_TEST(out[1] == 1);
    BOOST_TEST(out[2] == 1);

    double diffProb = 0.5;
    std::complex<double> indexOfRefractionRatio{1.2, 0.};
    lambert = lambert_scattering_distribution(diffProb, indexOfRefractionRatio, 0.);
    double theta_i = std::uniform_real_distribution(0., M_PI/2)(rng);
    auto [specular_reflections, specular_transmissions, diffuse_reflections, diffuse_transmissions] = integrateScatteringDistribution(lambert, theta_i, rng);

    double reflProb = specularReflectionProbability(std::cos(theta_i), indexOfRefractionRatio);
    BOOST_TEST(specular_reflections == reflProb*(1 - diffProb), boost::test_tools::tolerance(0.005));
    BOOST_TEST(specular_transmissions == (1 - reflProb)*(1 - diffProb), boost::test_tools::tolerance(0.005));
    BOOST_TEST(diffuse_reflections == reflProb*diffProb, boost::test_tools::tolerance(0.005));
    BOOST_TEST(diffuse_transmissions == (1 - reflProb)*diffProb, boost::test_tools::tolerance(0.005));
}


BOOST_AUTO_TEST_CASE(microroughnessScatteringDistributionTest){
    std::random_device rd;
    std::mt19937 rng(rd());

    // parameters for Cu surface from Steyerl's microroughness paper (DOI:10.1007/BF01380066)
    double kCu = 0.00894/1e-10; // critical wave number
    double FermiCu = static_cast<double>(hbar*hbar*kCu*kCu/2./m_n/ele_e/ele_e);
    double bCu = 35e-10;
    double wCu = 250e-10;
    std::array<double, 3> normal = {0., 0., -1.};

    // parameters for incident neutron
    double wavelength = 300e-10;
    double theta_i = 73.5/180.*M_PI;
    double vabs = 2.*M_PI*hbar/m_n/ele_e/wavelength;
    double E = 0.5*m_n*vabs*vabs;
    std::array<double, 3> v = {vabs*sin(theta_i), 0., vabs*cos(theta_i)};

    // lines of constant diffuse reflection probabilities on copper in theta-phi plane,
    // extracted from Fig. 2 in Steyerl's microroughness paper (DOI:10.1007/BF01380066) using WebPlotDigitizer
    std::map<double, std::vector<std::pair<double, double> > > steyerlCopperReflectionData = {
        {0.0005, {{65.665, 44.107}, {77.432, 38.304}, {86.072, 27.490}, {88.573, 12.460}, {88.512, 359.44}, {87.751, 342.54}, {85.252, 330.72}, {77.299, 320.26}, {65.424, 314.54}, {52.683, 317.02}, {40.018, 317.73}, {26.161, 326.73}, {19.707, 359.86}, {30.638, 37.198}, {48.195, 42.103}}},
        {0.005, {{65.331, 35.519}, {78.811, 27.358}, {85.445, 14.884}, {86.546, 1.8522}, {85.757, 347.28}, {82.212, 336.17}, {64.792, 323.53}, {49.401, 328.01}, {33.618, 337.52}, {28.415, 359.69}, {37.325, 25.888}, {53.162, 31.649}}},
        {0.05, {{64.831, 334.14}, {54.423, 342.46}, {42.927, 359.56}, {48.403, 12.349}, {65.131, 24.720}, {74.116, 18.524}, {79.809, 9.3596}, {81.342, 359.45}, {79.503, 349.12}, {73.894, 340.63}}},
        {0.25, {{64.925, 12.537}, {70.860, 5.8185}, {72.309, 358.08}, {68.241, 349.91}, {64.589, 348.26}, {62.906, 358.69}}}
    };

    std::map<double, std::vector<std::pair<double, double> > > steyerlCopperTransmissionData = {
        {0.0005, {{22.951, 359.78}, {31.661, 31.397}, {49.094, 38.944}, {63.088, 38.478}, {77.196, 34.875}, {85.020, 27.877}, {88.153, 12.757}, {88.512, 359.44}, {87.861, 345.73}, {85.582, 333.07}, {73.933, 322.75}, {58.343, 320.05}, {45.141, 320.35}, {31.700, 326.53}}},
        {0.005, {{34.418, 1.6963}, {42.021, 21.572}, {55.352, 26.862}, {71.838, 25.515}, {82.857, 17.530}, {84.947, 10.403}, {85.781, 359.09}, {84.552, 347.23}, {77.954, 336.48}, {60.926, 331.80}, {46.729, 334.11}, {37.302, 344.43}}},
        {0.046, {{65.928, 4.7918}, {69.920, 0.60844}, {67.458, 355.91}, {61.103, 356.37}, {61.019, 3.5955}}}
    };

    microroughness_scattering_distribution MRdist(E, {0., 0.}, {FermiCu, 0.}, bCu, wCu, 0.);
    for (auto isoLine: steyerlCopperReflectionData){
        for (auto dataPoint: isoLine.second){
            BOOST_TEST_CONTEXT("Parameters: theta = " << dataPoint.first << ", phi = " << dataPoint.second << ", MRprob = " << isoLine.first){
                auto [MRreflection, MRtransmission] = MRdist._distribution(theta_i, dataPoint.first*M_PI/180., dataPoint.second*M_PI/180.);
                if (MRreflection > 0.01) BOOST_CHECK_CLOSE_FRACTION(MRreflection, isoLine.first, 0.2);
                else BOOST_CHECK_SMALL(MRreflection - isoLine.first, 0.005);
            }
        }
    }
    for (auto isoLine: steyerlCopperTransmissionData){
        for (auto dataPoint: isoLine.second){
            BOOST_TEST_CONTEXT("Parameters: theta = " << dataPoint.first << ", phi = " << dataPoint.second << ", MRprob = " << isoLine.first){
                auto [MRreflection, MRtransmission] = MRdist._distribution(theta_i, dataPoint.first*M_PI/180., dataPoint.second*M_PI/180.);
                if (MRtransmission > 0.01) BOOST_CHECK_CLOSE_FRACTION(MRtransmission, isoLine.first, 0.2);
                else BOOST_CHECK_SMALL(MRtransmission - isoLine.first, 0.005);
            }
        }
    }

    // total diffuse reflection probabilities on copper at given wavelengths in Angstrom,
    // extracted from Fig. 3 in Steyerl's microroughness paper (DOI:10.1007/BF01380066) using WebPlotDigitizer
    std::vector<std::pair<double, double> > steyerlTotalReflectionData = {{78.415,0.039293 }, {111.64,0.046966 }, {152.03,0.054637 }, {197.79,0.062128 }, {240.85,0.068547 }, {291.07,0.07532  }, {336.81,0.081023 }, {394.20, 0.087795}, {464.13,0.095455 }, {534.06,0.10276  }, {594.11,0.10774  },
                                                                          {638.91,0.11076  }, {696.23,0.11271  }, {727.42,0.10162  }, {756.83,0.092312 }, {788.94,0.083364 }, {838.03,0.072265 }, {884.47,0.063312 }, {924.66,0.056506 }, {963.08,0.050773 }};

    auto start = std::chrono::high_resolution_clock::now();
    for (auto dataPoint: steyerlTotalReflectionData){
        double k = 2.*M_PI/dataPoint.first/1e-10;
        double E = std::pow(hbar/ele_e*k, 2)/2/m_n;
        microroughness_scattering_distribution mr(E, {0., 0.}, {FermiCu, 0.}, bCu, wCu, 0.);
        if (k > kCu) theta_i = acos(2./3.*kCu/k);
        else theta_i = 48./180.*M_PI;
        auto [specular_reflections, specular_transmissions, diffuse_reflections, diffuse_transmissions] = integrateScatteringDistribution(mr, theta_i, rng);
        BOOST_TEST_CONTEXT("Parameters: wavelength = " << dataPoint.first << "A, reflection = " << diffuse_reflections << ", transmission = " << diffuse_transmissions){
            BOOST_CHECK_CLOSE_FRACTION(diffuse_reflections + diffuse_transmissions, dataPoint.second, 0.05);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Calculated " << 2*steyerlTotalReflectionData.size() << " microroughness integrals in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";

/*
    FermiNi = 230.e-9;
    bNi = 35e-10;
    wNi = 250e-10;
    double kNi = sqrt(2.*m_n*FermiNi)*ele_e/hbar;

    if (k > kNi) theta_i = acos(2./3.*kNi/k);
    else theta_i = 48./180.*pi;
    v = {vabs*sin(theta_i), 0., vabs*cos(theta_i)};
    double nickelReflection = 0, nickelTransmission = 0;
    nickelReflection = MR::MRProb(false, &v[0], &normal[0], FermiNi, bNi, wNi);
    nickelTransmission = MR::MRProb(true, &v[0], &normal[0], FermiNi, bNi, wNi);
*/
}