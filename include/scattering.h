/**
 * Collection of scattering methods
*/

#ifndef SCATTERING_H_
#define SCATTERING_H_

#include "mc.h"

#include "vectormath.h"

#include <boost/math/tools/minima.hpp>


/**
 * Calculate probability of specular reflection for a wave incident on a surface
 * 
 * @param costheta Cosine of incident polar angle
 * @param indexOfRefractionRatio Ratio n1/n2 of complex indices of refraction for current medium n1 and new medium n2
 * @param Returns Probability of specular reflection
*/
template<typename T>
T specularReflectionProbability(const T costheta, std::complex<T> indexOfRefractionRatio){
    std::complex<T> ncostheta = indexOfRefractionRatio*costheta;
    return std::norm((ncostheta - std::sqrt(1. - std::pow(indexOfRefractionRatio, 2) + std::pow(ncostheta, 2)))/
                     (ncostheta + std::sqrt(1. - std::pow(indexOfRefractionRatio, 2) + std::pow(ncostheta, 2))));
}


/**
 * Scattering distribution of wave incident on a surface following Lambert's law (outgoing polar angle cosine distributed around surface normal)
*/
template<typename T>
class lambert_scattering_distribution{
private:
    T _diffuseProbability; ///< Total probability of diffuse reflection
    std::complex<T> _indexOfRefractionRatio; ///< Ratio n1/n2 of indices of refraction for current medium n1 and new medium n2
    T _lossPerBounce; ///< Constant loss per bounce subtracted from reflection probability
    sincos_distribution<T> _polarSampler; ///< Distribution of polar angle (P(theta) = cos(theta) sin(theta) dtheta)
    std::uniform_real_distribution<T> _azimuthSampler; ///< Distribution of azimuth angle (constant)

    /**
     * Initialize member variables
     * 
     * @param diffuseProbability Total probability of diffuse reflection
     * @param indexOfRefractionRatio Ratio n1/n2 of complex indices of refraction for current medium n1 and new medium n2
     * @param lossPerBounce Constant loss per bounce subtracted from reflection probability
    */
    void _init(const T diffuseProbability, const std::complex<T> indexOfRefractionRatio, const T lossPerBounce){
        if (diffuseProbability < 0 or diffuseProbability > 1 or lossPerBounce < 0 or lossPerBounce > 1)
            throw std::runtime_error("Invalid probabilities in lambert_scattering_distribution!");
        _diffuseProbability = diffuseProbability;
        _indexOfRefractionRatio = indexOfRefractionRatio;
        _lossPerBounce = lossPerBounce;
        _polarSampler = sincos_distribution<T>(0, M_PI/2);
        _azimuthSampler = std::uniform_real_distribution<T>(0, 2*M_PI);
    }
public:
    /**
     * Construct distribution based on indices of refraction
     * 
     * @param diffuseProbability Total probability of diffuse reflection
     * @param indexOfRefractionRatio Ratio n1/n2 of complex indices of refraction for current medium n1 and new medium n2
     * @param lossPerBounce Constant loss per bounce subtracted from reflection probability
    */
    lambert_scattering_distribution(const T diffuseProbability, const std::complex<T> indexOfRefractionRatio, const T lossPerBounce){
        _init(diffuseProbability, indexOfRefractionRatio, lossPerBounce);
    }

    /**
     * Construct distribution for neutron scattering on Fermi potentials
     * 
     * @param Energy Kinetic energy of incident neutron
     * @param FermiPotential1 Fermi potential of current medium
     * @param FermiPotential2 Fermi potential of new medium
     * @param diffuseProbability Total probability of diffuse reflection
     * @param lossPerBounce Constant loss per bounce subtracted from reflection probability
    */
    lambert_scattering_distribution(const T Energy, const std::complex<T> &FermiPotential1, const std::complex<T> &FermiPotential2, const T diffuseProbability, const T lossPerBounce){
    	std::complex<T> k1 = std::sqrt(std::complex(Energy, std::imag(FermiPotential1))); // wavenumber in first solid
        std::complex<T> k2 = std::sqrt(Energy - (FermiPotential2 - std::real(FermiPotential1)));
        _init(diffuseProbability, k1/k2, lossPerBounce);
    }

    /**
     * Sample scattering distribution for wave incident from direction (sin(incidentPolarAngle), 0, cos(incidentPolarAngle))
     * 
     * @param incidentPolarAngle Angle of incidence with respect to surface normal [0..pi/2]
     * @param rng Random number generator
     * 
     * @returns Tuple of polar angle [0..pi] and azimuth [0..2pi] of scattered direction. A polar angle > pi/2 indicates transmission
    */
    template<class Random>
    std::tuple<T, T> operator()(const T incidentPolarAngle, Random &rng){
        if (incidentPolarAngle < 0 or incidentPolarAngle > M_PI/2)
            throw std::runtime_error("Invalid incident angle in lambert_scattering_distribution!");
        T reflectionProbability = specularReflectionProbability(std::cos(incidentPolarAngle), _indexOfRefractionRatio) - _lossPerBounce;
        if (std::generate_canonical<T, std::numeric_limits<T>::digits>(rng) < _diffuseProbability){
            T theta = _polarSampler(rng);
            T phi = _azimuthSampler(rng);
            if (std::generate_canonical<T, std::numeric_limits<T>::digits>(rng) < reflectionProbability){
                return std::make_tuple(theta, phi);
            }
            else{
                return std::make_tuple(M_PI - theta, phi);
            }
        }
        else{
            if (std::generate_canonical<T, std::numeric_limits<T>::digits>(rng) < reflectionProbability){
                return std::tuple<T, T>(incidentPolarAngle, 0.);
            }
            else{
                T sinSpecularRefractionAngle = std::real(_indexOfRefractionRatio)*std::sin(incidentPolarAngle);
                if (sinSpecularRefractionAngle <= 1){
                    return std::tuple<T, T>(M_PI - std::asin(sinSpecularRefractionAngle), 0.);
                }
                else{
                    return std::tuple<T, T>(M_PI, 0.);
                }
            }
        }
    }

    /**
     * Return ratio n1/n2 of complex indices of refraction for current medium n1 and new medium n2
    */
    std::complex<T> indexOfRefractionRatio() const{
        return _indexOfRefractionRatio;
    }
};


/**
 * Calculate scattered velocity vector from velocity vector incident on a surface and a given scattering distribution.
 * 
 * Calculates the incident polar angle, samples the scattering angle distribution, and calculates the scattered vector, all with respect to the surface normal
 * 
 * @param incidentVelocity Velocity incident on surface
 * @param surfaceNormal Normal vector of surface
 * @param scatteringDistribution Distribution of polar and azimuth scattering angles for given incident polar angle
 * @param rng Random number generator
 * 
 * @return Scattered velocity vector
*/
template<class Vector, class ScatterDistribution, class Random>
Vector scattered_vector(const Vector &incidentVelocity, const Vector &surfaceNormal, ScatterDistribution &scatteringDistribution, Random &rng){
    auto vabs = std::hypot(incidentVelocity[0], incidentVelocity[1], incidentVelocity[2]);
    auto nabs = std::hypot(surfaceNormal[0], surfaceNormal[1], surfaceNormal[2]);
    auto vnormal = (incidentVelocity[0]*surfaceNormal[0] + incidentVelocity[1]*surfaceNormal[1] + incidentVelocity[2]*surfaceNormal[2])/nabs;
    auto theta_i = std::acos(std::abs(vnormal)/vabs);
    auto [theta, phi] = scatteringDistribution(theta_i, rng);
    if (vnormal > 0){
        theta = M_PI - theta;
    }
    Vector out{std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)};
    RotateVector(out, surfaceNormal, incidentVelocity);
    auto outnormal = (out[0]*surfaceNormal[0] + out[1]*surfaceNormal[1] + out[2]*surfaceNormal[2])/nabs;
    if (vnormal*outnormal > 0){
        vabs = vabs/std::real(scatteringDistribution.indexOfRefractionRatio());
    }
    return {out[0]*vabs, out[1]*vabs, out[2]*vabs};
};


/**
 * Neutron scattering distribution from a surface based on the microroughness model 
 * See A. Steyerl, Effect of surface roughness on the total reflexion and transmission of slow neutrons, Z. Phys. A 254, 169–188 (1972)
 * https://doi.org/10.1007/BF01380066
*/
template<typename T>
class microroughness_scattering_distribution{
private:
    std::complex<T> _ki2; ///< squared incident wave number
    std::complex<T> _kc2; ///< squared critical wave number
    std::complex<T> _kt2; ///< squared transmitted wave number
    std::complex<T> _indexOfRefractionRatio; ///< Ratio n1/n2 of complex indices of refraction for current medium n1 and new medium n2
    T _RMSroughness; ///< RMS roughness amplitude of surface
    T _correlationLength2; ///< Correlation length of surface roughness autocorrelation function
    T _constantFactor; ///< Precomputed constant multiplication factor
    T _lossPerBounce; ///< Constant loss per bounce subtracted from reflection probability

    /**
     * Calculate absolute square of incident or reflected amplitude
     * 
     * @param costheta Cosine of incident/reflected polar angle
     * 
     * @returns Absolute square of amplitude of incident/reflected wave
    */
    T _normS(const T costheta) const {
        return std::norm(2./(1. + std::sqrt(1. - _kc2/_ki2/costheta/costheta)));
    }

    /**
     * Calculate absolute square of transmitted amplitude
     * 
     * @param costheta Cosine of transmitted polar angle
     * 
     * @returns Absolute square of amplitude of transmitted wave
    */
    T _normSprime(const T costheta) const {
        return std::norm(2./(1. + std::sqrt(1. + _kc2/_kt2/costheta/costheta)));
    }

    /**
     * Calculate fourier transform of momentum transfer F(mu) in surface plane
     * 
     * @param theta_i Incident polar angle
     * @param ko2 Squared scattered (reflected or transmitted) wave number
     * @param theta Scattered polar angle
     * @param phi Scattered azimuth
     * 
     * @returns Fourier transform of momentum transfer F(mu) in surface plane
    */
    T _Fmu(const T theta_i, const T ko2, const T theta, const T phi) const {
        T q2 = std::real(_ki2)*std::pow(std::sin(theta_i), 2) + ko2*std::pow(std::sin(theta), 2) - 2*std::sqrt(std::real(_ki2)*ko2)*std::sin(theta_i)*std::sin(theta)*std::cos(phi);
        return std::exp(-_correlationLength2/2*q2);
    }

public:
    /**
     * Calculate distribution of scattering angles for neutron incident from direction (sin(theta_i), 0, cos(theta_i))
     * 
     * @param theta_i Incident polar angle
     * @param theta Scattered polar angle
     * @param phi Scattered azimuth
     * 
     * @returns Tuple of probabilities for diffuse reflection and transmission
    */
    std::tuple<T, T> _distribution(const T theta_i, const T theta, const T phi) const {
        if (theta < 0 || theta > M_PI/2){
            return std::tuple<T, T>(0, 0);
        }
        T costheta_i = std::cos(theta_i);
        T costheta = std::cos(theta);
        T Iplus = _constantFactor/costheta_i * _normS(costheta_i) * _normS(costheta) * _Fmu(theta_i, std::real(_ki2), theta, phi);
        T Iminus = 0;
        if (std::real(_kt2) > 0 and costheta*costheta > -std::real(_kc2)/std::real(_kt2)){
            Iminus = _constantFactor/costheta_i * _normS(costheta_i) * _normSprime(costheta) / std::real(_indexOfRefractionRatio) * _Fmu(theta_i, std::real(_kt2), theta, phi);
        }
        return std::make_tuple(Iplus, Iminus);
    }

    /**
     * Construct scattering distribution for a neutron incident on a surface with representing a step in Fermi potential
     * 
     * @param Energy Kinetic energy of incident neutron
     * @param FermiPotential1 Fermi potential of current medium
     * @param FermiPotential2 Fermi potential of new medium
     * @param RMSroughness RMS roughness amplitude of surface
     * @param correlationLength Correlation length of surface roughness autocorrelation function
     * @param lossPerBounce Constant loss per bounce subtracted from reflection probability
    */
    microroughness_scattering_distribution(const T Energy, const std::complex<T> &FermiPotential1, const std::complex<T> &FermiPotential2, const T RMSroughness, const T correlationLength, const T lossPerBounce){
        _ki2 = static_cast<T>(2*m_n*std::pow(ele_e/hbar, 2))*std::complex(Energy, std::imag(FermiPotential1));
        std::complex<T> potentialStep = FermiPotential2 - std::real(FermiPotential1);
        _kc2 = static_cast<T>(2*m_n*std::pow(ele_e/hbar, 2))*potentialStep;
        _kt2 = _ki2 - _kc2;
        _indexOfRefractionRatio = std::sqrt(_ki2/_kt2);
        _RMSroughness = RMSroughness;
        _correlationLength2 = std::pow(correlationLength, 2);
        _constantFactor = std::pow(std::real(_kc2)*_RMSroughness, 2)*_correlationLength2/8/M_PI;
        _lossPerBounce = lossPerBounce;
    }

    /**
     * Sample scattering distribution for neutron incident from direction (sin(incidentPolarAngle), 0, cos(incidentPolarAngle))
     * 
     * @param incidentPolarAngle Angle of incidence with respect to surface normal [0..pi/2]
     * @param rng Random number generator
     * 
     * @returns Tuple of polar angle [0..pi] and azimuth [0..2pi] of scattered direction. A polar angle > pi/2 indicates transmission
    */
    template<class Random>
    std::tuple<T, T> operator()(const T incidentPolarAngle, Random &rng){
        auto [theta_min, I_min] = boost::math::tools::brent_find_minima(
            [incidentPolarAngle, this](const T theta){
                auto [Iplus, Iminus] = _distribution(incidentPolarAngle, theta, 0.); 
                return -Iplus - Iminus; 
            }, 0., M_PI/2, 10);
        unsigned int N_samples = std::ceil(-I_min*2*M_PI);
        T sample_max = static_cast<T>(N_samples)/2/M_PI;

        for (unsigned int i = 0; i < N_samples; ++i){
            T theta = std::acos(std::generate_canonical<T, std::numeric_limits<T>::digits>(rng));
            T phi = std::generate_canonical<T, std::numeric_limits<T>::digits>(rng)*2*M_PI;
            T u = std::generate_canonical<T, std::numeric_limits<T>::digits>(rng)*sample_max;
            auto [Iplus, Iminus] = _distribution(incidentPolarAngle, theta, phi);
            if (u < Iplus)
                return std::make_tuple(theta, phi);
            else if (u < Iplus + Iminus)
                return std::make_tuple(M_PI - theta, phi);
        }

        if (std::generate_canonical<T, std::numeric_limits<T>::digits>(rng) < _specularReflectionProbability(incidentPolarAngle)){
            return std::tuple<T, T>(incidentPolarAngle, 0.);
        }
        else if (std::real(_kt2) < 0){
            return std::tuple<T, T>(M_PI, 0.); // absorption during total reflection
        }
        else{
            return std::tuple<T, T>(M_PI - std::asin(std::real(_indexOfRefractionRatio)*std::sin(incidentPolarAngle)), 0.);
        }
    }

    /**
     * Return ratio n1/n2 of complex indices of refraction calculated from Fermi potentials and neutron energy
    */
    std::complex<T> indexOfRefractionRatio() const{
        return _indexOfRefractionRatio;
    }

    /**
     * Calculate probability of specular reflection with microroughness correction
     * 
     * @param incidentPolarAngle Polar angle of incident direction
     * 
     * @returns Probability of specular reflection
    */
    T _specularReflectionProbability(const T incidentPolarAngle) const{
        return 1. - (1. + std::real(_kc2)*_RMSroughness*_RMSroughness)*(1. - specularReflectionProbability(std::cos(incidentPolarAngle), _indexOfRefractionRatio)) - _lossPerBounce;
    }
};


/**
 * Scattering distribution taking into account distribution of microfacets of a surface in addition to a diffuse scattering distribution
 * See E. Heitz, E. d'Eon: "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
 * https://dx.doi.org/10.1111/cgf.12417
*/
template<typename T, class ScatterDistribution>
class microfacet_scattering_distribution{
private:
    T _beckmannWidth; ///< Width of microfacet (Beckmann) distribution
    ScatterDistribution _scattering_distribution; ///< Diffuse scattering distribution
    beckmann_visible_x_slope_distribution<T> _beckmann_x_sampler; ///< Distribution of surface slope parallel to incident direction
    std::normal_distribution<T> _beckmann_y_sampler; ///< Distribution of surface slope orthogonal to incident direction

    /**
     * Calculate Smith masking function G1 for given incident or scattered polar angle
     * 
     * @param polarAngle Incident or scattered polar angle
     * 
     * @returns Probability that given direction is not masked
    */
    T smithShadowing(const T polarAngle){
        if (polarAngle == 0 or polarAngle == M_PI){
            return 1.;
        }
        T a = static_cast<T>(1)/_beckmannWidth/std::abs(std::tan(polarAngle));
        T erf_a = std::erf(a);
        T exp_a2 = std::exp(-a*a);
        return static_cast<T>(2)/(1 + erf_a + exp_a2*M_2_SQRTPI/2/a);
    }
public:
    /**
     * Construct scattering distribution based on given diffuse scattering distribution and width of microfacet distribution
     * 
     * @param scattering_distribution Distribution of polar and azimuth scattering angles for given incident polar angle
     * @param beckmannWidth Width of microfacet (Beckmann) distribution
    */
    microfacet_scattering_distribution(const ScatterDistribution &scattering_distribution, const T beckmannWidth){
        _scattering_distribution = scattering_distribution;
        _beckmann_x_sampler = beckmann_visible_x_slope_distribution(0., beckmannWidth);
        _beckmann_y_sampler = std::normal_distribution<T>(0, M_SQRT1_2*beckmannWidth);
        _beckmannWidth = beckmannWidth;
    }

    /**
     * Sample scattering distribution for wave incident from direction (sin(incidentPolarAngle), 0, cos(incidentPolarAngle)),
     * taking into account microfacet distribution, diffuse scattering distribution, and masking/shadowing.
     * See E. Heitz, E. d'Eon: "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
     * https://dx.doi.org/10.1111/cgf.12417
     * 
     * @param incidentPolarAngle Angle of incidence with respect to surface normal [0..pi/2]
     * @param rng Random number generator
     * 
     * @returns Tuple of polar angle [0..pi] and azimuth [0..2pi] of scattered direction. A polar angle > pi/2 indicates transmission
    */
    template<class Random>
    std::tuple<T, T> operator()(const T incidentPolarAngle, Random &rng){
        for(;;){
            T slope_x = _beckmann_x_sampler(rng, incidentPolarAngle);
            T slope_y = _beckmann_y_sampler(rng);
            auto [theta, phi] = _scattering_distribution(incidentPolarAngle, rng);
            std::array<T, 3> o{std::sin(theta)*std::cos(phi),
                               std::sin(theta)*std::sin(phi),
                               std::cos(theta)};
            std::array<T, 3> n{-slope_x, -slope_y, 1.};
            std::array<T, 3> i{std::sin(incidentPolarAngle), 0, -std::cos(incidentPolarAngle)};
            RotateVector(o, n, i);
            T theta_o = std::acos(o[2]);
            if (theta < M_PI/2 xor theta_o < M_PI/2){
                continue; // don't allow scattering angle to be rotated into the other hemisphere
            }
            T G1 = smithShadowing(theta_o);
            if (std::generate_canonical<T>(rng) < G1){ // if outgoing direction is not masked/shadowed return scattered angles
                return std::make_tuple(theta, std::atan2(o[1], o[0]));
            }
        }
    }
};


#endif /*SCATTERING_H_*/