/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

#include "particle.h"
#include "microroughness.h"

extern const char* NAME_NEUTRON; ///< name of TNeutron class

/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Create neutron
	 *
	 * Wraps basic TParticle constructor
	 *
	 * @param number Particle number
	 * @param t Starting time
	 * @param x Initial x coordinate
	 * @param y Initial y coordinate
	 * @param z Initial z coordinate
	 * @param E Initial kinetic energy
	 * @param phi Initial azimuth of velocity vector
	 * @param theta Initial polar angle of velocity vector
	 * @param polarisation polarisation
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TNeutron(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation,
			TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield);
	
protected:
	/**
	 * Calculate neutron-specific potential in material
	 *
	 * @param mat Material to calculate potential of
	 * @param y Particle state used to calculate potential
	 *
	 * @return Return material potential [eV]
	 */
	double MaterialPotential(const double t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const int pol, const material &mat) const{
		return mat.FermiReal*1e-9 - mat.InternalBField*GetMagneticMoment()*pol/ele_e;
	}

	/**
	 * Calculate energy step from one material to another, including spin-dependent energy due to magnetization
	 *
	 * @param leaving Material that particle is leaving
	 * @param entering Material that particle is entering
	 * @param y Particle state containing spin polarization
	 *
	 * @return Returns step in Fermi potential [eV]
	 */
	double CalcPotentialStep(const material &leaving, const material &entering, const TStep &stepper) const{
		double t1 = stepper.GetStartTime();
		return MaterialPotential(stepper.GetTime(), stepper.GetPosition(), stepper.GetVelocity(), stepper.GetPolarization(), entering)
			 - MaterialPotential(t1, stepper.GetPosition(t1), stepper.GetVelocity(t1), stepper.GetPolarization(t1), leaving);
	}
	
	/**
	 * Check for reflection/transmission/absorption on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Diffuse reflection can be done according to Lambert model or Micro Roughness model.
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(TStep &stepper, const double normal[3],
			const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts neutron velocity.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Transmit(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering) const;

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitMR(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const;

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitLambert(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects neutron specularly.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Reflect(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectMR(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectLambert(TStep &stepper,
			const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc) const;

	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(TStep &stepper,
			const solid &currentsolid, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Simulate neutron beta decay.
	 *
	 * Calculates velocities of neutron decay products proton and electron.
	 *
	 * Reaction(s):
	 *
	 *     n0  ->  p+  +  R-
	 *
	 *     R-  ->  e-  +  nue
	 *
	 * Procedure:
	 * (1) Dice the energy of the decay proton according to globals.h#ProtonBetaSpectrum in the rest frame of the neutron and
	 *     calculate the proton momentum via the energy momentum relation.
	 * (2) Dice isotropic orientation of the decay proton.
	 * (3) Calculate 4-momentum of rest R- via 4-momentum conservation.
	 * (4) Get fixed electron energy from two body decay of R-.
	 * (5) Dice isotropic electron orientation in the rest frame of R-.
	 * (6) Lorentz boost electron 4-momentum into moving frame of R-.
	 * (7) Calculate neutrino 4-momentum via 4-momentum conservation.
	 * (8) Boost all 4-momentums into moving neutron frame.
	 *
	 * Cross-check:
	 * (9) Print neutrino 4-momentum invariant mass (4-momentum square, should be zero).
	 *
	 * @param t Decay time
	 * @param y Particle state before decay
	 * @param mc TMCGenerator class to generate random numbers
	 * @param geom TGeometry class
	 * @param field TFieldManager class
	 * @param secondaries Electron and proton from neutron decay are added to this list
	 */
	void Decay(const TStep &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;


	/**
	 * Potential energy, other than electromagnetic or gravitational, which is added to total energy of neutron
	 *
	 * @param t Time
	 * @param y Position vector
	 * @param field TFieldManager for electromagnetic potential
	 * @param sld Solid in which the particle is currently.
	 *
	 * @return Returns potential energy plus Fermi-Potential of solid
	 */
	double GetPotentialEnergy(const double t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const int pol, const TFieldManager &field, const solid &sld) const override;
};


#endif // NEUTRON_H_
