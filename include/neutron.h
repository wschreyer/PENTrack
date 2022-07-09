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
	 * @param spinprojection component of semi-classical (Bloch) spin vector parallel to magnetic field
	 * @param amc Random number generator
	 * @param geometry Experiment geometry
	 * @param afield Optional fields (can be NULL)
	 */
	TNeutron(const int number, const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation, const double spinprojection,
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
	double MaterialPotential(const material &mat, const state_type &y) const{
		return mat.FermiReal*1e-9 - mat.InternalBField*GetMagneticMoment()*y[7]/ele_e;
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
	double CalcPotentialStep(const material &leaving, const material &entering, const state_type &y) const{
		return MaterialPotential(entering, y) - MaterialPotential(leaving, y);
	}
	
	/**
	 * Check for reflection/transmission/absorption on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Diffuse reflection can be done according to Lambert model or Micro Roughness model.
	 *
	 * For parameter doc see TParticle::OnHit
	 */
	void OnHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const double normal[3],
			const solid &leaving, const solid &entering, TMCGenerator &mc, stopID &ID, std::vector<TParticle*> &secondaries) const;


	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts neutron velocity.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Transmit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const double Estep) const;

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const;

	/**
	 * Transmit neutron through surface.
	 *
	 * Refracts or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void TransmitLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects neutron specularly.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void Reflect(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3]) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Micro Roughness model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectMR(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const double Estep, const material &mat, TMCGenerator &mc) const;

	/**
	 * Reflect neutron from surface.
	 *
	 * Reflects or scatters the neutron according to Lambert model.
	 * For parameter documentation see TNeutron::OnHit.
	 */
	void ReflectLambert(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const double normal[3], const material &mat, TMCGenerator &mc) const;

	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * For parameter doc see TParticle::OnStep
	 */
	void OnStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
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
	void Decay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, std::vector<TParticle*> &secondaries) const;


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
	value_type GetPotentialEnergy(const value_type t, const state_type &y, const TFieldManager &field, const solid &sld) const;
};


#endif // NEUTRON_H_
