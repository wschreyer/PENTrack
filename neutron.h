/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_


#include <complex>

#include "globals.h"
#include "particle.h"
#include "proton.h"
#include "electron.h"

/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Constructor, create neutron, set start values randomly according to particle.in.
	 *
	 * Sets start time TParticle::tstart according to [SOURCE] in geometry.in.
	 * Sets start TParticle::polarisation according to particle.in.
	 * Sets start energy according to TMCGenerator::NeutronSpectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 * For volume sources the total energy TParticle::Hstart is diced by TMCGenerator::NeutronSpectrum
	 * and the position search is repeated until the kinetic energy TParticle::Estart is >0.
	 * Additionally the kinetic energy spectrum is weighted by sqrt(Estart) (see Golub/Richardson/Lamoreaux p. 82).
	 * For surface sources the kinetic energy TParticle::Estart is diced and all positions are allowed.
	 *
	 * @param number Particle number
	 * @param ageometry Geomtry in which the particle will be simulated
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TField used to calculate energies
	 */
	TNeutron(int number, TGeometry &ageometry, TSource &src, TMCGenerator &mcgen, TField *afield)
			: TParticle(NEUTRON, 0, m_n, mu_nSI){
		tstart = mcgen.UniformDist(0,src.ActiveTime);

		polarisation = mcgen.DicePolarisation(type);
		Hstart = mcgen.NeutronSpectrum();

		printf("Dice starting position for E_neutron = %LG neV ",Hstart*1e9);
		fflush(stdout);
		long double phi, theta;
		long double B[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
		for (int nroll = 0; nroll <= MAX_DICE_ROLL; nroll++){ // try to create particle only MAX_DICE_ROLL times
			if (nroll % 1000000 == 0){
				printf("."); // print progress
				fflush(stdout);
			}
			src.RandomPointInSourceVolume(mcgen, tstart, Hstart, ystart[0], ystart[1], ystart[2], phi, theta);
			ageometry.GetSolids(tstart, ystart, currentsolids);
			if (src.sourcemode == "volume" || src.sourcemode == "customvol"){ // dice H for volume source
				Estart = Hstart - Epot(tstart, &ystart[0], afield);
				if (mcgen.UniformDist(0,1) > sqrt(Estart/Hstart))
					Estart = -1; // correct weighting of phase space (see Golub/Richardson/Lamoreaux p. 82)
			}
			else{ // dice E for surface source
				Estart = Hstart;
				Hstart = Estart + Epot(tstart, &ystart[0], afield);
			}
			if (Estart >= 0){
				printf(" Found! %i dice(s) rolled\n", nroll+1);
				break;
			}
			else if (nroll >= MAX_DICE_ROLL){
				ID = ID_INITIAL_NOT_FOUND;
				printf(" ABORT: Failed %i times to find a compatible spot!! NO particle will be simulated!!\n\n", MAX_DICE_ROLL);
				return;
			}
		}

		InitE(number, tstart, mcgen.LifeTime(type), ystart[0], ystart[1], ystart[2],
			Estart, phi, theta, polarisation, mcgen.MaxTrajLength(type), ageometry, afield);
	};

protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities.
	 * Each reflection has a certain chance of being diffuse.
	 * Diffuse reflection angles are cosine-distributed around the normal vector of the surface.
	 * Logs surface hits if reflectlog is set in config.in.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering
	 * @param resetintegration Tell integrator to restart integration at x2 because e.g. y2 was changed
	 * @return Returns true if particle was reflected
	 */
	bool Reflect(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, long double normal[3], const solid *leaving, const solid *entering, bool &resetintegration){
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		bool log = false, refl = false;
		long double prob = mc->UniformDist(0,1);

		long double Estep = entering->mat.FermiReal*1e-9 - leaving->mat.FermiReal*1e-9;
//		cout << "Leaving " << leaving->ID << " Entering " << entering->ID << " Enormal = " << Enormal << " Estep = " << Estep;
		if (Enormal > Estep){ // transmission only possible if E > Estep
			long double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
			long double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
			long double transprob = 4*k1*k2/(k1 + k2)/(k1 + k2); // transmission probability
//			cout << " TransProb = " << transprob << '\n';
			if (prob < transprob){

				for (int i = 0; i < 3; i++)
					y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)

				resetintegration = true;
				if (reflectlog & 2)
					log = true;
			}
			else // no transmission -> reflection
				refl = true;
		}
		else{
			long double k1 = sqrt(Enormal); // wavenumber in first solid (only real part)
			complex<long double> iEstep(Estep, -entering->mat.FermiImag*1e-9); // potential step using imaginary potential V - i*W of second solid
			complex<long double> k2 = sqrt(Enormal - iEstep); // wavenumber in second solid (including imaginary part)
			long double reflprob = pow(abs((k1 - k2)/(k1 + k2)), 2); // reflection probability
//			cout << " ReflProb = " << reflprob << '\n';
			if (prob > reflprob){ // particle is not reflected
				ID = entering->ID; // -> absorption on reflection
				resetintegration = true;
			}
			else // no absorption -> reflection
				refl = true;
		}

		long double phi_r = 0, theta_r = 0;
		long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		if (refl){
			//particle was neither transmitted nor absorbed, so it has to be reflected
			prob = mc->UniformDist(0,1);
			if (prob >= entering->mat.DiffProb){
				//************** specular reflection **************
//				printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
				if(reflectlog & 1)
					log = true;
				x2 = x1;
				y2 = y1;
				y2[3] -= 2*vnormal*normal[0]; // reflect velocity
				y2[4] -= 2*vnormal*normal[1];
				y2[5] -= 2*vnormal*normal[2];
			}
			else{
				//************** diffuse reflection ************
				phi_r = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
				theta_r = mc->SinCosDist(0, 0.5*pi);
				if (vnormal > 0) theta_r += pi; // if normal points out of volume rotate by 180 degrees
				if(reflectlog & 1)
					log = true;
				x2 = x1;
				y2 = y1;
				y2[3] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
				y2[4] = vabs*sin(phi_r)*sin(theta_r);
				y2[5] = vabs*cos(theta_r);
				RotateVector(&y2[3],normal); // rotate coordinate system so that new z-axis lies on normal
//				printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
			}
		}

		const solid *currentsolid;
		if (refl)
			currentsolid = leaving;
		else
			currentsolid = entering;

		if (log){
			long double H = 0.5*m_n*vabs*vabs + m_n*gravconst*y1[2] + currentsolid->mat.FermiReal; // total neutron energy
			if (field){ // if there is a field, add magnetic potential to total neutron energy
				long double B[4][4];
				field->BFeld(y1[0], y1[1], y1[2], x1, B);
				H += -polarisation*mu_nSI/ele_e*B[3][0];
			}
			int pol = 3;
			fprintf(REFLECTLOG, "%i %i %i %i %i %i "
								"%.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG "
								"%.10LG %.10LG %.10LG %.10LG "
								"%.10LG %.10LG\n",
								jobnumber, particlenumber, type, polarisation, refl, currentsolid->ID,
								x1,y1[0],y1[1],y1[2],y1[3],y1[4],y1[5],
								normal[0],normal[1],normal[2],H,
								phi_r,theta_r);
			fflush(REFLECTLOG);
		}

		return refl;
	};


	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * Caclulates adiabacity condition (TNeutron::frac).
	 * Estimates spin flip probability according to Vladimirskii (TNeutron::vlad)
	 * and by doing a bruteforce integration of the Bloch equation describing spin precession.
	 * Additionally it can print out a spatial neutron distribution matrix.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle was absorbed
	 */
	bool Absorb(long double x1, VecDoub_I &y1, long double &x2, VecDoub_IO &y2, const solid *currentsolid){
		bool result = false;
		if (currentsolid->mat.FermiImag > 0){
			long double prob = mc->UniformDist(0,1);
			complex<long double> E(0.5*m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]), currentsolid->mat.FermiImag*1e-9); // E + i*W
			complex<long double> k = sqrt(2*m_n*E)*ele_e/hbar; // wave vector
			long double l = sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2)); // travelled length
			if (prob > exp(-imag(k)*l)){ // exponential probability decay
				x2 = x1 + mc->UniformDist(0,1)*(x2 - x1); // if absorbed, chose a random time between x1 and x2
				for (int i = 0; i < 6; i++)
					y2[i] = stepper->dense_out(i, x2, stepper->hdid); // set y2 to position at random time
				ID = currentsolid->ID;
				result = true; // stop integration
			}
		}

		// do special calculations for neutrons (spinflipcheck, snapshots, etc)
		if (neutdist == 1)
			fillndist(x1, y1, x2, y2); // write spatial neutron distribution

		if (field){
			long double B[4][4];
			field->BFeld(y1[0],y1[1],y1[2],x1,B);
/*
			if (B[3][0] > 0){
				// spin flip properties according to Vladimirsky and thumbrule
				vlad = vladimirsky(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				frac = thumbrule(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				vladtotal *= 1-vlad;
				if (vlad > 1e-99)
					vladmax = max(vladmax,log10(vlad));
				if (frac > 1e-99)
					thumbmax = max(thumbmax,log10(frac));
			}
*/
			long double B2[4][4];
			field->BFeld(y2[0],y2[1],y2[2],x2,B2);
			long double sp = BruteForceIntegration(x1,y1,B,x2,y2,B2); // integrate spinflip probability
//			if (1-sp > 1e-30) logBF = log10(1-sp);
//			else logBF = -99;
//			BFsurvprob *= sp;
			// flip the spin with a probability of 1-BFsurvprob
			if (flipspin && mc->UniformDist(0,1) < 1-sp)
			{
				polarisation *= -1;
				Nspinflip++;
				printf("\n The spin has flipped! Number of flips: %i\n",Nspinflip);
				result = true;
			}
		}

		return result;
	};


	/**
	 * Neutron decay.
	 *
	 * Create decay proton and electron and add them to secondaries list
	 */
	void Decay(){
		long double v_p[3], v_e[3];
		TParticle *p;
		mc->NeutronDecay(&yend[3], v_p, v_e);
		p = new TProton(particlenumber, tend, mc->LifeTime(PROTON),
							yend[0], yend[1], yend[2], v_p[0], v_p[1], v_p[2],
							mc->DicePolarisation(PROTON), mc->MaxTrajLength(PROTON), *geom, field);
		secondaries.push_back(p);
		p = new TElectron(particlenumber, tend, mc->LifeTime(ELECTRON),
							yend[0], yend[1], yend[2], v_e[0], v_e[1], v_e[2],
							mc->DicePolarisation(ELECTRON), mc->MaxTrajLength(ELECTRON), *geom, field);
		secondaries.push_back(p);
	};


	/**
	 * Potential energy, other than electromagnetic or gravitational, which is added to total energy of neutron
	 *
	 * @return Returns additional potential energy due to Fermi-Potential of solid
	 */
	long double Epot(const long double t, const long double y[3], TField *field){
		return TParticle::Epot(t, y, field) + currentsolids.rbegin()->mat.FermiReal*1e-9;
	}


};


#endif // NEUTRON_H_
