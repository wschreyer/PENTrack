//
// Created by wolfgang on 11/22/19.
//

#include <sstream>
#include <random>

#include "tracking.h"

using namespace std;

TTracker::TTracker(TConfig& config){
    logger = CreateLogger(config);
}

void TTracker::IntegrateParticle(std::unique_ptr<TParticle>& p, const double tmax, std::map<std::string, std::string> &particleconf,
        TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field){
    double tau = 0;
    istringstream(particleconf["tau"]) >> tau;
    if (tau > 0){
        exponential_distribution<double> expdist(1./tau);
        tau = expdist(mc);
    }
    else
        istringstream(particleconf["tmax"]) >> tau;

    double maxtraj;
    istringstream(particleconf["lmax"]) >> maxtraj;

//	cout << "Particle no.: " << particlenumber << " particle type: " << name << '\n';
//	cout << "x: " << yend[0] << "m y: " << yend[1] << "m z: " << yend[2]
//		 << "m E: " << GetFinalKineticEnergy() << "eV t: " << tend << "s tau: " << tau << "s lmax: " << maxtraj << "m\n";

    if (p->GetFinalTime() > tmax)
        throw runtime_error("Tried to start trajectory simulation past the maximum simulation time. Check time settings.");
    // set initial values for integrator
    TStep stepper(p->GetFinalTime(), p->GetFinalState());

    logger->PrintTrack(p, stepper, p->GetFinalSpin(), p->GetFinalSolid(), field);

    bool flipspin = false;
    istringstream(particleconf["flipspin"]) >> flipspin;

    bool spininterpolatefields = false;
    istringstream(particleconf["interpolatefields"]) >> spininterpolatefields;

    double SpinBmax = 0;
    vector<double> SpinTimes;
    istringstream(particleconf["Bmax"]) >> SpinBmax;
    istringstream SpinTimess(particleconf["spintimes"]);
    do{
        double t;
        SpinTimess >> t;
        if (SpinTimess)
            SpinTimes.push_back(t);
    }while(SpinTimess.good());
    state_type spin = p->GetFinalSpin();

    currentsolids = geom.GetSolids(p->GetFinalTime(), &p->GetFinalState()[0]);
    p->SetStopID(ID_UNKNOWN);

    while (p->GetStopID() == ID_UNKNOWN){ // integrate as long as nothing happened to particle
        stopID ID = stepper.next([&](const state_type &y, state_type &dydx, const value_type &x){ p->derivs(y, dydx, x, field); }, geom);    //std::bind(&TParticle::derivs, p.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, field), geom);
        p->SetStopID(ID);
        if (ID != ID_UNKNOWN)
            break;

        if (stepper.IsStepOverCollision()){
            DoHit(p, stepper, mc, geom);
        }
        else{
            if (tau > 0 && stepper.GetProperTime() >= tau){ // if proper time is larger than lifetime
                stepper.SetStepEndToMatchComponent(6, tau);
                p->SetStopID(ID_DECAYED);
            }
            if (stepper.GetTime() >= tmax){	//If stepsize overshot max simulation time
                stepper.SetStepEnd(tmax);
                p->SetStopID(ID_NOT_FINISH);
            }
            if (stepper.GetPathLength() >= maxtraj){
                stepper.SetStepEndToMatchComponent(8, maxtraj);
                p->SetStopID(ID_NOT_FINISH);
            }

            DoStep(p, stepper, GetCurrentsolid(), mc, field);
        }

        // take snapshots at certain times
        logger->PrintSnapshot(p, stepper, spin, geom, field);

        IntegrateSpin(p, spin, stepper, SpinTimes, field, spininterpolatefields, SpinBmax, mc, flipspin); // calculate spin precession and spin-flip probability

        logger->PrintTrack(p, stepper, spin, GetCurrentsolid(), field);
    }

    if (p->GetStopID() == ID_DECAYED){
        p->DoDecay(stepper.GetTime(), stepper.GetState(), mc, geom, field);
    }

    p->SetFinalState(stepper.GetTime(), stepper.GetState(), spin, GetCurrentsolid());
    logger->Print(p, stepper.GetTime(), stepper, spin, geom, field);


//	cout << "x: " << yend[0];
//	cout << " y: " << yend[1];
//	cout << " z: " << yend[2];
//	cout << " E: " << GetFinalKineticEnergy();
//	cout << " Code: " << ID;
//	cout << " t: " << tend;
//	cout << " l: " << yend[8];
//	cout << " hits: " << Nhit;
//	cout << " spinflips: " << Nspinflip << '\n';
//	cout << "Computation took " << Nstep << " steps\n";
//	cout << "Done!!\n\n";
//			cout.flush();
}


void TTracker::DoStep(const std::unique_ptr<TParticle>& p, TStep &stepper, const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field) {
    p->DoStep(stepper, currentsolid, mc, field);
}

void TTracker::DoHit(const std::unique_ptr<TParticle>& p, TStep &stepper, TMCGenerator &mc, const TGeometry &geom) {
    bool traversed = true;

    multimap<TCollision, bool> colls = stepper.GetCollisions();

    vector<pair<solid, bool> > newsolids = currentsolids;
    for (auto coll: colls){
//    cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << coll.first.ID << endl;
        solid sld = geom.GetSolid(coll.first.ID);
        auto foundsld = find_if(newsolids.begin(), newsolids.end(), [&sld](const std::pair<solid, bool> s){ return s.first.ID == sld.ID; });
        if (coll.first.distnormal < 0){ // if entering solid
            if (foundsld != newsolids.end()){ // if solid has been entered before (self-intersecting surface)
//	cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << sld.name << endl;
//                if (coll.first.s > 0) // if collision happened right at the start of the step it is likely that the hit solid was already added to the list in the previous step and we will ignore this one
                    newsolids.push_back(make_pair(sld, foundsld->second)); // add additional entry to list, with ignore state as on first entry
            }
            else
                newsolids.push_back(make_pair(sld, coll.second)); // add solid to list
        }
        else if (coll.first.distnormal > 0){ // if leaving solid
            if (foundsld == newsolids.end()){ // if solid was not entered before something went wrong
//                if (coll.first.s > 0){ // if collision happened right at the start of the step it is likely that the hit solid was already removed from the list in the previous step and this is not an error
//	  cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << sld.name << endl;
//          throw runtime_error((boost::format("Particle inside '%1%' which it did not enter before!") % sld.name).str());
                    cout << "Particle inside solid " << sld.name << " which it did not enter before. Stopping it!\n";
                    p->SetStopID(ID_GEOMETRY_ERROR);
                    exit(-1);
                    return;
                }
//            }
            else
                newsolids.erase(foundsld); // remove solid from list
        }
        else{
//      throw runtime_error("Particle crossed surface with parallel track!");
            cout << "Particle crossed surface with parallel track. Stopping it!\n";
            p->SetStopID(ID_GEOMETRY_ERROR);
            return;
        }
    }

    solid leaving = GetCurrentsolid(); // particle can only leave highest-priority solid
    solid entering = geom.defaultsolid;
    for (auto sld: newsolids){
        if (!sld.second && sld.first.ID > entering.ID)
            entering = sld.first;
    }
//  cout << "Leaving " << leaving.name << ", entering " << entering.name << '\n';
    if (leaving.ID != entering.ID){ // if the particle actually traversed a material interface
        auto coll = find_if(colls.begin(), colls.end(), [&leaving, &entering](const pair<TCollision, bool> &c){ return !c.second && (c.first.ID == leaving.ID || c.first.ID == entering.ID); });
        if (coll == colls.end())
            throw std::runtime_error("Did not find collision corresponding to entering/leaving solid!");
        value_type prevtime = stepper.GetTime();
        auto prevpos = stepper.GetPosition();
        p->DoHit(stepper, coll->first.normal, leaving, entering, mc); // do particle specific things
        if (stepper.GetTime() == prevtime && stepper.GetPosition() == prevpos){
            traversed = true;
        }
        else if (stepper.GetTime() == stepper.GetStartTime() && stepper.GetPosition() == stepper.GetPosition(stepper.GetStartTime())){
            traversed = false;
        }
        else{
            throw std::runtime_error("OnHit routine returned inconsistent position. That should not happen!");
        }

        logger->PrintHit(p, stepper, coll->first.normal, leaving, entering); // print collision to file if requested
    }

    if (traversed){
        currentsolids = newsolids; // if surface was traversed (even if it was  physically ignored) replace current solids with list of new solids
    }
}

const solid& TTracker::GetCurrentsolid() const{
    auto sld = max_element(currentsolids.begin(), currentsolids.end(), [](const pair<solid, bool> &s1, const pair<solid, bool> &s2){ return s1.second || (!s2.second && s1.first.ID < s2.first.ID); });
    return sld->first;
}


void TTracker::IntegrateSpin(const std::unique_ptr<TParticle>& p, state_type &spin, const TStep &stepper,
        const std::vector<double> &times, const TFieldManager &field,
        const bool interpolatefields, const double Bmax, TMCGenerator &mc, const bool flipspin) const{
    value_type x1 = stepper.GetStartTime();
    value_type x2 = stepper.GetTime();
    if (p->GetGyromagneticRatio() == 0 || x1 == x2)
        return;

    state_type y1 = stepper.GetState(x1);
    state_type y2 = stepper.GetState(x2);
    double B1[3], B2[3], polarisation;
    field.BField(y1[0], y1[1], y1[2], x1, B1);
    field.BField(y2[0], y2[1], y2[2], x2, B2);
    double Babs1 = sqrt(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2]);
    double Babs2 = sqrt(B2[0]*B2[0] + B2[1]*B2[1] + B2[2]*B2[2]);

    if (Babs2 == 0){ // if there's no magnetic field or particle is leaving magnetic field, do nothing
        std::fill(spin.begin(), spin.end(), 0); // set all spin components to zero
        return;
    }
    else if (Babs1 == 0 && Babs2 > 0){ // if particle enters magnetic field
//		std::cout << "Entering magnetic field\n";
        p->DoPolarize(x2, y2, 0., flipspin, mc); // if spin flips are allowed, choose random polarization
        spin[0] = y2[7]*B2[0]/Babs2; // set spin (anti-)parallel to field
        spin[1] = y2[7]*B2[1]/Babs2;
        spin[2] = y2[7]*B2[2]/Babs2;
        return;
    }
    else{ // if particle already is in magnetic field, calculate spin-projection on magnetic field
        polarisation = (spin[0]*B1[0] + spin[1]*B1[1] + spin[2]*B1[2])/Babs1/sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
    }

    bool integrate1 = false, integrate2 = false;;
    for (unsigned int i = 0; i < times.size(); i += 2){
        integrate1 |= (x1 >= times[i] && x1 < times[i+1]);
        integrate2 |= (x2 >= times[i] && x2 < times[i+1]);
    }

    if ((integrate1 || integrate2) && (Babs1 < Bmax || Babs2 < Bmax)){ // do spin integration only, if time is in specified range and field is smaller than Bmax
//		if ((!integrate1 && integrate2) || (Babs1 > Bmax && Babs2 < Bmax))
//			std::cout << x1 << "s " << y1[7] - polarisation << " ";

        vector<alglib::spline1dinterpolant> omega_int;
        if (interpolatefields){
            alglib::real_1d_array ts; // set up values for interpolation of precession axis
            vector<alglib::real_1d_array> omega(3);
            omega_int.resize(3);
            const int int_points = 10;
            ts.setlength(int_points + 1);
            for (int i = 0; i < 3; i++)
                omega[i].setlength(int_points + 1);

            for (int i = 0; i <= int_points; i++){ // calculate precession axis at several points along trajectory step
                double t = x1 + i*(x2 - x1)/int_points;
                ts[i] = t;
                p->SpinPrecessionAxis(t, stepper, field, omega[0][i], omega[1][i], omega[2][i]);
            }

            for (int i = 0; i < 3; i++)
                alglib::spline1dbuildcubic(ts, omega[i], omega_int[i]); // interpolate all three components of precession axis
        }


        dense_stepper_type spinstepper = boost::numeric::odeint::make_dense_output(1e-12, 1e-12, stepper_type());
        spinstepper.initialize(spin, x1, std::abs(pi/p->GetGyromagneticRatio()/Babs1)); // initialize integrator with step size = half rotation
        logger->PrintSpin(p, spinstepper, stepper, field);
        unsigned int steps = 0;
        while (true){
            // take an integration step, SpinDerivs contains right-hand side of equation of motion
            spinstepper.do_step([&](const state_type &y, state_type &dydx, const value_type &x){ p->SpinDerivs(y, dydx, x, stepper, field, omega_int); }); //std::bind(&TParticle::SpinDerivs, p.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, stepper, field, omega_int));
            steps++;
            double t = spinstepper.current_time();
            if (t > x2){ // if stepper overshot, calculate end point and stop
                t = x2;
                spinstepper.calc_state(t, spin);
            }
            else
                spin = spinstepper.current_state();

            logger->PrintSpin(p, spinstepper, stepper, field);

            if (t >= x2)
                break;
        }

        // calculate new spin projection
        polarisation = (spin[0]*B2[0] + spin[1]*B2[1] + spin[2]*B2[2])/Babs2/sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
    }
    else if ((Babs1 < Bmax || Babs2 < Bmax)){ // if time outside selected ranges, parallel-transport spin along magnetic field
        if (polarisation*polarisation >= 1){ // catch rounding errors
            spin[0] = 0;
            spin[1] = 0;
        }
        else{
            std::uniform_real_distribution<double> phidist(0, 2.*pi);
            double spinaz = phidist(mc); // random azimuth
            spin[0] = sqrt(1 - polarisation*polarisation)*sin(spinaz);
            spin[1] = sqrt(1 - polarisation*polarisation)*cos(spinaz);
        }
        spin[2] = polarisation;
        RotateVector(&spin[0], B2); // rotate spin vector onto magnetic field vector
    }


    if (Babs2 > Bmax){ // if magnetic field grows above Bmax, collapse spin state to one of the two polarisation states
        p->DoPolarize(x2, y2, polarisation, flipspin, mc);
        spin[0] = B2[0]*y2[7]/Babs2;
        spin[1] = B2[1]*y2[7]/Babs2;
        spin[2] = B2[2]*y2[7]/Babs2;
    }
}
