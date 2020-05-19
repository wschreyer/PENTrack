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

//	cout << "Particle no.: " << p->GetParticleNumber() << " particle type: " << p->GetName() << '\n';
//	cout << "x: " << p->GetFinalState()[0] << "m y: " << p->GetFinalState()[1] << "m z: " << p->GetFinalState()[2]
//		 << "m E: " << p->GetFinalKineticEnergy() << "eV t: " << p->GetFinalTime() << "s tau: " << tau << "s lmax: " << maxtraj << "m\n";

    if (p->GetFinalTime() > tmax)
        throw runtime_error("Tried to start trajectory simulation past the maximum simulation time. Check time settings.");
    // set initial values for integrator
    TStep stepper(p->GetFinalTime(), p->GetFinalPosition(), p->GetFinalVelocity(), p->GetFinalProperTime(), p->GetFinalPolarization(), p->GetFinalPathlength());

    logger->PrintTrack(p, stepper, p->GetFinalSolid(), field);

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

    currentsolids = geom.GetSolids(p->GetFinalTime(), &p->GetFinalPosition()[0]);
    p->SetStopID(ID_UNKNOWN);

    while (p->GetStopID() == ID_UNKNOWN){ // integrate as long as nothing happened to particle
        if (quit.load())
            return;
        stopID ID = stepper.next(
            [&](const double &t, const array<double, 3>& pos, const array<double, 3>& v, const double& polarization){
                return p->EquationOfMotion(t, pos, v, polarization, field);
            }, geom);
        p->SetStopID(ID);
        if (ID != ID_UNKNOWN)
            break;

        if (stepper.IsStepOverCollision()){
            DoHit(p, stepper, mc, geom);
        }
        else{
            if (tau > 0 && stepper.GetProperTime() >= tau){ // if proper time is larger than lifetime
                stepper.SetStepEndToMatch([&](const double &t){ return stepper.GetProperTime(t); }, tau);
                p->SetStopID(ID_DECAYED);
            }
            if (stepper.GetTime() >= tmax){	//If stepsize overshot max simulation time
                stepper.SetStepEnd(tmax);
                p->SetStopID(ID_NOT_FINISH);
            }
            if (stepper.GetPathLength() >= maxtraj){
                stepper.SetStepEndToMatch([&](const double &t){ return stepper.GetPathLength(t); }, maxtraj);
                p->SetStopID(ID_NOT_FINISH);
            }

            DoStep(p, stepper, GetCurrentsolid(), mc, field);
        }

        // take snapshots at certain times
        logger->PrintSnapshot(p, stepper, geom, field);

        bool spincollapse = stepper.IntegrateSpin([&](const double &t, const TStep &stepper){
                                                    return p->SpinPrecessionAxis(t, stepper, field); 
                                                }, SpinTimes, abs(p->GetGyromagneticRatio()*SpinBmax), spininterpolatefields,
                                                [&](const double &t, const array<double, 3> &spin, const TStep &trajectory_stepper){
                                                    logger->PrintSpin(p, t, spin, trajectory_stepper, field);
                                                });
        if (spincollapse)
            p->DoPolarize(stepper, field, flipspin, mc);

        logger->PrintTrack(p, stepper, GetCurrentsolid(), field);
    }

    if (p->GetStopID() == ID_DECAYED){
        p->DoDecay(stepper, mc, geom, field);
    }

    p->SetFinalState(stepper.GetTime(), stepper.GetPosition(), stepper.GetVelocity(), stepper.GetProperTime(), stepper.GetPolarization(), stepper.GetPathLength(),
                    stepper.GetSpin(), stepper.GetSpinIntegrationTime(), stepper.GetSpinPhase(), GetCurrentsolid());
    logger->Print(p, stepper.GetTime(), stepper, geom, field);


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
        double prevtime = stepper.GetTime();
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


