//
// Created by wolfgang on 11/22/19.
//

#include <sstream>
#include <random>
#include <boost/format.hpp>

#include "tracking.h"

using namespace std;

TTracker::TTracker(TConfig& config, std::shared_ptr<TLogger>& alogger){  
  logger = alogger.get();
}

void TTracker::IntegrateParticle(std::unique_ptr<TParticle>& p, const double tmax, std::map<std::string, std::string> &particleconf, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field){

  vector<std::pair<solid, bool> > currentsolids;
  
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

  // set initial values for integrator
  value_type x = p->GetFinalTime();
  if (x > tmax)
    throw runtime_error("Tried to start trajectory simulation past the maximum simulation time. Check time settings.");
  state_type y = p->GetFinalState();

  bool resetintegration = false;

  // std::cout << "print track init ? " << std::endl;
  logger->PrintTrack(p, x, y, x, y, p->GetFinalSpin(), p->GetFinalSolid(), field, true); // add sly
  // std::cout << "track printed init ?? " << std::endl;


  
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

  // dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type()); // original value
  dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(1e-7, 1e-7, stepper_type()); // to DoTo sly
  stepper.initialize(y, x, 10.*MAX_TRACK_DEVIATION/sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])); // initialize stepper with fixed spatial length

  //	progress_display progress(100, cout, ' ' + to_string(particlenumber) + ' ');

  currentsolids = geom.GetSolids(x, &y[0]);
  p->SetStopID(ID_UNKNOWN);

  while (p->GetStopID() == ID_UNKNOWN){ // integrate as long as nothing happened to particle
    if (resetintegration){
      stepper.initialize(y, x, stepper.current_time_step()); // (re-)start integration with last step size
    }
    value_type x1 = x; // save point before next step
    state_type y1 = y;

    try{
      stepper.do_step(std::bind(&TParticle::derivs, p.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, &field));
      x = stepper.current_time();
      y = stepper.current_state();
    }
    catch(...){ // catch Exceptions thrown by odeint
      p->SetStopID(ID_ODEINT_ERROR);
    }

    if (tau > 0 && y[6] > tau){ // if proper time is larger than lifetime
      x = x1 + (x - x1)*(tau - y1[6])/(y[6] - y1[6]); // interpolate decay time in lab frame
      stepper.calc_state(x, y);
    }
    if (x > tmax){	//If stepsize overshot max simulation time
      x = tmax;
      stepper.calc_state(x, y);
    }

    while (x1 < x){ // split integration step in pieces (x1,y1->x2,y2) to reduce chord length, go through all pieces
      if (quit.load())
	return;

      double l2 = pow(y[8] - y1[8], 2); // actual length of step squared
      double d2 = pow(y[0] - y1[0], 2) + pow(y[1] - y1[1], 2) + pow(y[2] - y1[2], 2); // length of straight line between start and end point of step squared
      double dev2 = 0.25*(l2 - d2); // max. possible squared deviation of real path from straight line
      value_type x2 = x;
      state_type y2 = y;
      if (dev2 > MAX_TRACK_DEVIATION*MAX_TRACK_DEVIATION){ // if deviation is larger than MAX_TRACK_DEVIATION
	//				cout << "split " << x - x1 << " " << sqrt(l2) << " " << sqrt(d2) << " " << sqrt(dev2) << "\n";
	x2 = x1 + (x - x1)/ceil(sqrt(dev2)/MAX_TRACK_DEVIATION); // split step to reduce deviation
	stepper.calc_state(x2, y2);
	assert(x2 <= x);
      }
      //			l2 = pow(y2[8] - y1[8], 2);
      //			d2 = pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2);
      //			cout << x2 - x1 << " " << sqrt(l2) << " " << sqrt(d2) << " " << 0.5*sqrt(l2 - d2) << "\n";

      resetintegration = CheckHit(p, x1, y1, x2, y2, stepper, mc, geom, field, currentsolids); // check if particle hit a material boundary or was absorbed between y1 and y2
      if (resetintegration){
	x = x2; // if particle path was changed: reset integration end point
	y = y2;
      }

      x1 = x2;
      y1 = y2;
    }

    // take snapshots at certain times
    logger->PrintSnapshot(p, stepper.previous_time(), stepper.previous_state(), x, y, spin, stepper, geom, field);

    IntegrateSpin(p, spin, stepper, x, y, SpinTimes, field, spininterpolatefields, SpinBmax, mc, flipspin); // calculate spin precession and spin-flip probability

    logger->PrintTrack(p, stepper.previous_time(), stepper.previous_state(), x, y, spin, GetCurrentsolid(currentsolids), field, false);

    //		progress += 100*max(y[6]/tau, max((x - tstart)/(tmax - tstart), y[8]/maxtraj)) - progress.count();

    if (p->GetStopID() == ID_UNKNOWN && y[6] >= tau) // proper time >= tau?
      p->SetStopID(ID_DECAYED);
    else if (p->GetStopID() == ID_UNKNOWN && (x >= tmax || y[8] >= maxtraj)) // time > tmax or trajectory length > max length?
      p->SetStopID(ID_NOT_FINISH);
  }

  //	cout << "Done" << endl;

  if (p->GetStopID() == ID_DECAYED){ // if particle reached its lifetime call TParticle::Decay
    //		cout << "Decayed!\n";
    p->DoDecay(x, y, mc, geom, field);
  }

  p->DoPolarize(x, y, y[7], flipspin, mc);  //added by Niki //project the polarisation aligned or anti-aligned to the B-field at the end of the Simulation
  p->SetFinalState(x, y, spin, GetCurrentsolid(currentsolids));

  // std::cout << "print track? " << std::endl;
  // logger->PrintTrack(p, x, y, x, y, p->GetFinalSpin(), p->GetFinalSolid(), field, true); // add sly
  logger->PrintTrack(p, x, y, x, y, spin, GetCurrentsolid(currentsolids), field, true); // add sly
  // std::cout << "track printed ?? " << std::endl;

  logger->Print(p, x, y, spin, geom, field);


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


bool TTracker::CheckHit(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
			const dense_stepper_type &stepper, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field, vector<std::pair<solid, bool> > &currentsolids){
  if (!geom.CheckSegment(&y1[0], &y2[0])){ // check if start point is inside bounding box of the simulation geometry
    //    printf("\nParticle has hit outer boundaries: Stopping it! t=%g x=%g y=%g z=%g\n",x2,y2[0],y2[1],y2[2]);
    p->SetStopID(ID_HIT_BOUNDARIES);
    return true;
  }
  if (x2 == x1)
    return false;

  solid currentsolid = GetCurrentsolid(currentsolids);

  multimap<TCollision, bool> colls;
  bool collfound = false;
  try{
    collfound = geom.GetCollisions(x1, &y1[0], x2, &y2[0], colls);
  }
  catch(...){
    p->SetStopID(ID_CGAL_ERROR);
    return true;
  }

  if (collfound){	// if there is a collision with a wall
    value_type xc1 = x1, xc2 = x2;
    //    for (auto c: colls)
    //      cout << x1 << " " << x2 - x1 << " " << c.first.distnormal << " " << c.first.s << " " << c.first.ID << endl;
    state_type yc1 = y1, yc2 = y2;
    if (iterate_collision(xc1, yc1, xc2, yc2, colls.begin()->first, stepper, geom)){
      if (xc1 > x1 && DoStep(p, x1, y1, xc1, yc1, stepper, currentsolid, mc, field)){
	x2 = xc1;
	y2 = yc1;
	return true;
      }

      if (DoHit(p, xc1, yc1, xc2, yc2, stepper, mc, geom, currentsolids)){
	x2 = xc2;
	y2 = yc2;
	return true;
      }

      if (CheckHit(p, xc2, yc2, x2, y2, stepper, mc, geom, field, currentsolids)){
	return true;
      }
    }

  }
  else{ // if there was no collision: just check for absorption in solid with highest priority
    return DoStep(p, x1, y1, x2, y2, stepper, currentsolid, mc, field); // sly currentsolid
  }
  return false;
}

bool TTracker::iterate_collision(value_type &x1, state_type &y1, value_type &x2, state_type &y2,
				 const TCollision &coll, const dense_stepper_type &stepper, const TGeometry &geom, unsigned int iteration){
  if (pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2) < REFLECT_TOLERANCE*REFLECT_TOLERANCE){
    return true; // successfully iterated collision point
  }
  if (x2 - x1 < 4*(x1 + x2)*numeric_limits<value_type>::epsilon()){
    cout << "Collision point iteration limited by numerical precision.\n";
    return true;
  }
  if (iteration >= 100){
    cout << "Collision point iteration reached max. iterations. " << x1 << " " << x2 - x1 << " " << coll.distnormal << " " << coll.s << "\n";
    return true;
  }

  //  value_type xc = x1 + (x2 - x1)*coll.s;
  //  if (xc == x1 || xc == x2)
  value_type xc = x1 + (x2 - x1)*0.5;
  state_type yc(STATE_VARIABLES);
  stepper.calc_state(xc, yc);
  multimap<TCollision, bool> colls;
  if (geom.GetCollisions(x1, &y1[0], xc, &yc[0], colls)){ // if collision in first segment, further iterate
    //    cout << "1 " << x1 << " " << xc1 - x1 << endl;
    if (iterate_collision(x1, y1, xc, yc, colls.begin()->first, stepper, geom, iteration + 1)){
      x2 = xc;
      y2 = yc;
      return true; // if successfully iterated
    }
  }
  if (geom.GetCollisions(xc, &yc[0], x2, &y2[0], colls)){ // if collision in second segment, further iterate
    //    cout << "2 " << xc1 << " " << xc2 - xc1 << endl;
    if (iterate_collision(xc, yc, x2, y2, colls.begin()->first, stepper, geom, iteration + 1)){
      x1 = xc;
      y1 = yc;
      return true; // if successfully iterated
    }
  }
  return false; // if no iteration successful, actual trajectory curves past surface
}


bool TTracker::DoStep(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		      const dense_stepper_type &stepper, const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field) {
  value_type x2temp = x2;
  state_type y2temp = y2;
  p->DoStep(x1, y1, x2, y2, stepper, currentsolid, mc, field);
  if (x2temp == x2 && y2temp == y2)
    return false;
  else{
    // std::cout << "Do particle " << p->GetParticleNumber() << " has altered path" << std::endl;
    return true;
  }
}

bool TTracker::DoHit(const std::unique_ptr<TParticle>& p, const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
		     const dense_stepper_type &stepper, TMCGenerator &mc, const TGeometry &geom, vector<pair<solid, bool> > &currentsolids) {
  bool trajectoryaltered = false, traversed = true;

  multimap<TCollision, bool> colls;
  if (!geom.GetCollisions(x1, &y1[0], x2, &y2[0], colls))
    throw std::runtime_error("Called DoHit for a trajectory segment that does not contain a collision!");

  vector<pair<solid, bool> > newsolids = currentsolids;
  for (auto coll: colls){
    //    cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << coll.first.ID << endl;
    solid sld = geom.GetSolid(coll.first.ID);
    auto foundsld = find_if(newsolids.begin(), newsolids.end(), [&sld](const std::pair<solid, bool> s){ return s.first.ID == sld.ID; });
    if (coll.first.distnormal < 0){ // if entering solid
      if (foundsld != newsolids.end()){ // if solid has been entered before (self-intersecting surface)
	//	cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << sld.name << endl;
	if (coll.first.s > 0) // if collision happened right at the start of the step it is likely that the hit solid was already added to the list in the previous step and we will ignore this one
	  newsolids.push_back(make_pair(sld, foundsld->second)); // add additional entry to list, with ignore state as on first entry
      }
      else
	newsolids.push_back(make_pair(sld, coll.second)); // add solid to list
    }
    else if (coll.first.distnormal > 0){ // if leaving solid
      if (foundsld == newsolids.end()){ // if solid was not entered before something went wrong
	if (coll.first.s > 0){ // if collision happened right at the start of the step it is likely that the hit solid was already removed from the list in the previous step and this is not an error
	  //	  cout << x1 << " " << x2 - x1 << " " << coll.first.distnormal << " " << coll.first.s << " " << sld.name << endl;
	  //          throw runtime_error((boost::format("Particle inside '%1%' which it did not enter before!") % sld.name).str());
	  cout << "Particle inside solid " << sld.name << " which it did not enter before. Stopping it!\n";
	  p->SetStopID(ID_GEOMETRY_ERROR);
	  return true;
	}
      }
      else
	newsolids.erase(foundsld); // remove solid from list
    }
    else{
      //      throw runtime_error("Particle crossed surface with parallel track!");
      cout << "Particle crossed surface with parallel track. Stopping it!\n";
      p->SetStopID(ID_GEOMETRY_ERROR);
      return true;
    }
  }

  solid leaving = GetCurrentsolid(currentsolids); // particle can only leave highest-priority solid
  solid entering = geom.defaultsolid;
  for (auto sld: newsolids){
    if (!sld.second && sld.first.ID > entering.ID)
      entering = sld.first;
  }
  //  cout << "Leaving " << leaving.name << ", entering " << entering.name << '\n';
  if (leaving.ID != entering.ID){ // if the particle actually traversed a material interface
    auto coll = find_if(colls.begin(), colls.end(), [&leaving, &entering](const pair<TCollision, bool> &c){ return c.first.ID == leaving.ID or (c.first.ID == entering.ID and not c.second); });
    if (coll == colls.end())
      throw std::runtime_error((boost::format("Did not find collision going from %5% to %6%! t=%1%s, x=%2%, y=%3%, z=%4%") % x1 % y1[0] % y1[1] % y1[2] % leaving.ID % entering.ID).str());
    value_type x2temp = x2;
    state_type y2temp = y2;
    p->DoHit(x1, y1, x2, y2, coll->first.normal, leaving, entering, mc); // do particle specific things
    if (x2temp == x2 && y2temp == y2){ // if end point of step was not modified
      trajectoryaltered = false;
      traversed = true;
    }
    else{
      trajectoryaltered = true;
      if (x2 == x2temp && y2[0] == y2temp[0] && y2[1] == y2temp[1] && y2[2] == y2temp[2]){
	traversed = true;
      }
      else if (x2 == x1 && y2[0] == y1[0] && y2[1] == y1[1] && y2[2] == y1[2]){
	traversed = false;
      }
      else{
	throw std::runtime_error("OnHit routine returned inconsistent position. That should not happen!");
      }
    }

    logger->PrintHit(p, x1, y1, y2, coll->first.normal, leaving, entering); // print collision to file if requested
  }

  if (traversed){
    currentsolids = newsolids; // if surface was traversed (even if it was  physically ignored) replace current solids with list of new solids
  }

  if (trajectoryaltered || p->GetStopID() != ID_UNKNOWN)
    return true;

  return false;
}

const solid& TTracker::GetCurrentsolid(vector<pair<solid, bool> > currentsolids) const{
  auto sld = max_element(currentsolids.begin(), currentsolids.end(), [](const pair<solid, bool> &s1, const pair<solid, bool> &s2){ return s1.second || (!s2.second && s1.first.ID < s2.first.ID); });
  return sld->first;
}


// const solid& TTracker::GetCurrentsolid() const{
//   auto sld = max_element(currentsolids.begin(), currentsolids.end(), [](const pair<solid, bool> &s1, const pair<solid, bool> &s2){ return s1.second || (!s2.second && s1.first.ID < s2.first.ID); });
//   return sld->first;
// }


void TTracker::IntegrateSpin(const std::unique_ptr<TParticle>& p, state_type &spin, const dense_stepper_type &stepper,  const double x2, state_type &y2, const std::vector<double> &times, const TFieldManager &field, const bool interpolatefields, const double Bmax, TMCGenerator &mc, const bool flipspin) const{
  value_type x1 = stepper.previous_time();
  if (p->GetGyromagneticRatio() == 0 || x1 == x2)
    return;

  state_type y1 = stepper.previous_state();
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


    dense_stepper_type spinstepper = boost::numeric::odeint::make_dense_output(1e-7, 1e-7, stepper_type());
    // dense_stepper_type spinstepper = boost::numeric::odeint::make_dense_output(1e-12, 1e-12, stepper_type()); // original
    spinstepper.initialize(spin, x1, std::abs(pi/p->GetGyromagneticRatio()/Babs1)); // initialize integrator with step size = half rotation
    // spinstepper.initialize(spin, x1, std::abs(pi/p->GetGyromagneticRatio()/Babs1/36)); // initialize integrator with step size = half rotation ToDo sly
    logger->PrintSpin(p, x1, spinstepper, stepper, field);
    unsigned int steps = 0;
    while (true){
      if (quit.load())
	return;

      // take an integration step, SpinDerivs contains right-hand side of equation of motion
      spinstepper.do_step(std::bind(&TParticle::SpinDerivs, p.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, stepper, &field, omega_int));
      steps++;
      double t = spinstepper.current_time();
      if (t > x2){ // if stepper overshot, calculate end point and stop
	t = x2;
	spinstepper.calc_state(t, spin);
      }
      else
	spin = spinstepper.current_state();

      logger->PrintSpin(p, t, spinstepper, stepper, field);

      if (t >= x2)
	break;

      // std::cout << "time t = " << t << " Polarisation proj = " << polarisation << " Sx = " << spin[0] << " Sy = " << spin[1] << " Sz = " << spin[2] << std::endl;
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


  // if (Babs2 > Bmax){ // original // if magnetic field grows above Bmax, collapse spin state to one of the two polarisation states
  if (flipspin){ // add sly
    p->DoPolarization(x2, y2, polarisation, flipspin, mc); //added by Niki
    if (Babs2 > Bmax){
      p->DoPolarize(x2, y2, polarisation, flipspin, mc);
      spin[0] = B2[0]*y2[7]/Babs2;
      spin[1] = B2[1]*y2[7]/Babs2;
      spin[2] = B2[2]*y2[7]/Babs2;
    }
  }
}
