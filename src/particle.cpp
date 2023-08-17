/**
 * \file
 * Particle base class definition.
 */

#include <vector>
#include <algorithm>
#include <iomanip>

#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

#include "particle.h"

using namespace std;

double TParticle::GetInitialTotalEnergy(const TGeometry &geom, const TFieldManager &field) const{
  return GetKineticEnergy(&ystart[3]) + GetPotentialEnergy(tstart, ystart, field, geom.GetSolid(tstart, &ystart[0]));
}


double TParticle::GetFinalTotalEnergy(const TGeometry &geom, const TFieldManager &field) const{
  return GetKineticEnergy(&yend[3]) + GetPotentialEnergy(tend, yend, field, geom.GetSolid(tend, &yend[0]));
}


double TParticle::GetInitialKineticEnergy() const{
  return GetKineticEnergy(&ystart[3]);
}


double TParticle::GetFinalKineticEnergy() const{
  return GetKineticEnergy(&yend[3]);
}

TParticle::TParticle(const char *aname, const  double qq, const long double mm, const long double mumu, const long double agamma, const int number,
		     const double t, const double x, const double y, const double z, const double E, const double phi, const double theta, const int polarisation,
		     const double spinprojection, TMCGenerator &amc, const TGeometry &geometry, const TFieldManager &afield)
  : name(aname), q(qq), m(mm), mu(mumu), gamma(agamma), particlenumber(number), ID(ID_UNKNOWN),
    tstart(t), tend(t), Hmax(0), Nhit(0), Nspinflip(0), noflipprob(1), Nstep(0){

  // std::cout << "TParticle" << std::endl;

  // for small velocities Ekin/m is very small and the relativstic claculation beta^2 = 1 - 1/gamma^2 gives large round-off errors
  // the round-off error can be estimated as 2*epsilon
  // if a series expansion to order O(Eoverm^5) has a smaller error then the round-off error, we will use the series expansion
  double Eoverm = E/m/c_0/c_0; // gamma - 1
  double beta2;
  //cout << pow(Eoverm, 6) << " " << ((((6.0*Eoverm - 5.0)*Eoverm + 4.0)*Eoverm - 3.0)*Eoverm + 2.0)*Eoverm << " " << 1.0 - 1.0/(Eoverm + 1)/(Eoverm + 1) << "\n";
  if (pow(Eoverm, 6) < 2*numeric_limits<double>::epsilon()) // if error in series expansion smaller than round-off error
    beta2 = ((((6.0*Eoverm - 5.0)*Eoverm + 4.0)*Eoverm - 3.0)*Eoverm + 2.0)*Eoverm; // use series expansion
  else
    beta2 = 1.0 - 1.0/(Eoverm + 1)/(Eoverm + 1); // relativstic beta^2
  double vstart = c_0*sqrt(beta2);

  if (polarisation < -1 || polarisation > 1)
    throw std::runtime_error("Polarisation has to be between -1 and 1");

  ystart.resize(STATE_VARIABLES, 0);
  ystart[0] = x; // position
  ystart[1] = y;
  ystart[2] = z;
  ystart[3] = vstart*cos(phi)*sin(theta); // velocity
  ystart[4] = vstart*sin(phi)*sin(theta);
  ystart[5] = vstart*cos(theta);
  ystart[6] = 0; // proper time
  ystart[7] = polarisation; // spin polarization used for trajectory tracking in magnetic field
  ystart[8] = 0;
  yend = ystart;

  spinstart.resize(SPIN_STATE_VARIABLES, 0);

  double B[3];
  afield.BField(x, y, z, t, B);
  double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  if (Babs > 0){
    std::uniform_real_distribution<double> phidist(0., 2.*pi);
    double spinaz = phidist(amc); // random azimuth angle of spin
    spinstart[0] = sqrt(1 - spinprojection*spinprojection)*sin(spinaz); // set initial spin vector
    spinstart[1] = sqrt(1 - spinprojection*spinprojection)*cos(spinaz);
    spinstart[2] = spinprojection;
    RotateVector(&spinstart[0], B); // rotate initial spin vector such that its z-component is parallel to magnetic field

    spinstart[3] = 0; // initial integration time is zero
    spinstart[4] = 0; // initial phase is zero
  }

  spinend = spinstart;

  solidend = solidstart = geometry.GetSolid(t, &ystart[0]); // set to solid with highest priority
  Hmax = GetInitialTotalEnergy(geometry, afield);
}



void TParticle::derivs(const state_type &y, state_type &dydx, const value_type x, const TFieldManager *field) const{
  double B[3], dBidxj[3][3], E[3], V; // magnetic/electric field and electric potential in lab frame
  if (q != 0 || (mu != 0 && y[7] != 0)) // if particle has charge or magnetic moment, calculate magnetic field
    field->BField(y[0],y[1],y[2], x, B, dBidxj);
  if (q != 0) // if particle has charge calculate electric field
    field->EField(y[0],y[1],y[2], x, V, E);


  // double reltrace = (dBidxj[0][0] + dBidxj[1][1] + dBidxj[2][2])/sqrt(dBidxj[0][0]*dBidxj[0][0] + dBidxj[1][1]*dBidxj[1][1] + dBidxj[2][2]*dBidxj[2][2]);
  // double relrotxy = (dBidxj[0][1] - dBidxj[1][0])/(dBidxj[0][1] + dBidxj[1][0]);
  // double relrotxz = (dBidxj[0][2] - dBidxj[2][0])/(dBidxj[0][2] + dBidxj[2][0]);
  // double relrotyz = (dBidxj[1][2] - dBidxj[2][1])/(dBidxj[1][2] + dBidxj[2][1]);

  // double trace = dBidxj[0][0] + dBidxj[1][1] + dBidxj[2][2];
  // double rotxy = dBidxj[0][1] - dBidxj[1][0];
  // double rotxz = dBidxj[0][2] - dBidxj[2][0];
  // double rotyz = dBidxj[1][2] - dBidxj[2][1];

  // if(abs(reltrace) > 0.01){ //} || abs(relrotxy > 0.01 || relrotxz > 0.01 || relrotyz > 0.01){
  // // if(abs(trace) > 1){ // || abs(rotxy) > 0.01 || abs(rotxz) > 0.01 || abs(rotyz) > 0.01){
  //   std::cout << " \ndBidxj, x=" << y[0] << ", r=" << sqrt(y[1]*y[1] + y[2]*y[2]) << ", y=" << y[1] << ", z=" << y[2] << std::endl;
  //   for(int i=0; i<3;++i){
  //     for(int j=0; j<3; ++j){
  // 	std::cout << dBidxj[i][j] << ", ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << " reltrace = " << reltrace << ", relrotxy = " << relrotxy << ", relrotxz = " << relrotxz << ", relrotyz" << relrotyz << std::endl;
  //   std::cout << " trace = " << trace << ", rotxy = " << rotxy << ", rotxz = " << rotxz << ", rotyz" << rotyz << std::endl;
  // }

  // for(int i=0; i<3; ++i){
  //   dBidxj[i][i] = dBidxj[i][i] - trace/3;
  //     for(int j=i+1; j<3; ++j){
  // 	dBidxj[i][j] = (dBidxj[i][j] + dBidxj[j][i])/2;
  // 	dBidxj[j][i] = dBidxj[i][j];
  //     }
  // }
  
  EquationOfMotion(y, dydx, x, B, dBidxj, E);

}

void TParticle::EquationOfMotion(const state_type &y, state_type &dydx, const value_type x, const double B[3], const double dBidxj[3][3], const double E[3]) const{
  dydx[0] = y[3]; // time derivatives of position = velocity
  dydx[1] = y[4];
  dydx[2] = y[5];

  value_type F[3] = {0,0,0}; // Force in lab frame
  F[2] += -gravconst*m*ele_e; // add gravitation to force
  if (q != 0){
    F[0] += q*(E[0] + y[4]*B[2] - y[5]*B[1]); // add Lorentz-force
    F[1] += q*(E[1] + y[5]*B[0] - y[3]*B[2]);
    F[2] += q*(E[2] + y[3]*B[1] - y[4]*B[0]);
  }
  if (mu != 0 && y[7] != 0 && (B[0] != 0 || B[1] != 0 || B[2] != 0)){
    double Babs = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    double dBdxi[3] = {	(B[0]*dBidxj[0][0] + B[1]*dBidxj[1][0] + B[2]*dBidxj[2][0])/Babs,
			(B[0]*dBidxj[0][1] + B[1]*dBidxj[1][1] + B[2]*dBidxj[2][1])/Babs,
			(B[0]*dBidxj[0][2] + B[1]*dBidxj[1][2] + B[2]*dBidxj[2][2])/Babs}; // derivatives of |B|
    F[0] += y[7]*mu*dBdxi[0]; // add force on magnetic dipole moment
    F[1] += y[7]*mu*dBdxi[1];
    F[2] += y[7]*mu*dBdxi[2];
  }
  double v2 = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  value_type inversegamma = sqrt(1 - v2/(c_0*c_0)); // relativstic factor 1/gamma
  dydx[3] = inversegamma/m/ele_e*(F[0] - (y[3]*y[3]*F[0] + y[3]*y[4]*F[1] + y[3]*y[5]*F[2])/c_0/c_0); // general relativstic equation of motion
  dydx[4] = inversegamma/m/ele_e*(F[1] - (y[4]*y[3]*F[0] + y[4]*y[4]*F[1] + y[4]*y[5]*F[2])/c_0/c_0); // dv/dt = 1/gamma/m*(F - v * v^T * F / c^2)
  dydx[5] = inversegamma/m/ele_e*(F[2] - (y[5]*y[3]*F[0] + y[5]*y[4]*F[1] + y[5]*y[5]*F[2])/c_0/c_0);

  dydx[6] = inversegamma; // derivative of proper time is 1/gamma
  dydx[7] = 0; // polarisaton does not change
  dydx[8] = sqrt(v2); // derivative of path length is abs(velocity)
}





void TParticle::SpinPrecessionAxis(const double t, const dense_stepper_type &stepper, const TFieldManager &field, double &Omegax, double &Omegay, double &Omegaz) const{
  double B[3], dBidxj[3][3], V, E[3];
  state_type y(STATE_VARIABLES), dydt(STATE_VARIABLES);
  stepper.calc_state(t, y); // calculate particle state at time t
  field.BField(y[0], y[1], y[2], t, B, dBidxj);
  field.EField(y[0], y[1], y[2], t, V, E);
  EquationOfMotion(y, dydt, t, B, dBidxj, E); // calculate velocity and acceleration required for vxE effect and Thomas precession
  SpinPrecessionAxis(t, B, E, dydt, Omegax, Omegay, Omegaz); // calculate precession axis
}

void TParticle::SpinPrecessionAxis(const double t, const double B[3], const double E[3], const state_type &dydt, double &Omegax, double &Omegay, double &Omegaz) const{
  double v2 = dydt[0]*dydt[0] + dydt[1]*dydt[1] + dydt[2]*dydt[2];
  double gamma_rel = 1./sqrt(1. - v2/c_0/c_0);
  double Bdotv = B[0]*dydt[0] + B[1]*dydt[1] + B[2]*dydt[2];
  double Bparallel[3];
  for (int j = 0; j < 3; j++)
    Bparallel[j] = Bdotv*dydt[j]/v2; // magnetic-field component parallel to velocity

  // spin precession axis due to relativistically distorted magnetic field, omega_B = -gyro/gamma * ( (1 - gamma)*(v.B)*v/v^2 + gamma*B - gamma*(v x E)/c^2 )
  double OmegaB[3];
  OmegaB[0] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[0] + gamma_rel*B[0] - gamma_rel*(dydt[1]*E[2] - dydt[2]*E[1])/c_0/c_0);
  OmegaB[1] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[1] + gamma_rel*B[1] - gamma_rel*(dydt[2]*E[0] - dydt[0]*E[2])/c_0/c_0);
  OmegaB[2] = -gamma/gamma_rel * ((1 - gamma_rel)*Bparallel[2] + gamma_rel*B[2] - gamma_rel*(dydt[0]*E[1] - dydt[1]*E[0])/c_0/c_0);

  // Thomas precession in lab frame, omega_T = gamma^2/(gamma + 1)/c^2*(dv/dt x v)
  double OmegaT[3] = {0,0,0};
  OmegaT[0] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[4]*dydt[2] - dydt[5]*dydt[1])/c_0/c_0;
  OmegaT[1] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[5]*dydt[0] - dydt[3]*dydt[2])/c_0/c_0;
  OmegaT[2] = gamma_rel*gamma_rel/(gamma_rel + 1)*(dydt[3]*dydt[1] - dydt[4]*dydt[0])/c_0/c_0;

  // Total spin precession is sum of magnetic-field precession and Thomas precession
  Omegax = OmegaB[0] + OmegaT[0];
  Omegay = OmegaB[1] + OmegaT[1];
  Omegaz = OmegaB[2] + OmegaT[2];
}


void TParticle::SpinDerivs(const state_type &y, state_type &dydx, const value_type x, const dense_stepper_type &stepper, const TFieldManager *field, const std::vector<alglib::spline1dinterpolant> &omega) const{
  double omegax, omegay, omegaz;
  if (omega.size() == 3){ // if interpolator exists, use it
    std::cout << "inside Spinderivs interpolator !" << std::endl;
    omegax = alglib::spline1dcalc(omega[0], x);
    omegay = alglib::spline1dcalc(omega[1], x);
    omegaz = alglib::spline1dcalc(omega[2], x);
  }
  else
    SpinPrecessionAxis(x, stepper, *field, omegax, omegay, omegaz); // else calculate precession axis directly

  dydx[0] = omegay*y[2] - omegaz*y[1]; // dS/dt = W x S
  dydx[1] = omegaz*y[0] - omegax*y[2];
  dydx[2] = omegax*y[1] - omegay*y[0];
  dydx[3] = 1.; // integrate time
  dydx[4] = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz); // integrate precession phase
}


void TParticle::DoStep(const value_type x1, const state_type &y1, value_type &x2, state_type &y2, const dense_stepper_type &stepper,
                       const solid &currentsolid, TMCGenerator &mc, const TFieldManager &field){
  state_type y2temp = y2;
  vector<TParticle*> secs;
  OnStep(x1, y1, x2, y2, stepper, currentsolid, mc, ID, secs);
  for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
  Hmax = max(GetKineticEnergy(&y2[3]) + GetPotentialEnergy(x2, y2, field, currentsolid), Hmax);
  if (y2temp[7] != y2[7])
    Nspinflip++;
  Nstep++;
}

void TParticle::DoHit(const value_type x1, const state_type &y1, value_type &x2, state_type &y2,
                      const double normal[3], const solid &leaving, const solid &entering, TMCGenerator &mc){
  state_type y2temp = y2;
  vector<TParticle*> secs;
  OnHit(x1, y1, x2, y2, normal, leaving, entering, mc, ID, secs); // do particle specific things
  for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
  if (y2temp[7] != y2[7])
    Nspinflip++;
  Nhit++;
}

void TParticle::DoDecay(const double t, const state_type &y, TMCGenerator &mc, const TGeometry &geom, const TFieldManager &field){
  vector<TParticle*> secs;
  Decay(t, y, mc, geom, field, secs);
  for (auto s: secs) secondaries.push_back(unique_ptr<TParticle>(s));
}

void TParticle::DoPolarization(const double t, state_type &y, const double polarization, const bool flipspin, TMCGenerator &mc){  //added by Niki
  y[7] = polarization ;
}

void TParticle::DoPolarize(const double t, state_type &y, const double polarization, const bool flipspin, TMCGenerator &mc){
  double flipprob = 0.5*(1 - y[7]*polarization);
  noflipprob *= 1. - flipprob;
  if (flipspin){
    double prevpol = y[7];
    y[7] = polarization_distribution<double>(polarization)(mc);
    if (y[7] != prevpol)
      ++Nspinflip;
  }
}


double TParticle::GetKineticEnergy(const value_type v[3]) const{
  double beta2 = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/c_0/c_0;
  double gammaminusone;

  // for small velocities the relativistic gamma factor = 1/sqrt(1 - beta^2) is very close to one, givin large round-off errors when calculatin E = m*c^2*(gamma - 1)
  // calculating gamma - 1 has a round-off error of epsilon
  // if a series expansion to order O(beta^8) has a smaller error than epsilon, we will use the series expansion

  //cout << pow(beta2, 5) << " " << (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2 << " " << 1.0/sqrt(1.0 - beta2) - 1.0 << "\n";
  if (pow(beta2, 5) < numeric_limits<value_type>::epsilon()) // if error in series expansion O(beta^10) is smaller than rounding error in gamma factor
    gammaminusone = (0.5 + (3.0/8.0 + (5.0/16.0 + 35.0/128.0*beta2)*beta2)*beta2)*beta2; // use series expansion for energy calculation with small beta
  else
    gammaminusone = 1.0/sqrt(1.0 - beta2) - 1.0; // use relativistic formula for larger beta
  return c_0*c_0*m*gammaminusone;
}


double TParticle::GetPotentialEnergy(const value_type t, const state_type &y, const TFieldManager &field, const solid &sld) const{
  value_type result = 0;
  if (q != 0 || mu != 0){
    double B[3], E[3], V;
    if (mu != 0){
      field.BField(y[0],y[1],y[2],t,B);
      result += -y[7]*mu/ele_e*sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    }
    if (q != 0){
      field.EField(y[0],y[1],y[2],t,V,E);
      result += q/ele_e*V;
    }
  }
  result += m*gravconst*y[2];
  return result;
}

void TParticle::SetFinalState(const value_type& x, const state_type& y, const state_type& spin, const solid& sld) {
  tend = x;
  yend = y;
  spinend = spin;
  solidend = sld;
}
