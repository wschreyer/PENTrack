#include "globals.h"

#include "interpolation.h"

#include "stepper.h"

using namespace std;

const static int STATE_VARIABLES = 8;
const static int SPIN_STATE_VARIABLES = 5;

static const double MAX_TRACK_DEVIATION = 0.001; ///< max deviation of actual trajectory from straight line between start and end points of a step used for geometry-intersection test. If deviation is larger, the step will be split
static const int MAX_COLLISION_ITERATIONS = 100;

TStep::TStep(const double &t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const double &propertime, const int &polarization, const double &pathlength){
    xstart = xend = t;
    state_type y = {pos[0], pos[1], pos[2], v[0], v[1], v[2], propertime, pathlength};
    if (abs(polarization) != 1)
        throw std::runtime_error("Tried to set invalid polarization. Polarization can only be +/-1");
    polstart = polend = polarization;
    stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type());
    stepper.initialize(y, t, 10.*MAX_TRACK_DEVIATION/sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]));
    yend = y;
    next_step_over_collision = false;
    step_over_collision = false;

    spinstepper = boost::numeric::odeint::make_dense_output(1e-12, 1e-12, stepper_type());
    spinstart.resize(SPIN_STATE_VARIABLES, 0.);
    spinend = spinstart;
    noflipprob = 1.;
}

TStep::state_type TStep::GetState(const value_type &x) const {
    if (x < xstart or x > xend)
        throw std::runtime_error((boost::format("Requested particle state outside of step range! %.20f %.20f %.20f") % xstart % x % xend).str());
    if (x == stepper.current_time())
        return stepper.current_state();
    else if (x == stepper.previous_time())
        return stepper.previous_state();
    state_type y(STATE_VARIABLES);
    stepper.calc_state(x, y);
    return y;
}


stopID TStep::next(const std::function<std::array<double, 3>(const double&, const std::array<double, 3>&, const std::array<double, 3>&, const double&)>& EquationOfMotion, const TGeometry &geometry){
    if (yend != GetState(xend) or polend != polstart){
        stepper.initialize(yend, xend, stepper.current_time_step());
    }

    step_over_collision = false;
    state_type ystart;
    polstart = polend;
    spinstart = spinend;
    spinend = {0., 0., 0., spinstart[3], spinstart[4]};
    if (stepper.current_time() == xend){
        try{
            stepper.do_step([&](const state_type& y, state_type &dydt, const double &x){
                auto a = EquationOfMotion(x, {y[0], y[1], y[2]}, {y[3], y[4], y[5]}, polstart);
                double v2 = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
                value_type inversegamma = sqrt(1 - v2/(c_0*c_0));
                dydt = {y[3], y[4], y[5], a[0], a[1], a[2], inversegamma, sqrt(v2)};
            });
        }catch(...){
            return ID_ODEINT_ERROR;
        }
        xstart = stepper.previous_time();
        ystart = stepper.previous_state();
        xend = stepper.current_time();
        yend = stepper.current_state();
    }
    else if (next_step_over_collision){
        xstart = xend;
        ystart = yend;
        xend = next_step;
        stepper.calc_state(xend, yend);
        step_over_collision = true;
        next_step_over_collision = false;
        return ID_UNKNOWN;
    }
    else{
        xstart = xend;
        ystart = yend;
        xend = stepper.current_time();
        yend = stepper.current_state();
    }

    collisions.clear();
    if (not geometry.CheckSegment(&ystart[0], &yend[0]))
        return ID_HIT_BOUNDARIES;

    double l2 = pow(yend[7] - ystart[7], 2); // actual length of step squared
    double d2 = pow(yend[0] - ystart[0], 2) + pow(yend[1] - ystart[1], 2) + pow(yend[2] - ystart[2], 2); // length of straight line between start and end point of step squared
    double dev2 = 0.25*(l2 - d2); // max. possible squared deviation of real path from straight line
    int nsubsteps = dev2 > 0 ? static_cast<int>(ceil(sqrt(dev2)/MAX_TRACK_DEVIATION)) : 1;
    for (int i = 0; i < nsubsteps; ++i){
        double x1 = xstart + i*(xend - xstart)/nsubsteps;
        double x2 = (i + 1 == nsubsteps) ? xend : xstart + (i+1)*(xend-xstart)/nsubsteps;
        if (iterate_collision(x1, GetState(x1), x2, GetState(x2), geometry, 0))
            break;
    }
    return ID_UNKNOWN;
}


bool TStep::iterate_collision(const value_type &x1, const state_type &y1, const value_type &x2, const state_type &y2, const TGeometry &geometry, int iteration){
    multimap<TCollision, bool> colls;
    if (not geometry.GetCollisions(x1, &y1[0], x2, &y2[0], colls))
        return false;

    bool found = false;
    if (pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2) < REFLECT_TOLERANCE*REFLECT_TOLERANCE){
        found = true; // successfully iterated collision point
    }
    if (x2 - x1 < 4*(x1 + x2)*numeric_limits<value_type>::epsilon()){
        cout << "Collision point iteration limited by numerical precision. " << x1 << " " << x2 - x1 << " " << colls.begin()->first.distnormal << " " << colls.begin()->first.s << "\n";
        found = true;
    }
    if (iteration == MAX_COLLISION_ITERATIONS){
        cout << "Collision point iteration reached max. iterations. " << x1 << " " << x2 - x1 << " " << colls.begin()->first.distnormal << " " << colls.begin()->first.s << "\n";
        found = true;
    }
    if (found){
        collisions = colls;
        xend = x1;
        yend = y1;
        next_step = x2;
        next_step_over_collision = true;
        return true;
    }

    double xc = 0.5*(x1 + x2);
    state_type yc = GetState(xc);
    if (iterate_collision(x1, y1, xc, yc, geometry, iteration + 1)){ // check if collision in first segment
        return true;
    }
    if (iterate_collision(xc, yc, x2, y2, geometry, iteration + 1)){ // check if collision in second segment
        return true;
    }

    return false;
}

void TStep::SetStepEndToMatch(const std::function<double(const double& t)> &func, const double &match){
    value_type f1 = func(xstart);
    value_type f2 = func(xend);
    if (match < f1 or match > f2)
        throw std::runtime_error("Tried to find state outside of step range");
    value_type t = xstart + (xend - xstart)*(match - f1)/(f2 - f1);
    SetStepEnd(t);
}



double TStep::GetSpinPolarization(const std::array<double, 3>& axis) const{
    double a = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    double s = sqrt(spinend[0]*spinend[0] + spinend[1]*spinend[1] + spinend[2]*spinend[2]);
    return (spinend[0]*axis[0] + spinend[1]*axis[1] + spinend[2]*axis[2])/a/s;
}

bool TStep::IntegrateSpin(const std::function<std::array<double, 3>(const double& t, const TStep &stepper)>& SpinPrecessionAxis,
                            const std::vector<double>& times, const double& maxfreq, const bool& interpolate,
                            const std::function<void(const double&, const std::array<double, 3>&, const TStep&)> &Observer){
    bool integrate1 = false, integrate2 = false;
    for (unsigned int i = 0; i < times.size(); i += 2){
        integrate1 |= (xstart >= times[i] && xstart < times[i+1]);
        integrate2 |= (xend >= times[i] && xend < times[i+1]);
    }
    auto Omega1 = SpinPrecessionAxis(xstart, *this);
    double w1 = sqrt(Omega1[0]*Omega1[0] + Omega1[1]*Omega1[1] + Omega1[2]*Omega1[2]);
    auto Omega2 = SpinPrecessionAxis(xend, *this);
    double w2 = sqrt(Omega2[0]*Omega2[0] + Omega2[1]*Omega2[1] + Omega2[2]*Omega2[2]);
    integrate1 &= w1 < maxfreq;
    integrate2 &= w2 < maxfreq;
    if (not integrate1 and not integrate2)
        return false;
    cout << " tend " << xend << " " << integrate2 << " ";
    copy(times.begin(), times.end(), ostream_iterator<double>(cout, " "));
    if ((not integrate1 and integrate2) or (spinstart[0] == 0. and spinstart[1] == 0. and spinstart[2] == 0.)){
        spinstart = {polstart*Omega1[0]/w1, polstart*Omega1[1]/w2, polstart*Omega1[2]/w2, spinstart[3], spinstart[4]};
        cout << "\nStarting spin " << xstart << " ";
        spinstepper.initialize(spinstart, xstart, pi/w1);
    }
    spinend = spinstart;
    if (not interpolate){
        boost::numeric::odeint::integrate_adaptive(spinstepper,
            [&](const state_type &y, state_type &dydx, const value_type &x){
                auto W = SpinPrecessionAxis(x, *this);
                dydx[0] = W[1]*y[2] - W[2]*y[1];
                dydx[1] = W[2]*y[0] - W[0]*y[2];
                dydx[2] = W[0]*y[1] - W[1]*y[0];
                dydx[3] = 1.;
                dydx[4] = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);
            },
            spinend, xstart, xend, static_cast<double>(pi/w1)/*,
            [&](const state_type &y, const value_type &x){
                logger->PrintSpin(p, spinstepper, stepper, field);
            }*/);
    }
    else{
        array<alglib::spline1dinterpolant, 3> omega_int;
        alglib::real_1d_array ts; // set up values for interpolation of precession axis
        vector<alglib::real_1d_array> omega(3);
        const int int_points = 10;
        ts.setlength(int_points + 1);
        for (int i = 0; i < 3; i++){
            omega[i].setlength(int_points + 1);
            omega[i][0] = Omega1[i];
            omega[i][int_points] = Omega2[i];
        }

        for (int i = 1; i < int_points; i++){ // calculate precession axis at several points along trajectory step
            ts[i] = xstart + i*(xend - xstart)/int_points;
            auto W = SpinPrecessionAxis(ts[i], *this);
            omega[0][i] = W[0];
            omega[1][i] = W[1];
            omega[2][i] = W[2];
        }

        for (int i = 0; i < 3; i++)
            alglib::spline1dbuildcubic(ts, omega[i], omega_int[i]); // interpolate all three components of precession axis

        boost::numeric::odeint::integrate_adaptive(spinstepper, 
            [&](const state_type &y, state_type &dydx, const value_type &x){
                array<double, 3> W = {alglib::spline1dcalc(omega_int[0], x), alglib::spline1dcalc(omega_int[1], x), alglib::spline1dcalc(omega_int[2], x)};
                dydx[0] = W[1]*y[2] - W[2]*y[1];
                dydx[1] = W[2]*y[0] - W[0]*y[2];
                dydx[2] = W[0]*y[1] - W[1]*y[0];
                dydx[3] = 1.;
                dydx[4] = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);
            },
            spinend, xstart, xend, static_cast<double>(pi/w1),
            [&](const state_type &y, const value_type &x){
                Observer(x, {y[0], y[1], y[2]}, *this);
            });
    }
    return not integrate2;
}