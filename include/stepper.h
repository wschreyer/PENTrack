#ifndef STEPPER_H_
#define STEPPER_H_

#include <array>
#include <functional>

#include <boost/numeric/odeint.hpp>

#include "geometry.h"
#include "globals.h"


class TStep{
private:
    typedef double value_type; ///< data type used for trajectory integration
    typedef std::vector<value_type> state_type; ///< type representing current particle state (position, velocity, proper time, and polarization)
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, value_type> stepper_type; ///< basic integration stepper (5th-order Runge-Kutta)
    typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type; ///< integration step length controller
    typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type; ///< integration step interpolator
    //	typedef boost::numeric::odeint::bulirsch_stoer_dense_out<state_type, value_type> dense_stepper_type2; ///< alternative stepper type (Bulirsch-Stoer)

    dense_stepper_type stepper;
    value_type xstart;
    value_type xend;
    state_type yend;
    int polstart;
    int polend;
    bool next_step_over_collision;
    bool step_over_collision;
    value_type next_step;
    std::multimap<TCollision, bool> collisions;
    bool iterate_collision(const value_type &x1, const state_type &y1, const value_type &x2, const state_type &y2, const TGeometry &geometry, int iteration);
    state_type GetState() const { return yend; };
    state_type GetState(const value_type &x) const;
    dense_stepper_type& GetStepper(){ return stepper; };
    void SetStepEnd(const value_type &x, const state_type &y){
        xend = x;
        yend = y;
        next_step_over_collision = false;
        next_step = x;
    };

    dense_stepper_type spinstepper;
    state_type spinstart;
    state_type spinend;
    double noflipprob;
public:
    TStep(const double &t, const std::array<double, 3> &pos, const std::array<double, 3> &v, const double &propertime, const int &polarization, const double &pathlength);
    stopID next(const std::function<std::array<double, 3>(const double&, const std::array<double, 3>&, const std::array<double, 3>&, const double&)>& EquationOfMotion, const TGeometry &geometry);
    double GetStartTime() const { return xstart; };
    double GetTime() const { return xend; };
    std::array<double, 3> GetPosition() const { return {yend[0], yend[1], yend[2]}; };
    std::array<double, 3> GetPosition(const double &x) const{
        state_type y = GetState(x);
        return {y[0], y[1], y[2]};
    };
    std::array<double, 3> GetVelocity() const { return {yend[3], yend[4], yend[5]}; };
    std::array<double, 3> GetVelocity(const double &x) const {
        state_type y = GetState(x);
        return {y[3], y[4], y[5]};
    };
    double GetProperTime() const { return yend[6]; };
    double GetProperTime(const double &x) const { return GetState(x)[6]; };
    int GetPolarization() const { return polend; };
    int GetPolarization(const double &x) const { return x < xend ? polstart : polend; };
    double GetPathLength() const { return yend[7]; };
    double GetPathLength(const double &x) const { return GetState(x)[7]; };
    bool IsStepOverCollision() const { return step_over_collision; };
    std::multimap<TCollision, bool> GetCollisions() const { return collisions; };
    void SetStepEndToMatch(const std::function<double(const double& t)> &func, const double &match);
    void SetStepEnd(const double &x){
        state_type y = GetState(x);
        SetStepEnd(x, y);
    };
    void SetVelocity(const std::array<double, 3> &v){
        yend[3] = v[0];
        yend[4] = v[1];
        yend[5] = v[2];
    };
    void SetPolarization(const int &pol){
        if (abs(pol) != 1)
            throw std::runtime_error("Tried to set invalid polarization. Polarization can only be +/-1");
        polend = pol;
    };

    std::array<double, 3> GetSpin() const { return {spinend[0], spinend[1], spinend[2]}; };
    double GetSpinIntegrationTime() const { return spinend[3]; };
    double GetSpinPhase() const { return spinstart[4]; };
    double GetSpinPolarization(const std::array<double, 3>& axis) const;
    bool IntegrateSpin(const std::function<std::array<double, 3>(const double& t, const TStep &stepper)>& SpinPrecessionAxis,
                        const std::vector<double>& times, const double& maxfreq, const bool& interpolate,
                        const std::function<void(const double&, const std::array<double, 3>&, const TStep&)> &Observer);
};

#endif //STEPPER_H_