#ifndef STEPPER_H_
#define STEPPER_H_

#include <array>
#include <functional>

#include <boost/numeric/odeint.hpp>

#include "geometry.h"
#include "globals.h"

typedef double value_type; ///< data type used for trajectory integration
typedef std::vector<value_type> state_type; ///< type representing current particle state (position, velocity, proper time, and polarization)
typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, value_type> stepper_type; ///< basic integration stepper (5th-order Runge-Kutta)
typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type; ///< integration step length controller
typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type; ///< integration step interpolator
//	typedef boost::numeric::odeint::bulirsch_stoer_dense_out<state_type, value_type> dense_stepper_type2; ///< alternative stepper type (Bulirsch-Stoer)

const static int STATE_VARIABLES = 9;
const static int SPIN_STATE_VARIABLES = 5;


class TStep{
private:
    dense_stepper_type stepper;
    value_type xstart;
    value_type xend;
    state_type yend;
    bool next_step_over_collision;
    bool step_over_collision;
    value_type next_step;
    std::multimap<TCollision, bool> collisions;
    void reset(const value_type &x, const state_type &y, const value_type &dx = 0.);
    bool iterate_collision(const value_type &x1, const state_type &y1, const value_type &x2, const state_type &y2, const TGeometry &geometry, int iteration);
public:
    TStep(const value_type &x, const state_type &y);
    stopID next(std::function<void(const state_type&, state_type&, const value_type&)> derivs, const TGeometry &geometry);
    value_type GetStartTime() const { return xstart; };
    value_type GetTime() const { return xend; };
    state_type GetState() const { return yend; };
    state_type GetState(const value_type &x) const;
    dense_stepper_type& GetStepper(){ return stepper; };
    std::array<value_type, 3> GetPosition() const { return {yend[0], yend[1], yend[2]}; };
    std::array<value_type, 3> GetPosition(const value_type &x) const{
        state_type y = GetState(x);
        return {y[0], y[1], y[2]};
    };
    std::array<value_type, 3> GetVelocity() const { return {yend[3], yend[4], yend[5]}; };
    std::array<value_type, 3> GetVelocity(const value_type &x) const {
        state_type y = GetState(x);
        return {y[3], y[4], y[5]};
    };
    value_type GetProperTime() const { return yend[6]; };
    value_type GetProperTime(const value_type &x) const { return GetState(x)[6]; };
    value_type GetPolarization() const { return yend[7]; };
    value_type GetPolarization(const value_type &x) const { return GetState(x)[7]; };
    value_type GetPathLength() const { return yend[8]; };
    value_type GetPathLength(const value_type &x) const { return GetState(x)[8]; };
    bool IsStepOverCollision() const { return step_over_collision; };
    std::multimap<TCollision, bool> GetCollisions() const { return collisions; };
    void SetStepEndToMatchComponent(const int componentindex, const value_type &component){
        value_type x = xstart + (xend - xstart)*(component - GetState(xstart)[componentindex])/(GetState()[componentindex] - GetState(xstart)[componentindex]);
        SetStepEnd(x);
    };
    void SetStepEnd(const value_type &x){
        state_type y = GetState(x);
        SetStepEnd(x, y);
    };
    void SetStepEnd(const value_type &x, const state_type &y){
        xend = x;
        yend = y;
        next_step_over_collision = false;
        next_step = x;
    };
};

#endif //STEPPER_H_