#include "globals.h"

#include "stepper.h"

using namespace std;

static const double MAX_TRACK_DEVIATION = 0.001; ///< max deviation of actual trajectory from straight line between start and end points of a step used for geometry-intersection test. If deviation is larger, the step will be split
static const int MAX_COLLISION_ITERATIONS = 100;

TStep::TStep(const value_type &x, const state_type &y){
    xstart = x;
    reset(x, y, 10.*MAX_TRACK_DEVIATION/sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])); // initialize stepper with fixed spatial length
}

state_type TStep::GetState(const value_type &x) const {
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


void TStep::reset(const value_type &x, const state_type& y, const value_type &dx){
    stepper.initialize(y, x, dx==0 ? stepper.current_time_step() : dx);
    xend = x;
    yend = y;
    next_step_over_collision = false;
    step_over_collision = false;
    collisions.clear();
}

stopID TStep::next(std::function<void(const state_type&, state_type&, const value_type&)> derivs, const TGeometry &geometry){
    step_over_collision = false;
    state_type ystart;
    if (yend != GetState(xend)){
        stepper.initialize(yend, xend, stepper.current_time_step());
    }
    if (stepper.current_time() == xend){
        try{
            stepper.do_step(derivs);
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

    double l2 = pow(yend[8] - ystart[8], 2); // actual length of step squared
    double d2 = pow(yend[0] - ystart[0], 2) + pow(yend[1] - ystart[1], 2) + pow(yend[2] - ystart[2], 2); // length of straight line between start and end point of step squared
    double dev2 = 0.25*(l2 - d2); // max. possible squared deviation of real path from straight line
    int nsubsteps = dev2 > 0 ? static_cast<int>(ceil(sqrt(dev2)/MAX_TRACK_DEVIATION)) : 1;
    for (int i = 0; i < nsubsteps; ++i){
        double x1 = xstart + i*(xend - xstart)/nsubsteps;
        double x2 = i + 1 == nsubsteps ? xend : xstart + (i+1)*(xend-xstart)/nsubsteps;
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
