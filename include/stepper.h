#ifndef INCLUDE_STEPPER_H_
#define INCLUDE_STEPPER_H_

#include <array>
#include <boost/numeric/odeint.hpp>

#include "globals.h"
#include "cgaltypes.h"

static const double MAX_TRACK_DEVIATION = 0.001;

class TTrajectoryState{
public:
	typedef std::array<double, 8> state_type; ///< type representing current particle state (position, velocity, proper time, and trajectory length)
private:
	double _time;
	state_type _state;
	int _polarization;
public:
	TTrajectoryState(const double time, const state_type &state, const TTrajectoryState &previous):
		_time(time), _state(state), _polarization(previous._polarization){ }
	double time() const { return _time; }
	state_type state() const { return _state; }
	template<class equation_of_motion> double guess_time_step(const equation_of_motion &eom) const { return 10*MAX_TRACK_DEVIATION/sqrt(velocity().squared_length()); }
	CPoint position() const { return CPoint(_state[0], _state[1], _state[2]); }
	CVector velocity() const { return CVector(_state[3], _state[4], _state[5]); }
	double beta() const { return sqrt(velocity().squared_length()); }
	double proper_time() const { return _state[6]; }
	double polarization() const { return _polarization; }
	double trajectory_length() const { return _state[7]; }

	void velocity(const CVector &v){ std::copy(v.cartesian_begin(), v.cartesian_end(), _state.begin() + 3); }
	void flip_polarization(const int polarization){ _polarization *= -1; }
	template<class equation_of_motion> void state_derivs(const equation_of_motion &eom, const state_type &state, state_type &derivs, const double time);
};

bool operator==(const TTrajectoryState &lhs, const TTrajectoryState &rhs);
bool operator!=(const TTrajectoryState &lhs, const TTrajectoryState &rhs);


class TSpinState{
public:
	typedef std::array<double, 3> state_type; ///< type representing current spin state (spin vector in lab frame, total time, total phase of rotation)
private:
	double _time;
	state_type _state;
public:
	TSpinState(const double time, const state_type &state, const TSpinState &previous): _time(time), _state(state){ };
	state_type state() const { return _state; }
	double time() const { return _time; }
	template<class equation_of_motion> double guess_time_step(const equation_of_motion &eom) const { return M_PI_2/sqrt(eom(*this).squared_length()); }
	CVector spin() const { return CVector(_state[0], _state[1], _state[2]); }
};


template<class TState>
class TStep{
private:
	typedef boost::numeric::odeint::runge_kutta_dopri5<typename TState::state_type, double> stepper_type; ///< basic integration stepper (5th-order Runge-Kutta)
	typedef boost::numeric::odeint::controlled_runge_kutta<stepper_type> controlled_stepper_type; ///< integration step length controller
	typedef boost::numeric::odeint::dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type; ///< integration step interpolator

	dense_stepper_type _stepper;
	TState _start;
	TState _end;
	TStep(const TState &start, const TState &end, const dense_stepper_type &stepper): _stepper(stepper), _start(start), _end(end){ };
public:
	template<class equation_of_motion> TStep(const TState &state, const equation_of_motion &eom);
	double start_time() const { return _start.time(); }
	double end_time() const { return _end.time(); }
	TState state(const double time) const;
	template<class equation_of_motion> void next_step(const TState &s, const equation_of_motion &eom);
	template<class equation_of_motion> void next_step(const equation_of_motion &eom){ next_step(_end, eom); }
	TStep<TState> sub_step(const double starttime, const double endtime);
};

using TTrajectoryStep = TStep<TTrajectoryState>;
using TSpinStep = TStep<TSpinState>;


#endif /* INCLUDE_STEPPER_H_ */
