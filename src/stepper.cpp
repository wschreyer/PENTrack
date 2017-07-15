#include "stepper.h"

template<class equation_of_motion>
void TTrajectoryState::state_derivs(const equation_of_motion &eom, const state_type &state, state_type &derivs, const double time){
	TTrajectoryState s(time, state, _polarization);
	derivs[0] = s.velocity().x(); // derivatives of position = velocity
	derivs[1] = s.velocity().y();
	derivs[2] = s.velocity().z();
	CVector a = eom(s); // calculate acceleration with equation of motion
	derivs[3] = a.x(); // derivatives of velocity = acceleration
	derivs[4] = a.y();
	derivs[5] = a.z();
	double v2 = s.velocity().squared_length();
	derivs[6] = sqrt(1. - v2/c_0/c_0); // derivative of proper time = 1/gamma
	derivs[7] = sqrt(v2); // derivative of trajectory length = velocity
}


bool operator==(const TTrajectoryState &lhs, const TTrajectoryState &rhs){
	return rhs.time() == lhs.time() && rhs.state() == lhs.state() && rhs.polarization() == lhs.polarization();
}

bool operator!=(const TTrajectoryState &lhs, const TTrajectoryState &rhs){
	return !(rhs == lhs);
}


template<class TState> template<class equation_of_motion>
TStep<TState>::TStep(const TState &state, const equation_of_motion &eom){
	_start = state;
	_stepper = boost::numeric::odeint::make_dense_output(1e-9, 1e-9, stepper_type());
	_stepper.initialize(_start.state(), _start.time(), _start.guess_time_step(eom));
	next_step(eom);
}


template<class TState>
inline TState TStep<TState>::state(const double time) const {
	if (time < start_time() || time > end_time())
		throw std::runtime_error("Tried to calculate particle state outside of step range");
	else if (time == start_time())
		return _start;
	else if (time == end_time())
		return TParticleState(time, _stepper.current_state(), _start.polarization());
	else{
		typename TState::state_type state;
		_stepper.calc_state(time, state);
		return TState(time, state, _start);
	}
}

template<class TState> template <class equation_of_motion>
void TStep<TState>::next_step(const TState &s, const equation_of_motion &eom){
	if (s != TState(_stepper.current_time, _stepper.current_state, _start))
		_stepper.initialize(_end.state(), _end.time(), _stepper.current_time_step());
	_start = s;
	_stepper.do_step(std::bind(&TState::state_derivs, _start, eom, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
	_end = TState(_stepper.current_time(), _stepper.current_state(), _start);
}

template<class TState>
TStep<TState> TStep<TState>::sub_step(const double starttime, const double endtime){
	if (starttime < start_time() || starttime > end_time() || endtime < start_time() || endtime > end_time())
		throw std::runtime_error("Tried to divide trajectory step outside of step range");
	return TStep<TState>(state(starttime), state(endtime), _stepper);
}
