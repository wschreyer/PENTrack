#include "bruteforce.h"
#include "globals.h"
#include <math.h>

TBFIntegrator::TBFIntegrator(double agamma, std::string aparticlename, std::map<std::string, std::string> &conf, std::ofstream &spinout)
				: gamma(agamma), particlename(aparticlename), Bmax(0), BFBminmem(std::numeric_limits<double>::infinity()),
				  spinlog(false), spinloginterval(5e-7), intsteps(0), fspinout(spinout), starttime(0), spinlogdegrees(1E-5), t1(0),
				  t2(0), initialAngle(0), phaseAngle(-4), newPhaseAngle(0), numRotations(0), deltaPhi(0), 
				  prevDeltaPhi(0), larmFreq(0), blochPolar(0) {
	std::istringstream(conf["BFmaxB"]) >> Bmax;
	std::istringstream BFtimess(conf["BFtimes"]);
	do{
		double t;
		BFtimess >> t;
		if (BFtimess.good())
			BFtimes.push_back(t);
	}while(BFtimess.good());
	std::istringstream(conf["spinlog"]) >> spinlog;
	std::istringstream(conf["spinloginterval"]) >> spinloginterval;
	std::istringstream(conf["spinlogdegrees"]) >> spinlogdegrees;
}

void TBFIntegrator::Binterp(value_type t, value_type B[3]){
	value_type x = (t-t1)/(t2 - t1);
	B[0] = ((cx[0]*x + cx[1])*x + cx[2])*x + cx[3];
	B[1] = ((cy[0]*x + cy[1])*x + cy[2])*x + cy[3];
	B[2] = ((cz[0]*x + cz[1])*x + cz[2])*x + cz[3];
}

void TBFIntegrator::operator()(state_type y, state_type &dydx, value_type x){
	value_type B[3];
	Binterp(x,B);
	dydx[0] = -gamma * (y[1] * B[2] - y[2] * B[1]);
	dydx[1] = -gamma * (y[2] * B[0] - y[0] * B[2]);
	dydx[2] = -gamma * (y[0] * B[1] - y[1] * B[0]);
}

void TBFIntegrator::operator()(const state_type &y, value_type x){
	if (!spinlog)
		return;
	if (!fspinout.is_open()){
		std::ostringstream BFoutfile1;
		BFoutfile1 << outpath << "/" << std::setw(12) << std::setfill('0') << jobnumber << std::setw(0) << particlename << "spin.out";
		std::cout << "Creating " << BFoutfile1.str() << '\n';
		fspinout.open(BFoutfile1.str().c_str());
		if(!fspinout.is_open())
		{
			std::cout << "Could not open " << BFoutfile1.str() << '\n';
			exit(-1);
		}
		fspinout.precision(10);
		fspinout << "t Babs Polar logPolar Ix Iy Iz Bx By Bz larmFreq\n";
	}

	value_type B[3];
	Binterp(x, B);
	value_type BFBws = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	value_type BFpol = (y[0]*B[0] + y[1]*B[1] + y[2]*B[2])/BFBws;
	value_type BFlogpol = 0;
	if (BFpol<0.5)
		BFlogpol = log10(0.5-BFpol);
	else if (BFpol==0.5)
		BFlogpol = 0.0;
	
	//calculate the larmor frequency
	blochPolar=BFpol;
	newPhaseAngle=atan2(y[1], y[0]); //calculate phase angle from the current iteration
	
	fspinout << std::setprecision(std::numeric_limits<double>::digits); //to obtain maximum of larmFreq in log file	
	
	//when the phase angle has completed a revolution (i.e. the previous and new phase angles around the initial angle)
	if ( (initialAngle < phaseAngle && initialAngle >= newPhaseAngle && gamma < 0 ) || ( initialAngle > phaseAngle && initialAngle <= newPhaseAngle && gamma > 0 ) && phaseAngle != -4 ) 
		numRotations+=1;
	
	deltaPhi = 2.0*pi*numRotations+fabs(initialAngle-newPhaseAngle); 
	
	if ( deltaPhi < prevDeltaPhi ) //when the newPhaseAngle has come up "behind" the initial angle
                deltaPhi = 2.0*pi*numRotations+(2*pi-fabs(initialAngle-newPhaseAngle));	
	 
	larmFreq = (deltaPhi/x)/(2*pi); //convert from angular to regular frequency
	phaseAngle = newPhaseAngle; //update the current angle to the previous angle for the next iteration
	prevDeltaPhi = deltaPhi; //update value of previousDeltaPhi
	
	if ( deltaPhi-prevOut >= spinlogdegrees*pi/180 ) { //only output after 10 full rotation
		fspinout << x << " " << BFBws << " " << BFpol << " " << BFlogpol << " "
      	        	<< 2*y[0] << " " << 2*y[1] << " " << 2*y[2] << " "
	               	<< B[0]/BFBws << " " << B[1]/BFBws << " " << B[2]/BFBws  
			<< " " << larmFreq << '\n';
		prevOut=deltaPhi;
	}
}

long double TBFIntegrator::Integrate(double x1, double y1[6], double B1[4][4],
					double x2, double y2[6], double B2[4][4]){
	if (gamma == 0)
		return 1;
	
	bool BruteForce1 = false, BruteForce2 = false;;
	for (unsigned int i = 0; i < BFtimes.size(); i += 2){
		BruteForce1 |= (x1 >= BFtimes[i] && x1 < BFtimes[i+1]);
		BruteForce2 |= (x2 >= BFtimes[i] && x2 < BFtimes[i+1]);
	}
	if (BruteForce1 || BruteForce2){
		BFBminmem = std::min(BFBminmem, static_cast<double>(std::min(B1[3][0],B2[3][0]))); // save lowest value for info

		// check if this value is worth for Bloch integration
		if (B1[3][0] < Bmax || B2[3][0] < Bmax){
			if (I_n.empty()){
				I_n.resize(3, 0);
				if (B1[3][0] > 0){
				//	I_n[0] = B1[0][0]/B1[3][0]*0.5;
				//	I_n[1] = B1[1][0]/B1[3][0]*0.5;
				//	I_n[2] = B1[2][0]/B1[3][0]*0.5;
				
					I_n[0] = 0; 
					I_n[1] = 0.5;
					I_n[2] = 0;
				}
				else
					I_n[2] = 0.5;
				
				initialAngle=atan2(I_n[1], I_n[0]);
				prevOut=0;
				starttime = x1;
				std::cout << "\nBF starttime, " << x1 << " ";
//					stepper = boost::numeric::odeint::make_dense_output((value_type)1e-12, (value_type)1e-12, stepper_type());
//					stepper = boost::numeric::odeint::make_controlled(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12), stepper_type());
//					stepper = stepper_type(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12));
			}

			// calculate temporal derivative of field from spatial derivative and particle's speed dBi/dt = dBi/dxj * dxj/dt
			t1 = x1;
			t2 = x2;
			value_type dBxdt1 = B1[0][1]*y1[3] + B1[0][2]*y1[4] + B1[0][3]*y1[5];
			value_type dBydt1 = B1[1][1]*y1[3] + B1[1][2]*y1[4] + B1[1][3]*y1[5];
			value_type dBzdt1 = B1[2][1]*y1[3] + B1[2][2]*y1[4] + B1[2][3]*y1[5];

			value_type dBxdt2 = B2[0][1]*y2[3] + B2[0][2]*y2[4] + B2[0][3]*y2[5];
			value_type dBydt2 = B2[1][1]*y2[3] + B2[1][2]*y2[4] + B2[1][3]*y2[5];
			value_type dBzdt2 = B2[2][1]*y2[3] + B2[2][2]*y2[4] + B2[2][3]*y2[5];

			// calculate coefficients of cubic splines Bi(x) = ci[0]*x^3 + ci[1]*x^2 + ci[2]*x + ci[3]   (i = x,y,z)
			// with boundary conditions Bi(x1) = Bi1, Bi(x2) = Bi2, dBidt(x1) = dBidt1, dBidt(x2) = dBidt2
			value_type h = x2-x1;
			cx[0] = 2*B1[0][0] - 2*B2[0][0] + dBxdt1*h + dBxdt2*h;
			cx[1] = 3*B2[0][0] - 3*B1[0][0] - 2*dBxdt1*h - dBxdt2*h;
			cx[2] = dBxdt1*h;
			cx[3] = B1[0][0];
			cy[0] = 2*B1[1][0] - 2*B2[1][0] + dBydt1*h + dBydt2*h;
			cy[1] = 3*B2[1][0] - 3*B1[1][0] - 2*dBydt1*h - dBydt2*h;
			cy[2] = dBydt1*h;
			cy[3] = B1[1][0];
			cz[0] = 2*B1[2][0] - 2*B2[2][0] + dBzdt1*h + dBzdt2*h;
			cz[1] = 3*B2[2][0] - 3*B1[2][0] - 2*dBzdt1*h - dBzdt2*h;
			cz[2] = dBzdt1*h;
			cz[3] = B1[2][0];

			if (spinlog){
				std::vector<value_type> times;
				for (value_type x = x1; x < x2; x += spinloginterval)
					times.push_back(x);
				times.push_back(x2);
				// use dense output stepper if spin trajectory should be logged
				dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12), stepper_type());
				// integrate(ODEsystem functor, initial state, start time, end time, initial time step, observer functor)
				intsteps += boost::numeric::odeint::integrate_times(
						stepper, boost::ref(*this), I_n, times.begin(),
						times.end(), static_cast<value_type>(1e-9), boost::ref(*this));
			}
			else{
				// use simpler error controlling stepper when output not needed
				controlled_stepper_type stepper = boost::numeric::odeint::make_controlled(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12), stepper_type());
				// integrate(ODEsystem functor, initial state, start time, end time, initial time step, observer functor)
				intsteps += boost::numeric::odeint::integrate_adaptive(stepper, boost::ref(*this), I_n, x1, x2, static_cast<value_type>(1e-9));
			}

			if (B2[3][0] > Bmax || !BruteForce2){
				// output of polarisation after BF int completed
				// calculate polarisation at end of step BFpol = (I_n*B/|B|) in [-1/2,1/2]
				value_type BFpol = (I_n[0]*B2[0][0]
									 + I_n[1]*B2[1][0]
									 + I_n[2]*B2[2][0])/B2[3][0];

				std::cout << "BF dt " << x2 - starttime << ", BFflipprop " << 1 - (BFpol + 0.5) << ", intsteps taken " << intsteps << ", Bmin " << BFBminmem << " ";

				BFBminmem = std::numeric_limits<double>::infinity(); // reset values when done
				intsteps = 0;
				I_n.clear();

				return BFpol + 0.5;
			}
		}
	}

	return 1;
}
