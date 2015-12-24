#include "bruteforce.h"
#include "globals.h"

TBFIntegrator::TBFIntegrator(double agamma, std::string aparticlename, std::map<std::string, std::string> &conf, std::ofstream &spinout)
				: gamma(agamma), particlename(aparticlename), Bmax(0), BFBminmem(std::numeric_limits<double>::infinity()),
				  spinlog(false), spinloginterval(5e-7), nextspinlog(0), intsteps(0), fspinout(spinout), starttime(0), 
				  wLstarttime(0), wL(0), blochPolar(0), startpol(1){
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
	std::istringstream(conf["startpol"]) >> startpol;
	
	if ( fabs(startpol) > 1 ) //if user specified polarization doesn't make sense
		startpol = 1;
}

void TBFIntegrator::Binterp(value_type t, value_type B[3]){
	for (int i = 0; i < 3; i++){
		B[i] = alglib::spline1dcalc(Binterpolant[i], t);
	}
}

void TBFIntegrator::operator()(state_type y, state_type &dydx, value_type x){
	value_type B[3];
	Binterp(x,B);
	dydx[0] = -gamma * (y[1] * B[2] - y[2] * B[1]);
	dydx[1] = -gamma * (y[2] * B[0] - y[0] * B[2]);
	dydx[2] = -gamma * (y[0] * B[1] - y[1] * B[0]);
}

void TBFIntegrator::LogSpin(value_type x1, const state_type &y1, const double pv1[6], const double E1[3], const double dE1[3][3],
   			    value_type x2, const state_type &y2, const double pv2[6], const double E2[3], const double dE2[3][3]){
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
		
		//need the maximum accuracy in spinoutlog for the larmor frequency to see any difference
		fspinout << std::setprecision(std::numeric_limits<double>::digits);
		fspinout << "t Polar logPolar Ix Iy Iz Bx By Bz wL "
			    "x1 y1 z1 v1x v1y v1z E1x E1y E1z dE1xdx dE1xdy dE1xdz dE1ydx dE1ydy dE1ydz " 
			    "dE1zdx dE1zdy dE1zdz x2 y2 z2 v2x v2y v2z E2x E2y E2z dE2xdx dE2xdy dE2xdz "
			    "dE2ydx dE2ydy dE2ydz dE2zdx dE2zdy dE2zdz\n" ;
	}
	
	value_type B[3];
	Binterp(x2, B);
	value_type BFpol = (y2[0]*B[0] + y2[1]*B[1] + y2[2]*B[2])/sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	value_type BFlogpol = 0;
	if (BFpol<0.5)
		BFlogpol = log10(0.5-BFpol);
	else if (BFpol==0.5)
		BFlogpol = 0.0;
	
	wL = LarmorFreq(x1, y1, x2, y2);
		
	fspinout << x2 << " " << BFpol << " " << BFlogpol << " "
		<< 2*y2[0] << " " << 2*y2[1] << " " << 2*y2[2] << " "
		<< B[0] << " " << B[1] << " " << B[2] << " " << wL << " ";
		
	//write out position, velocity, E field, and E field derivatives at x1 (the x2 from TBFIntegrator::Integrate)	
	for (int i = 0; i < 6; i++)
		fspinout << pv1[i] << " "; 
	for ( int i = 0; i < 3; i++) 
		fspinout << E1[i] << " "; 
	for ( int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
			fspinout << dE1[i][j] << " "; 
	
	//write out position, velocity, E field and E field derivatives at x2 (the x2 from TBFIntegrator::Integrate)
	for (int i = 0; i < 6; i++)
		fspinout << pv2[i] << " "; 
	for ( int i = 0; i < 3; i++) 
		fspinout << E2[i] << " "; 
	for ( int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
			fspinout << dE2[i][j] << " "; 
	
	fspinout << "\n";	
}


long double TBFIntegrator::Integrate(double x1, double y1[6], double dy1dx[6], double B1[4][4], double E1[3], double dE1[3][3], 
					double x2, double y2[6], double dy2dx[6], double B2[4][4], double E2[3], double dE2[3][3]){
	if (gamma == 0 || x1 == x2 ) 
		return 1;
		
	bool BruteForce1 = false, BruteForce2 = false;;
	for (unsigned int i = 0; i < BFtimes.size(); i += 2){
		BruteForce1 |= (x1 >= BFtimes[i] && x1 < BFtimes[i+1]);
		BruteForce2 |= (x2 >= BFtimes[i] && x2 < BFtimes[i+1]);
	}
	if (BruteForce1 || BruteForce2){
		BFBminmem = std::min(BFBminmem, std::min(B1[3][0],B2[3][0])); // save lowest value for info

		// check if this value is worth for Bloch integration
		if (B1[3][0] < Bmax || B2[3][0] < Bmax){
			if (I_n.empty()){
				I_n.resize(3, 0);
				if (B1[3][0] > 0){
					state_type test(3);
					//1) obtain the spin with the startpol polarization as if the Bfield were along the z-axis, with the y-component being 0
					//a) Find angle that spin vector needs to be rotated (assuming startpol represents z component of spin, i.e. probability of finding spin up or down)
					long double theta = acosl(startpol);
					//b) Define spin vector rotated away from z-axis by angle theta towards x-axis
					I_n[0] = 0.5*sinl(theta);
					I_n[1] = 0;
					I_n[2] = 0.5*cosl(theta);
					
					//2) Rotate the spin vector into coordinate system where the B0 field defines the z-axis and choose the starting point of the spin to be 
					//either along the x or y axes projected onto the plane defined by B0 depending on whichever axes has smaller dot product with the B0 field			
					const double xaxis[3] = { 1, 0, 0 };
					const double yaxis[3] = { 0, 1, 0 };
					startBField[0] = B1[0][0]/B1[3][0];
					startBField[1] = B1[1][0]/B1[3][0];
					startBField[2] = B1[2][0]/B1[3][0];
										
//					const double startBField[3] = { 0, 1, 0 };
//					std::cout << "startBField: " << startBField[0] << ", startBField: " << startBField[1] << ", startBField: " << startBField[2] << std::endl; 
//					for (int i = 0; i < 3; i++)
//						std::cout << "Istart[" << i << "]: " << I_n[i] << std::endl;
					
					//to ensure that the x-axis of the coordinate transform is not chosen to be along the B0 field
					if ( fabs(startBField[0]) <= fabs(startBField[1]) ) // if the Bx is smaller than By the x-axis will represent the x-axis of the coordinate transform
						RotateVector(&I_n[0], startBField, xaxis);
					else 
						RotateVector(&I_n[0], startBField, yaxis);
					
//					for (int i = 0; i < 3; i++)
//						std::cout << "Ifinal[" << i << "]: " << I_n[i] << std::endl;

				}
				else
					I_n[2] = 0.5;
				starttime = x1; wLstarttime = x1;
				std::cout << "\nBF starttime, " << x1 << " ";
			}
			
			
			//calculate temporal derivatives of electric field from spatial derivatives (needed for vxE derivative which is needed for B spline)
			double dE1dt[3], dE2dt[3];
			for (int i = 0; i < 3; i++ ) {
				dE1dt[i] = y1[3]*dE1[i][0] + y1[4]*dE1[i][1] + y1[5]*dE1[i][2];
				dE2dt[i] = y2[3]*dE2[i][0] + y2[4]*dE2[i][1] + y2[5]*dE2[i][2];
			} 
				
			
			// set up cubic spline for each magnetic field component for fast field interpolation between x1 and x2
			alglib::real_1d_array t, B, dBdt;
			t.setlength(2);
			B.setlength(2);
			dBdt.setlength(2);
			t[0] = x1, t[1] = x2;
			//go through all the coordinates and calculate B1i, dB1idt, B2i and dB2idt where 1 and 2 represent the B field at x1, x2
			//and i represents the coordinates x, y, z. The results of each interpolation are stored in the Binterpolant array of spline1dcubic objects
			for (int i = 0; i < 3; i++){
				//Bi_corr = Bi_default + Bi_from_vxE
				//dBidt = dBi_defaultdt + dBi_from_vxEdt
				B[0] = B1[i][0] + (y1[3 + (i + 1) % 3]*E1[(i + 2) % 3] - y1[3 + (i + 2) % 3]*E1[(i + 1) % 3]) / (c_0*c_0); //Adding vxE effect to Bx, By, Bz
				dBdt[0] = y1[3]*B1[i][1] + y1[4]*B1[i][2] + y1[5]*B1[i][3]
						+ ( (dy1dx[3 + (i + 1) % 3]*E1[(i + 2) % 3] - dy1dx[3 + (i + 2) % 3]*E1[(i + 1) % 3])
						+   (y1[3 + (i + 1) % 3]*dE1dt[(i + 2) % 3] - y1[3 + (i + 2) % 3]*dE1dt[(i + 1) % 3]) ) / (c_0*c_0); //Contribution to temporal Bfield derivative from vxE part
				
				B[1] = B2[i][0] + (y2[3 + (i + 1) % 3]*E2[(i + 2) % 3] - y2[3 + (i + 2) % 3]*E2[(i + 1) % 3]) / (c_0*c_0);
				dBdt[1] = y2[3]*B2[i][1] + y2[4]*B2[i][2] + y2[5]*B2[i][3]
						+ ( (dy2dx[3 + (i + 1) % 3]*E2[(i + 2) % 3] - dy2dx[3 + (i + 2) % 3]*E2[(i + 1) % 3]) 
						+   (y2[3 + (i + 1) % 3]*dE2dt[(i + 2) % 3] - y2[3 + (i + 2) % 3]*dE2dt[(i + 1) % 3]) ) / (c_0*c_0);

				/**Parameters to spline1dbuildcubic :
				* @param spline nodes (i.e. independent var)
				* @param function vals ( dependent vals)
				* @param boundLType - boundary conditions on the left boundary (1 = first derivative)
				* @param left boundary condition ( first derivative since boundLType is 1 for us)
				* @param boundRType - right boundary condition type (1 = first derivative condition), 
				* @param boundR - boundary condition on the right boundary ( first derivative is used since boundRType = 1 )
				* @param splineinterpolant the interpolated B field is calculated and stored in Binterpolant because it is passed by reference.
				*
				* In this case we have three different spline1dinterpolant's all stored together in the Binterpolant array
				**/
				try { // this call is known to produce errors sometimes becuase t[0] and t[1] are sometimes the same value
					alglib::spline1dbuildcubic(t, B, 2, 1, dBdt[0], 1, dBdt[1], Binterpolant[i]); // create cubic spline with known derivatives as boundary conditions
				} catch ( alglib::ap_error e ) {
//					std::cout << "t1: " << t[0] << ", t2: " << t[1] << ", difference: " << t[0] - t[1] << std::endl;
					std::cout << "Error message from 1dsplinebuild: " << e.msg.c_str() << std::endl;
				}
			}

			// set up integrator and spin logging
			if (spinlog && nextspinlog < x1) // set spin log time to x1 if smaller
				nextspinlog = x1;
			dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12), stepper_type());
			stepper.initialize(I_n, x1, std::abs(1./gamma/B1[3][0]/8/pi)); // initialize stepper with step length = quarter of one precession
			while(true){
				stepper.do_step(boost::ref(*this)); // do step
				intsteps++;
				
				//calculate larmor precession frequency and the projection of spin onto magnetic field (s_z). Required for EDM measurement simulation.
				wL = LarmorFreq(stepper.previous_time(), stepper.previous_state(), stepper.current_time(), stepper.current_state());
					
				value_type curB[3]; 
				Binterp(stepper.current_time(), curB);
				blochPolar = (I_n[0]*curB[0] + I_n[1]*curB[1] + I_n[2]*curB[2])/sqrt(curB[0]*curB[0] + curB[1]*curB[1] + curB[2]*curB[2]);
				
				double prevspinlog = stepper.previous_time();  
				state_type prevspinstate = stepper.previous_state();
				 
				while (spinlog && nextspinlog <= x2 && nextspinlog <= stepper.current_time()){ // log spin if step ended after nextspinlog
					stepper.calc_state(nextspinlog, I_n);
					LogSpin(prevspinlog, prevspinstate, y1, E1, dE1, nextspinlog, I_n, y2, E2, dE2);
					prevspinlog = nextspinlog;
					prevspinstate = I_n;
					nextspinlog += spinloginterval;
				}
				if (stepper.current_time() >= x2){ // if stepper reached/overshot x2
					stepper.calc_state(x2, I_n); // get final state
					break;
				}
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

double TBFIntegrator::LarmorFreq(value_type x1, const state_type &y1, value_type x2, const state_type &y2){
	double B1[3] = { startBField[0], startBField[1], startBField[2] }; 
	double B2[3] = { startBField[0], startBField[1], startBField[2] };
//	double B1[3], B2[3];
//	Binterp((x1+x2)/2, B1); // calculate field at t1
//	Binterp((x1+x2)/2, B2); // calculate field at t2
//	std::cout << "B1x: " << B1[0] << ", B1y: " << B1[1] << ", B1z: " << B1[2] << std::endl;
//	std::cout << "B2x: " << B2[0] << ", B2y: " << B2[1] << ", B2z: " << B2[2] << std::endl;
	state_type I1perp(3), I2perp(3);
	double IB1 = y1[0]*B1[0] + y1[1]*B1[1] + y1[2]*B1[2];
	double IB2 = y2[0]*B2[0] + y2[1]*B2[1] + y2[2]*B2[2];
	double B1abs2 = B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2];
	double B2abs2 = B2[0]*B2[0] + B2[1]*B2[1] + B2[2]*B2[2];
	for (int i = 0; i < 3; i++){
		I1perp[i] = y1[i] - IB1*B1[i]/B1abs2; // calculate projection of spin onto plane orthogonal to magnetic field
		I2perp[i] = y2[i] - IB2*B2[i]/B2abs2; // Iperp = I - Iparallel = I - (I*B)*B/|B|^2
	}
	double I1perpabs = sqrt(I1perp[0]*I1perp[0] + I1perp[1]*I1perp[1] + I1perp[2]*I1perp[2]);
	double I2perpabs = sqrt(I2perp[0]*I2perp[0] + I2perp[1]*I2perp[1] + I2perp[2]*I2perp[2]);
	double deltaphi = acos((I1perp[0]*I2perp[0] + I1perp[1]*I2perp[1] + I1perp[2]*I2perp[2])/I1perpabs/I2perpabs); // calculate angle between I1perp und I2perp
	// the final wL is the weighted average of the value obtained from the previous steps and the current step
	// this method is equivalent to obtaining a cumulative deltaPhi since start time and dividing by the total time since passed
	if ( boost::math::isfinite(deltaphi/(x2-x1)) ) {
		//if the previous calculation of wL produced an error (nan or inf), then wL should be reinitialized so that the weighted average is not done with -1,
		//the wLstarttime must also be reset if wL is reset to ensure the average is done correctly
		if ( wL == -1) {
			wL = 0;
			wLstarttime = x1;
		}
		return wL*(x1-wLstarttime)/(x2-wLstarttime) + (deltaphi/(x2 - x1)/2/pi)*(x2-x1)/(x2-wLstarttime);
	}
	else //wL is set to default value of -1 when an error occured
		return -1;
}

double TBFIntegrator::getLarmorFreq () { return wL; }  
double TBFIntegrator::getBlochPolar () { return blochPolar; }

