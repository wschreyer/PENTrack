#include "bruteforce.h"
#include "globals.h"

TBFIntegrator::TBFIntegrator(double agamma, std::string aparticlename, std::map<std::string, std::string> &conf, std::ofstream &spinout)
				: gamma(agamma), particlename(aparticlename), Bmax(0), BFBminmem(std::numeric_limits<double>::infinity()),
				  spinlog(false), spinloginterval(5e-7), nextspinlog(0), intsteps(0), fspinout(spinout), starttime(0){
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

void TBFIntegrator::LogSpin(const state_type &y, value_type x){
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
		fspinout << "t Polar logPolar Ix Iy Iz Bx By Bz\n";
	}

	value_type B[3];
	Binterp(x, B);
	value_type BFpol = (y[0]*B[0] + y[1]*B[1] + y[2]*B[2])/sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	value_type BFlogpol = 0;
	if (BFpol<0.5)
		BFlogpol = log10(0.5-BFpol);
	else if (BFpol==0.5)
		BFlogpol = 0.0;
		
	fspinout << x << " " << BFpol << " " << BFlogpol << " "
		<< 2*y[0] << " " << 2*y[1] << " " << 2*y[2] << " "
		<< B[0] << " " << B[1] << " " << B[2] << '\n';
}


long double TBFIntegrator::Integrate(double x1, double y1[6], double dy1dx[6], double B1[4][4], double E1[3],
					double x2, double y2[6], double dy2dx[6], double B2[4][4], double E2[3]){
	if (gamma == 0)
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
					I_n[0] = B1[0][0]/B1[3][0]*0.5;
					I_n[1] = B1[1][0]/B1[3][0]*0.5;
					I_n[2] = B1[2][0]/B1[3][0]*0.5;
				}
				else
					I_n[2] = 0.5;
				starttime = x1;
				std::cout << "\nBF starttime, " << x1 << " ";
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
				B[0] = B1[i][0] + (y1[3 + (i + 1) % 3]*E1[(i + 2) % 3] - y1[3 + (i + 2) % 3]*E1[(i + 1) % 3]) / (c_0*c_0); //Adding vxE effect to Bx, By, Bz, Note dEdt=0 (for now)
				dBdt[0] = y1[3]*B1[i][1] + y1[4]*B1[i][2] + y1[5]*B1[i][3]
						+ (dy1dx[3 + (i + 1) % 3]*E1[(i + 2) % 3] - dy1dx[3 + (i + 2) % 3]*E1[(i + 1) % 3]) / (c_0*c_0); //Contribution to temporal Bfield derivative from vxE part, dEdt=0 (for now)
				
				B[1] = B2[i][0] + (y2[3 + (i + 1) % 3]*E2[(i + 2) % 3] - y2[3 + (i + 2) % 3]*E2[(i + 1) % 3]) / (c_0*c_0);
				dBdt[1] = y2[3]*B2[i][1] + y2[4]*B2[i][2] + y2[5]*B2[i][3]
						+ (dy2dx[3 + (i + 1) % 3]*E2[(i + 2) % 3] - dy2dx[3 + (i + 2) % 3]*E2[(i + 1) % 3]) / (c_0*c_0);

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
				alglib::spline1dbuildcubic(t, B, 2, 1, dBdt[0], 1, dBdt[1], Binterpolant[i]); // create cubic spline with known derivatives as boundary conditions
			}

			// set up integrator and spin logging
			if (spinlog && nextspinlog < x1) // set spin log time to x1 if smaller
				nextspinlog = x1;
			dense_stepper_type stepper = boost::numeric::odeint::make_dense_output(static_cast<value_type>(1e-12), static_cast<value_type>(1e-12), stepper_type());
			stepper.initialize(I_n, x1, std::abs(1./gamma/B1[3][0]/8/pi)); // initialize stepper with step length = quarter of one precession
			while(true){
				stepper.do_step(boost::ref(*this)); // do step
				intsteps++;
				double wL = LarmorFreq(stepper.previous_time(), stepper.previous_state(), stepper.current_time(), stepper.current_state());
				while (spinlog && nextspinlog <= x2 && nextspinlog <= stepper.current_time()){ // log spin if step ended after nextspinlog
					stepper.calc_state(nextspinlog, I_n);
					LogSpin(I_n, nextspinlog);
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
	double B1[3], B2[3];
	Binterp(x1, B1); // calculate field at t1
	Binterp(x2, B2); // calculate field at t2
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
	return deltaphi/(x2 - x1)/2/pi;
}

