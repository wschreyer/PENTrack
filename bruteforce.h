/**
 * \file
 * Do "brute force" integration of the Bloch equation.
 */

#ifndef BRUTEFORCE_H_
#define BRUTEFORCE_H_


/**
 * Bloch equation calculator.
 *
 * This class interpolates the magnetic field components between two track points using a cubic spline and calculates the Bloch equation dS/dt = -gamma S x B.
 */
struct TBFderivs{
	long double xx1; ///< field interpolation start time
	long double xx2; ///< field interpolation end time
	long double cx[4]; ///< cubic spline coefficients for magnetic field x-component
	long double cy[4]; ///< cubic spline coefficients for magnetic field y-component
	long double cz[4]; ///< cubic spline coefficients for magnetic field z-component
	long double gamma; ///< gyromagnetic ratio used in Bloch equation

	/**
	 * Constructor.
	 *
	 * @param agamma Gyromagnetic ratio of the particle whose spin is to be tracked.
	 * @param x1 Time of first track point.
	 * @param y1 State vector of first track point.
	 * @param B1 Magnetic field components at first track point.
	 * @param x2 Time of second track point.
	 * @param y2 State vector of sectond track point.
	 * @param B2 Magnetic field components at second track point.
	 */
	TBFderivs(long double agamma, long double x1, long double y1[6], long double B1[4][4], long double x2, long double y2[6], long double B2[4][4]): gamma(agamma){
		xx1 = x1;
		xx2 = x2;

		// calculate temporal derivative of field from spatial derivative and particle's speed dBi/dt = dBi/dxj * dxj/dt
		long double dBxdt1 = B1[0][1]*y1[3] + B1[0][2]*y1[4] + B1[0][3]*y1[5];
		long double dBydt1 = B1[1][1]*y1[3] + B1[1][2]*y1[4] + B1[1][3]*y1[5];
		long double dBzdt1 = B1[2][1]*y1[3] + B1[2][2]*y1[4] + B1[2][3]*y1[5];

		long double dBxdt2 = B2[0][1]*y2[3] + B2[0][2]*y2[4] + B2[0][3]*y2[5];
		long double dBydt2 = B2[1][1]*y2[3] + B2[1][2]*y2[4] + B2[1][3]*y2[5];
		long double dBzdt2 = B2[2][1]*y2[3] + B2[2][2]*y2[4] + B2[2][3]*y2[5];

		// calculate coefficients of cubic splines Bi(x) = ci[0]*x^3 + ci[1]*x^2 + ci[2]*x + ci[3]   (i = x,y,z)
		// with boundary conditions Bi(x1) = Bi1, Bi(x2) = Bi2, dBidt(x1) = dBidt1, dBidt(x2) = dBidt2
		long double h = x2-x1;
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
	};

	/**
	 * Do cubic spline interpolation of magnetic field components with coefficients determined in TBFderivs::TBFderivs
	 *
	 * @param x Time
	 * @param B Magnetic field components
	 */
	void Binterp(long double x, long double B[3]){
		x = (x-xx1)/(xx2-xx1);
		B[0] = ((cx[0]*x + cx[1])*x + cx[2])*x + cx[3];
		B[1] = ((cy[0]*x + cy[1])*x + cy[2])*x + cy[3];
		B[2] = ((cz[0]*x + cz[1])*x + cz[2])*x + cz[3];
	};

	/**
	 *  Bloch equation integrator calls BFderivs(x,y,dydx) to get derivatives
	 *
	 *  @param x Time
	 *  @param y Spin vector
	 *  @param dydx Returns temporal derivative of spin vector
	 */
	void operator()(Doub x, VecDoub_I &y, VecDoub_O &dydx){
		long double B[3];
		Binterp(x,B);
		dydx[0] = -gamma * (y[1] * B[2] - y[2] * B[1]);
		dydx[1] = -gamma * (y[2] * B[0] - y[0] * B[2]);
		dydx[2] = -gamma * (y[0] * B[1] - y[1] * B[0]);
	};
};


/**
 * Bloch equation integrator.
 *
 * Create this class to do "brute force" tracking of your particle's spin in magnetic fields along its track.
 * It starts tracking the spin by integrating the Bloch equation when the absolute magnetic field drops below TBFIntegrator::Bmax and stops it when the field rises above this value again.
 * Then it calculates the spin flip probability after such a low-field-pass.
 */
struct TBFIntegrator{
	long double gamma; ///< Particle's gyromagnetic ration
	long double Bmax; ///< Spin tracking is only done when absolut magnetic field drops below this value.
	vector<double> BFtimes; ///< Pairs of absolute time in between which spin tracking shall be done.
	long double BFBminmem; ///< Stores minimum field during one spin track for information
	VecDoub *I_n; ///< Spin vector
	bool spinlog; ///< Should the tracking be logged to file?
	double spinloginterval; ///< Time interval between log file entries.
	string particlename; ///< Name of particle whose spin is to be tracked, needed for logging.
	long int intsteps; ///< Count integrator steps during spin tracking for information.

	/**
	 * Constructor.
	 *
	 * Set initial values and read options from config file.
	 *
	 * @param agamma Gyromagnetic ration of particle whose spin is to be tracked.
	 * @param aparticlename Particle name.
	 * @param conf Option map containing particle specific spin tracking options.
	 */
	TBFIntegrator(long double agamma, string aparticlename, map<string, string> conf): gamma(agamma), particlename(aparticlename),
		Bmax(0), BFBminmem(numeric_limits<long double>::infinity()), I_n(NULL), spinlog(false), spinloginterval(5e-7), intsteps(0){
		istringstream(conf["BFmaxB"]) >> Bmax;
		istringstream BFtimess(conf["BFtimes"]);
		do{
			double t;
			BFtimess >> t;
			if (BFtimess.good())
				BFtimes.push_back(t);
		}while(BFtimess.good());
		istringstream(conf["spinlog"]) >> spinlog;
		istringstream(conf["spinloginterval"]) >> spinloginterval;
	}

	/**
	 * Track spin between two particle track points.
	 *
	 * Checks if spin should be tracked between given track points. Returns spin flip probability when tracking is finished.
	 *
	 * @param x1 Time at first track point.
	 * @param y1 State vector at first track point.
	 * @param x2 Time at second track point.
	 * @param y2 State vector at second track point.
	 * @param field TFieldManager to calculate magnetic field at track points.
	 * @param spinout File stream to which spin tracking should be logged.
	 *
	 * @return Probability, that NO spin flip occured (usually close to 1).
	 */
	long double Integrate(long double x1, long double y1[6], long double x2, long double y2[6], TFieldManager *field, ofstream *&spinout){
		if (gamma == 0)
			return 1;

		bool BruteForce1 = false, BruteForce2 = false;;
		for (unsigned int i = 0; i < BFtimes.size(); i += 2){
			BruteForce1 |= (x1 >= BFtimes[i] && x1 < BFtimes[i+1]);
			BruteForce2 |= (x2 >= BFtimes[i] && x2 < BFtimes[i+1]);
		}
		if (BruteForce1 || BruteForce2){
			long double B1[4][4], B2[4][4];
			field->BField(y1[0], y1[1], y1[2], x1, B1);
			field->BField(y2[0], y2[1], y2[2], x2, B2);
			BFBminmem = min(BFBminmem,min(B1[3][0],B2[3][0])); // save lowest value for info

			// check if this value is worth for Bloch integration
			if (B1[3][0] < Bmax || B2[3][0] < Bmax){
				if (I_n == NULL){
					I_n = new VecDoub(3,0.0); // fill spin vector with values
					if (B1[3][0] > 0){
						(*I_n)[0] = B1[0][0]/B1[3][0]*0.5;
						(*I_n)[1] = B1[1][0]/B1[3][0]*0.5;
						(*I_n)[2] = B1[2][0]/B1[3][0]*0.5;
					}
					else
						(*I_n)[2] = 0.5;
					cout << "\nBF starttime " << x1 << " ";
				}

				TBFderivs BFderivs(gamma, x1, y1, B1, x2, y2, B2); // create derivs-struct
				int nsave = 0;
				if (spinlog)
					nsave = (int)((x2-x1)/spinloginterval); // save integrations steps in BFdxsav intervals
				Output out(nsave); // create output object
				// create integrator<stepper<derivobject> >, params: spin vector, start time, end time, atol, rtol, first step size, min step size, output object, derivobject
				Odeint<StepperDopr853<TBFderivs> > ode(*I_n, x1, x2, 1e-13, 0, 1e-9, 0, out, BFderivs);
				ode.integrate(); // integrate
				intsteps += ode.nok + ode.nbad; // add up integration steps
				if (spinlog)
					PrintBFStep(spinout, out, BFderivs); // print integrations steps

				if (B2[3][0] > Bmax || !BruteForce2){
					// output of polarisation after BF int completed
					// calculate polarisation at end of step BFpol = (I_n*B/|B|) in [-1/2,1/2]
					long double BFpol = ((*I_n)[0]*B2[0][0] + (*I_n)[1]*B2[1][0] + (*I_n)[2]*B2[2][0])/B2[3][0];

					cout << " BF endtime " << x2 << ", BFflipprop " << 1 - (BFpol + 0.5) << ", intsteps taken " << intsteps << ", Bmin " << BFBminmem << " ";

					BFBminmem = numeric_limits<long double>::infinity(); // reset values when done
					intsteps = 0;
					delete I_n;
					I_n = NULL;

					return BFpol + 0.5;
				}
			}
		}

		return 1;
	}

	/**
	 * Print spin track to file.
	 *
	 * @param spinout Filestream to which should be logged.
	 * @param out Output struct of Odeint integrator.
	 * @param BFderivs Magnetic field interpolator to print field values.
	 */
	void PrintBFStep(ofstream *&spinout, Output &out, TBFderivs &BFderivs){
		for (int i = 0; i < out.count; i++){
			if (!spinout){
				ostringstream BFoutfile1;
				BFoutfile1 << outpath << "/" << setw(12) << setfill('0') << jobnumber << setw(0) << particlename << "spin.out";
				cout << "Creating " << BFoutfile1.str() << '\n';
				spinout = new ofstream(BFoutfile1.str().c_str());
				if(!spinout || !spinout->is_open())
				{
					cout << "Could not open " << BFoutfile1.str() << '\n';
					exit(-1);
				}
				*spinout << "t Babs Polar logPolar Ix Iy Iz Bx By Bz\n";
			}

			long double B[3];
			BFderivs.Binterp(out.xsave[i],B);
			long double BFBws =  sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
			long double BFpol = (out.ysave[0][i]*B[0] + out.ysave[1][i]*B[1] + out.ysave[2][i]*B[2])/BFBws;
			long double BFlogpol = 0;
			if (BFpol<0.5)
				BFlogpol = log10(0.5-BFpol);
			else if (BFpol==0.5)
				BFlogpol = 0.0;
			*spinout << out.xsave[i] << " " << BFBws << " " << BFpol << " " << BFlogpol << " "
					<< 2*out.ysave[0][i] << " " << 2*out.ysave[1][i] << " " << 2*out.ysave[2][i] << " "
					<< B[0]/BFBws << " " << B[1]/BFBws << " " << B[2]/BFBws << '\n';
		}
	}

};


#endif // BRUTEFORCE_H_
