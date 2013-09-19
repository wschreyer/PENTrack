PENTrack
========

PENTrack - a simulation tool for ultra-cold neutrons, protons and electrons

To compile the code, you need to copy some files from Numerical Recipes (http://www.nr.com, only used v3.02 so far) into the /nr/ directory:
nr3.h (main header, for best compatibility define the "Doub" type as "long double" instead of "double")
interp_1d.h (cubic spline interpolation)
odeint.h (ODE integration main header)
stepper.h (ODE stepper)
stepperdopr853.h (8th order Runge-Kutta method)
stepperdopr853.c (8th order Runge-Kutta method)

You can modify the simulation on four different levels:
1. Modify the /in/* configuration files and use the implemented particles, sources, spectra, etc. (more info in the corresponding *.in files)
2. Modify the code in main.c and combine TParticle-, TGeometry-, TField-, TSource-classes etc. into your own simulation (look at the Doxygen documentation)
3. Implement your own particles by inheriting from the TParticle class and fill the Reflect-, Absorb- and Decay-Routines with the corresponding physics
4. Make low level changes to the existing classes
