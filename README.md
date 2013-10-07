PENTrack
========

PENTrack - a simulation tool for ultra-cold neutrons, protons and electrons

Compilation
-----------

To compile the code, you need to copy some files from Numerical Recipes (http://www.nr.com, only used v3.02 so far) into the /nr/ directory:
* nr3.h (main header, for best compatibility define the "Doub" type as "long double" instead of "double")
* interp_1d.h (cubic spline interpolation)
* odeint.h (ODE integration main header)
* stepper.h (ODE stepper)
* stepperdopr853.h (8th order Runge-Kutta method)
* stepperdopr853.c (8th order Runge-Kutta method)

Writing your own simulation
---------------------------

You can modify the simulation on four different levels:

1. Modify the /in/* configuration files and use the implemented particles, sources, spectra, etc. (more info in the corresponding *.in files)
2. Modify the code in main.c and combine TParticle-, TGeometry-, TField-, TSource-classes etc. into your own simulation (you can generate a Doxygen documentation)
3. Implement your own particles by inheriting from the TParticle class and fill the Reflect-, Absorb- and Decay-Routines with the corresponding physics
4. Make low level changes to the existing classes

Defining your experiment
------------------------

* Geometry can be imported using binary STL files which describe surfaces as triangle meshes (Check your coordinate system! "Up" with respect to gravity is the z direction).
* Magnetic and electric fields (axisymmetric 2D and 3D) can be included from text-based field maps (right now, only "Vectorfields OPERA" maps in cgs units are supported). You can also define analytic fields from straight, finite conductors.
* Particle sources can be defined using STL files or manual parameter ranges, however source spectra and velocity distributions are hardcoded into the TMCGenerator class at the moment

Run the simulation
------------------

Type "make" to compile the code, then run the executable. Some information will be shown during runtime. Log files (start- and end-values, tracks and snapshots of the particles) will be written to the /out/ directory, depending on the options chosen in /in/config.in.

Physics
-------

All particles use the same relativistic equation of motion, including gravity, Lorentz force and magnetic force on their magnetic moment. UCN interaction with matter is described via the Fermi potential formalism and the Lambert model for diffuse reflection. Protons and electrons do not have any interaction so far, they are just stopped when hitting a wall. UCN spin tracking is also included.

