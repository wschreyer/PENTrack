PENTrack
========

PENTrack - a simulation tool for ultra-cold neutrons, protons and electrons

Compilation
-----------

Always use the master branch. The dev branch contains code in active development and might not work at all. They will be merged from time to time.

To compile the code, you need to copy some files from [Numerical Recipes v3.00-v3.04] (http://www.nr.com) into the /nr/ directory:
* nr3.h (main header, one minor modification needed, patch file is provided)
* ran.h (random number generator)
* interp_1d.h (cubic spline interpolation)
* interp_linear.h (linear interpolation)
* interp_2d.h (bicubic interpolation, some modifications needed, patch file is provided)
* odeint.h (ODE integration main header, one minor midification needed, patch file is provided)
* stepper.h (ODE integration step control
* stepperdopr853.h (8th order Runge Kutta method)

Patch files can be applied by executing:  
`patch originalheader.h patchfile.diff`

You may have to change the nr file format from DOS/MAC to UNIX text file format by executing:  
`dos2unix originalheader.h`

Numerical Recipes forces us to put almost all code into header files, which makes it less easy to read but it also avoids duplicate work in code- and header-files.

Writing your own simulation
---------------------------

You can modify the simulation on four different levels:

1. Modify the /in/*.in configuration files and use the implemented particles, sources, spectra, etc. (more info in the corresponding *.in files)
2. Modify the code in main.c and combine TParticle-, TGeometry-, TField-, TSource-classes etc. into your own simulation (you can generate a Doxygen documentation by typing `doxygen doxygen.config`)
3. Implement your own particles by inheriting from the TParticle class and fill the Reflect-, Absorb- and Decay-Routines with the corresponding physics
4. Make low level changes to the existing classes

Defining your experiment
------------------------

* Geometry can be imported using [binary STL files] (http://en.wikipedia.org/wiki/STL_%28file_format%29) which describe surfaces as triangle meshes (Check your coordinate system! "Up" with respect to gravity is the z direction).
* Magnetic and electric fields (axisymmetric 2D and 3D) can be included from text-based field maps (right now, only "Vectorfields OPERA" maps in cgs units are supported). You can also define analytic fields from straight, finite conductors.
* Particle sources can be defined using STL files or manual parameter ranges, however source spectra and velocity distributions are hardcoded into the TMCGenerator class at the moment

Run the simulation
------------------

Type "make" to compile the code, then run the executable. Some information will be shown during runtime. Log files (start- and end-values, tracks and snapshots of the particles) will be written to the /out/ directory, depending on the options chosen in /in/config.in.

Physics
-------

All particles use the same relativistic equation of motion, including gravity, Lorentz force and magnetic force on their magnetic moment. UCN interaction with matter is described via the Fermi potential formalism and the Lambert model for diffuse reflection. Protons and electrons do not have any interaction so far, they are just stopped when hitting a wall. UCN spin tracking by bruteforce integration of the Bloch equation is also included.

