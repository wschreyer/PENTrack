PENTrack
========

PENTrack - a simulation tool for ultra-cold neutrons, protons and electrons

If you just want to do simulations, you should check out the stable releases, which have been tested. If you want the latest and greatest features and probably some bugs, or even want to contribute, you should check out the latest revision from master branch.

External libraries
-----------

### Numerical recipes

To compile the code, you need to copy some files from [Numerical Recipes] (http://www.nr.com) (v3.00 to v3.04 should work fine) into the /nr/ directory:
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

### libtricubic

[Lekien and Marsden] (http://dx.doi.org/10.1002/nme.1296) developed a tricubic interpolation method in three dimensions. It is included in the repository.

Writing your own simulation
---------------------------

You can modify the simulation on four different levels:

1. Modify the /in/*.in configuration files and use the implemented particles, sources, spectra, etc. (more info below and in the corresponding *.in files)
2. Modify the code in main.c and combine TParticle-, TGeometry-, TField-, TSource-classes etc. into your own simulation (you can generate a Doxygen documentation by typing `doxygen doxygen.config`)
3. Implement your own particles by inheriting from the TParticle class and fill the Reflect-, Absorb- and Decay-Routines with the corresponding physics
4. Make low level changes to the existing classes

Defining your experiment
------------------------

### Geometry

Geometry can be imported using [binary STL files] (http://en.wikipedia.org/wiki/STL_%28file_format%29).
STL files use a list of triangles to describe 3D-surfaces. They can be created with most CAD software, e.g. Solidworks via "File - Save As..." using the STL file type. Solidworks makes sure that the surface normals always point outside a solid body. This is used to check if points lie inside a solid object.
Note the "Options..." button in the Save dialog, there you can change format (only binary supported), resolution, coordinate system etc.

PENTrack expects the STL files to be in unit Meters.

Gravity acts in negative z-direction, so choose your coordinate system accordingly.

Do not choose a too high resolution. Usually, setting the deviation tolerance to a value small enough to represent every detail in your geometry and leaving angle tolerance at the most coarse setting is good enough. Low angle tolerance quickly increases triangle count and degeneracy.

If you want to export several parts of a Solidworks assembly you can do the following:

1. Select the part(s) to be exported and rightclick.
2. Select "Invert selection" and then use "Suppress" on all other parts.
3. Now you can save that single part as STL (make sure the option "Do not translate STL output data to positive space" is checked and to use the same coordinate system for every part or else the different parts will not fit together).
4. You can check the positioning of the parts with e.g. [MeshLab](http://meshlab.sourceforge.net/), SolidView, Minimagics, Solidworks...

### Fields

Magnetic and electric fields (rotationally symmetric 2D and 3D) can be included from text-based field maps
(right now, only "[Vectorfields OPERA](https://www.cobham.com/about-cobham/aerospace-and-security/about-us/antenna-systems/specialist-technical-services-and-software/products-and-services/design-simulation-software/opera.aspx)" maps in cgs units are supported).
You can also define analytic fields from straight, finite conductors.

### Particle sources

Particle sources can be defined using STL files or manual parameter ranges, however particle spectra and velocity distributions are hardcoded into the TMCGenerator class at the moment

Run the simulation
------------------

Type "make" to compile the code, then run the executable. Some information will be shown during runtime. Log files (start- and end-values, tracks and snapshots of the particles) will be written to the /out/ directory, depending on the options chosen in the *.in files.

Physics
-------

All particles use the same relativistic equation of motion, including gravity, Lorentz force and magnetic force on their magnetic moment. UCN interaction with matter is described via the Fermi potential formalism and the Lambert model for diffuse reflection. Protons and electrons do not have any interaction so far, they are just stopped when hitting a wall. UCN spin tracking by bruteforce integration of the Bloch equation is also included.

