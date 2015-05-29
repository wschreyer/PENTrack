PENTrack
========

PENTrack - a simulation tool for ultra-cold neutrons, protons and electrons

If you just want to do simulations, you should check out the stable releases, which have been tested. If you want the latest and greatest features and probably some bugs, or even want to contribute, you should check out the latest revision from master branch.


External libraries
-----------

### CGAL

The [Computational Geometry Algorithms Library](http://www.cgal.org/) is used to detect collisions of particle tracks with the experiment geometry defined by triangle meshes using AABB trees.
Some Linux distributions (e.g. Ubuntu, Debian) include the libcgal-dev package, for all others, it has to be downloaded and installed manually from the website. In the latter case, you may have to adjust the CGAL_INCLUDE, CGAL_LIB and CGAL_SHAREDLIB paths in the Makefile.

CGAL v4.1 - v4.5 have been tested.

### Boost

The [Boost C++ libraries](https://www.boost.org/) are a prerequisite for the CGAL library. Additionally, the simulation uses a 64bit Mersenne Twister pseudo-random number generator and a Runge-Kutta integrator included in Boost 1.53.0 and newer.

Boost 1.53.0 - 1.57.0 have been tested. 1.56.0 and newer seem to require a C++11 capable compiler (i.e. GCC 4.9.0+) and possibly the option -std=c++11 in the CFLAGS in the Makefile.

### muparser

[muparser](http://muparser.beltoforion.de/) is a fast formula parser and is used to interpret energy distributions etc. given by the user in particle.in.

It is included in the repository.

### ALGLIB

[ALGLIB](http://www.alglib.net) is used to do 1D and 2D interpolation for field calculations.
It also provides numerical optimization and integration routines required for the MicroRoughness model for UCN interaction with matter.

It is included in the repository.

### libtricubic

[Lekien and Marsden] (http://dx.doi.org/10.1002/nme.1296) developed a tricubic interpolation method in three dimensions. It is included in the repository.


Run the simulation
------------------

Type "make" to compile the code, then run the executable. Some information will be shown during runtime. Log files (start- and end-values, tracks and snapshots of the particles) will be written to the /out/ directory, depending on the options chosen in the *.in files.
Three command line parameters can be passed to the executable: a job number (default: 0) which is appended to all log file names, a path from where the *.in files should be read (default: in/) and a path to which the out file will be written (default: out/).


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

Particle sources can be defined using STL files or manual parameter ranges. Particle spectra and velocity distributions can also be conviniently defined in the particle.in file.


Physics
-------

All particles use the same relativistic equation of motion, including gravity, Lorentz force and magnetic force on their magnetic moment. UCN interaction with matter is described with the Fermi potential formalism and either the [Lambert model](https://en.wikipedia.org/wiki/Lambert%27s_cosine_law) or the MicroRoughness model (see see [Z. Physik 254, 169--188 (1972)](http://link.springer.com/article/10.1007%2FBF01380066) and [Eur. Phys. J. A 44, 23-29 (2010)](http://ucn.web.psi.ch/papers/EPJA_44_2010_23.pdf)) for diffuse reflection and includes spin flips on wall bounce. Protons and electrons do not have any interaction so far, they are just stopped when hitting a wall. Spin tracking by bruteforce integration of the Bloch equation is also included.



Writing your own simulation
---------------------------

You can modify the simulation on four different levels:

1. Modify the /in/*.in configuration files and use the implemented particles, sources, fields, etc. (more info below and in the corresponding *.in files)
2. Modify the code in main.c and combine TParticle-, TGeometry-, TField-, TSource-classes etc. into your own simulation (you can generate a Doxygen documentation by typing `doxygen doxygen.config`)
3. Implement your own particles, sources or fields by inheriting from the corresponding base classes and fill the virtual routines with the corresponding physics
4. Make low level changes to the existing classes


Output
-------

Output files are separated by particle type, (e.g. electron, neutron and proton) and type of output (endlog, tracklog, ...). Output files are only created if particles of the specific type are simulated and can also be completely disabled for each particle type individually by adding corresponding variables in 'particle.in'. All output files are tables with space separated columns; the first line contains the column name.

Types of output: endlog, tracklog, hitlog, snapshotlog, spinlog.

### Endlog

The endlog keeps track of the starting and end parameters of the particles simulated. In the endlog, you get the following parameters:

- jobnumber: corresponds to job number of the PENTrack run (command line parameter of PENTrack). Only used when running multiple instances of PENTrack simultaneously.
- particle: number which corresponds to particle being simulated 
- tstart: time at which particle is created in the simulation [s]
- xstart, ystart, zstart: coordinate from which the particle starts [m]
- vxstart, vystart,vzstart: velocity with which the particle starts [m/s]
- polstart: initial polarization of the particle (-1,1)
- Hstart: initial total energy of particle [eV]
- Estart: initial kinetic energy of the particle [eV]
- Bstart: magnetic field at starting point [T]
- Ustart: electric potential at starting point [V]
- solidstart: number of geometry part the particle started in, see geometry.in
- tend: time at which particle simulation is stopped [s]
- xend, yend, zend: coordinate at which the particle ends [m]
- vxend, vyend, vzend: velocity with which the particle ends [m/s]
- polend: final polarization of the particle (-1,1)
- Hend: final total energy of particle [eV]
- Eend: final kinetic energy of the particle [eV]
- Bend: magnetic field at stopping point [T]
- Uend: electric potential at stopping point [V]
- solidend: number of geometry part the particle stopped in, defined in geometry.in
- stopID: code which identifies why the particle was stopped (defined in globals.h)
  - 0: not categorized
  - -1: did not finish (reached max. simulation time)
  - -2: hit outer poundaries
  - -3: produced error during trajectory integration
  - -4: decayed
  - -5: found no initial position
  - -6: produced error during geometry collision detection
  - -7: produced error during tracking of crossed material boundaries
  - 1: absorbed in bulk material (see solidend)
  - 2: absorbed on total reflection on surface (see solidend)
- NSpinflip: number of spin flips that the particle underwent during simulation
- spinflipprob: probability that the particle has undergone a spinflip, calculated by bruteforce integration of the Bloch equation
- Nhit: number of times particle hit a geometry surface (e.g. wall of guidetube)
- Nstep: number of steps that it took to simulate particle
- trajlength: the total length of the particle trajectory from creation to finish [m]
- Hmax: the maximum total energy that the particle had during trajectory [eV]

### Snapshotlog

Switching on snapshotlog in "particle.in" will output the particle parameters at additional 'snapshot times' (also defined in "particle.in") in the snapshotlog. It contains the same data fields as the endlog.

### Tracklog

The tracklog contains the complete trajectory of the particle, with following parameters:

- jobnumber: corresponds to job number of the PENTrack run (command line parameter of PENTrack). Only used when running multiple instances of PENTrack simultaneously.
- particle: number which corresponds to particle being simulated 
- polarization: the polarization of the particle (-1,1)
- t: time [s]
- x,y,z coordinates [m]
- vx,vy,vz: velocity [m/s]
- H: total energy [eV]
- E: kinetic energy [eV]
- Bx: x-component of Magnetic Field at coordinates [T]
- dBxdx, dBxdy, dBxdz: gradient of x-component of Magnetic field at coordinates [T/m]
- By : y-component of magnetic field at coordinates[T]
- dBydx, dBydy, dBydz : gradient of y-component of Magnetic field at coordinates [T/m]
- Bz: y-component of magnetic field at coordinates[T]
- dBzdx, dBzdy, dBzdz: gradient of z-component of Magnetic field at coordinates [T/m]
- Ex, Ey, Ez: X, Y, and Z component of electric field at coordinates [V/m]
- V: electric potential at coordinates [V]

### Hitlog

This log contains all the points at which particle hits a surface of the experiment geometry. This includes both reflections and transmissions.

In the hitlog, you get the following parameters:

- jobnumber: number which corresponds to job number of the PENTrack run (command line parameter of PENTrack). Only used when running multiple instances of PENTrack simultaneously.
- particle: number which corresponds to particle being simulated
- t: time at which particle hit the surface [s]
- x, y, z: coordinate of the hit [m]
- v1x, v1y,v1z: velocity of the particle before it interacts with the new geometry/surface [m/s]
- pol1: polarization of particle before it interacts with the new geometry/surface (-1,1)
- v2x, v2y,v2z: velocity of the particle after it interacted with the new geometry/surface [m/s]
- pol2: polarization of particle after it interacted with the new geometry/surface (-1,1)
- nx,ny,nz: normal to the surface of the geometry that the particle hits
- solid1: ID number of the geometry that the particle starts in
- solid2: ID number of the geometry that the particle hits

### Spinlog

If bruteforce integration of particle spin precession is active, it will be logged into the spinlog. In the spinlog, you get the following parameters:

- t: time [s]
- Babs: absolute magnetic field [T]
- Polar: projection of Bloch spin vector onto magnetic field direction vector (Bloch vector has dimensionless length 0.5, so this quantity is also dimensionless)
- logPolar: logarithm of above value [dimensionless]
- Ix, Iy, Iz: x, y and z components of the Bloch vector (times 2) [dimensionless]
- Bx, By, Bz: x, y and z components of magnetic field direction vector (B[i]/Babs) [dimensionless]

### Writing Output files to ROOT readable files

merge_all.c: A [ROOT](http://root.cern.ch) script that writes all out-files (or those, whose filename matches some pattern) into a single ROOT file containing trees for each log- and particle-type.

