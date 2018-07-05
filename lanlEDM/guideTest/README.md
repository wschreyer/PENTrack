guideTest
=========================

Neutrons bounce around in a cylindrical guide tube with perfect
specular reflection. Both ends of the tube are capped, with one end
being a polyethelyne absorber.

Cylinder is 5 m long with inner radius 0.29m (goes from x = [-2.5,2.5]).
Neutrons are created in a spherical source (diameter = .05m) at the
left (x = -2.45) end of the tube

There is an exponentially decreasing field along x, and a linear field in the z direction,
and a static field in the z direction

For more, see config.in

Run the test with ./RunTest.sh, edit parameters in the config file

Visualize output by running graphstuff.py with the appropriate command line flags
(python graphstuff.py --help)
