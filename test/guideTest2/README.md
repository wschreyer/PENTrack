guideTest2
=========================

Neutrons bounce around in a cylindrical guide tube with perfect
specular reflection. Both ends of the tube are capped, with one end
being a polyethelyne absorber.

Cylinder is 2 m long with inner diameter 0.06m (goes from x = [-2,0]).
Neutrons are created in a spherical source (diameter = .02m) at the
left (x = -1.975) end of the tube

The tube is at a height of z = -0.215m

There is a magnetic field, generated in COMSOL, of a segmented solenoid
with roughly 1 uT of field in the positive z direction (inside).
The solenoid is ~1.2 m in height and ~1 m in diameter

There is also a exponentially decreasing magnetic field in the x direction.
The B_x component starts at around ~8 uT at x = -2m,
and drops to ~1 mG at x = -0.5m (where the pipe enters the solenoid)

For more, see config.in

Run the test with ./RunTest.sh, edit parameters in the config file

Visualize depolarization of neutrons by running graphstuff.py
