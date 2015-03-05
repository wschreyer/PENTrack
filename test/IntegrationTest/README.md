Integration test
================

This test simulates ultra-cold neutrons on an analytically calculatable trajectory and compares that with the numerical result.

The neutrons start at a defined point inside a cubic shell. The shell is permeated by two overlapping fields, a 2D tabulated field and a 3D tabulated field, both with a linear vertical gradient of dBz/dz = 0.5 T/m.
So, the UCN feel a negative vertical acceleration due to gravity and magnetic gradient and travel on a parabolic trajectory.
Once they hit the shell, the trajectory end point is subtracted from the analytical solution and the difference is displayed in an histogram.

Run RunTest.sh to run the test.
