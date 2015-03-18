2D field test
=============

This test loads a 2D tabulated field whose z-component depends on cos(sqrt(x*x + y*y)*cos(z).
It then runs PENTrack to print a cut through the field and compares the result with an analytical calculation of the field.
The difference is then shown in a 3D scatter plot. It should display the typical under-/overshoot of bicubic interpolation at the edges and should be close to zero in the center.

Run RunTest.sh to run the test.
