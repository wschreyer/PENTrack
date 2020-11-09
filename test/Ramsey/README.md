Ramsey cycles
=======================================

This test tracks spins of ultracold neutrons during a Ramsey cycle.
The UCN-storage chamber is 100 mm high and has a radius of 200 mm. A vertical magnetic field, with a strength of 1 µT and a vertical gradient of 0.001 µT/m, and an electric field, pointing downwards with a strength of 1 MV/m, are applied.
UCN are randomly created in the chamber at t = 0. After about 0.3 s, a pi/2 spin flip is induced by a horizontally oscillating magnetic field with an amplitude of 0.01 µT and a frequency of 183.019 rad/s, applied for 1.7144 s. Then, the UCN are left to precess freely for 50 s until another pi/2 flip is induced. The simulation ends after 54 s.

Ramsey_up.in contains the configuration with electric field "up", Ramsey_down.in containes the configuration with electric field "down"

RunTest.sh scans the frequency of the pi/2-flipping field from 183.2 rad/s to 183.3 rad/s, the results can be compared to the expected Ramsey fringe pattern.
