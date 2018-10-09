lanlEDM/comsolFieldTest
========================

Reads in a 3D arrow plot from COMSOL containing the magnetic field from a helmholtz coil
https://www.comsol.com/model/magnetic-field-of-a-helmholtz-coil-15

Outputs a cut of the B field (simtype 4)
Can be any text file with 6 columns, x y z Bx By Bz
(Lines beginning with % or # will be skipped. Columns may be delineated by space, comma, or tab)
Units are assumed to be in meters and Tesla, though units can be edited in comsolField3D.h

container.stl is an empty cylinder with OD 1m and height 0.5 meter, centered at (0,0,0).
It has walls 1cm thick. For simtype = 1 we put a neutron in this field and watch it precess 


Run the test with ./RunTest.sh


graphstuff.py will plot the magnetic field in the y-z plane from the simtype=4 output, the magnetic field in the y-z plane from the original file comsolField.txt, and the spin trajectory from the simtype=1 output. The latter should be a circle in the z-x plane. 
