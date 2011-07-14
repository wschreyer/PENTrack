/**
 * \file
 * Different estimations of spin flip probability.
 */

/**
 * Probability for NO neutron spin flip after Rabi (from Matora paper)
 *
 * @param Br Radial magnetic field component
 * @param Bz Axial magnetic field component
 * @param dBrdr Derivative of radial component with respect to radial coordinate
 * @param dBrdz Derivative of radial component with respect to axial coordinate
 * @param dBzdr Derivative of axial component with respect to radial coordinate
 * @param dBzdz Derivative of axial component with respect to axial coordinate
 * @param vr_n Radial velocity of neutron
 * @param vz_n Axial velocity of neutron
 * @param t Time
 *
 * @return Returns probability that neutron undergoes NO spin flip
 */
extern long double rabiplus(double Br,double Bz,double dBrdr,double dBrdz,double dBzdr,double dBzdz,double vr_n,double vz_n,double t);


/**
 * Probability for neutron spin flip after Rabi (from Matora paper)
 *
 * @param Br Radial magnetic field component
 * @param Bz Axial magnetic field component
 * @param dBrdr Derivative of radial component with respect to radial coordinate
 * @param dBrdz Derivative of radial component with respect to axial coordinate
 * @param dBzdr Derivative of axial component with respect to radial coordinate
 * @param dBzdz Derivative of axial component with respect to axial coordinate
 * @param vr_n Radial velocity of neutron
 * @param vz_n Axial velocity of neutron
 * @param t Time
 *
 * @return Returns probability that neutron undergoes spin flip
 */
extern long double rabimin(double Br,double Bz,double dBrdr,double dBrdz,double dBzdr,double dBzdz,double vr_n,double vz_n,double t);


/**
 * Probability for neutron spin flip after Vladimirsky (Soviet Physics JETP, Volume 12, Number 4, April 1961)
 *
 * @param Bx Magnetic field x component
 * @param By Magnetic field y component
 * @param Bz Magnetic field z component
 * @param dBxdx Derivative of x component with respect to x
 * @param dBxdy Derivative of x component with respect to y
 * @param dBxdz Derivative of x component with respect to z
 * @param dBydx Derivative of y component with respect to x
 * @param dBydy Derivative of y component with respect to y
 * @param dBydz Derivative of y component with respect to z
 * @param dBzdx Derivative of z component with respect to x
 * @param dBzdy Derivative of z component with respect to y
 * @param dBzdz Derivative of z component with respect to z
 * @param vx Velocity of neutron in x direction
 * @param vy Velocity of neutron in y direction
 * @param vz Velocity of neutron in z direction
 *
 * @return Returns probability that neutron undergoes spin flip
 */
extern long double vladimirsky(long double Bx,long double By, long double Bz,long double dBxdx, long double dBxdy, long double dBxdz, long double dBydx, long double dBydy, long double dBydz, long double dBzdx, long double dBzdy, long double dBzdz, long double Bws, long double vx, long double vy, long double vz);


/**
 * Calculate adiabacity criterion dBdt*hbar/2/mu/B^2 (should be <<1 for adiabacity)
 *
 * @param Bx Magnetic field x component
 * @param By Magnetic field y component
 * @param Bz Magnetic field z component
 * @param dBxdx Derivative of x component with respect to x
 * @param dBxdy Derivative of x component with respect to y
 * @param dBxdz Derivative of x component with respect to z
 * @param dBydx Derivative of y component with respect to x
 * @param dBydy Derivative of y component with respect to y
 * @param dBydz Derivative of y component with respect to z
 * @param dBzdx Derivative of z component with respect to x
 * @param dBzdy Derivative of z component with respect to y
 * @param dBzdz Derivative of z component with respect to z
 * @param vx Velocity of neutron in x direction
 * @param vy Velocity of neutron in y direction
 * @param vz Velocity of neutron in z direction
 *
 * @return Returns adiabacity criterion dBdt*hbar/2/mu/B^2
 */
extern long double thumbrule(long double Bx, long double By, long double Bz,long double dBxdx, long double dBxdy, long double dBxdz, long double dBydx, long double dBydy, long double dBydz, long double dBzdx, long double dBzdy, long double dBzdz, long double Bws, long double vx, long double vy, long double vz);

