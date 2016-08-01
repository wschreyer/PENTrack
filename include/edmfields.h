/**
 * \file
 * The EDM fields used in the Ramsey cycle. 
 */

#ifndef STATICEDMFIELD_H_
#define STATICEDMFIELD_H_

#include "field.h"

/**
 * Magnetic field definition for the static B_0 field experienced by neutrons throughout the Ramsey cycle.
 *
 * Constructor requires:
 *
 * @param xoff,	the x-coordinate offset from the origin
 * @param yoff,	the y-coordinate offset from the origin
 * @param zoff, the z-coordinate offset from the origin
 * @param ang1, the polar angle theta
 * @param ang2, the azimuthal angle phi
 * @param abz,	the z component of the B0 field at the origin in Tesla
 * @param adB0zdz, the z derivative of z component of the B0 field: dB0zdz at the origin in Tesla
 * @param AC, A boolean to toggle between AC or DC field
 * @param frq, 	the frequency of the AC field
 * @param tstart1, the start time of AC
 * @param tend1, the end time of AC
 * @param pshift, the phase shift for the AC field
 * @param BoundaryWidth, Distance from the edges where folding begins
 * @param xmax, maximum x value for the field to permiate
 * @param xmin, minimum x value for the field to permiate
 * @param ymax, maximum y value for the field to permiate
 * @param ymin, minimum y value for the field to permiate
 * @param zmax, maximum z value for the field to permiate
 * @param zmin, minimum z value for the field to permiate
 *
 **/

struct TEDMStaticB0GradZField: public TField{
	double edmB0xoff;
	double edmB0yoff;
	double edmB0zoff;
	double pol_ang1;
	double azm_ang2;
	double edmB0z0;
	double edmdB0z0dz;
	bool ac;
	double f;
	double on1;
	double off1;
	double phase;
	double BoundaryWidth;
	double xmax;
	double xmin;
	double ymax;
	double ymin;
	double zmax;
	double zmin;
	double Rot1[3][3];
	double Rot2[3][3];
	double Rot3[3][3];
	double dB[9];
	double Bd[6][3];
	
	
	TEDMStaticB0GradZField(double xoff, double yoff, double zoff, double ang1, double ang2, double abz, double adB0zdz, bool AC, double frq, double tstart1, double tend1, double pshift, double bW, double _xmax, double _xmin, double _ymax, double _ymin, double _zmax, double _zmin, std::string Bscale);
	
	/**
	 * The calculation for the B0 components are derived from a first order approximation of
	 * an axially symmetric static magnetic field. See for example: http://physics.princeton.edu/~mcdonald/examples/axial.pdf
	 *
	 * @param x, the x-coordinate in the field's coordinate system
	 * @param y, the x-coordinate in the field's coordinate system
	 * @param z, the x-coordinate in the field's coordinate system
	 * @param t, the time (not used since the field is presumed to be static)
	 * @param B, Magnetic field component matrix to which the values are added
	**/
	void BField(double x, double y, double z, double t, double B[4][4]);
	
	/**
	 * Adds no electric field. 
	 * Required because we are inheriting from abstract class TField which has virtual void EField function.
	 *
	 * @param x, Cartesian x coordinate
	 * @param y, Cartesian y coordinate
	 * @param z, Cartesian z coordinate
	 * @param t, Time
	 * @param V, Electric potential
	 * @param Ei, Electric field components
	 * @param dEidxj, Spatial derivatives of electric field components
	 **/
	void EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3] = NULL) {};
	
	/**
	 * Folds the BField component (Bx, By, or Bz) and its spatial derivative components 
	 *
	 * @param x, Cartesian x coordinate
	 * @param y, Cartesian y coordinate
	 * @param z, Cartesian z coordinate
	 * @param Bxi, The Bfield component (Bx,By, or Bz) to be folded
	 * @param dBScaled, The three derivatives of each Bfield component; Only the ones pertaining to Bxi are updated in a call
	 * @param i, Integer to specify which Bxi component you passed; x = 0, y = 1, z = 2
	 **/
	void FieldSmthr(double x, double y, double z, double *Bxi, double *dBScaled, int i);
	
	/**
	 * Computes how much scale the field and its derivatives in the x y and z direction based off of how close they are to their respective boundaries
	 *
	 * @param x, Cartesian x coordinate
	 * @param y, Cartesian y coordinate
	 * @param z, Cartesian z coordinate
	 * @param compFactors, First three elements to compress in the x,y,z directions respectively. The next 3 are to indicate upper boundary (a value of 1) or a lower boundary (a value of 0)
	 **/
	void CompressionFactor(double x, double y, double z, double *compFactors);
	
	/**
	 * Smooth function used to scale the field at the edges
	 *
	 * @param x, function parameter (x = 0..1)
	 *
	 * @return Returns number between 0 and 1, smoothly rising with x
	 **/
	double SmthrStp(double x);
	
	/**
	 * Derivative of SmthStp
	 *
	 * @param x, function parameter (x = 0..1)
	 *
	 * @return Returns derivative of SmthrStp at parameter x
	 **/

	double SmthrStpDer(double x);		
};

/**
 * The static homogenous field that is parallel to the B0 field.
 * 
 */
struct TEDMStaticEField: public TField {
	double exMag, eyMag, ezMag; 	

	/**
 	* Constructor requires just the magnitude of the field. 
 	*/
	TEDMStaticEField(double aexMag, double aeyMag, double aezMag, std::string Escale);

	/**
	 * Adds no magnetic field.
	 *
	 * @param x, the x-coordinate in the field's coordinate system
	 * @param y, the x-coordinate in the field's coordinate system
	 * @param z, the x-coordinate in the field's coordinate system
	 * @param t, the time (not used since the field is presumed to be static)
	 * @param B, Magnetic field component matrix to which the values are added
	**/
	void BField(double x, double y, double z, double t, double B[4][4]) {};
	
	/**
	 * The static E field should be parallel to the B0 field in Ramsey cycle.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param V Electric potential
	 * @param Ei Electric field components
	 * @param dEidxj Spatial derivatives of electric field components
	 */
	void EField(double x, double y, double z, double t, double &V, double Ei[3], double dEidxj[3][3] = NULL);

};


#endif /*STATICEDMFIELD_H_*/
