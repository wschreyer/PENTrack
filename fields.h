/**
 * \file
 * Manage all the fields.
 */

#ifndef FIELDS_H_
#define FIELDS_H_

#include <vector>

#include "field.h"
#include "field_2d.h"
#include "field_3d.h"
#include "conductor.h"
#include "globals.h"
#include "edmfields.h"

/**
 * Contains list of all fields (2D/3D-maps, conductors, ...).
 */
struct TFieldManager{
	protected: 
		//required for AddvCrossE function
		typedef double value_type;
		typedef std::vector<value_type> state_type;
	public:
		std::vector<TField*> fields; ///< list of fields
		int FieldOscillation; ///< If =1 field oscillation is turned on
		double OscillationFraction; ///< Field oscillation amplitude
		double OscillationFrequency; ///< Field oscillation frequency
		
		/**
		 * Constructor.
		 *
		 * Reads [FIELDS] section of configuration file and loads all field maps/conductors given there
		 *
		 * @param conf TConfig map containing field options
		 * @param aFieldOscillation Turn on field oscillations
		 * @param aOscillationFraction Amplitude of field oscillation
		 * @param aOscillationFrequency Frequency of field oscillation
		 */
		TFieldManager(TConfig &conf, int aFieldOscillation = 0, double aOscillationFraction = 0, double aOscillationFrequency = 0);

		
		/**
		 * Destructor, delete all fields.
		 */
		~TFieldManager();


		/**
		 * Calculate magnetic field at a given position and time.
		 *
		 * Chooses the right map for this position and adds racetrack fields.
		 * field matrix B:
		 *	Bx,		dBxdx,	dBxdy,	dBxdz;
		 *	By,		dBydx,	dBydy,	dBydz;
		 *	Bz,		dBzdx,	dBzdy,	dBzdz;
		 *	Babs,	dBdx,	dBdy,	dBdz;
		 *
		 * @param x Cartesian x coordinate
		 * @param y Cartesian y coordinate
		 * @param z Cartesian z coordinate
		 * @param t Time
		 * @param B Returns magnetic field component matrix
		 */
		void BField (double x, double y, double z, double t, double B[4][4]);

		
		/**
		 * Calculate electric field and potential at a given position.
		 *
		 * Chooses the right map for this position and returns interpolated values.
		 *
		 * @param x Cartesian x coordinate
		 * @param y Cartesian y coordinate
		 * @param z Cartesian z coordinate
		 * @param t Time
		 * @param V Return electric potential (!=0 only if a map with potential was loaded)
		 * @param Ei Returns electric field vector
		 */
		void EField(double x, double y, double z, double t, double &V, double Ei[3]);
		
		/**
 		* Add the motional BField induced by the v x E effect.
 		* 
 		* @param y, the 6D vector containing the particle's position and velocity
 		* @param B, the matrix containing the components of the magnetic field 
 		* @param Ei[3] the 3D components of the electric field
 		*/
		void AddvCrossE (state_type y, double B[4][4], double Ei[3] ); 
	private:

		/**
		 * Scale magnetic field due to field oscillation.
		 *
		 * @param t Time
		 */
		double BFieldScale(double t);
		
};

#endif // FIELDS_H_
