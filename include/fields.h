/**
 * \file
 * Manage all the fields.
 */

#ifndef FIELDS_H_
#define FIELDS_H_

#include <vector>

#include "field.h"
#include "config.h"
#include "globals.h"

/**
 * Contains list of all fields (2D/3D-maps, conductors, ...).
 */
class TFieldManager{
private:
    std::vector< TFieldContainer > fields; ///< list of fields
		
public:
	TFieldManager(const TFieldManager &f) = delete; ///< TFieldManager is not copyable
	TFieldManager& operator=(const TFieldManager &f) = delete; ///< TFieldManager is not copyable

	/**
	 * Constructor.
	 *
	 * Reads [FIELDS] section of configuration file and loads all field maps/conductors given there
	 *
	 * @param conf TConfig map containing field options
	 */
	explicit TFieldManager(TConfig &conf);


	/**
	 * Calculate superposition of all loaded magnetic fields at a given position and time.
	 *
	 * @param x Cartesian x coordinate
	 * @param y Cartesian y coordinate
	 * @param z Cartesian z coordinate
	 * @param t Time
	 * @param B Returns magnetic x, y, and z components of magnetic field
	 * @param dBidxj Returns spatial derivatives of each magnetic-field component (optional)
	 */
	void BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3] = nullptr) const;


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
	void EField(const double x, const double y, const double z, const double t,
			double &V, double Ei[3]) const;
};

#endif // FIELDS_H_
