/**
 * \file
 * Contains everything needed to read/write PENTrack-style config files:
 * 
 * [SECTION]
 * variable parameter1 parameter2 parameter3 [...]
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <map>
#include <iosfwd>

#include "exprtk.hpp"

/**
 * TConfig class is read from in file with structure
 *
 * [SECTION]
 * variable parameter1 parameter2 parameter3 [...]
 * variable2 parameter
 *
 * [SECTION2]
 * ...
 *
 * and maps it into a a nested map<map<string, string> > structure
 */
class TConfig{
private:
    typedef std::map<std::string, std::map<std::string, std::string> > map_type;
    map_type _map; ///< internally stored map structure
public:
	/**
	 * Create empty TConfig
	 */
	TConfig(){};

	/**
	 * Create config from section/variable map
	 */
	explicit TConfig(const std::map<std::string, std::map<std::string, std::string> > &map): _map(map){};

	/**
	 * Create TConfig and read sections and variables from file
	 *
	 * @param inpath Path of file to read from
	 */
	TConfig(const std::string &inpath){
		ReadFromFile(inpath);
	}

	/**
	 * Read sections and variables from file
	 *
	 * @param inpath Path of file to read from
	 */
	void ReadFromFile(const std::string &inpath);

	/**
	 * Access section, maps to operator[] of internally stored map structure
	 *
	 * @param section Section to access
	 *
	 * @return Returns map<string, string> of variables and parameters
	 */
    std::map<std::string, std::string>& operator[](const std::string &section){
		auto it = _map.lower_bound(section);
		if (it == _map.end() or it->first != section)
			throw std::runtime_error("Section " + section + " not found in config file");
		return it->second;
    }

	/**
	 * Get iterator to first section
	 *
	 * @return Returns begin iterator of internally stored map
	 */
	map_type::const_iterator begin(){
		return _map.begin();
	}

	/**
	 * Get iterator to end of section list
	 *
	 * @return Returns end iterator of internally stored map
	 */
	map_type::const_iterator end(){
		return _map.end();
	}
	
	/**
	 * Convert old-style configuration to new style
	 * Loads geometry.in and particle.in from the same path the config file is found in and combines them
	 *
	 * @param configpath Path of config file
	 */
	void convert(const std::string &configpath);

	/**
	 * Output operator, writes TConfig to file
	 *
	 * @param str Stream to write to
	 * @param conf TConfig that should be written
	 *
	 * @return Returns reference to output stream
	 */
	friend std::ostream& operator<<(std::ostream &str, const TConfig &conf);
};

/**
 * Find and evaluate formulate with given name using the list of variables
 * 
 * @param config TConfig containing configuration variables
 * @param formulaname Name of the formula to be evaluated
 * @param variables Map of variables names and values used to evaluate formula
 * 
 * @return Returns value of formula
 */
double EvalFormula(TConfig &config, const std::string formulaname, const std::map<std::string, double> &variables);


#endif /* CONFIG_H_ */
