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

class TConfig{
private:
	typedef std::map<std::string, std::map<std::string, std::string> > map_type;
	map_type _map;
public:
	TConfig(){};
	TConfig(const std::string &inpath){
		ReadFromFile(inpath);
	}
	void ReadFromFile(const std::string &inpath);

	std::map<std::string, std::string>& operator[](const std::string &section){
		return _map[section];
	}
	map_type::const_iterator begin(){
		return _map.begin();
	}
	map_type::const_iterator end(){
		return _map.end();
	}
	
	friend std::ostream& operator<<(std::ostream &str, const TConfig &conf);
};



#endif /* CONFIG_H_ */
