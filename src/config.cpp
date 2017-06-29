#include "config.h"

#include <iostream>
#include <fstream>

#include <boost/format.hpp>

//read variables from *.in file into map
void TConfig::ReadFromFile(const std::string &inpath){
	std::cout << "Reading config from " << inpath << "\n";
	std::ifstream infile(inpath);
	if (!infile)
		throw std::runtime_error((boost::format("Could not open %s!") % inpath).str());
	
	char c;
	std::string rest,section,key;
	while (infile && (infile >> std::ws) && (c = infile.peek())){
		if (c == '[' && infile.ignore()){
			if (infile.peek() == '/'){
				section = "";
			}
			else{
				std::getline(infile, section, ']');
//				std::cout << "\nsection: " << section.c_str() << '\n';
			}
			std::getline(infile,rest);
		}
		else if (c == '#')
			std::getline(infile,rest);
		else if (section != ""){
			infile >> key;
			std::getline(infile,rest);
			if (infile){
				std::string::size_type l = rest.find('#');
				if (_map.count(section) > 0 && _map[section].count(key) > 0)
					throw std::runtime_error((boost::format("Variable %s in config file is ambiguous!") % key).str());
				if (l == std::string::npos)
					_map[section][key] = rest;
				else
					_map[section][key] = rest.substr(0,l);
//				std::cout << key << " " << vars[section][key] << '\n';
			}
		}
		else
			std::getline(infile,rest);
	}
}


std::ostream& operator<<(std::ostream &str, const TConfig &conf){
	for (auto section : conf._map){
		str << "[" << section.first << "]\n";
		for (auto var : section.second){
			str << var.first << ": " << var.second << "\n";
		}
	}
	return str;
}
