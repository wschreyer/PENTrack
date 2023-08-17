#include "config.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

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
				_map.insert(std::make_pair(section, std::map<std::string, std::string>()));
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

void TConfig::convert(const std::string &configpath){
	// let user decide what to do
	if (_map.find("GEOMETRY") == _map.end() || _map.find("FIELDS") == _map.end() ||
			_map.find("SOURCE") == _map.end() || _map.find("PARTICLES") == _map.end()){
		std::cout << "The config file seems to be incomplete.\nAre you using an old-style configuration with geometry.in and particle.in files?\nIf so, should I try to convert those to a new-style configuration? [y/N]\n";
		std::string answer;
		std::getline(std::cin, answer);
		if (answer != "y" && answer != "Y")
			exit(-1);
	}
	else
		return;

	boost::filesystem::path parentpath = boost::filesystem::path(configpath).parent_path();
	boost::filesystem::path geometryin = parentpath / "geometry.in";
	boost::filesystem::path particlein = parentpath / "particle.in";
	TConfig gconf(geometryin.native());
	std::swap(_map["GLOBAL"], _map["global"]); // rename global section
	_map.erase("global");
	_map["MATERIALS"] = gconf["MATERIALS"]; // copy materials, geometry and fields from geometry.in
	_map["GEOMETRY"] = gconf["GEOMETRY"];
	_map["FIELDS"] = gconf["FIELDS"];

	// read source parameters
	std::string sourcemode;
	std::istringstream(gconf["SOURCE"].begin()->first) >> sourcemode;
	_map["SOURCE"]["sourcemode"] = sourcemode;
	std::string particle, STLfile, ActiveTime, Enormal, PhaseSpaceWeighting;
	std::array<std::string, 6> p;
	std::istringstream sourceparams(gconf["SOURCE"].begin()->second);
	if (sourcemode == "STLsurface")
		sourceparams >> particle >> STLfile >> ActiveTime >> Enormal;
	else if (sourcemode == "STLvolume")
		sourceparams >> particle >> STLfile >> ActiveTime >> PhaseSpaceWeighting;
	else if (sourcemode == "cylvolume" || sourcemode == "boxvolume")
		sourceparams >> particle >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> ActiveTime >> PhaseSpaceWeighting;
	else if (sourcemode == "cylsurface")
		sourceparams >> particle >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> ActiveTime >> Enormal;

	TConfig pconf(particlein.native());

	_map["SOURCE"]["sourcemode"] = sourcemode; // reformat [SOURCE] parameters from geometry.in
	_map["SOURCE"]["STLfile"] = STLfile;
	_map["SOURCE"]["particle"] = particle;
	_map["SOURCE"]["ActiveTime"] = ActiveTime;
	_map["SOURCE"]["Enormal"] = Enormal;
	_map["SOURCE"]["PhaseSpaceWeighting"] = PhaseSpaceWeighting;
	_map["SOURCE"]["parameters"] = p[0] + ' ' + p[1] + ' ' + p[2] + ' ' + p[3] + ' ' + p[4] + ' ' + p[5];

	std::vector<std::string> sourcevars = { "Emin", "Emax", "spectrum", "phi_v_min", "phi_v_max", "phi_v",
											"theta_v_min", "theta_v_max", "theta_v", "polarization" , "pulseWidth", "pulseGap"};
	std::vector<std::string> anglevars = { "phi_v_min", "phi_v_max", "theta_v_min", "theta_v_max" };
	for (auto var : sourcevars){
		if (pconf[particle].count(var) > 0)
			_map["SOURCE"][var] = pconf[particle][var]; // copy spectrum and angular distributions from particle.in to [SOURCE]
		else
			_map["SOURCE"][var] = pconf["all"][var];

		if (std::find(anglevars.begin(), anglevars.end(), var) != anglevars.end())
			_map["SOURCE"][var] = std::to_string(std::stof(_map["SOURCE"][var])*180./M_PI);
	}

	for (auto section : pconf){ // copy every variable from every section of particle.in, except for those moved to [SOURCE]
		std::copy_if(section.second.begin(), section.second.end(), std::inserter(_map[section.first], _map[section.first].end()),
					[&sourcevars](const std::pair<std::string, std::string> &var){
						return std::find(sourcevars.begin(), sourcevars.end(), var.first) == sourcevars.end();
					}
		);
	}
	std::swap(_map["PARTICLES"], _map["all"]); // rename "all" section to "PARTICLES"
	_map.erase("all");

	// set up a list of all known variables
	std::vector<std::string> knownvars =
		{"simtype", "simcount", "simtime", "materials_file", "secondaries", "BCutPlane", "MRSolidAngleDRP", "MRThetaIEnergy",
		"sourcemode", "STLfile", "parameters", "particle", "ActiveTime", "Enormal", "PhaseSpaceWeighting",
		"tau", "tmax", "lmax", "endlog", "tracklog", "trackloginterval", "hitlog", "snapshotlog", "snapshots", "spinlog", "spinloginterval",
		"spintimes", "Bmax", "flipspin", "interpolatefields"};
	std::vector<std::pair<std::string, bool> > knownvarsfound;
	std::transform(knownvars.begin(), knownvars.end(), std::back_inserter(knownvarsfound),
			[](const std::string &var){ return std::make_pair(var, false); });
	std::transform(sourcevars.begin(), sourcevars.end(), std::back_inserter(knownvarsfound),
			[](const std::string &var){ return std::make_pair(var, false); });

	std::cout << "The old-style config might use paths for STL files and field tables that are relative to the executable.\nE.g. " << (++_map["GEOMETRY"].begin())->second << '\n';
	std::cout << "The new-style config uses paths relative to the config file. Do you want me to replace certain paths?\nE.g. replace \"in/\" with \"\"\n";
	std::string search, replace;
	std::cout << "Search: ";
	std::getline(std::cin, search);
	std::cout<< "Replace with: ";
	std::getline(std::cin, replace);
	std::cout << '\n';

	for (auto &section : _map){ // go through all variables and check if they are in the list of known variables and replace paths
		if (section.first == "MATERIALS") // ignore user-defined variables in materials section
			continue;
		for (auto &var : section.second){
			auto foundpos = var.second.find(search);
			if (foundpos != std::string::npos)
				var.second = var.second.replace(foundpos, search.length(), replace);

			bool isnumber = std::all_of(var.first.begin(), var.first.end(), [](const char c){ return std::isdigit(c); });
			auto found = std::find_if(knownvarsfound.begin(), knownvarsfound.end(),
										[&var](const std::pair<std::string, bool> &v){ return v.first == var.first; });

			if (found != knownvarsfound.end() || ((section.first == "GEOMETRY" || section.first == "FIELDS") && isnumber))
				found->second = true; // if variable is in list of known variables or is a number in the GEOMETRY and FIELD sections, mark variable as existing
			else
				std::cout << "Section " << section.first << " contains unknown variable " << var.first << ". It might have become deprecated or incompatible. See the example config for the correct format.\n"; // else, warn that variable is unknown
		}
	}
	for (auto var : knownvarsfound){ // warn if any known variables were not found in new config
		if (!var.second)
			std::cout << "Variable " << var.first << " was not found in new config and will be loaded with a default value. You might want to add it.\n";
	}

	// rename old config files
	boost::filesystem::rename(configpath, parentpath / "config.old");
	boost::filesystem::rename(geometryin, parentpath / "geometry.old");
	boost::filesystem::rename(particlein, parentpath / "particle.old");

	// write new config
	std::ofstream confout(configpath);
	confout << *this << std::endl;

	std::cout << "\nDone. Press ENTER to continue simulation.\n";
	std::getline(std::cin, search);
}


std::ostream& operator<<(std::ostream &str, const TConfig &conf){
	for (auto section : conf._map){
		str << "\n[" << section.first << "]\n";
		for (auto var : section.second){
			str << var.first << " " << var.second << "\n";
		}
	}
	return str;
}

double EvalFormula(TConfig &config, const std::string formulaname, const std::map<std::string, double> &variables){
    auto formula = config["FORMULAS"].lower_bound(formulaname);
    if (formula == config["FORMULAS"].end() or formula->first != formulaname)
        throw std::runtime_error("Formula " + formulaname + " not found in config file");
    exprtk::symbol_table<double> symbols;
    for (auto var: variables){
        if (not symbols.add_constant(var.first, var.second)){
			throw std::runtime_error("Error parsing variable " + var.first);
		}
    }
    exprtk::parser<double> parser;
    exprtk::expression<double> expr;
    expr.register_symbol_table(symbols);
    if (not parser.compile(formula->second, expr))
        throw std::runtime_error("Could not evaluate formula " + formula->first + ": " + parser.error());
    return expr.value();
}
