#include "geometry.h"

#include <iostream>
#include <algorithm>

#include "globals.h"

using namespace std;

std::istream& operator>>(std::istream &str, material &mat){
	str >> mat.name;
    str >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb >> mat.SpinflipProb >> mat.RMSRoughness >> mat.CorrelLength;
	if (!str)
		throw std::runtime_error((boost::format("Could not read material %s!") % mat.name).str());
	if (mat.DiffProb < 0 or mat.DiffProb > 1)
	    throw std::range_error("You set a diffuse-reflection probability outside range 0..1 for material " + mat.name +"!");
    if (mat.SpinflipProb < 0 or mat.SpinflipProb > 1)
        throw std::range_error("You set a spin-flip probability outside range 0..1 for material " + mat.name +"!");
    str >> mat.InternalBField;
    if (!str){
        std::cout << "No internal magnetic field set for material " << mat.name << ". Assuming 0T.\n";
        mat.InternalBField = 0.;
    }
    str >> mat.ModifiedLambertProb;
	if (not str){
		std::cout << "No modified-Lambert reflection set for material " << mat.name << ". Assuming 0.\n";
	    mat.ModifiedLambertProb = 0.;
	}
    if (mat.SpinflipProb < 0 or mat.SpinflipProb > 1){
        throw std::range_error("You set a modified-Lambert probability outside range 0..1 for material " + mat.name +"!");
    }

	str >> mat.LossPerBounce;
	if (not str){
	    std::cout << "No loss per bounce set for material " << mat.name << ". Assuming 0.\n";
	    mat.LossPerBounce = 0.;
	}
    if (mat.LossPerBounce < 0 or mat.LossPerBounce > 1)
        throw std::range_error("You set a loss-per-bounce probability outside range 0..1 for material " + mat.name +"!");

    int diffmodels = 0;
    if (mat.DiffProb != 0) ++diffmodels;
    if (mat.RMSRoughness != 0 || mat.CorrelLength != 0) ++diffmodels;
    if (mat.ModifiedLambertProb != 0) ++diffmodels;
    if (diffmodels > 1){
        throw std::runtime_error((boost::format("You have set parameters for more than one diffuse-reflection model for material %s! "\
												"Now I don't know which to use.") % mat.name).str());
    }
    if (diffmodels == 0) {
        std::cout << "No diffuse reflection set for material " << mat.name << "\n";
    }

	str >> mat.microfacetDistributionWidth;
	if (not str){
		std::cout << "No microfacet distribution width defined for material " << mat.name << ". Assuming 0\n";
		mat.microfacetDistributionWidth = 0;
	}
	if (mat.microfacetDistributionWidth < 0){
		throw std::runtime_error("You set a microfacet distribution width smaller than 0 for " + mat.name + "!");
	}
	str >> mat.microfacetDistributionWidthExponent;
	if (not str){
		std::cout << "No microfacet distribution exponent defined for material " << mat.name << ". Assuming 0\n";
		mat.microfacetDistributionWidthExponent = 0;
	}

	return str;
}
std::ostream& operator<<(std::ostream &str, const material &mat){
	str << mat.name << " " << mat.FermiReal << " " << mat.FermiImag << " " << mat.DiffProb << " " << mat.SpinflipProb << " ";
	str << mat.RMSRoughness << " " << mat.CorrelLength << " " << mat.InternalBField << " " << mat.ModifiedLambertProb << " ";
	str << mat.LossPerBounce << " " << mat.microfacetDistributionWidth << " " << mat.microfacetDistributionWidthExponent << "\n";
	if (!str)
		throw std::runtime_error((boost::format("Could not write material %s!") % mat.name).str());
	return str;
}

TGeometry::TGeometry(TConfig &geometryin){
	boost::filesystem::path matpath;
	std::istringstream(geometryin["GLOBAL"]["materials_file"]) >> matpath; // check if there is a materials file linked in the config file
	if (matpath.empty()){ // if there is no materials file given load materials from config file
		std::cout << "Loading materials from " << configpath << "\n";
	}
	else{
		matpath = boost::filesystem::absolute(matpath, configpath.parent_path()); // make path absolute, relative paths are assumed to be relative to the config file's path
		std::cout << "Loading materials from " << matpath << "\n";
		geometryin.ReadFromFile(matpath.native());
	}

	std::map<std::string, material> materials;
	for (auto materialParameters: geometryin["MATERIALS"]){
		material mat;
		std::istringstream(materialParameters.first + " " + materialParameters.second) >> mat;
		materials[materialParameters.first] = mat;
	}

	std::map<std::string, TPathInterpolator<std::vector<double> > > paths;
	for (auto pathParameters: geometryin["TRANSFORMATIONS"]){
		paths.emplace(pathParameters.first, constructPathInterpolator(pathParameters.second));
	}

	for (auto sldparams : geometryin["GEOMETRY"]){
		unsigned ID;
		std::istringstream(sldparams.first) >> ID;
		if (ID == 1){
			std::string ignored;
			std::string materialName;
			std::istringstream(sldparams.second) >> ignored >> materialName;
			defaultMaterial = materials.at(materialName);
		}
		else{
			solids.push_back(constructBody(ID, sldparams.second, materials, paths));
		}
	}

	if (std::unique(solids.begin(), solids.end(), [](auto &s1, auto &s2){ return s1.ID == s2.ID; }) != solids.end()) // check if IDs of each solid are unique
		throw std::runtime_error("You defined solids with identical ID! IDs have to be unique!");

}

bool TGeometry::CheckSegment(const double tstart, const std::array<double, 3> &start, const double tend, const std::array<double, 3> &end) const{
	for (auto &body: solids){
		if (body.path and inMovingBoundingBox(body.mesh, tstart, start, tend, end, *body.path))
			return true;
		else if (not body.path and inBoundingBox(body.mesh, start, end))
			return true;
	}
	return false;
}


bool TGeometry::GetCollisions(const double tstart, const std::array<double, 3> &start, const double tend, const std::array<double, 3> &end, std::multimap<TCollision, bool> &colls) const{
	for (auto &body: solids){
		if (body.path){
			auto collisions = findCollisionsWithMovingMesh(body.mesh, tstart, start, tend, end, *body.path);
			for (auto &coll: collisions){
				colls.emplace(TCollision{(coll.collisionTime - tstart)/(tend - tstart), coll.surfaceNormal, coll.surfaceVelocity, body.ID}, false);
			}
		}
		else{
			auto collisions = findCollisionsWithMesh(body.mesh, start, end);
			for (auto &coll: collisions){
				auto s = boost::qvm::mag(coll.first - start)/boost::qvm::mag(end - start);
				colls.emplace(TCollision{tstart * s*(tend - tstart), coll.second, {0, 0, 0}, body.ID}, false);
			}
		}
	}
	return !colls.empty();
}


std::vector<std::pair<solid, bool> > TGeometry::GetSolids(const double t, const std::array<double, 3> &p) const{
	std::vector<std::pair<solid, bool> > solidList{{{1, defaultMaterial}, false}};
	std::array<double, 3> pv{p[0], p[1], p[2]};
	for (auto &body: solids){
		if ((body.path and inMovingSolid(body.mesh, t, pv, *body.path)) or (not body.path and inSolid(body.mesh, pv))){
			solidList.push_back({{body.ID, body.material}, false});
		}
	}
	return solidList;
}
