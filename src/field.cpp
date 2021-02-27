#include <cmath>

#include "field.h"

using namespace std;

double TFieldScaler::scalingFactor(const double t) const{
    *tvar = t;
    return scaler.value();
}

void TFieldScaler::scaleScalarField(const double t, double &F, double dFdxi[3]) const{
    double scaling = scalingFactor(t);
    F *= scaling;
    if (dFdxi != nullptr){
        dFdxi[0] *= scaling;
        dFdxi[1] *= scaling;
        dFdxi[2] *= scaling;
    }
}


void TFieldScaler::scaleVectorField(const double t, double F[3], double dFidxj[3][3]) const{
    for (int i = 0; i < 3; ++i){
        if (dFidxj == nullptr){
            scaleScalarField(t, F[i], nullptr);
        }
        else{
            scaleScalarField(t, F[i], dFidxj[i]);
        }
    }
}


TFieldScaler::TFieldScaler(const std::string &scalingFormula){
    tvar = unique_ptr<double>(new double(0.0));
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("t",*tvar);
    symbol_table.add_constants();
    scaler.register_symbol_table(symbol_table);
    exprtk::parser<double> parser;
    if (not parser.compile(scalingFormula, scaler)){
        throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing formula '" + scalingFormula + "': " + parser.get_error(0).diagnostic);
    }
}


double TFieldBoundaryBox::smthrStp(const double x) const{
    return ((6*x - 15)*x + 10)*x*x*x;
}


double TFieldBoundaryBox::smthrStpDer(const double x) const{
    return ((30*x - 60)*x + 30)*x*x;
}

bool TFieldBoundaryBox::hasBounds() const{
    return xmin < xmax and ymin < ymax and zmin < zmax;
}

bool TFieldBoundaryBox::inBounds(const double x, const double y, const double z) const{
    if (not hasBounds()){
        return true;
    }
    else{
        return x >= xmin and x < xmax and y >= ymin and y < ymax and z >= zmin and z < zmax;
    }
}


void TFieldBoundaryBox::scaleScalarFieldAtBounds(const double x, const double y, const double z, double &F, double dFdxi[3]) const{
    if (not hasBounds() or (F == 0 and dFdxi == nullptr) or (F == 0 and dFdxi[0] == 0 and dFdxi[1] == 0 and dFdxi[2] == 0)){ // skip if no boundary is set or field is zero
        return;
    }
    else if (not inBounds(x, y, z)){ // set field to zero if coordinates are outside bounding box
        F = 0.;
        if (dFdxi != nullptr){
            dFdxi[0] = dFdxi[1] = dFdxi[2] = 0.;
        }
    }
    else if (boundaryWidth > 0){
        // F'(x,y,z) = F(x,y,z)*f(x)*f(y)*f(z)
        // dF'/dx = dF/dx*f(x)*f(y)*f(z) + F*df(x)/dx*f(y)*f(z) --> similar for dF'/dy and dF'/dz

        double Fscale = 1.; // f(x)*f(y)*f(z)
        std::array<double, 3> dFadd = {0., 0., 0.}; // F*df(x_i)/dx_i / f(x_i)

        std::array<double, 3> distanceFromLoBoundary = {(x - xmin)/boundaryWidth, (y - ymin)/boundaryWidth, (z - zmin)/boundaryWidth};
        std::array<double, 3> distanceFromHiBoundary = {(xmax - x)/boundaryWidth, (ymax - y)/boundaryWidth, (zmax - z)/boundaryWidth}; // calculate distance to edges in units of BoundaryWidth
        for (int i = 0; i < 3; ++i){
            if (0 <= distanceFromLoBoundary[i] and distanceFromLoBoundary[i] <= 1){
                Fscale *= smthrStp(distanceFromLoBoundary[i]);
                if (dFdxi != nullptr) dFadd[i] += F*smthrStpDer(distanceFromLoBoundary[i])/boundaryWidth/smthrStp(distanceFromLoBoundary[i]);
            }
            else if (0 <= distanceFromHiBoundary[i] and distanceFromHiBoundary[i] <= 1){
                Fscale *= smthrStp(distanceFromHiBoundary[i]);
                if (dFdxi != nullptr) dFadd[i] += -F*smthrStpDer(distanceFromHiBoundary[i])/boundaryWidth/smthrStp(distanceFromHiBoundary[i]);
            }
        }
        if (Fscale != 1.){
            F *= Fscale; // scale field value
            if (dFdxi != nullptr){
                for (int i = 0; i < 3; i++){
                    dFdxi[i] = dFdxi[i]*Fscale + dFadd[i]*Fscale; // scale derivatives according to product rule
                }
            }
        }
    }
    else{ // coordinates are in bounds but field doesn't need to be scaled
        return;
    }
}


void TFieldBoundaryBox::scaleVectorFieldAtBounds(const double x, const double y, const double z, double F[3], double dFidxj[3][3]) const{
    for (int i = 0; i < 3; ++i){
        if (dFidxj == nullptr){
            scaleScalarFieldAtBounds(x, y, z, F[i], nullptr);
        }
        else{
            scaleScalarFieldAtBounds(x, y, z, F[i], dFidxj[i]);
        }
        
    }
}


TFieldBoundaryBox::TFieldBoundaryBox(const double _xmax, const double _xmin, const double _ymax, const double _ymin, const double _zmax, const double _zmin, const double _boundaryWidth):
                    xmin(_xmin), xmax(_xmax), ymin(_ymin), ymax(_ymax), zmin(_zmin), zmax(_zmax), boundaryWidth(_boundaryWidth) {
    if (boundaryWidth < 0 or 
        boundaryWidth > (xmax - xmin)/2 or 
        boundaryWidth > (ymax - ymin)/2 or 
        boundaryWidth > (zmax - zmin)/2){
        throw std::runtime_error("You defined a field with an invalid boundary. Check field parameters.");
    }

}

void TFieldContainer::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
    if (not boundary->inBounds(x, y, z) or BScaler.scalingFactor(t) == 0.){
        for (int i = 0; i < 3; ++i){
            B[i] = 0.;
            if (dBidxj != nullptr){
                for (int j = 0; j < 3; ++j){
                    dBidxj[i][j] = 0.;
                }
            }
        }
    }
    else{
        field->BField(x, y, z, t, B, dBidxj);
        BScaler.scaleVectorField(t, B, dBidxj);
        boundary->scaleVectorFieldAtBounds(x, y, z, B, dBidxj);
    }
}

void TFieldContainer::EField(const double x, const double y, const double z, const double t, double &V, double Ei[3]) const{
    if (not boundary->inBounds(x, y, z) or EScaler.scalingFactor(t) == 0.){
        V = 0.;
        for (int i = 0; i < 3; ++i){
            Ei[i] = 0.;
        }
    }
    else{
        field->EField(x, y, z, t, V, Ei);
        EScaler.scaleScalarField(t, V, Ei);
        boundary->scaleScalarFieldAtBounds(x, y, z, V, Ei);
    }
}
