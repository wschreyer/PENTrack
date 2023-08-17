// #ifndef PYTHONFIELDS_H_
// #define PYTHONFIELDS_H_

// #define PY_SSIZE_T_CLEAN
// #include <Python.h>

#include "field.h"
#include <vector>
#include <iostream>
#include "boost/format.hpp"
#include "pythonFields.h"



TPythonField::TPythonField(const std::string ft, const bool tp){

  temporal = tp;
  std::cout << "temporal = " << temporal << std::endl;
  // PyRun_SimpleString("import sys"); 
  // PyRun_SimpleString("sys.path.append(\".\")");

  // PyObject* magpymodule = PyImport_ImportModule("magpylib");
  // pBFieldFunc = PyObject_GetAttrString(magpymodule, "getB");

  std::cout << " pythonField function loading" << std::endl;

  boost::python::object bmagpymodule = boost::python::import("magpylib");
  bpBFieldFunc = bmagpymodule.attr("getB");
  // boost::python::object bmagpymodule = boost::python::import("pythonField");
  // bpBFieldFunc = bmagpymodule.attr("BField");

  std::cout << " python magnetic source building" << std::endl;
  
  // char* ftc = const_cast<char*>(ft.c_str());
  boost::python::str bftc(ft);
  boost::python::object bmagnetmodule = boost::python::import(bftc);
  // boost::python::object buildSourceFunc = bmagnetmodule.attr("buildSource");
  buildSourceFunc = bmagnetmodule.attr("buildSource");

  try
    {
      bpSourceObject = buildSourceFunc(0);
    }
  catch (const std::exception& e)
    {
      // Catch and handle C++ exceptions
      std::cout << "Caught C++ exception: " << e.what() << std::endl;
    }
  catch (...)
    {
      // Catch and handle any other C++ exceptions
      std::cout << "Caught unknown C++ exception" << std::endl;
    }
  
  
  // bpSourceObject = boost::python::call_method<boost::python::object>(, NULL);
  
  // PyObject* magnetmodule = PyImport_ImportModule(ftc);
  // PyObject* pMagnetFunc = PyObject_GetAttrString(magnetmodule, "BuildMagnet");
  // pMagnetObject = PyObject_CallObject(pMagnetFunc, NULL);
  // Py_DECREF(magnetmodule);
  // Py_DECREF(pMagnetFunc);
      
}




void TPythonField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {

  const double h = 1e-8;
  int dimarg = (dBidxj != nullptr) ? 7 : 1;
  double xyz[3];
  boost::python::list bpList;
  
  for(int i=0; i<dimarg; ++i){

    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    
    switch(i) {
    case 0:
      break;
    case 1:
      xyz[0] -= h;
      break;
    case 2:
      xyz[0] += h;
      break;
    case 3:
      xyz[1] -= h;
      break;
    case 4:
      xyz[1] += h;
      break;
    case 5:
      xyz[2] -= h;
      break;
    case 6:
      xyz[2] += h;
      break;
    default:
      printf("out of case");
    }

    boost::python::list bpxyzList;

    for(int j=0; j<3; ++j){
      bpxyzList.append(1000 * xyz[j]);
    }
    bpList.append(bpxyzList);
  }

    
  double Bs[dimarg][3];
  boost::python::object bnpArray;
  boost::python::object bpSourceObjectlocal;

  if (temporal)
    bpSourceObjectlocal = buildSourceFunc(t);
  else
    bpSourceObjectlocal = bpSourceObject;

  bnpArray = bpBFieldFunc(bpSourceObjectlocal, bpList);

  
  if (dBidxj != nullptr){
    for (int i=0; i<dimarg; ++i) {
      for (int j=0; j<3; ++j) { 
	Bs[i][j] = 0.001 * boost::python::extract<double>(bnpArray[i][j]);
      }
    }
  }
  else{
    for (int j=0; j<3; ++j) {
      Bs[0][j] = 0.001 * boost::python::extract<double>(bnpArray[j]);
    }
  }

  for (int i=0; i<3; ++i) {
    B[i] = Bs[0][i];
  }
  
  if (dBidxj != nullptr){

    double trace_3;
    double dBi_dxj[3][3];

    for(int i=0; i<3; ++i){
      dBi_dxj[i][0] = (Bs[2][i] - Bs[1][i])/(2*h);
      dBi_dxj[i][1] = (Bs[4][i] - Bs[3][i])/(2*h);
      dBi_dxj[i][2] = (Bs[6][i] - Bs[5][i])/(2*h); 
    }
    
    trace_3 = 0; //(dBi_dxj[0][0] + dBi_dxj[1][1] + dBi_dxj[2][2])/3;
    // std::cout << trace_3 << std::endl;
    
    dBidxj[0][0] = dBi_dxj[0][0] - trace_3;
    dBidxj[1][1] = dBi_dxj[1][1] - trace_3;
    dBidxj[2][2] = dBi_dxj[2][2] - trace_3;

    dBidxj[1][0] = (dBi_dxj[0][1] + dBi_dxj[1][0])/2;
    dBidxj[2][0] = (dBi_dxj[0][2] + dBi_dxj[2][0])/2;
    dBidxj[2][1] = (dBi_dxj[1][2] + dBi_dxj[2][1])/2;

    dBidxj[0][1] = dBidxj[1][0];
    dBidxj[1][2] = dBidxj[2][1];
    dBidxj[0][2] = dBidxj[2][0];
    
  }    
}


// void TPythonField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {

//   const double h = 1e-8;
//   int dimarg = (dBidxj != nullptr) ? 7 : 1;
//   double xyz[3];
//   boost::python::list bpList;
  
//   for(int i=0; i<dimarg; ++i){

//     xyz[0] = x;
//     xyz[1] = y;
//     xyz[2] = z;
    
//     switch(i) {
//     case 0:
//       break;
//     case 1:
//       xyz[0] -= h;
//       break;
//     case 2:
//       xyz[0] += h;
//       break;
//     case 3:
//       xyz[1] -= h;
//       break;
//     case 4:
//       xyz[1] += h;
//       break;
//     case 5:
//       xyz[2] -= h;
//       break;
//     case 6:
//       xyz[2] += h;
//       break;
//     default:
//       printf("out of case");
//     }

//     boost::python::list bpxyzList;

//     for(int j=0; j<3; ++j){
//       bpxyzList.append(1000 * xyz[j]);
//     }
//     bpList.append(bpxyzList);
//   }

//   boost::python::object bnpArray;
//   boost::python::object bpSourceObjectlocal;

//   if (temporal)
//     bpSourceObjectlocal = buildSourceFunc(t);
//   else
//     bpSourceObjectlocal = bpSourceObject;

//   bnpArray = bpBFieldFunc(bpSourceObjectlocal, bpList);
  
//   if (dBidxj != nullptr){
//     double Bs[7][3];
//     for (int i=0; i<7; ++i) {
//       for (int j=0; j<3; ++j) { 
//         Bs[i][j] = 0.001 * boost::python::extract<double>(bnpArray[i][j]);
//       }
//     }

//     for (int i=0; i<3; ++i)
//       B[i] = Bs[0][i];
  

//     double dBi_dxj[3][3];
//     double trace_3 = 0.0;

//     for(int i=0; i<3; ++i){
//       dBi_dxj[i][0] = (Bs[2][i] - Bs[1][i])/(2*h);
//       dBi_dxj[i][1] = (Bs[4][i] - Bs[3][i])/(2*h);
//       dBi_dxj[i][2] = (Bs[6][i] - Bs[5][i])/(2*h); 
//       trace_3 += dBi_dxj[i][i];
//     }
    
//     trace_3 /= 3.0;

//     for(int i=0; i<3; ++i){
//       dBidxj[i][0] = dBi_dxj[i][0] - trace_3;
//       dBidxj[i][1] = (dBi_dxj[i][1] + dBi_dxj[1][i])/2.0;
//       dBidxj[i][2] = (dBi_dxj[i][2] + dBi_dxj[2][i])/2.0;
//       dBidxj[0][i] = dBidxj[i][0];
//       dBidxj[1][i] = dBidxj[i][1];
//       dBidxj[2][i] = dBidxj[i][2];
//     }
//   }
//   else{
//     for (int j=0; j<3; ++j) {
//       B[j] = 0.001 * boost::python::extract<double>(bnpArray[j]);
//     }
//   }
// }




    

  //   double reltrace = (dBidxj[0][0] + dBidxj[1][1] + dBidxj[2][2])/sqrt(dBidxj[0][0]*dBidxj[0][0] + dBidxj[1][1]*dBidxj[1][1] + dBidxj[2][2]*dBidxj[2][2]);
  //   double relrotxy = (dBidxj[0][1] - dBidxj[1][0])/(dBidxj[0][1] + dBidxj[1][0]);
  //   double relrotxz = (dBidxj[0][2] - dBidxj[2][0])/(dBidxj[0][2] + dBidxj[2][0]);
  //   double relrotyz = (dBidxj[1][2] - dBidxj[2][1])/(dBidxj[1][2] + dBidxj[2][1]);
    
  //   double trace = dBidxj[0][0] + dBidxj[1][1] + dBidxj[2][2];
  //   double rotxy = dBidxj[0][1] - dBidxj[1][0];
  //   double rotxz = dBidxj[0][2] - dBidxj[2][0];
  //   double rotyz = dBidxj[1][2] - dBidxj[2][1];

  //   if(abs(reltrace) > 0.01){ //} || abs(relrotxy > 0.01 || relrotxz > 0.01 || relrotyz > 0.01){
  //     // if(abs(trace) > 1){ // || abs(rotxy) > 0.01 || abs(rotxz) > 0.01 || abs(rotyz) > 0.01){
  //     std::cout << " \ndBidxj, x=" <<x << ", r=" << sqrt(y*y + z*z) << ", y=" << y << ", z=" << z << std::endl;
  //     for(int i=0; i<3;++i){
  // 	for(int j=0; j<3; ++j){
  // 	  std::cout << dBidxj[i][j] << ", ";
  // 	}
  // 	std::cout << std::endl;
  //     }
  //     std::cout << " reltrace = " << reltrace << ", relrotxy = " << relrotxy << ", relrotxz = " << relrotxz << ", relrotyz" << relrotyz << std::endl;
  //     std::cout << " trace = " << trace << ", rotxy = " << rotxy << ", rotxz = " << rotxz << ", rotyz" << rotyz << std::endl


//////////////////////


  // Py_DECREF(npArray);
  // Py_DECREF(pList);
  // PyTuple_SET_ITEM(pArgs, 0, pxyzList);
  // Py_DECREF(pxyzList);
  // Py_DECREF(pArgs);
