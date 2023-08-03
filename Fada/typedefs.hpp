//
//  typedefs.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef typedefs_h
#define typedefs_h

// #define _LONG_LONG

#include  <map>
#include  <utility>
#include  <vector>
#include  <armadillo>
#include  <cassert>

typedef arma::Col<double> armavec;
typedef arma::Mat<double> armamat;
typedef arma::Col<int> armaicvec;
typedef arma::Row<int> armairvec;
typedef arma::Mat<int> armaimat;
typedef double (*function2d)(double, double);
typedef double (*function3d)(double, double, double);

typedef std::map<std::string,std::shared_ptr<armavec>> PointDataMap;

inline void _not_written_(std::string msg="")
{
    std::cerr << "--NotWritten--" << msg << "\n";
    assert(0);
    exit(1);    
}

#endif
