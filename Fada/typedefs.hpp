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

#include  <armadillo>

typedef arma::Col<double> armavec;
typedef arma::Mat<double> armamat;
typedef arma::Col<int> armaicvec;
typedef arma::Row<int> armairvec;
typedef arma::Mat<int> armaimat;
typedef double (*function2d)(double, double);
typedef double (*function3d)(double, double, double);


#endif /* typedefs_h */
