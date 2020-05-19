//
//  typedefs.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef typedefs_h
#define typedefs_h

//#define ARMA_USE_SUPERLU 1
//#define ARMA_SUPERLU_INCLUDE_DIR /usr/local/include/superlu/

#include  <armadillo>
#include  "array.hpp"
#include  "vector.hpp"


typedef arma::Col<double> armavec;
typedef arma::Mat<double> armamat;
typedef arma::Col<int> armaicvec;
typedef arma::Row<int> armairvec;
typedef arma::Mat<int> armaimat;
typedef double (*function2d)(double, double);
typedef double (*function3d)(double, double, double);
typedef Array<Vector> VectorMG;


#endif /* typedefs_h */
