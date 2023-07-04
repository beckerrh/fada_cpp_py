//
//  smootherumf.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef smootherumf_hpp
#define smootherumf_hpp

#include  "smootherinterface.hpp"
#include  "umfmatrix.hpp"
#include  <string>

/*-------------------------------------------------*/
class SmootherUmf
{
protected:
  std::shared_ptr<MatrixInterface const> _matrix;
  UmfMatrix _umfmat;
  void set_matrix(std::shared_ptr<MatrixInterface const> matrix);

public:
  SmootherUmf() {}
  SmootherUmf(const SmootherUmf& smoother) {}
  SmootherUmf(std::shared_ptr<MatrixInterface const> matrix, std::string type) {set_matrix(matrix);}

  void solve(armavec& out, const armavec& in) const;
};

#endif /* smootherumf_hpp */
