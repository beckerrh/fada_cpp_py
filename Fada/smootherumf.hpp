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
class SmootherUmf : public SmootherInterface
{
protected:
  std::shared_ptr<MatrixInterface> _matrix;
  UmfMatrix _umfmat;

public:
  SmootherUmf() {}
  SmootherUmf(const SmootherUmf& smoother) {}

  void set_matrix(std::shared_ptr<MatrixInterface> matrix);
  void solve(Vector& out, const Vector& in) const;
};

#endif /* smootherumf_hpp */
