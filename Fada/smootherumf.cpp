//
//  smootherumf.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "smootherumf.hpp"
#include  "matrixinterface.hpp"
#include  "vector.hpp"

/*-------------------------------------------------*/
void SmootherUmf::set_matrix(std::shared_ptr<MatrixInterface> matrix)
{
  _matrix = matrix;
  matrix->get_sparse_matrix(_umfmat.getSparseMatrix());
  _umfmat.computeLu();
}
/*-------------------------------------------------*/
void SmootherUmf::solve(Vector& out, const Vector& in) const
{
  _umfmat.solve(out.data(), in.data());
//    u.fill(0);
//    for(int iter=0;iter<100;iter++)
//    {
//      d = f;
//      _umfmat.getSparseMatrix().dot(d.data(), u.data(), -1.0);
//  //    residual(l, d, u, f);
//      double res = d.norm();
//      printf("-- %3d %10.2e\n",iter, res);
//      if(res<1e-14) break;
//      _mgmatrix(l)->jacobi(w,d);
//      u.add(1.0, w);
//  //    _mgupdatesmooth(l)->addUpdate(w, u, d);
//    }
//    _spmat.solve(w.data(), f.data());
//    w.add(-1.0,u);
//    assert(w.norm()<1e-12);
}
