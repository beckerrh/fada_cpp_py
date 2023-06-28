//
//  smootherumf.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "smootherumf.hpp"
#include  "matrixinterface.hpp"
#include  "vectorinterface.hpp"
#include  "feminterface.hpp"

/*-------------------------------------------------*/
void SmootherUmf::set_matrix(std::shared_ptr<MatrixInterface const> matrix)
{
  _umfmat.init(matrix);
  // _matrix = matrix;
  // std::shared_ptr<FemAndMatrixInterface const> stencil = std::dynamic_pointer_cast<FemAndMatrixInterface const>(matrix);
  // if(stencil)
  // {
  //   std::shared_ptr<FemAndMatrixInterface const> stencil = std::dynamic_pointer_cast<FemAndMatrixInterface const>(matrix);
  //   arma::umat locations;
  //   armavec values;
  //   stencil->get_locations_values(locations, values);
  //   typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
  //   _umfmat.init(std::shared_ptr<MatrixInterface>(new SparseMatrixDef(locations, values)));
  // }
  // else
  // {
  //   _umfmat.init(matrix);
  // }
  // typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
  // std::shared_ptr<SparseMatrixDef const> spmatrix = std::dynamic_pointer_cast<SparseMatrixDef const>(matrix);
  // if(spmatrix)
  // {
  //   _umfmat.init(spmatrix);
  // }
  // else
  // {
  //   std::shared_ptr<FemAndMatrixInterface const> stencil = std::dynamic_pointer_cast<FemAndMatrixInterface const>(matrix);
  //   arma::umat locations;
  //   armavec values;
  //   stencil->get_locations_values(locations, values);
  //   _umfmat.init(std::shared_ptr<MatrixInterface>(new SparseMatrixDef(locations, values)));
  // }
  _umfmat.computeLu();
}
/*-------------------------------------------------*/
// void SmootherUmf::solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
// {
//   _umfmat.solve(out->data(), in->data());
// }
void SmootherUmf::solve(armavec& out, const armavec& in) const
{
  _umfmat.solve(out, in);
}
