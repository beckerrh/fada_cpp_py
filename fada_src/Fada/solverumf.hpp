//
//  solverumf.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef solverumf_hpp
#define solverumf_hpp

#include  "coarsesolverinterface.hpp"
#include  "smootherinterface.hpp"
#include  "umfmatrix.hpp"
#include  <string>

/*-------------------------------------------------*/
class SolverUmf
{
protected:
    std::shared_ptr<MatrixInterface const> _matrix;
    UmfMatrix _umfmat;
    void set_matrix(std::shared_ptr<MatrixInterface const> matrix);

public:
    SolverUmf(std::shared_ptr<MatrixInterface const> matrix, std::string type) {set_matrix(matrix);}

    void solve(armavec& out, const armavec& in) const;
    void update(std::shared_ptr<MatrixInterface const> matrix);
};

/*-------------------------------------------------*/
class SolverUmfSystem : public CoarseSolverInterface
{
protected:
    std::shared_ptr<MatrixInterface const> _matrix;
    UmfMatrix _umfmat;
    void set_matrix(std::shared_ptr<MatrixInterface const> matrix);
    armaicvec _ofs;
    mutable armavec _in, _out;
    void from(std::shared_ptr<VectorInterface const> in) const;
    void to(std::shared_ptr<VectorInterface> out) const;

public:
    SolverUmfSystem(std::shared_ptr<MatrixInterface const> matrix, std::string type=""): CoarseSolverInterface() {set_matrix(matrix);}

    void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
    void update(std::shared_ptr<MatrixInterface const> matrix);
};

#endif