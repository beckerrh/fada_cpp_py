//
//  solverumf.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverumf.hpp"
#include  "matrixinterface.hpp"
#include  "vectorinterface.hpp"
#include  "sparsematrix.hpp"
#include  "systemmatrix.hpp"
#include  "systemvector.hpp"

/*-------------------------------------------------*/
void SolverUmf::update(std::shared_ptr<MatrixInterface const> matrix) 
{
    set_matrix(matrix);
}

/*-------------------------------------------------*/
void SolverUmf::set_matrix(std::shared_ptr<MatrixInterface const> matrix)
{
    _umfmat.init(matrix);
    _umfmat.computeLu();
}

/*-------------------------------------------------*/
void SolverUmf::solve(armavec& out, const armavec& in) const
{
    _umfmat.solve(out, in);
}


/*-------------------------------------------------*/
void SolverUmfSystem::update(std::shared_ptr<MatrixInterface const> matrix) 
{
    set_matrix(matrix);
}

/*-------------------------------------------------*/
void SolverUmfSystem::set_matrix(std::shared_ptr<MatrixInterface const> matrix)
{
    auto p = std::dynamic_pointer_cast<SystemMatrix const>(matrix);
    assert(p);
    _ofs = p->get_ofs();
    _in.resize(_ofs.back());
    _out.resize(_ofs.back());
    auto psm = std::make_shared<SparseMatrix>(p);
    // std::cerr << "SparseMatrix\n";
        // psm->save(std::cerr);  
    _umfmat.init(psm);
    _umfmat.computeLu();
}
/*-------------------------------------------------*/
void SolverUmfSystem::from(std::shared_ptr<VectorInterface const> in) const
{
    auto pin = std::dynamic_pointer_cast<SystemVector const>(in);
    assert(pin);
    for(int i=0;i<_ofs.n_elem-1;i++)
    {
        auto pini = std::dynamic_pointer_cast<armavec const>(pin->get(i));
        _in(arma::span(_ofs[i],_ofs[i+1]-1)) = *pini;
    }
}
void SolverUmfSystem::to(std::shared_ptr<VectorInterface> out) const
{
    auto pout = std::dynamic_pointer_cast<SystemVector>(out);
    assert(pout);
    for(int i=0;i<_ofs.n_elem-1;i++)
    {
        auto pouti = std::dynamic_pointer_cast<armavec>(pout->get(i));
        *pouti = _out(arma::span(_ofs[i],_ofs[i+1]-1));
    }    
}

/*-------------------------------------------------*/
void SolverUmfSystem::solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
{
    from(in);
    _umfmat.solve(_out, _in);
    to(out);
}
