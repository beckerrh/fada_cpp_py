#include  <sstream>
#include  "smoother_chorin.hpp"
#include  "stokesvector.hpp"
#include  "../typedefs.hpp"
#include  "../gridinterface.hpp"
#include  "../matrixinterface.hpp"
#include  "../vectorinterface.hpp"

/*-------------------------------------------------*/
Smoother_Chorin::Smoother_Chorin(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) : SmootherInterface()
{
    auto smatrix = std::dynamic_pointer_cast<SystemMatrix>(matrix);
    assert(smatrix);
    _dim = grid->dim();
    _A.resize(_dim);
    _B.resize(_dim);
    _SV.resize(_dim);
    _D = std::make_shared<SparseMatrix>();
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "A" << i; 
        _A[i] = std::dynamic_pointer_cast<SparseMatrix>(smatrix->get_matrices().at(ss.str()));
        _SV[i] = std::make_shared<SmootherSimple<SparseMatrix>>("GS", _A[i]);
        std::stringstream ss2;
        ss2 << "B" << i; 
        _B[i] = std::dynamic_pointer_cast<SparseMatrix>(smatrix->get_matrices().at(ss2.str()));
        _D->addBBT(*_B[i], 1.0, _A[i]->getDiag());  
    }
    _SP = std::make_shared<SmootherSimple<SparseMatrix>>("GS",_D);
    _helpp.resize(_B[0]->ncols());
}


/*-------------------------------------------------*/
void Smoother_Chorin::presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
{
    auto sout = std::dynamic_pointer_cast<StokesVector>(out);
    assert(sout);
    auto sin = std::dynamic_pointer_cast<StokesVector const>(in);
    assert(sin);
    // _helpp.fill(0.0);
    _helpp = *std::dynamic_pointer_cast<armavec const>(sin->get_p());
    for(int i=0;i<_dim;i++)
    {
        _SV[i]->presmooth(*sout->get_v(i), *sin->get_v(i));
        _B[i]->dot(_helpp, *sout->get_v(i), -1.0);
    }
    _SP->presmooth(*sout->get_p(), _helpp);
    for(int i=0;i<_dim;i++)
    {
        _B[i]->Tdot(*sout->get_v(i), *sout->get_p(),1.0);
    }
}
/*-------------------------------------------------*/
void Smoother_Chorin::postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
{
    presmooth(out, in);
}
/*-------------------------------------------------*/
void Smoother_Chorin::update(std::shared_ptr<MatrixInterface const> matrix)
{
    _not_written_();
}
