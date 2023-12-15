//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef smoother_chorin_hpp
#define smoother_chorin_hpp


#include  "smootherinterface.hpp"
#include  "smoothersimple.hpp"
#include  "sparsematrix.hpp"

class GridInterface;

/*-------------------------------------------------*/
class Smoother_Chorin : public SmootherInterface
{
protected:
    int _dim;
    mutable armavec _helpp;
    std::vector<std::shared_ptr<SparseMatrix const>> _A, _B;
    std::shared_ptr<SparseMatrix> _D;
    std::vector<std::shared_ptr<SmootherSimple<SparseMatrix>>> _SV;
    std::shared_ptr<SmootherSimple<SparseMatrix>> _SP;

public:
    ~Smoother_Chorin()
    {
    }
    Smoother_Chorin(const Smoother_Chorin& smoother) : SmootherInterface(smoother) {}
    Smoother_Chorin(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix);

    void presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
    void postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
    void update(std::shared_ptr<MatrixInterface const> matrix);
};
#endif  
