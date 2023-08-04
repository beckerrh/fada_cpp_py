//
//  StokesVector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef StokesVector_h
#define StokesVector_h

#include  <armadillo>
#include  "typedefs.hpp"
#include  "gridvector.hpp"
#include  "vectorinterface.hpp"

/*-------------------------------------------------*/
class StokesVector : public SystemVector
{
protected:
    int _dim;
public:
    StokesVector(const StokesVector& svector):SystemVector(svector), _dim(svector._dim) {}
    StokesVector(int dim): SystemVector(dim+1), _dim(dim) {}
    StokesVector& operator=(const StokesVector& v)
    {
        _not_written_();
        return *this;
    }
    std::shared_ptr<GridVector> get_p() {return std::dynamic_pointer_cast<GridVector>(get(_dim));}
    std::shared_ptr<GridVector const> get_p() const {return std::dynamic_pointer_cast<GridVector const>(get(_dim));}
    std::shared_ptr<GridVector> get_v(int i) {return std::dynamic_pointer_cast<GridVector>(get(i));}
    std::shared_ptr<GridVector const> get_v(int i) const {return std::dynamic_pointer_cast<GridVector const>(get(i));}
};

#endif
