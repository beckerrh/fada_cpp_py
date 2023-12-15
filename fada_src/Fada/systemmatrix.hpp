//
//  matrixinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef systemmatrix_h
#define systemmatrix_h

#include  "matrixinterface.hpp"

/*-------------------------------------------------*/
class SystemMatrix : public  virtual MatrixInterface
{
protected:
    armaicvec _ofs;
    std::map<std::string,std::shared_ptr<MatrixInterface>> _matrices;
    std::map<std::pair<int,int>,std::pair<std::string, std::string>> _patterns;
public:
    SystemMatrix():_matrices(), _patterns() {}
    SystemMatrix(const SystemMatrix& matrix):_matrices(matrix._matrices), _patterns(matrix._patterns) {}
    SystemMatrix(std::map<std::string,std::shared_ptr<MatrixInterface>> matrices, std::map<std::pair<int,int>,std::pair<std::string, std::string>> patterns);
    const std::map<std::string,std::shared_ptr<MatrixInterface>>&  get_matrices() const {return _matrices;}
    const std::map<std::pair<int,int>,std::pair<std::string, std::string>>& get_patterns() const {return _patterns;}
    void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1.0) const;
    void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1.0) const;
    int nrows() const {return _ofs.back();}
    int ncols() const {return _ofs.back();}
    int nelem() const {_not_written_();return 0;}
    const armaicvec& get_ofs() const {return _ofs;}
};


#endif
