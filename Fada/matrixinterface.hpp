//
//  matrixinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef matrixinterface_h
#define matrixinterface_h

#include  "typedefs.hpp"
#include  "vectorinterface.hpp"

class GridInterface;
/*-------------------------------------------------*/
class MatrixInterface
{
public:
    virtual ~MatrixInterface() {}
    MatrixInterface() {}
    MatrixInterface(const MatrixInterface& matrix) {}

    virtual int nrows() const=0;
    virtual int ncols() const=0;
    virtual int nelem() const=0;

    virtual void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const=0;
    virtual void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const=0;

    virtual void jacobi       (std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{_not_written_();}
    virtual void gauss_seidel1(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{_not_written_();}
    virtual void gauss_seidel2(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{_not_written_();}
    virtual void set_elements(const arma::umat& locations, const armavec& values, bool compute_diag_flag=true) {_not_written_();}
    virtual void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{_not_written_();}
};

/*-------------------------------------------------*/
template<typename MATRIX, class VECTOR>
class Matrix : public  virtual MATRIX, public  virtual MatrixInterface
{
protected:
    const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
    VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}

public:
    Matrix<MATRIX, VECTOR>() : MATRIX(), MatrixInterface() {}
    Matrix<MATRIX, VECTOR>(const MATRIX& matrix) : MATRIX(matrix), MatrixInterface() {}
    Matrix<MATRIX, VECTOR>(const arma::umat& locations, const armavec& values, bool compute_diag=true): MATRIX(locations,values,compute_diag), MatrixInterface(){}

    MATRIX& get() { return static_cast<MATRIX&>(*this); }
    MATRIX const& get() const { return static_cast<MATRIX const&>(*this); }
    void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().dot(getVector(out),getVector(in), d);}
    void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().Tdot(getVector(out),getVector(in), d);}
    void set_elements(const arma::umat& locations, const armavec& values, bool compute_diag_flag=true) {get().set_elements(locations, values,compute_diag_flag);}
    void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{get().save(out, datatype);}
    int nrows() const{return get().nrows();};
    int ncols() const{return get().ncols();};
    int nelem() const{return get().nelem();};
};

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
    void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const
    {
        auto pout = std::dynamic_pointer_cast<SystemVector>(out);
        assert(pout);
        auto pin = std::dynamic_pointer_cast<SystemVector const>(in);
        assert(pin);
        for(auto p:_patterns)
        {
            if(p.second.second=="T")
            {
                _matrices.at(p.second.first)->Tdot(pout->get(p.first.first), pin->get(p.first.second), d);                
            }
            else if(p.second.second=="MT")
            {
                _matrices.at(p.second.first)->Tdot(pout->get(p.first.first), pin->get(p.first.second), -d);                
            }
            else
            {
                _matrices.at(p.second.first)->dot(pout->get(p.first.first), pin->get(p.first.second), d);                
            }
        }
    }
    void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const
    {
        auto pout = std::dynamic_pointer_cast<SystemVector>(out);
        assert(pout);
        auto pin = std::dynamic_pointer_cast<SystemVector const>(in);
        assert(pin);
        for(auto p:_patterns)
        {
            if(p.second.second=="T")
            {
                _matrices.at(p.second.first)->dot(pout->get(p.first.first), pin->get(p.first.second), d);                
            }
            else if(p.second.second=="MT")
            {
                _matrices.at(p.second.first)->dot(pout->get(p.first.first), pin->get(p.first.second), -d);                
            }
            else
            {
                _matrices.at(p.second.first)->Tdot(pout->get(p.first.first), pin->get(p.first.second), d);                
            }
        }
    }
    int nrows() const {return _ofs.back();}
    int ncols() const {return _ofs.back();}
    int nelem() const {_not_written_();return 0;}
    const armaicvec& get_ofs() const {return _ofs;}
};


#endif
