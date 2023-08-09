    //
//  GridVector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef GridVector_h
#define GridVector_h

#include  <armadillo>
#include  "boundary_conditions.hpp"
#include  "typedefs.hpp"
#include  "vectorinterface.hpp"

/*-------------------------------------------------*/
class GridVector : public armavec
{
protected:
    bool _mean;
    BoundaryConditionsBool _bc;
    armaicvec _n, _ofs;
    void set_ofs()
    {
        if(_mean)
        {
            for(auto p:_bc)
            {
                if(p[0] or p[1])
                {
                    _not_written_("cannot have dir");
                }
            }
        }
        _ofs.set_size(_n.n_elem);
        _ofs.fill(1);
        for(int i=0;i<_n.n_elem;i++)
        {
            for(int j=i+1;j<_n.n_elem;j++)
            {
                _ofs[i] *= _n[j];
            }
        }
    }

public:
    GridVector() : armavec(), _n(), _ofs(), _bc(nullptr) {}
    GridVector(const armaicvec& n, std::shared_ptr<BoundaryConditions const> bc=nullptr, bool mean=false) : armavec(arma::prod(n)), _n(n), _bc(bc), _mean(mean)
    {
        set_ofs();
    }
    const armavec& get_arma() const {return static_cast<const armavec&>(*this);}
    armavec& get_arma() {return static_cast<armavec&>(*this);}
    
    GridVector& operator=(const GridVector& v)
    {
        armavec::operator=(v);
        assert(arma::all(_n==v.n()));
        _bc = v._bc;
        _n = v._n;
        _ofs = v._ofs;
        _mean = v._mean;
        return *this;
    }
    GridVector& operator=(const armavec& v)
    {
        armavec::operator=(v);
        return *this;
    }
    void set_size(const armaicvec& n)
    {
        //    std::cerr << "GridVector::set_size() n = " << n.t();
        _n = n;
        set_ofs();
        armavec::set_size(arma::prod(n));
    }
    int get_size() const{ return get_arma().n_elem;}
    void set_size(const GridVector& u)
    {
        set_size(u.n());
    }
    void after_prolongate();
    void after_restrict();
    void after_smooth();
    void project_mean();
    void boundary_zero();
    int dim() const {return (int) _n.n_elem;}
    const armaicvec& n() const {return _n;}
    const armaicvec& ofs() const {return _ofs;}
    double& at(int ix, int iy)
    {
        return (*this)[_ofs[0]*ix+iy];
    }
    const double& at(int ix, int iy) const
    {
        return (*this)[_ofs[0]*ix+iy];
    }
    double& at(int ix, int iy, int iz)
    {
        return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
    }
    const double& at(int ix, int iy, int iz) const
    {
        return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
    }
    void scale(double d)
    {
        *this *= d;
    }
    double norm(double p=2) const
    {
        // const armavec& tarma = static_cast<const armavec&>(*this);
        // return arma::norm(tarma,p);
        return arma::norm(get_arma(), p);
    }
    void add(double d, const GridVector& v)
    {
        // armavec& tarma = static_cast<armavec&>(*this);
        // const armavec& varma = static_cast<const armavec&>(v);
        // tarma += d*varma;
        get_arma() += d * v.get_arma();
    }
    void equal(const GridVector& v)
    {
        armavec::operator=(v);
    }
    double dot(const GridVector& v) const
    {
        // const armavec& tarma = static_cast<const armavec&>(*this);
        // const armavec& varma = static_cast<const armavec&>(v);
        // return arma::dot(tarma, varma);
        return arma::dot(get_arma(), v.get_arma());
    }
    void savehdf5(const std::string& filename) const
    {
        save(arma::hdf5_name(filename));
    }
};
std::ostream& operator<<(std::ostream& os, const GridVector& v);

#endif /* GridVector_h */
