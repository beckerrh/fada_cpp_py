//
//  stencil2d.hpp
//  Fada
//
//  Created by Roland Becker on 01/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef stencil_hpp
#define stencil_hpp

#include  "../gridvector.hpp"
#include  "../seamvector.hpp"

class MatrixInterface;
/*-------------------------------------------------*/
template<int DIM, int N>
class Stencil
{
protected:
  arma::vec::fixed<N> _coef;
  mutable SeamVector _seam;
  int _nx, _ny, _nz;
  std::string _smoother;
  void _boundary(GridVector& out) const
  {
  // std::cerr << "_boundary() _nx="<<_nx << "_ny=" << _ny <<"\n";
  //  std::cerr << "Stencil2d9() _coef="<<_coef.t();
  if(DIM==2)
  {
    for (int ix = 0; ix < _nx; ix++)
    {
      out.at(ix, 0)       = 0.0;
      out.at(ix, _ny - 1) = 0.0;
    }
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(0, iy)       = 0.0;
      out.at(_nx - 1, iy) = 0.0;
    }
  }
  else
  {
    for (int ix = 0; ix < _nx; ix++)
    {
      for (int iy = 0; iy < _ny; iy++)
      {
        out.at(ix, iy, 0)       = 0.0;
        out.at(ix, iy, _nz - 1) = 0.0;
      }
    }
    for (int ix = 0; ix < _nx; ix++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, 0, iz)       = 0.0;
        out.at(ix, _ny - 1, iz) = 0.0;
      }
    }
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(0, iy, iz)       = 0.0;
        out.at(_nx - 1, iy, iz) = 0.0;
      }
    }
  }
}

public:
  Stencil<DIM,N>() : _coef(), _seam() {}
  Stencil<DIM,N>(const Stencil<DIM,N>& stencil) : _coef(stencil._coef), _seam(stencil._seam), _nx(stencil._nx), _ny(stencil._ny), _nz(stencil._nz), _smoother(stencil._smoother) {}
  Stencil<DIM,N>(const armaicvec& n, const armavec& coef, std::string smoother) : _coef(coef), _smoother(smoother)
  {
    _seam.set_size(n + 2);
    assert(coef.n_elem==N);
    // _coef = coef;
    if(DIM==2)
    {
      assert(n.n_elem == 2);
      _nx   = n[0];
      _ny   = n[1];
    }
    else
    {
      assert(n.n_elem == 3);
      _nx   = n[0];
      _ny   = n[1];
      _nz   = n[2];
    }
  }
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{out<<_coef;}
  void set_grid(const armaicvec& n, const armavec& dx){assert(0); exit(1);}
  void presmooth(GridVector& out, const GridVector& in) const
  {
    if(_smoother=="Jac")
    {
      this->jacobi(out, in);
    }
    else if(_smoother=="GS")
    {
      this->gauss_seidel2(out, in);
    }
    else if(_smoother=="GS1")
    {
      this->gauss_seidel1(out, in);
    }
    else if(_smoother=="GS2")
    {
      this->gauss_seidel2(out, in);
    }
    else
    {
      throw std::runtime_error("Stencil<DIM,N>: unknown smoother " + _smoother);
    }
  }
  void postsmooth(GridVector& out, const GridVector& in) const
  {
    if(_smoother=="Jac")
    {
      this->jacobi(out, in);
    }
    else if(_smoother=="GS")
    {
      this->gauss_seidel1(out, in);
    }
    else if(_smoother=="GS1")
    {
      this->gauss_seidel1(out, in);
    }
    else if(_smoother=="GS2")
    {
      this->gauss_seidel2(out, in);
    }
    else
    {
      throw std::runtime_error("Stencil<DIM,N>: unknown smoother " + _smoother);
    }
  }
  virtual void jacobi       (GridVector& out, const GridVector& in) const=0;
  virtual void gauss_seidel1(GridVector& out, const GridVector& in) const=0;
  virtual void gauss_seidel2(GridVector& out, const GridVector& in) const=0;
  void update(std::shared_ptr<MatrixInterface const> matrix) {_not_written_();}
  void Tdot(GridVector& out, const GridVector& in, double d) const{_not_written_();}
};

/*-------------------------------------------------*/
class Stencil2d9 : public Stencil<2,9>
{
public:
  Stencil2d9() : Stencil<2,9>() {}
  Stencil2d9(const Stencil2d9& stencil) : Stencil<2,9>(stencil) {}
  Stencil2d9(const armaicvec& n, const armavec& coef, std::string smoother) : Stencil<2,9>(n, coef, smoother) {}

  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

/*-------------------------------------------------*/
class Stencil2d5 : public Stencil<2,5>
{
public:
  Stencil2d5() : Stencil<2,5>() {}
  Stencil2d5(const Stencil2d5& stencil) : Stencil<2,5>(stencil) {}
  Stencil2d5(const armaicvec& n, const armavec& coef, std::string smoother) : Stencil<2,5>(n, coef, smoother) {}

  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

/*-------------------------------------------------*/
class Stencil3d27 : public Stencil<3,27>
{
public:
  Stencil3d27() : Stencil<3,27>() {}
  Stencil3d27(const Stencil3d27& stencil) : Stencil<3,27>(stencil) {}
  Stencil3d27(const armaicvec& n, const armavec& coef, std::string smoother) : Stencil<3,27>(n, coef, smoother) {}

  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

/*-------------------------------------------------*/
class Stencil3d7 : public Stencil<3,7>
{
public:
  Stencil3d7() : Stencil<3,7>() {}
  Stencil3d7(const Stencil3d7& stencil) : Stencil<3,7>(stencil) {}
  Stencil3d7(const armaicvec& n, const armavec& coef, std::string smoother) : Stencil<3,7>(n, coef, smoother) {}

  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

#endif
