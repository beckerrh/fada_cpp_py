//
//  transferq1.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferq1_h
#define transferq1_h

#include  "gridvector.hpp"
#include  "seamvector.hpp"
#include  "boundary_conditions.hpp"
#include  "typedefs.hpp"


/*-------------------------------------------------*/
class TransferBase
{
protected:
  // std::shared_ptr <BoundaryConditions const> _boundaryconditions;
  mutable SeamVector _seam;
  int _nx, _ny, _nz, _dim;
  // void _boundary(GridVector& u) const;

public:
  ~TransferBase()
  {
  }
  TransferBase()
  {
  }

  TransferBase(const TransferBase& transfer) : _seam(transfer._seam), _nx(transfer._nx), _ny(transfer._ny), _nz(transfer._nz)
  {
  }

  TransferBase(const armaicvec& n, const armavec& dx)
  {
    set_grid(n, dx);
  }

  void set_grid(const armaicvec& n, const armavec& dx);
};

/*-------------------------------------------------*/
class TransferQ12d : public TransferBase
{
protected:
public:
  ~TransferQ12d()
  {
  }

  TransferQ12d(const TransferQ12d& transfer) : TransferBase(transfer)
  {
  }

  TransferQ12d(const armaicvec& n, const armavec& dx) : TransferBase(n, dx)
  {
  }
  TransferQ12d(const arma::umat& locations, const armavec& values) : TransferBase()
  {
      _not_written_();
  }

  void restrict (GridVector & out, const GridVector& in) const;
  void prolongate(GridVector& out, const GridVector& in) const;
};


/*-------------------------------------------------*/
class TransferQ13d : public TransferBase
{
protected:
public:
  ~TransferQ13d()
  {
  }

  TransferQ13d(const TransferQ13d& transfer) : TransferBase(transfer)
  {
  }

  TransferQ13d(const armaicvec& n, const armavec& dx) : TransferBase(n, dx)
  {
  }
  TransferQ13d(const arma::umat& locations, const armavec& values) : TransferBase()
  {
      _not_written_();
  }

  void restrict (GridVector & out, const GridVector& in) const;
  void prolongate(GridVector& out, const GridVector& in) const;
};

#endif
