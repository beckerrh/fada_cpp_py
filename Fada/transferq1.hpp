
//
//  transferq1.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferq1_h
#define transferq1_h

#include  "transferinterface.hpp"
#include  "typedefs.hpp"
#include  "nodevector.hpp"

/*-------------------------------------------------*/
class TransferQ12d
{
protected:
  int _nx, _ny;
  void _boundary(NodeVector& u) const;

public:
  ~TransferQ12d();
  TransferQ12d() {}
  TransferQ12d(const TransferQ12d& transfer) {}
  TransferQ12d(const armaicvec& n, const armavec& dx)
  {
    set_grid(n, dx);
  }

  void set_grid(const armaicvec& n, const armavec& dx);
  void restrict(NodeVector& out, const NodeVector& in) const;
  void prolongate(NodeVector& out, const NodeVector& in) const;
};


/*-------------------------------------------------*/
class TransferQ13d
{
protected:
  int _nx, _ny, _nz;
  void _boundary(NodeVector& u) const;
  
public:
  ~TransferQ13d();
  TransferQ13d() {}
  TransferQ13d(const TransferQ13d& transfer) {}
  TransferQ13d(const armaicvec& n, const armavec& dx)
  {
    set_grid(n, dx);
  }

  void set_grid(const armaicvec& n, const armavec& dx);
  void restrict(NodeVector& out, const NodeVector& in) const;
  void prolongate(NodeVector& out, const NodeVector& in) const;
};

#endif /* transferq1_h */
