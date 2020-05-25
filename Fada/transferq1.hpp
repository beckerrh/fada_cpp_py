
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

/*-------------------------------------------------*/
class TransferQ12d : public TransferInterface
{
protected:
  int _nx, _ny;
  void _boundary(Vector& u) const;

public:
  ~TransferQ12d();
  TransferQ12d() : TransferInterface() {}
  TransferQ12d(const TransferQ12d& transfer) : TransferInterface(transfer) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void restrict(Vector& out, const Vector& in) const;
  void prolongate(Vector& out, const Vector& in) const;
};


/*-------------------------------------------------*/
class TransferQ13d : public TransferInterface
{
protected:
  int _nx, _ny, _nz;
  void _boundary(Vector& u) const;
  
public:
  ~TransferQ13d();
  TransferQ13d() : TransferInterface() {}
  TransferQ13d(const TransferQ13d& transfer) : TransferInterface(transfer) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void restrict(Vector& out, const Vector& in) const;
  void prolongate(Vector& out, const Vector& in) const;
};

#endif /* transferq1_h */
