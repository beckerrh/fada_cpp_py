//
//  transferinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferinterface_h
#define transferinterface_h

#include  <memory>
#include  "typedefs.hpp"

class BoundaryConditions;
class VectorInterface;
/*-------------------------------------------------*/
class TransferInterface
{
public:
  virtual ~TransferInterface()
  {
  }

  TransferInterface()
  {
  }

  TransferInterface(const TransferInterface& transfer)
  {
  }
  virtual void restrict (std::shared_ptr <VectorInterface> out, std::shared_ptr <VectorInterface const> in) const = 0;
  virtual void prolongate(std::shared_ptr <VectorInterface> out, std::shared_ptr <VectorInterface const> in) const = 0;
};

/*-------------------------------------------------*/
template <typename TRANSFER, class VECTOR>
class Transfer : public TRANSFER, public TransferInterface
{
protected:
  TRANSFER& get()
  {
    return(static_cast <TRANSFER&>(*this));
  }

  TRANSFER const& get() const
  {
    return(static_cast <TRANSFER const&>(*this));
  }

  const VECTOR& getVector(std::shared_ptr <VectorInterface const> u) const
  {
    return(static_cast <const VECTOR&>(*u));
  }

  VECTOR& getVector(std::shared_ptr <VectorInterface> u) const
  {
    return(static_cast <VECTOR&>(*u));
  }

public:
  Transfer <TRANSFER, VECTOR>() : TRANSFER(), TransferInterface()
  {
  }

  Transfer <TRANSFER, VECTOR>(const armaicvec& n, const armavec& dx) : TRANSFER(n, dx), TransferInterface()
  {
  }
  Transfer <TRANSFER, VECTOR>(const arma::umat& locations, const armavec& values) : TRANSFER(locations, values), TransferInterface()
  {
  }
  

  // void set_grid(const armaicvec& n, const armavec& dx)
  // {
  //   get().set_grid(n, dx);
  // }
  //
  void restrict (std::shared_ptr <VectorInterface> out, std::shared_ptr <VectorInterface const> in) const { get().restrict (getVector(out), getVector(in)); }
  void prolongate(std::shared_ptr <VectorInterface> out, std::shared_ptr <VectorInterface const> in) const
  {
    get().prolongate(getVector(out), getVector(in));
  }
};

#endif /* transferinterface_h */
