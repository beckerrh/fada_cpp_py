#ifndef __vecteurmg_h
#define __vecteurmg_h

#include "array.hpp"
#include "vecteur.hpp"

/*-------------------------------------------------*/
class VecteurMG
{
private:
  Array<Vecteur>  v_val;

public:

  void set_size(const VecteurMG& u);

  Array<Vecteur>&       val()       {return v_val;}
  const Array<Vecteur>& val() const {return v_val;}
  Vecteur&              operator()(int l)       {return v_val(l);}
  const Vecteur&        operator()(int l) const {return v_val(l);}
};

/*-------------------------------------------------*/
inline void VecteurMG::set_size(const VecteurMG& u)
{
  int lev = u.val().n();
  v_val.set_size(lev);
  for(int l=0;l<lev;l++)
    {
//      (*this)(l).set_size(u(l).nx(),u(l).ny());
      (*this)(l).set_size(u(l).n());
    }
}

#endif
