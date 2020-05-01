#ifndef __vecteurmg_h
#define __vecteurmg_h

#include "array.h"
#include "vecteur.h"

/**************************************************/

class VecteurMG
{
private:
  Array<Vecteur>  v_val;

public:

  void reinit(const VecteurMG& u);

  Array<Vecteur>&       val()       {return v_val;}
  const Array<Vecteur>& val() const {return v_val;}
  Vecteur&              operator()(int l)       {return v_val(l);}
  const Vecteur&        operator()(int l) const {return v_val(l);}
};

/**************************************************/

inline void VecteurMG::reinit(const VecteurMG& u)
{
  int lev = u.val().n();
  v_val.reinit(lev);
  for(int l=0;l<lev;l++)
    {
      (*this)(l).reinit(u(l).nx(),u(l).ny());
    }
}

#endif
