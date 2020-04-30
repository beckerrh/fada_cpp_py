#ifndef __operateur_h
#define __operateur_h

#include "vecteur.h"
#include "vecteurmg.h"
#include "info.h"

class Operateur
{
private:

  int            o_smoother;
  int            o_levels, o_n0, o_m0;
  Array<int>     o_n, o_m;

  Array<Vecteur>          ocgmem;
  Array<VecteurMG>        omgmem;
  VecteurMG               umg;
 
public:

  Operateur(int s, int l, int n, int m);

  void reinit();
  int minlevel() const   { return 0;}
  int levels()   const   { return o_levels;}
  int n(int l)   const   { return o_n(l);}
  int m(int l)   const   { return o_m(l);}
  int n()        const   { return o_n(levels()-1);}
  int m()        const   { return o_m(levels()-1);}
  Array<Vecteur>& cgmem()  { return ocgmem; }

  int smoother() const   { return o_smoother;}

  void reinit (VecteurMG& v) const;
//  void cg2mg  (VecteurMG&, const Vecteur&);
//  void mg2cg  (Vecteur&, const VecteurMG&) const;


  void vmult(Vecteur& out, const Vecteur& in) const;
  void solve(Vecteur& out, const Vecteur& in, Info& info);


  void jacobi            (int l, VecteurMG&, double);
  void gauss_seidel_pre  (int l, VecteurMG&);
  void gauss_seidel_post (int l, VecteurMG&);
  void smooth_pre        (int l, VecteurMG&);
  void smooth_post       (int l, VecteurMG&);
  void smooth            (int l, VecteurMG&);
//  void smooth            (VecteurMG&);
  void restrict          (int l, VecteurMG&);
//  void restrict          (VecteurMG&);
  void prolongate        (int l, VecteurMG&);
//  void prolongate        (VecteurMG&);
  void residual          (int l, VecteurMG&, VecteurMG&, VecteurMG&);
};

#endif
