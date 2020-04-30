#ifndef __vecteur_h
#define __vecteur_h

#include "array.h"
#include <stdio.h>

/**************************************************/

class Vecteur
{
private:
  Array<double>  v_val;
  int            v_m,v_n; 

  int& pn()                        {return v_n;}
  int& pm()                        {return v_m;}
  Array<double>& val()             {return v_val;}
  const Array<double>& val() const {return v_val;}

public:

  Vecteur();
  Vecteur(int n, int m);
  Vecteur(const Vecteur&);

  int n() const {return v_n;}
  int m() const {return v_m;}

  void reinit(int nn, int mm)   { pn() = nn; pm() = mm; val().reinit(n()*m()); }
  void reinit(const Vecteur& u) { reinit(u.n(),u.m());}

  double&       operator()(int i, int j)        { return val()(i+n()*j); }
  const double& operator()(int i, int j) const  { return val()(i+n()*j); }
  const double& operator()(int i)        const  { return val()(i); }

  Vecteur& operator=(const Vecteur&);
  Vecteur& operator=(double);

  double operator* (const Vecteur& v) const;
  void   equ (double, const Vecteur& v);
  void   add (double, const Vecteur& v);
  void   sadd(double, double, const Vecteur& v);

  void   boundary(const Vecteur& v);
  void   boundary();
  void   right();
  void   null();

  void   output_plotmtv();
};

#endif
