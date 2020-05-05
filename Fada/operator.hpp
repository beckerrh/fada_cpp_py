//
//  operateur.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef __operateur_h
#define __operateur_h

#include  "typedefs.hpp"
#include  "vector.hpp"
#include  "vectormg.hpp"

/*-------------------------------------------------*/
class Operator
{
protected:
  armairvec _nall;
  armaimat  _n;
  Array<VectorMG> _mgmem;
  vector _u, _f;

  
  void restrict(int l, VectorMG& out, const VectorMG& in) const;
  void prolongate(int l, VectorMG& out, const VectorMG& in) const;
  void residual(int l, VectorMG& r, const VectorMG& u, const VectorMG& f) const;

  void jacobi       (vector& out, const vector& in) const;
  void gauss_seidel1(vector& out, const vector& in) const;
  void gauss_seidel2(vector& out, const vector& in) const;
  void smooth       (vector& out, const vector& in) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, VectorMG& Aw, double tol);
  void update(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, VectorMG& Aw, bool print=false) const;
  int solve(bool print=true);

public:
  std::string smoother;
  int maxiter, optmem;
  double tol_rel, tol_abs;

  void set_parameters();
  Operator();
  Operator(int nlevels, const armaicvec& n0);
  void set_size(int nlevels, const armaicvec& n0);

  int dim() const {return _n.n_rows;}
  int minlevel() const   { return 0;}
  int maxlevel() const   { return nlevels()-1;}
  int nlevels()  const   { return _n.n_cols;}
  arma::subview_col<int>  n(int l) const { return _n.col(l);}
  arma::subview_col<int>  n() const { return _n.col(nlevels()-1);}
  int nall()     const   { return _nall(nlevels()-1);}
  const vector& get_solution() const {return _u;}
  vector& get_solution() {return _u;}

  void set_size(VectorMG& v) const;

  void dot(vector& out, const vector& in, double d) const;
  int testsolve(bool print=true);
  void boundary(vector& v, const vector& w, double d=1.0) const;
  void boundary(vector& v) const;
  void right(vector& v) const;
};

#endif
