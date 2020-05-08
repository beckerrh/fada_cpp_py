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
#include  "array.hpp"
#include  "vector.hpp"
#include  "updaterinterface.hpp"

typedef Array<Vector> VectorMG;

/*-------------------------------------------------*/
class Operator
{
protected:
  armairvec _nall;
  armaimat  _n;
  Array<VectorMG> _mgmem;
  Vector _u, _f;
  mutable Array<UpdaterInterface*> _mgupdate;

  void vectormg2vector(int l, Vector& u, const VectorMG& umg) const;
  void vector2vectormg(int l, VectorMG& umg, const Vector& u) const;
  void restrict(int l, VectorMG& out, const VectorMG& in) const;
  void prolongate(int l, VectorMG& out, const VectorMG& in) const;
  void residual(int l, VectorMG& r, const VectorMG& u, const VectorMG& f) const;

  void jacobi       (int l, Vector& out, const Vector& in) const;
  void gauss_seidel1(int l, Vector& out, const Vector& in) const;
  void gauss_seidel2(int l, Vector& out, const Vector& in) const;
  void smooth       (int l, Vector& out, const Vector& in) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, VectorMG& Aw, double tol);
  void update(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, VectorMG& Aw, bool print=false) const;
  int solve(bool print=true);
  void _boundary(int l, Vector& out) const;

public:
  std::string smoother;
  int maxiter, optmem;
  double tol_rel, tol_abs;

  void set_parameters();
  ~Operator();
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
  const Vector& get_solution() const {return _u;}
  Vector& get_solution() {return _u;}

  void set_size(VectorMG& v) const;

  void dot(int l, Vector& out, const Vector& in, double d) const;
  int testsolve(bool print=true);
  void boundary(Vector& v, const Vector& w, double d=1.0) const;
  void boundary(Vector& v) const;
  void right(Vector& v) const;
};

#endif
