//
//  updater.c
//  Fada
//
//  Created by Roland Becker on 07/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "updater.hpp"
#include  "mgsolver.hpp"

/*--------------------------------------------------------------------------*/
void UpdaterConstant::setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat, int nvectors, const std::string& type, const std::string& solutiontype)
{
  _mat = mat;
  _mem = model.newVector(grid);
}

/*--------------------------------------------------------------------------*/
void UpdaterConstant::addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r, bool print)
{
  _mem->fill(0.0);
  _mat->dot(_mem, w, 1.0);
  u->add(_omega, *w);
  r->add(-_omega, *_mem);
}

/*--------------------------------------------------------------------------*/
void UpdaterSimple::addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r, bool print)
{
  _mem->fill(0.0);
  _mat->dot(_mem, w, 1.0);
  double d1 = w->dot(*r);
  double d2 = w->dot(*_mem);
  double omega = d1/d2;
  omega = fmax(fmin(omega, 10.0),0.1);
//  omega = fmax(fmin(omega, 2.0),0.2);
//  double omega = 1.0;
// std::cerr << "d1=" << d1 << " d2=" << d2 << " omega="<<omega<<"\n";
// omega = 1;
  u->add(omega, *w);
  r->add(-omega, *_mem);
}

/*--------------------------------------------------------------------------*/
Updater::~Updater(){}
Updater::Updater() : UpdaterInterface(), _mat(nullptr) {}
Updater::Updater(const Updater& updater) : UpdaterInterface()
{
  assert(0);
}
Updater& Updater::operator=(const Updater& updater)
{
  assert(0);
  return *this;
}
void Updater::setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat, int nvectors, const std::string& type, const std::string& solutiontype)
{
  _scale = false;
  _mat = mat;
  _nvectors = nvectors;
  assert(nvectors>0);
  _nextupdate = _nextproduct = _nmemory = _nextmemory = _niterafterrestar = 0;
  _type = type;
  _solutiontype = solutiontype;
  assert(_type == "cyc" or _type == "coef" or _type == "ortho" or _type == "restart");
  assert(_solutiontype == "ls" or _solutiontype == "gal" or _solutiontype == "gals");
  // std::cerr << "_type = " << _type << " _solutiontype " << _solutiontype << "\n";
  _H.set_size(_nvectors, _nvectors);
  _H.eye();
  _b.set_size(_nvectors);
  _b.fill(0);
  _x.set_size(_nvectors);
  _conditionmean = _conditionmax = 0.0;
  _mem.resize(2*_nvectors);
  for(int i=0;i<_mem.size();i++)
  {
    _mem[i] = model.newVector(grid);
  }
}

/*--------------------------------------------------------------------------*/
std::shared_ptr<VectorInterface> Updater::getV(int i)
{
  return _mem[2*i];
}
std::shared_ptr<VectorInterface> Updater::getAV(int i)
{
  return _mem[2*i+1];
}

/*--------------------------------------------------------------------------*/
void Updater::restart()
{
  _nmemory = 0;
  _nextmemory = 0;
}

/*--------------------------------------------------------------------------*/
void Updater::addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> res, bool print)
{
  // les impairs sont les produit avec la matrice des pairs !
//  std::cerr << _level << " addUpdate() _nextmemory="<< _nextmemory << " _nmemory="<< _nmemory<< " _nvectors="<< _nvectors << std::endl;
  _rnorm = res->norm();

  // getV(_nextmemory) = w;
  getV(_nextmemory)->equal(*w);
  if(_type == "ortho")
  {
    for(int i = 0; i < _nmemory; i++)
    {
      if(i == _nextmemory)
      {
        continue;
      }
      double scal = getV(i)->dot(*getV(_nextmemory));
      getV(_nextmemory)->add(-scal, *getV(i));
    }
  }
  double wnorm = w->norm();
  // cela ne change pas la solution, ni le résidu
  getV(_nextmemory)->scale(1.0/wnorm);
  getAV(_nextmemory)->fill(0.0);
//  _op->dot(_level, getAV(_nextmemory), getV(_nextmemory), 1.0);
  _mat->dot(getAV(_nextmemory), getV(_nextmemory), 1.0);

  int nmemory = _nmemory+1;
  _computeSmallSystem(_nextmemory, nmemory, res);
  _x = arma::solve( _H, _b, arma::solve_opts::fast);
//  _x = arma::pinv(_H)*_b;
  if(print) printf("'%4.2f' %10.3e %10.3e %10.3e\n", _x[0]/wnorm, w->norm(), getAV(_nextmemory)->norm()*wnorm, _rnorm);
  for(int i = 0; i < nmemory; i++)
  {
    if(_scale)
    {
      u->add(_rnorm*_x[i], *getV(i));
      res->add(-_rnorm*_x[i], *getAV(i));
    }
    else
    {
      u->add(_x[i], *getV(i));
      res->add(-_x[i], *getAV(i));
    }
  }

  if(_nmemory < _nvectors-1)
  {
    _nextmemory++;
    _nmemory++;
  }
  else
    // memoire pleine !
  {
    if( ( _type == "cyc" )or ( _type == "ortho") )
    {
      if(_nextmemory == _nvectors-1)
      {
        _nextmemory = 0;
      }
      else
      {
        _nextmemory++;
      }
    }
    else if(_type == "coef")
    {
      arma::uword ind;
      armavec s = arma::abs(_x);
      s.min(ind);
      _nextmemory = ind;
    }
    else if(_type == "restart")
    {
      getV(0) = getV(_nextmemory);
      getAV(0) =  getAV(_nextmemory);
      _nmemory = 0;
      _nextmemory = 0;
    }
    else
    {
      assert(0);
    }
  }
}

/*--------------------------------------------------------------------------*/
void Updater::_computeSmallSystem(int index, int nmemory, std::shared_ptr<VectorInterface const> r)
{
  _H.resize(nmemory, nmemory);
  _b.resize(nmemory);
  _x.resize(nmemory);
//  Vector& r = _mem(2*_nvectors);
//  const Vector& r = getV(_nvectors);

//  std::cerr << "index=" << index << " rnorm = " << r.norm() << "\n";

  //-------------------------------------------------
  if(_solutiontype == "gal")
  //-------------------------------------------------
  {
    for(int i = 0; i < nmemory; i++)
    {
      _H(index, i) = getAV(i)->dot(*getV(index));
      if(index != i)
      {
        _H(i, index) = getAV(index)->dot(*getV(i));
      }
      _b[i] = r->dot(*getV(i));
    }
  }
  //-------------------------------------------------
  else if(_solutiontype == "ls")
  //-------------------------------------------------
  {
    for(int i = 0; i < nmemory; i++)
    {
      double d = getAV(i)->dot(*getAV(index));
      _H(index, i) = d;
      if(index != i)
      {
        _H(i, index) = d;
      }
      _b[i] = r->dot(*getAV(i));
    }
  }
  //-------------------------------------------------
  else if(_solutiontype == "gals")
  //-------------------------------------------------
  {
    double alpha = 0.01;
    for(int i = 0; i < nmemory; i++)
    {
      //ls
      double d = getAV(i)->dot(*getAV(index));
      _H(index, i) = alpha*d;
      if(index != i)
      {
        _H(i, index) = alpha*d;
      }
      //gal
      _H(index, i) += getAV(i)->dot(*getV(index));
      if(index != i)
      {
        _H(i, index) += getAV(index)->dot(*getV(i));
      }
      //ls
      _b[i] = alpha* r->dot(*getAV(i));
      //gal
      _b[i] += r->dot(*getV(i));
    }
  }
  else
  {
    assert(0);
  }
  // std::cerr << "nmemory="<<nmemory<<" index="<<index<<"\n";
  // std::cerr << "_H="<<_H<<" _b="<<_b<<"\n";
  //  int nrang = arma::rank(_H);
  //  if(nrang != _H.n_rows)
  //  {
  //    std::cerr << getVisitor()->getClassName() << " rang(_H)="<<nrang<<" n="<<_H.n_rows<<" diff="<<_H.n_rows-nrang<<"\n";
  //  }
}
