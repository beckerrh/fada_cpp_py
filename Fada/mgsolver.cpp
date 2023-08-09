#include  <sstream>
#include  <math.h>
#include  <stdio.h>
#include  <armadillo>
#include  <algorithm>
#include  "mgsolver.hpp"
#include  "updater.hpp"
#include  "uniformgrid.hpp"

/*-------------------------------------------------*/
std::string MgSolver::toString() const
{
    std::stringstream ss;
    ss << "nlevels = " << _nlevels;
    return ss.str();
}
// /*-------------------------------------------------*/
// void MgSolver::set_parameters(int maxiter, double tol_rel, double tol_abs)
// {
//     _maxiter = maxiter;
//     _tol_rel = tol_rel;
//     _tol_abs = tol_abs;
// }
/*-------------------------------------------------*/
void MgSolver::set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<ModelInterface> model)
{
    _timer.enrol("smooth");
    _timer.enrol("transfer");
    _timer.enrol("residual");
    _timer.enrol("update");
    _timer.enrol("solvecoarse");
    _nlevels = mgrid->nlevels();
    // std::cerr << "mgrid" << mgrid->toString() << "\n";
    // std::cerr << "_nlevels=" << _nlevels << "\n";
    _mgrid = mgrid;
    _model = model;
    _mgmem.resize(4);
    for(int i=0;i<_mgmem.size();i++)
    {
        _set_size_vectormg(mgrid, _mgmem[i]);
    }
    _mgmatrix.resize(_nlevels);
    for(int l=0;l<_nlevels;l++)
    {
        _mgmatrix[l] = _model->newMatrix(mgrid->get(l));
    }
    _mgsmoother.resize(_nlevels-1);
    for(int l=0;l<_nlevels-1;l++)
    {
        _mgsmoother[l] = _model->newSmoother(mgrid->get(l), _mgmatrix[l]);
    }
    _mgcoarsesolver = _model->newCoarseSolver("direct", mgrid->get(_nlevels-1), _mgmatrix[_nlevels-1]);
    _mgtransfer.resize(_nlevels-1);
    for(int l=0;l<_nlevels-1;l++)
    {
        _mgtransfer[l] = _model->newTransfer(mgrid->get(l+1), mgrid->get_ref_factor());
    }
    _mgupdate.resize(_nlevels);
    _mgupdatesmooth.resize(_nlevels);
    for(int l=0;l<_nlevels;l++)
    {
        // if(updatemem==0)
        // {
        //     // _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterConstant(0.8));
        //     _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
        // }
        // else
        // {
        //     // _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new Updater);
        //     _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
        // }
        // _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterConstant(0.8));
        // _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
        // _mgupdatesmooth[l] = std::shared_ptr<UpdaterInterface>(new UpdaterConstant(0.8));
        // _mgupdatesmooth[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
        // _mgupdatesmooth[l] = std::shared_ptr<UpdaterInterface>(new Updater);
        std::stringstream ss;
        ss << "mg_" << l;
        std::stringstream ss2;
        ss2 << "sm_" << l;
        
        _mgupdate[l] = std::make_shared<Updater>(ss.str(), "cyc", "gal", 2);
        // _mgupdate[l] = std::make_shared<UpdaterSimple>(ss.str(), "gal");
        // _mgupdate[l] = std::make_shared<UpdaterConstant>(ss.str(), 1);

        _mgupdatesmooth[l] = std::make_shared<Updater>(ss2.str(), "cyc", "ls", 2);
        // _mgupdatesmooth[l] = std::make_shared<UpdaterSimple>(ss2.str(), "ls");
        // _mgupdatesmooth[l] = std::make_shared<UpdaterConstant>(ss2.str(), 1);
    }
    for(int l=0;l<_nlevels;l++)
    {
        _mgupdate[l]->setParameters(*_model, mgrid->get(l), _mgmatrix[l]);
        _mgupdatesmooth[l]->setParameters(*_model, mgrid->get(l), _mgmatrix[l]);
    }
    _niter_post = 1; 
    _niter_pre = 1;
}
/*-------------------------------------------------*/
void MgSolver::_set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const
{
    v.resize(_nlevels);
    for(int l=0;l<_nlevels;l++)
    {
        v[l] = _model->newVector(_mgrid->get(l));
    }
}
/*-------------------------------------------------*/
void MgSolver::update_coefficients(double dt)
{
    for(int l=0;l<_nlevels;l++)
    {
        _model->update_coefficients(_mgrid->get(l), _mgmatrix[l], dt);
    }
    for(int l=0;l<_nlevels-1;l++)
    {
        _mgsmoother[l]->update(_mgmatrix[l]);
    }
    _mgcoarsesolver->update(_mgmatrix[_nlevels-1]);
}


// /*-------------------------------------------------*/
// void MgSolver::residual(int l, std::shared_ptr<VectorInterface> r, std::shared_ptr<VectorInterface const> u, std::shared_ptr<VectorInterface const>  f) const
// {
//     // std::cerr << "\nu\n";
//     // u->save(std::cerr);
//
//     r->equal(f);
//     _mgmatrix[l]->dot(r, u, -1.0);
//     // r->boundary_zero();
//
//     // std::cerr << "\nr\n";
//     // r->save(std::cerr);
// }
/*-------------------------------------------------*/
int MgSolver::solve(bool print, MgSolver::IterationInfo info)
{
    VectorMG& umg = _mgmem[0];
    VectorMG& fmg = _mgmem[1];
    VectorMG& d   = _mgmem[2];
    VectorMG& w   = _mgmem[3];

    umg[0]->equal(fmg[0]);

    int iter, finest_level = 0;
    double res, resold, tol=0;
  
    // _maxiter=3;
  
  
    for(iter=0; iter<info.maxiter+1; iter++)
    {
        _timer.start("residual");
        // residual(finest_level, *d[finest_level], *umg[finest_level], *fmg[finest_level]);
        // residual(finest_level, d[finest_level], umg[finest_level], fmg[finest_level]);
        
        d[finest_level]->equal(fmg[finest_level]);
        _mgmatrix[finest_level]->dot(d[finest_level], umg[finest_level], -1.0);
        
        // boundary_zero// d[finest_level]->fill_bdry(0);
        _timer.stop("residual");
        resold = res;
        res = d[finest_level]->norm();
        if(iter==0)
        {
            resold = 2*res;
            tol = fmax(info.tol_abs, info.tol_rel*res);
            if(print) printf("-mg- ---tol= %10.3e ---\n", tol);
        }
        if(print) printf("-mg- %3d %10.3e\n", iter, res);
        if(res <= tol)
        {
            return iter;
        }
        else if(res > resold)
        {
            printf("-mg- *** no convergence  res = %10.3e resold = %10.3e ---\n", res, resold);
            exit(1);
        }
        mgstep(finest_level, umg, fmg, d, w, tol);
    }
    printf("-mg- *** no convergence  res = %10.3e iter = %4d ---\n", res, iter);
    exit(1);
    return -1;
}
/*-------------------------------------------------*/
void MgSolver::mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol)
{
    // std::cerr << "MgSolver::mgstep() l=" << l << "_nlevels="<<_nlevels<<"\n";
    if(l==_nlevels-1)
    {
        _timer.start("solvecoarse");
        // std::cerr << "\n BEFORE u\n";
        // u[l]->save(std::cerr);
        // std::cerr << "\n BEFORE f\n";
        // f[l]->save(std::cerr);
        // _mgsmoother[l]->solve(u[l], f[l]);

        // _mgcoarsesolver->solve(u[l], f[l]);
        // u[l]->after_smooth();
        _mgcoarsesolver->solve(w[l], d[l]);
        w[l]->after_smooth();
        _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);

        // std::cerr << "\n AFTER u\n";
        // u[l]->save(std::cerr);
        _timer.stop("solvecoarse");
    }
    else
    {
        for(int iter_sm=0;iter_sm<_niter_post;iter_sm++)
        {
            _timer.start("smooth");
            // _mgsmoother[l]->presmooth(w[l], d[l]);
            _mgsmoother[l]->postsmooth(w[l], d[l]);
            w[l]->after_smooth();
            _timer.stop("smooth");
            _timer.start("update");
            _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
            // _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
            _timer.stop("update");        
        }

        _timer.start("transfer");
        _mgtransfer[l]->restrict(d[l+1], d[l]);
        d[l+1]->after_restrict();
        _timer.stop("transfer");

        f[l+1]->equal(d[l+1]);
        u[l+1]->fill(0.0);
        mgstep(l+1, u, f, d, w, tol);

        _timer.start("transfer");
        _mgtransfer[l]->prolongate(w[l], u[l+1]);
        d[l+1]->after_prolongate();
        _timer.stop("transfer");
        // _timer.start("residual");
        // residual(l, *d[l], *u[l], *f[l]);
        // _timer.stop("residual");
        _timer.start("update");
        _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
        _timer.stop("update");
        for(int iter_sm=0;iter_sm<_niter_post;iter_sm++)
        {
            _timer.start("smooth");
            _mgsmoother[l]->postsmooth(w[l], d[l]);
            w[l]->after_smooth();
            _timer.stop("smooth");
            _timer.start("update");
            // _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
            _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
            _timer.stop("update");
        }
    }
}
