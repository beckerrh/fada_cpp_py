//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

// #include  <omp.h>
#include  <sstream>
#include  <iomanip>
#include  "solverstokes.hpp"
#include  "sparsematrix.hpp"
#include  "boundary_conditions.hpp"
#include  "uniformgrid.hpp"
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

// /*-------------------------------------------------*/
// bool SolverStokes::_directoryExists( std::string pzPath )const
// {
//     const fs::path p{pzPath.c_str()};
//     return fs::exists(p);
// }
/*-------------------------------------------------*/
void SolverStokes::save_for_visu() const
{
    const fs::path p{_datadir.c_str()};
    if(not fs::exists(p))
    {
        std::string command = "mkdir -p " + _datadir;
        system( command.c_str() );
    }    
    std::string filename(_datadir + "/solution.hdf");
    PointDataMap u_intp = _model_all->to_point_data(_mgsolver_all->getU(), _mggrid->get(0));
    armaicvec n(_dim+1);
    // auto sv = std::dynamic_pointer_cast<StokesVector>(_mgsolver_all->getU());
    // assert(sv);
    // n[0] = sv->get_p()->dim();
    // for(int i=0;i<_dim;i++)
    // {
    //     n[1+i] = sv->get_v(i)->dim();
    // }
    // arma::hdf5_name spec(filename, "n");
    // n.save(spec);
    int count=0;        
    for(auto p:u_intp)
    {
        // std::cerr << "### " << p.first << "\n";
        // p.second->save(std::cerr, arma::arma_ascii);
        n[count] = p.second->n_elem;
        count++;        
    }
    arma::hdf5_name spec(filename, "n");
    n.save(spec);
    for(auto p:u_intp)
    {
        arma::hdf5_name spec(filename, p.first, arma::hdf5_opts::append);
        p.second->save(spec);        
    }
    // std::stringstream ss;
    // ss << _datadir << "solution";
    // std::string filename(ss.str() + "_p.hdf");
    // arma::hdf5_name spec(filename);
    // PointDataMap p_intp = _model_p->to_point_data(_mgsolver_p->getU(), _mggrid->get(0));
    // p_intp.save(spec);
    // for(int i=0;i<_dim;i++)
    // {
    //     std::stringstream ss2;
    //     ss2  << "_v" << i;
    //     std::string filename(ss.str() + ss2.str());
    //     filename += ".hdf";
    //     arma::hdf5_name spec(filename);
    //     armavec v_intp = _model_p->to_point_data(_mgsolver_v[i]->getU(), _mggrid->get(0));
    //     v_intp.save(spec);
    // }
    _mggrid->get(0)->savehdf5(_datadir + "/grid.hdf");
    // auto ugrid = std::make_shared<UniformGrid>();
    // ugrid->loadhdf5(_datadir + "/grid.hdf");
    // std::cerr << "ugrid="<<ugrid->toString();
}

/*-------------------------------------------------*/
void SolverStokes::save(SolverStokes::PressureVector p, SolverStokes::VelocityVector v, const armaicvec& iters) const
{
    std::stringstream ss; 
    ss << _datadir << "/solution";
    for(auto p: iters)
    {
        ss << "_" << std::setfill('0') << std::setw(6) << p;                
    }
    std::string filename(ss.str() + "_p.hdf");
    arma::hdf5_name spec(filename);
    // armavec p_intp = _model_p->to_point_data(_mgsolver_p->getU(), _mggrid->get(0));
    // armavec p_intp = _model_p->to_point_data(*p, _mggrid->get(0));
    // p_intp.save(spec);
    p->save(spec);
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss2;
        ss2  << "_v" << i;
        std::string filename(ss.str() + ss2.str());
        filename += ".hdf";
        arma::hdf5_name spec(filename);
        // armavec v_intp = _model_p->to_point_data(_mgsolver_v[i]->getU(), _mggrid->get(0));
        // armavec v_intp = _model_p->to_point_data(*v[i]->getU(), _mggrid->get(0));
        // v_intp.save(spec);
        v[i]->save(spec);
    }
    // _mggrid->get(0)->savehdf5("grid.hdf");
}
/*-------------------------------------------------*/
void SolverStokes::load(std::vector<SolverStokes::VelocityVector>& vs, int& nend) const
{
    const fs::path p{_datadir+"/status"};
    if(not fs::exists(p))
    {
        nend=0;
        for(auto pv: vs)
        {
            for(int i=0;i<_dim;i++)
            {
                pv[i]->get_arma().fill(0.0);
            }            
        }        
        return;
    }
    std::ifstream file(_datadir+"/status");
    double end_time, q, dt;
    int niter_out, nt, ngrid;
    file >> end_time;
    file >> q;
    file >> dt;
    file >> niter_out;
    file >> nend;
    file >> nt;
    file >> ngrid;
    assert(nt==vs.size());
    assert(ngrid==vs[0][0]->dim());
    if(_end_time!=end_time or _q!=q or _dt!=dt or _niter_out!=niter_out or _nt!=nt)
    {
        _not_written_("wrong data");
    }
    file.close();      
    std::stringstream ss; 
    ss << _datadir << "/solution";
    ss << "_" << std::setfill('0') << std::setw(6) << nend;                
    ss << "_" << std::setfill('0') << std::setw(6) << 0;                
    for(int it=0;it<nt;it++)
    {
        for(int i=0;i<_dim;i++)
        {
            std::stringstream ss2;
            ss2 << "_" << std::setfill('0') << std::setw(6) << it;                
            ss2  << "_v" << i;
            std::string filename(ss.str() + ss2.str());
            filename += ".hdf";
            arma::hdf5_name spec(filename);
            vs[it][i]->get_arma().load(spec);
        }
    }
    nend++;
}


/*-------------------------------------------------*/
std::string SolverStokes::toString() const
{
    std::stringstream ss;
    ss << "model_p=" << _model_p->toString();
    ss << "mggrid=" << _mggrid->toString();
    return(ss.str());
}

/*-------------------------------------------------*/
void SolverStokes::make_coupled_solver()
{
    bool mgtimer=true, mgdebug=false;
    for(std::map<std::string,std::string>::const_iterator p = _parameters.begin(); p != _parameters.end(); p++)
    {
        if(p->first=="mgtimer")
        {
            mgtimer = p->second=="true";
        }
    }
    
    // _model_all = std::make_shared<Model<ModelStokes2d, StokesVector>>(_parameters, _boundaryconditions);
    _model_all = std::make_shared<ModelStokes2d>(_parameters, _application);    
    _mgsolver_all = std::make_shared<MgSolver>(mgtimer, mgdebug);    
    _mgsolver_all->set_sizes(_mggrid, _model_all);    
}

/*-------------------------------------------------*/
void SolverStokes::make_chorin_solver()
{
    bool mgtimer=true, mgdebug=false;
    for(std::map<std::string,std::string>::const_iterator p = _parameters.begin(); p != _parameters.end(); p++)
    {
        if(p->first=="mgtimer")
        {
            mgtimer = p->second=="true";
        }
    }
    
    _model_p = std::make_shared<Model<Q02d, Vector<GridVector>>>("p",_parameters, _application);    
    _mgsolver_p = std::make_shared<MgSolver>(mgtimer, mgdebug);
    _mgsolver_p->set_sizes(_mggrid, _model_p);    
    _model_v.resize(_dim);
    _mgsolver_v.resize(_dim);
    // for(int i=0;i<_dim;i++)
    // {
    //     _mgsolver_v[i] = std::make_shared<MgSolver>(mgtimer, mgdebug);
    //     // _mgsolver_v[i]->set_sizes(_mggrid, _model_p);
    // }
    _B.resize(_dim);
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "v" << i;        
        std::map<std::string,std::string> parameters2(_parameters.begin(), _parameters.end());
        parameters2["direction"] =  std::to_string(i);
        _model_v[i] =  std::make_shared<Model<Q1shifted2d, Vector<GridVector>>>(ss.str(),parameters2, _application); 
        _mgsolver_v[i] = std::make_shared<MgSolver>(mgtimer, mgdebug);
        _mgsolver_v[i]->set_sizes(_mggrid, _model_v[i]);
        auto b = getModelV(i).newMatrixDivergence(_mggrid->get(0));
        _B[i] = std::dynamic_pointer_cast<SparseMatrix const>(b);
        assert(_B[i]);
    }
}

/*-------------------------------------------------*/
void SolverStokes::set_data(const std::map<std::string,std::string>& parameters)
{
    int nlevels(8);
    _applicationname = "Sinus_per";
    _datadir = "";
    _dim=2;
    _dt=0.0;
    int n00 = 4;
    for(std::map<std::string,std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
    {
        if(p->first=="nlevels")
        {
            nlevels = std::stoi(p->second);
        }
        else if(p->first=="dim")
        {
            _dim = std::stoi(p->second);
        }
        else if(p->first=="application")
        {
            _applicationname = p->second.c_str();
        }
        else if(p->first=="datadir")
        {
            _datadir = p->second.c_str();
        }
        else if(p->first=="n0")
        {
            n00 = std::stoi(p->second);
        }
        else if(p->first=="dt")
        {
            _dt = std::stod(p->second);
        }
    }
    _parameters["dt"] = std::to_string(_dt);
    if(_datadir=="")
    {
        _datadir = "datadir_" +_applicationname;
    }
    try{
        _end_time = std::stod(_parameters.at("end_time"));
        _nt = std::stoi(_parameters.at("nt"));
        _niter_out = std::stoi(_parameters.at("niter_out"));
        _nsamples = std::stoi(_parameters.at("nsamples"));
        _q = std::stod(_parameters.at("q"));        
    }
    catch(...){
        std::cerr << ("no sde!");
        _end_time = 1.0; _nt=0; _niter_out=0; _nsamples=0; _q=0.0;
    }         
    _application = std::make_shared<StokesApplication>(_applicationname, _dim);
    armaicvec n0 = {n00, n00};
    int ref_factor = 2;
    _mggrid = std::make_shared<UniformMultiGrid>(nlevels, n0, nullptr, ref_factor);
}
/*-------------------------------------------------*/
SolverStokes::VelocityVector SolverStokes::newVelocityVector()
{
    SolverStokes::VelocityVector v(_dim);
    for(int i=0;i<_dim;i++)
    {
        v[i] = std::dynamic_pointer_cast <GridVector>(_model_v[i]->newVector(_mggrid->get(0)));
        assert(v[i]);
    }
    return v;
}
/*-------------------------------------------------*/
SolverStokes::PressureVector SolverStokes::newPressureVector()
{
    return std::dynamic_pointer_cast <GridVector>(_model_p->newVector(_mggrid->get(0)));   
}
/*-------------------------------------------------*/
StokesInfo SolverStokes::chorin_stationary(bool print)
{
    make_chorin_solver();    
    double dt = std::stod(_parameters.at("dt"));
    StokesInfo info;
    info.niter_mean_v=0.0;
    info.niter_mean_p=0.0;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(_mggrid->get(0));
    double h = ug->dx(0);
    SolverStokes::VelocityVector v = newVelocityVector(); 
    SolverStokes::VelocityVector vold = newVelocityVector(); 
    SolverStokes::VelocityVector fv = newVelocityVector(); 
    auto pold = std::dynamic_pointer_cast <GridVector>(_model_p->newVector(_mggrid->get(0)));   
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "v" << i;
        getModelV(i).rhs(*fv[i], _mggrid->get(0), _application->rhs(ss.str()));
        v[i]->fill(0.0);
    }
    int niter_out = 100;
    int niter_in = 100;
    bool print_mg=false;
    for(int iter_out=0; iter_out<niter_out; iter_out++)
    {
        for(int iter_in=0; iter_in<niter_in; iter_in++)
        {
            get_rhs_p().fill(0.0);
            for(int i=0;i<_dim;i++)
            {
                vold[i]->equal(*v[i]);
                get_rhs_v(i) = *(fv[i]);
                get_rhs_v(i).add(h*h/dt,*(v[i]));        
                int mgiter = _mgsolver_v[i]->solve(print_mg);
                info.niter_mean_v += mgiter;
                assert(mgiter<12);
                _B[i]->dot(get_rhs_p(), get_solution_v(i), -1.0/dt);
                v[i]->equal(get_solution_v(i));
            }
            pold->equal(get_solution_p());
            int mgiter = _mgsolver_p->solve(print_mg);        
            info.niter_mean_p += mgiter;        
            assert(mgiter<12);
            for(int i=0;i<_dim;i++)
            {
                _B[i]->Tdot(*(v[i]), get_solution_p(), dt/h/h);
            }
            pold->add(-1, get_solution_p());
            double diff_p = arma::norm(pold->get_arma());
            double diff_v = 0.0;
            for(int i=0;i<_dim;i++)
            {
                vold[i]->add(-1, *v[i]);
                diff_v += arma::norm(vold[i]->get_arma());
            }
            int niter = iter_in + iter_out*niter_in;
            if(print) printf("Iter: %6d Diff: %12.4e %12.4e\n", niter, diff_p, diff_v);
            if(diff_p+diff_v<1e-6*dt)
            {
                info.niter = niter;
                info.niter_mean_p /= (double) niter;
                info.niter_mean_v /= (double) _dim*niter;
                auto err = compute_error();
                info.err_p = err[0];
                info.err_v = 0.0;
                for(int i=0;i<_dim;i++)
                {
                    std::stringstream ss;
                    ss << "v" << i;
                    info.err_v += err.at(ss.str());
                }
                return info;
            }
        }
    }
    info.niter = -1;
    return info;
}
/*-------------------------------------------------*/
StokesInfoSde SolverStokes::chorin_sde(bool print)
{
    Timer T(true, false, std::set<std::string>({"mg", "update_mg", "error"}));
    make_chorin_solver();    
    std::stringstream ss;
    ss << "_" << _end_time << "_" << _q << "_" << _dt << "_" << _nt << "_" << _niter_out;
    _datadir += ss.str();
    const fs::path p{_datadir.c_str()};
    if(not fs::exists(p))
    {
        std::string command = "mkdir -p " + _datadir;
        system( command.c_str() );
    }    
    
    StokesInfoSde info(_nt-1);
    info.niter_mean_v=0.0;
    info.niter_mean_p=0.0;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(_mggrid->get(0));
    double h = ug->dx(0);
    std::cerr << "nall="<<ug->n_gridpoints()<<"\n";
    SolverStokes::VelocityVector v = newVelocityVector(); 
    SolverStokes::VelocityVector fv = newVelocityVector(); 
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "v" << i;
        getModelV(i).rhs(*fv[i], _mggrid->get(0), _application->rhs(ss.str()));
        v[i]->fill(0.0);
    }
    // int niter_out = 3;
    
    
    int niter_in = (int) ( _end_time/(double) _niter_out/_dt);
    std::cerr<< "niter_in=" <<niter_in<<"\n";
    std::cerr<< "_dt=" <<_dt<<"\n";
    std::cerr<< "_end_time/(double) _niter_out/_dt=" <<_end_time/(double) _niter_out/_dt<<"\n";
    
    bool print_mg=false;
    
    std::cerr<< "niter_out=" <<_niter_out<<"\n";
    std::cerr<< "q=" <<_q<<"\n";
    armaicvec n_inners(_nt);
    for(int i=0; i<_nt; i++) n_inners[i] = niter_in*pow(2,i);
    std::cerr << "n_inners" << n_inners.t();
    armavec dts(_nt);
    for(int i=0; i<_nt; i++) dts[i] = _dt/pow(2,i);
    for(int i=0; i<_nt-1; i++) info.dt[i] = dts[i];
    std::cerr << "dts" << dts;

    std::vector<SolverStokes::VelocityVector> v_out(_nt);
    std::vector<SolverStokes::PressureVector> p_out(_nt);
    SolverStokes::VelocityVector v_diff = newVelocityVector();
    SolverStokes::PressureVector p_diff = newPressureVector();
    for(int i=0; i<_nt; i++)
    {
        v_out[i] = newVelocityVector(); 
        p_out[i] = newPressureVector(); 
    }
    int isample_first=0;
    load(v_out, isample_first);
    std::cerr << "***************** isample_first = " << isample_first << "\n";
    
    info.err_v.fill(0.0);
    info.err_p.fill(0.0);
    arma::arma_rng::set_seed_random();
    
//     omp_set_dynamic(0);     // Explicitly disable dynamic teams
//     omp_set_num_threads(11); // Use 4 threads for all consecutive parallel regions
// #pragma omp parallel
//     printf("Hello from process: %d %d\n", omp_get_thread_num(), omp_get_num_threads());
// #pragma omp parallel for
    for(int isample=0;isample<_nsamples;isample++)
    {
        armavec white_noise = arma::randn(_niter_out*n_inners[_nt-1]);
        for(int iter_out=0; iter_out<_niter_out; iter_out++)
        {
            for(int it=0; it<_nt; it++)
            {
                double dt = dts[it];
                T.start("update_mg");        
                for(int i=0;i<_dim;i++)
                {
                    _mgsolver_v[i]->update_coefficients(dt);
                }
                T.stop("update_mg");        
                for(int i=0;i<_dim;i++)
                {
                    v[i]->equal(*v_out[it][i]);
                }
                for(int iter_in=0; iter_in<n_inners[it]; iter_in++)
                {
                    int index_wn = iter_out*n_inners[_nt-1] + pow(2,_nt-1-it)*iter_in;
                    // std::cerr<< "iter_out=" <<iter_out<<" it="<<it << " iter_in="<< iter_in << " index_wn="<<index_wn<<"\n";
                    double fscale = _q * white_noise[index_wn]/sqrt(dt);
                    get_rhs_p().fill(0.0);
                    for(int i=0;i<_dim;i++)
                    {
                        get_rhs_v(i) = fscale * fv[i]->get_arma();
                        // get_rhs_v(i) = fv[i]->get_arma();
                        get_rhs_v(i).add(h*h/dt,*(v[i]));
                        T.start("mg");
                        get_rhs_v(i).project_mean();        
                        int mgiter = _mgsolver_v[i]->solve(print_mg);
                        get_solution_v(i).project_mean();        
                        T.stop("mg");        
                        info.niter_mean_v += mgiter;
                        assert(mgiter<12);
                        _B[i]->dot(get_rhs_p(), get_solution_v(i), -1.0/dt);
                        v[i]->equal(get_solution_v(i));
                    }
                    T.start("mg");        
                    int mgiter = _mgsolver_p->solve(print_mg);
                    T.stop("mg");        
                    info.niter_mean_p += mgiter;        
                    assert(mgiter<12);
                    for(int i=0;i<_dim;i++)
                    {
                        _B[i]->Tdot(*(v[i]), get_solution_p(), dt/h/h);
                    }
                }
                p_out[it]->equal(get_solution_p());
                for(int i=0;i<_dim;i++)
                {
                    v_out[it][i]->equal(*v[i]);
                }
            }
            T.start("error");
            // iters[0] =  isample;
            // iters[1] =  iter_out;
            // iters[2] =  it;

            
            std::ofstream file(_datadir+"/status");
            file << _end_time << " " << _q << " " << _dt << " " << _niter_out << "\n";
            file << isample_first + isample << " " << _nt << " " << v[0]->get_arma().n_elem<<"\n";
            file.close();      
            for(int it=0; it<_nt; it++)
            {
                armaicvec iters({isample_first + isample, iter_out, it});
                save(p_out[it], v_out[it], iters);
            }
            for(int it=0; it<_nt-1; it++)
            {
                p_diff->equal(*p_out[it]);
                p_diff->add(-1, *p_out[_nt-1]);
                info.err_p[it] += h*arma::norm(p_diff->get_arma())/(double) _nsamples*_niter_out;
                // diff_p[it] += arma::dot(p_diff->get_arma(),p_diff->get_arma());
                // diff_p[it] = sqrt(diff_p[it])*h;
                for(int i=0;i<_dim;i++)
                {
                    v_diff[i]->equal(*v_out[it][i]);
                    v_diff[i]->add(-1, *v_out[_nt-1][i]);
                    info.err_v[it] += h*arma::norm(v_diff[i]->get_arma())/(double) _nsamples*_niter_out;
                    // diff_v[it] += arma::dot(v_diff[i]->get_arma(),v_diff[i]->get_arma());
                    // diff_v[it] = sqrt(diff_v[it])*h;
                }
            }
            T.stop("error");        
            if(print) std::cerr << "iter_out="<<iter_out << "\ninfo.err_p="<< info.err_p.t() << "info.err_v="<< info.err_v.t(); 
        }        
    }
    info.niter = -1;
    return info;
}

/*-------------------------------------------------*/
StokesInfo SolverStokes::solve_stationary(bool print)
{
    make_coupled_solver();    
    StokesInfo info;
    // _u.set_size(_mggrid->get(0)->n());
    // _u.fill(0);
    // _f.set_size(_u);
    _model_all->rhs(_mgsolver_all->getF(), _mggrid->get(0), _application);
    auto sv = std::dynamic_pointer_cast<StokesVector>(_mgsolver_all->getF());
    assert(sv);
    sv->get_p()->project_mean();
    for(int i=0;i<_dim;i++)
    {
        sv->get_v(i)->project_mean();
    }
    //     _model_v[i].Q1shiftedget_rhs_v(i), _mggrid->get(0), _application->rhs_v(i));
    //     getQ0().Q1shiftedget_rhs_p(), _mggrid->get(0), _application->rhs_p());
    // get_rhs_p().project_mean();
    // // get_Q1shifted).boundary_zero();
    // for(int i=0;i<_dim;i++)
    // {
    //     _model_v[i].Q1shiftedget_rhs_v(i), _mggrid->get(0), _application->rhs_v(i));
    //     get_rhs_v(i).project_mean();
    // }
    // std::cerr << "f_p " << p->get_p()->min() << " " << p->get_p()->max() << "\n";
    info.niter = _mgsolver_all->solve(print);
    sv = std::dynamic_pointer_cast<StokesVector>(_mgsolver_all->getU());
    assert(sv);
    sv->get_p()->project_mean();
    for(int i=0;i<_dim;i++)
    {
        sv->get_v(i)->project_mean();
    }
    // std::cerr << "p " << p->get_p()->min() << " " << p->get_p()->max() << "\n";
    auto err = compute_error();
    // std::cerr << "***" << err.t();
    info.err_p = err.at("p");
    info.err_v = err.at("v0");
    for(int i=1;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "v" << i;        
        info.err_v += err.at(ss.str());
    }
    //
    // std::cerr << "f_p " << get_rhs_p().min() << " " << get_rhs_p().max() << "\n";
    // int iter = _mgsolver_p->solve(print);
    // std::cerr << "p " << get_solution_p().min() << " " << get_solution_p().max() << "\n";
    // for(int i=0;i<_dim;i++)
    // {
    //     std::cerr << "f_v " << get_rhs_v(i).min() << " " << get_rhs_v(i).max() << "\n";
    //     int iter = _mgsolver_v[i]->solve(print);
    //     std::cerr << "p " << get_solution_v(i).min() << " " << get_solution_v(i).max() << "\n";
    // }
    return info;
}
/*-------------------------------------------------*/
std::map<std::string,double> SolverStokes::compute_error() const
{
    if(_mgsolver_all)
    {
        return _model_all->compute_error(_mgsolver_all->getU(), _mggrid->get(0), _application);
    }
    std::map<std::string,double> err;
    err["p"] = getQ0().compute_error(get_solution_p(), _mggrid->get(0), _application->solution("p")).at("p");
    for(int i=0;i<_dim;i++)
    {
        std::stringstream ss;
        ss << "v" << i;
        err[ss.str()] = getModelV(i).compute_error(get_solution_v(i), _mggrid->get(0), _application->solution(ss.str())).at(ss.str());
    }
    return err;
}
