//
//  vectorinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef systemvector_h
#define systemvector_h

#include  "typedefs.hpp"
#include  "vectorinterface.hpp"

/*-------------------------------------------------*/
class SystemVector : public  virtual VectorInterface
{
protected:
    std::vector<std::shared_ptr<VectorInterface>> _vectors;
public:
    SystemVector(int n): _vectors(n, nullptr) {}
    SystemVector(const SystemVector& vector):_vectors(vector._vectors) {}
    std::shared_ptr<VectorInterface>& get(int i) {return _vectors[i];}
    const std::shared_ptr<VectorInterface>& get(int i) const {return _vectors[i];}
    int n_vectors() const {return _vectors.size();}
    void set_size(const armaicvec& n)
    {
        for(auto p:_vectors) p->set_size(n);
    }
    int get_size() const{
        int n(0);
        for(auto p:_vectors) 
        {
            assert(p);
            n+=p->get_size();
        }
        return n;
    }
    void after_prolongate()
    {
        for(auto p:_vectors) p->after_prolongate();
    }
    void after_restrict()
    {
        for(auto p:_vectors) p->after_restrict();        
    }
    void after_smooth()
    {
        for(auto p:_vectors) p->after_smooth();                
    }
    void after_residual()
    {
        for(auto p:_vectors) p->after_residual();                
    }
    void fill(double d=0)
    {
        for(auto p:_vectors) p->fill(d);                        
    }
    double dot(std::shared_ptr<VectorInterface const> v)const
    {
        double d(0.0);
        auto pv = std::dynamic_pointer_cast<SystemVector const>(v);
        for(int i=0;i<_vectors.size();i++) d += get(i)->dot(pv->get(i));        
        return d;                        
    }
    double norm(double p=2)const
    {
        double d(0.0);
        for(auto q:_vectors) d += q->norm(p);
        return d;                                
    }
    void equal(std::shared_ptr<VectorInterface const> v)
    {
        auto pv = std::dynamic_pointer_cast<SystemVector const>(v);
        for(int i=0;i<_vectors.size();i++) get(i)->equal(pv->get(i));        
    }
    void add(double d, std::shared_ptr<VectorInterface const> v)
    {
        auto pv = std::dynamic_pointer_cast<SystemVector const>(v);
        for(int i=0;i<_vectors.size();i++) get(i)->add(d, pv->get(i));        
    }
    void scale(double d)
    {
        for(auto p:_vectors) p->scale(d);        
    }
    void save(std::ostream& os, arma::file_type ftype=arma::arma_binary) const
    {
        for(auto p:_vectors) p->save(os, ftype);                
    }
    void savehdf5(const std::string& filename) const
    {
        for(auto p:_vectors) p->savehdf5(filename);                
    }
};


#endif
