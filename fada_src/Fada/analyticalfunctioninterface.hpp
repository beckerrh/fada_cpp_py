//
//  vectorinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef analyticalfunctioninterface_h
#define analyticalfunctioninterface_h

#include  "typedefs.hpp"


/*-------------------------------------------------*/
class AnalyticalFunctionInterface
{
public:
  virtual ~AnalyticalFunctionInterface() {}
  AnalyticalFunctionInterface() {}
  AnalyticalFunctionInterface(const AnalyticalFunctionInterface& vector) {}

  virtual double operator()(double x, double y)const {_not_written_("AnalyticalFunctionInterface::operator()"); return 0;}
  virtual double operator()(double x, double y, double z)const {_not_written_("AnalyticalFunctionInterface::operator()"); return 0;}
};

typedef std::map<std::string, std::shared_ptr<AnalyticalFunctionInterface>> FunctionMap;
static std::ostream& operator<<(std::ostream& os, const FunctionMap& v)
{
    for(auto p:v)
    {
        os << p.first << " : " << p.second << "\n";
    }
    return os;
}

/*-------------------------------------------------*/
class ConstantFunction: public AnalyticalFunctionInterface
{
protected:
    double _d;
public:
    ConstantFunction(double d=1): AnalyticalFunctionInterface(), _d(d) {}
    ConstantFunction(const ConstantFunction& fct): AnalyticalFunctionInterface(fct), _d(fct._d) {}
    double get_constant() const{return _d;};
};

/*-------------------------------------------------*/
class RandomFunction: public AnalyticalFunctionInterface
{
protected:
    std::string _type;
public:
    RandomFunction(std::string type="uniform"): AnalyticalFunctionInterface(), _type(type) {}
    RandomFunction(const RandomFunction& fct): AnalyticalFunctionInterface(fct), _type(fct._type) {}
    std::string get_type() const{return _type;}
};


//
// /*-------------------------------------------------*/
// class SystemAnalyticalFunction : public AnalyticalFunctionInterface
// {
// protected:
//     std::vector<std::shared_ptr<AnalyticalFunctionInterface>> _vector;
// public:
//   SystemAnalyticalFunction(int n) : AnalyticalFunctionInterface(), _vector(n,nullptr) {}
//   SystemAnalyticalFunction(const SystemAnalyticalFunction& vector) : AnalyticalFunctionInterface(), _vector(vector._vector) {}
//   std::shared_ptr<AnalyticalFunctionInterface>& get(int i) {return _vector[i];}
//   const std::shared_ptr<AnalyticalFunctionInterface>& get(int i) const {return _vector[i];}
//   std::shared_ptr<AnalyticalFunctionInterface>& back() {return _vector.back();}
//   const std::shared_ptr<AnalyticalFunctionInterface>& back() const {return _vector.back();}
// };

#endif
