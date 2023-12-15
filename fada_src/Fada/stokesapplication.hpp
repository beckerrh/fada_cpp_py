//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  stokesapplication_hpp
#define  stokesapplication_hpp

#include  <sstream>
#include  <string>
#include  "applicationinterface.hpp"


/*-------------------------------------------------*/
class Sinus2DSolutionP: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DSolutionP(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2DRhsP: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DRhsP(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2DSolutionV0: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DSolutionV0(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2DRhsV0: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DRhsV0(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2DSolutionV1: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DSolutionV1(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2DRhsV1: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2DRhsV1(double A_=1, double B_=1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};

/*-------------------------------------------------*/
class SinusBis2DSolutionP: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DSolutionP(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
class SinusBis2DRhsP: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DRhsP(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
class SinusBis2DSolutionV0: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DSolutionV0(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
class SinusBis2DRhsV0: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DRhsV0(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
class SinusBis2DSolutionV1: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DSolutionV1(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
class SinusBis2DRhsV1: public AnalyticalFunctionInterface
{
    double A, B, C, D;
public:
    SinusBis2DRhsV1(double A_=1, double B_=1, double C_=1, double D_=1) : AnalyticalFunctionInterface(), A(A_), B(B_), C(C_), D(D_) {}
    double operator()(double x, double y) const;
};
/*-------------------------------------------------*/
class Linear2DSolutionP: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Linear2DRhsP: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Linear2DSolutionV0: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Linear2DRhsV0: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Linear2DSolutionV1: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Linear2DRhsV1: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};

/*-------------------------------------------------*/
class StokesApplication : public ApplicationInterface
{
protected:
    void set_application(int dim, std::string name);

public:
    StokesApplication(int dim, std::string name) : ApplicationInterface(name, dim, dim+1) 
    {
        set_application(dim, name);
    }
};


#endif
