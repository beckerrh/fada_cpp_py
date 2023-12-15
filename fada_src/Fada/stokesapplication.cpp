//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <sstream>
#include  "stokesapplication.hpp"
#include  "tokenize.hpp"


/*-------------------------------------------------*/
// #define  A 1
// #define  B 1
// #define  C 1
// #define  D 1
double SinusBis2DSolutionP::operator()(double x, double y) const {return cos(2*C*M_PI*x)*cos(2*D*M_PI*y);}
double SinusBis2DSolutionV0::operator()(double x, double y) const {return -2*B*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y);}
double SinusBis2DSolutionV1::operator()(double x, double y) const {return  2*A*M_PI*cos(A*2*M_PI*x)*sin(B*2*M_PI*y);}

// double SinusBis2DRhsP::operator()(double x, double y) const {return  (C*C+D*D)*4*M_PI*M_PI*cos(C*2*M_PI*x)*cos(D*2*M_PI*y);}
double SinusBis2DRhsP::operator()(double x, double y) const {return  0;}
double SinusBis2DRhsV0::operator()(double x, double y) const {return -2*B*M_PI*(A*A+B*B)*4*M_PI*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y)-2*M_PI*C*sin(2*C*M_PI*x)*cos(2*D*M_PI*y);}
double SinusBis2DRhsV1::operator()(double x, double y) const {return  2*A*M_PI*(A*A+B*B)*4*M_PI*M_PI*cos(A*2*M_PI*x)*sin(B*2*M_PI*y)-2*M_PI*D*cos(2*C*M_PI*x)*sin(2*D*M_PI*y);}

double Sinus2DSolutionP::operator()(double x, double y) const {return 0.0;}
double Sinus2DSolutionV0::operator()(double x, double y) const {return cos(B*2*M_PI*y);}
double Sinus2DSolutionV1::operator()(double x, double y) const {return cos(A*2*M_PI*x);}
// double Sinus2DSolutionV0::operator()(double x, double y) const {return y*(1-y);}
// double Sinus2DSolutionV1::operator()(double x, double y) const {return x*(1-x);}

double Sinus2DRhsP::operator()(double x, double y) const {return  0.0;}
double Sinus2DRhsV0::operator()(double x, double y) const {return B*B*4*M_PI*M_PI*cos(B*2*M_PI*y);}
double Sinus2DRhsV1::operator()(double x, double y) const {return A*A*4*M_PI*M_PI*cos(A*2*M_PI*x);}
// double Sinus2DRhsV0::operator()(double x, double y) const {return 2.0;}
// double Sinus2DRhsV1::operator()(double x, double y) const {return 2.0;}

double Linear2DSolutionP::operator()(double x, double y) const {return 0.0;}
double Linear2DSolutionV0::operator()(double x, double y) const {return 1.0;}
double Linear2DSolutionV1::operator()(double x, double y) const {return 0.0;}

double Linear2DRhsP::operator()(double x, double y) const {return  0.0;}
double Linear2DRhsV0::operator()(double x, double y) const {return 0.0;}
double Linear2DRhsV1::operator()(double x, double y) const {return 0.0;}

// /*-------------------------------------------------*/
// std::string StokesApplication::toString() const
// {
//     return _name;
// }

/*-------------------------------------------------*/
void StokesApplication::set_application(int dim, std::string name)
{
    // _sol = std::make_shared<SystemAnalyticalFunction>(dim+1);
    // _rhs = std::make_shared<SystemAnalyticalFunction>(dim+1);
    // _boundaryconditions = std::make_shared<BoundaryConditions>(dim);
    std::vector<std::string> tokens = tokenize(name,"_");
    assert(tokens.size()==2);
    if(tokens[0]=="Sinus")
    {
        _sol["p" ] = std::make_shared<Sinus2DSolutionP>();  
        _sol["v0"] = std::make_shared<Sinus2DSolutionV0>();  
        _sol["v1"] = std::make_shared<Sinus2DSolutionV1>();  
        _rhs["p" ] = std::make_shared<Sinus2DRhsP>();  
        _rhs["v0"] = std::make_shared<Sinus2DRhsV0>();  
        _rhs["v1"] = std::make_shared<Sinus2DRhsV1>();  
    }
    else if(tokens[0]=="SinusBis")
    {
        _sol["p" ] = std::make_shared<SinusBis2DSolutionP>();  
        _sol["v0"] = std::make_shared<SinusBis2DSolutionV0>();  
        _sol["v1"] = std::make_shared<SinusBis2DSolutionV1>();  
        _rhs["p" ] = std::make_shared<SinusBis2DRhsP>();  
        _rhs["v0"] = std::make_shared<SinusBis2DRhsV0>();  
        _rhs["v1"] = std::make_shared<SinusBis2DRhsV1>();  
    }
    else if(tokens[0]=="Linear")
    {
        _sol["p" ] = std::make_shared<Linear2DSolutionP>();  
        _sol["v0"] = std::make_shared<Linear2DSolutionV0>();  
        _sol["v1"] = std::make_shared<Linear2DSolutionV1>();  
        _rhs["p" ] = std::make_shared<Linear2DRhsP>();  
        _rhs["v0"] = std::make_shared<Linear2DRhsV0>();  
        _rhs["v1"] = std::make_shared<Linear2DRhsV1>();  
    }
    else
    {
        _not_written_("unknown app name="+name);
    }
    if(tokens[1]=="per")
    {
        for(auto& b: *_boundaryconditions)
        {
            b[0] = "per";
            b[1] = "per";
        }        
    }
    else if(tokens[1]=="dir")
    {
        for(auto& b: *_boundaryconditions)
        {
            b[0] = "dir";
            b[1] = "dir";
        }
        // std::cerr << "has_solution() = " << has_solution() << "\n";
        if(has_solution())
        {
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<2;j++)
                {
                    _boundaryconditions->get_bf(i,j)["p"] = _sol["p"];
                    for(int idim=0;idim<dim;idim++)
                    {
                        std::stringstream ss;
                        ss << "v" << idim;        
                        _boundaryconditions->get_bf(i,j)[ss.str()] = _sol[ss.str()];
            
                    }            
                }            
            }            
        }        
    }
    else
    {
        _not_written_("unknown name="+name);
    }
}
