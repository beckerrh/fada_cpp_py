//
//  testarma.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  <iostream>
#include  <armadillo>
#include  "Fada/vector.hpp"


/*-------------------------------------------------*/
int main(int argc, char** argv)
{
  armaicvec n;
  n << 4 << 2 << 3 << arma::endr;
  vector v(n), w(n);
  v.fill(2);
  w.randu();
  
  vector x(v);
  
  x = 3*v.arma() + w.arma();
  std::cerr << x << std::endl;
  std::cerr << "x(2,1,2) =" << x(2,1,2) << std::endl;
  std::cerr << arma::dot(v.arma(),w.arma()) << std::endl;
}
