//
//  vecteur.cpp
//  Fada
//
//  Created by Roland Becker on 01/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "vecteur.h"

double Vecteur::operator* (const Vecteur& v) const
{
//    return arma::dot(_val, v.val());
  double d = 0.;
  for(int i=0;i<_nx*_ny;i++)
  {
    d +=  val()(i)*v.val()(i);
  }
  return d;
}
