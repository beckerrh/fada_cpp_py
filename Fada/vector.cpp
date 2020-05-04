//
//  vector.c
//  Fada
//
//  Created by Roland Becker on 03/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  "vector.hpp"

std::ostream& operator<<(std::ostream& os, const vector& v)
{
  const armavec& tarma =static_cast<const armavec&>(v);
  os << tarma.t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
  return os;
}
