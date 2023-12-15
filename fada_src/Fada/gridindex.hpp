//
//  uniformmultigrid.hpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef gridindex_hpp
#define gridindex_hpp

#include  <armadillo>

/*-------------------------------------------------*/
class GridIndex
{
protected:
    armaicvec _ofs;
    void set_ofs(const armaicvec& n)
    {
        _ofs.set_size(n.n_elem);
        _ofs.fill(1);
        for(int i=0;i<n.n_elem;i++)
        {
            for(int j=i+1;j<n.n_elem;j++)
            {
                _ofs[i] *= n[j];
            }
        }
    }    
public:
  GridIndex(const armaicvec& n) {set_ofs(n);}
  inline int operator()(int ix, int iy)
  {
      return _ofs[0]*ix+iy;
  }
  inline int operator()(int ix, int iy, int iz)
  {
      return _ofs[0]*ix+_ofs[1]*iy+iz;
  }
};



#endif /* uniformmultigrid_hpp */
