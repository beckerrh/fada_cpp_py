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
    armaicvec _n, _ofs;
    void set_ofs()
    {
        _ofs.set_size(_n.n_elem);
        _ofs.fill(1);
        for(int i=0;i<_n.n_elem;i++)
        {
            for(int j=i+1;j<_n.n_elem;j++)
            {
                _ofs[i] *= _n[j];
            }
        }
    }    
public:
  GridIndex(const armaicvec& n) : _n(n) {set_ofs();}
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
