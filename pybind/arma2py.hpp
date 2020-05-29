//
//  arma2py.hpp
//  pyfada
//
//  Created by Roland Becker on 26/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef arma2py_h
#define arma2py_h

#include <armadillo>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/buffer_info.h>
#include <pybind11/detail/common.h>

/*-------------------------------------------------*/
template <typename T> arma::Col<T> arr2col(pybind11::array_t<T>& src)
{
  pybind11::buffer_info info = src.request();
  auto dims = src.ndim();
  if ((dims != 1) && (src.shape(1) != 1)) {
       throw std::runtime_error("arr2col(): Number of columns must <= 1");
   }
  auto ptr = static_cast<T*>(info.ptr);
  bool copy=false, strict=false;
  return arma::Col<T>(static_cast<T *>(info.ptr), src.size(), copy, strict);
}

/*-------------------------------------------------*/
template <typename T> pybind11::array_t<T> col2arr(arma::Col<T>& src)
{
  ssize_t tsize =  static_cast<ssize_t>(sizeof(T));
  ssize_t nrows = static_cast<ssize_t>(src.n_rows);
  T* data = src.memptr();

  pybind11::capsule free_when_done(data, [](void *f) {
       T *foo = reinterpret_cast<T *>(f);
       std::cerr << "col2arr (not) freeing memory @ " << f << "\n";
//      delete[] foo;
   });

  return pybind11::array_t<T>(
      {nrows, static_cast<ssize_t>(1)}, // shape
      {tsize, nrows * tsize}, // contiguous strides
      data, // the data pointer
      free_when_done // sans ça, la mémoire n'est pas délibérée....
     );
}


/*-------------------------------------------------*/
template <typename T> arma::Mat<T> arr2mat(pybind11::array_t<T>& src)
{
  pybind11::buffer_info info = src.request();
  auto dims = src.ndim();
  if (dims != 2)
  {
       throw std::runtime_error("arr2mat(): not a matrix");
   }
  auto ptr = static_cast<T*>(info.ptr);
  bool copy=false, strict=false;
  std::cerr << "arr2mat() " <<info.shape[0]<< " " <<  info.shape[1] << "\n";
  return arma::Mat<T>(static_cast<T *>(info.ptr), info.shape[0], info.shape[1], copy, strict);
}

/*-------------------------------------------------*/
template <typename T> pybind11::array_t<T> mat2arr(arma::Mat<T>& src)
{
  ssize_t tsize =  static_cast<ssize_t>(sizeof(T));
  ssize_t nrows = static_cast<ssize_t>(src.n_rows);
  ssize_t ncols = static_cast<ssize_t>(src.n_cols);
  T* data = src.memptr();
  arma::access::rw(src.mem) = 0;

  std::cerr << "mat2arr() " <<nrows<< " " <<  ncols << " data " << data << "\n";

  pybind11::capsule free_when_done(data, [](void *f) {
       T *foo = reinterpret_cast<T *>(f);
       std::cerr << "mat2arr (not) freeing memory @ " << f << "\n";
//      delete[] foo;
   });

  return pybind11::array_t<T>(
      {nrows, ncols}, // shape
      {tsize, nrows * tsize}, // F-style contiguous strides
      data, // the data pointer
      free_when_done // sans ça, la mémoire n'est pas délibérée....
     );
}

#endif /* arma2py_h */
