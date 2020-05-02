#ifndef __array_h
#define __array_h

/*-------------------------------------------------*/
template<class T>
class Array
{
private:
  int dim;
  T*  ptr;
  void copy(const Array<T>&);
public:
  Array(int n);
  Array();
  Array(const Array<T>&);
  ~Array();
  const T& operator()(int i) const;
  T& operator()(int i);
  int n() const;
  Array<T>& operator=(const Array<T>&);
  Array<T>& operator=(const T);
  void set_size(int n);
};

/*-------------------------------------------------*/
template<class T>
inline void Array<T>::copy(const Array<T>& a)
{
//  THROW1(dim!=a.dim, IntError(IntError::IllegalDimension, n, "Array"));
  T* p = ptr + dim;
  T* q = a.ptr + dim;
  while(p>ptr) *--p = *--q;
}

template<class T>
inline Array<T>::Array(int n)
{
  dim = n;
  ptr = new T[n];
  //THROWUNCOND(!ptr, Error( Error::NoMem,"Array::Array"));
 }

template<class T>
inline Array<T>::Array()
{
  dim = 0;
  ptr = 0;
}

template<class T>
inline Array<T>::Array(const Array<T>& a)
{
  dim = a.dim;
  ptr = new T[dim];
  //THROWUNCOND(!ptr, Error(Error::NoMem,"Array::Array"));
  copy(a);
}

template<class T>
inline Array<T>::~Array()
{
  delete[] ptr;
}

template<class T>
inline const T& Array<T>::operator()(int i) const
{
  return ptr[i];
}

template<class T>
inline T& Array<T>::operator()(int i)
{
  return ptr[i];
}

template<class T>
inline int Array<T>::n() const
{
  return dim;
}

template<class T>
inline Array<T>& Array<T>::operator=(const Array<T>& a)
{
  if(ptr!=a.ptr)
  {
    set_size(a.dim);
    copy(a);
  }
  return *this;
}

template<class T>
inline Array<T>& Array<T>::operator=(const T a)
{
  T* p = ptr + dim;
  while(p>ptr) *--p = a;
  return *this;
}

template<class T>
inline void Array<T>::set_size(int n)
{
  if(n!=dim)
  {
    if(ptr) delete[] ptr;
    dim = n;
    ptr = new T[dim];
   }
}

/*-------------------------------------------------*/
#endif
