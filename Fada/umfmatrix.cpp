//
//  umfmatrix.cpp
//  Fada
//
//  Created by Roland Becker on 19/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "umfmatrix.hpp"


#define UMFPACK_OK       0
#define UMFPACK_INFO    90
#define UMFPACK_CONTROL 20

#define UMFPACK_A       ( 0 )     /* Ax=b    */
#define UMFPACK_At      ( 1 )     /* A'x=b  */

/*-------------------------------------------------*/
#ifdef _LONG_LONG
extern "C" int umfpack_dl_symbolic
 (
  long long n,
  long long m,
  const long long Ap [],
  const long long Ai [],
  const double Ax [],
  void** Symbolic,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
  );
extern "C" int umfpack_dl_numeric
 (
  const long long Ap [],
  const long long Ai [],
  const double Ax [],
  void* Symbolic,
  void** Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
  );
extern "C" int umfpack_dl_solve
 (
  long long sys,
  const long long Ap [],
  const long long Ai [],
  const double Ax [],
  double X [],
  const double B [],
  void* Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
  );
extern "C" void umfpack_dl_free_symbolic(void** Symbolic);
extern "C" void umfpack_dl_free_numeric(void** Numeric);
extern "C" void umfpack_dl_report_status(const double Control [UMFPACK_CONTROL], long long status);
extern "C" void umfpack_dl_report_info(const double Control [UMFPACK_CONTROL], const double Info [UMFPACK_INFO]);
extern "C" int  umfpack_dl_report_numeric(const char name [], void* Numeric, const double Control [UMFPACK_CONTROL]);
extern "C" void umfpack_dl_defaults(const double Control [UMFPACK_CONTROL]);
#else
extern "C" int umfpack_di_symbolic
(
 int n,
 int m,
 const int Ap [],
 const int Ai [],
 const double Ax [],
 void** Symbolic,
 const double Control [UMFPACK_CONTROL],
 double Info [UMFPACK_INFO]
 );
extern "C" int umfpack_di_numeric
(
 const int Ap [],
 const int Ai [],
 const double Ax [],
 void* Symbolic,
 void** Numeric,
 const double Control [UMFPACK_CONTROL],
 double Info [UMFPACK_INFO]
 );
extern "C" int umfpack_di_solve
(
 int sys,
 const int Ap [],
 const int Ai [],
 const double Ax [],
 double X [],
 const double B [],
 void* Numeric,
 const double Control [UMFPACK_CONTROL],
 double Info [UMFPACK_INFO]
 );
extern "C" void umfpack_di_free_symbolic(void** Symbolic);
extern "C" void umfpack_di_free_numeric(void** Numeric);
extern "C" void umfpack_di_report_status(const double Control [UMFPACK_CONTROL], int status);
extern "C" void umfpack_di_report_info(const double Control [UMFPACK_CONTROL], const double Info [UMFPACK_INFO]);
extern "C" int  umfpack_di_report_numeric(const char name [], void* Numeric, const double Control [UMFPACK_CONTROL]);
extern "C" void umfpack_di_defaults(const double Control [UMFPACK_CONTROL]);
#endif
/*-------------------------------------------------*/
UmfMatrix::~UmfMatrix()
{
  if(Symbolic)
  {
    #ifdef _LONG_LONG
      umfpack_dl_free_symbolic (&Symbolic);
    #else
      umfpack_di_free_symbolic (&Symbolic);
    #endif
  }
  if(Numeric)
  {
    #ifdef _LONG_LONG
      umfpack_dl_free_numeric (&Numeric);
    #else
      umfpack_di_free_numeric (&Numeric);
    #endif
  }
  if(Control)
  {
    delete[] Control;
    Control = NULL;
  }
  if(Info)
  {
    delete[] Info;
    Info = NULL;
  }
}
/*-------------------------------------------------*/
UmfMatrix::UmfMatrix() : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
{
  Control = new double[UMFPACK_CONTROL];
  #ifdef _LONG_LONG
    umfpack_dl_defaults(Control);
  #else
    umfpack_di_defaults(Control);
  #endif
  Control[0] = 2;
}
/*-------------------------------------------------*/
UmfMatrix::UmfMatrix( const UmfMatrix& umfmatrixbase) : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
{
  Control = new double[UMFPACK_CONTROL];
  #ifdef _LONG_LONG
    umfpack_dl_defaults(Control);
  #else
    umfpack_di_defaults(Control);
  #endif
  Control[0] = 2;
}

/*---------------------------------------------------------*/
void UmfMatrix::save(std::ostream& os, arma::file_type datatype) const
{
  char name[] = "toto";
  #ifdef _LONG_LONG
    umfpack_dl_report_numeric(name, (void*) &Numeric, Control);
  #else
    umfpack_di_report_numeric(name, (void*) &Numeric, Control);
  #endif
}
/*---------------------------------------------------------*/
void UmfMatrix::init()
{
  n = _sp.nrows();
  mb = _sp.values().begin();
  #ifdef _LONG_LONG
    umfpack_dl_free_symbolic (&Symbolic);
    sb = (const long long*) _sp.rows().begin();
    cb = (const long long*) _sp.cols().begin();
    int status = umfpack_dl_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      std::cerr<<"*** ERROR UmfMatrix: umfpack_symbolic failed\n";
      _sp.save(std::cerr);
      assert(0);
      exit(1);
    }
  #else
    umfpack_di_free_symbolic (&Symbolic);
    sb = (const int*) _sp.rows().begin();
    cb = (const int*) _sp.cols().begin();
    int status = umfpack_di_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control, Info);
      umfpack_di_report_status(Control, status);
      std::cerr<<"*** ERROR UmfMatrix: umfpack_symbolic failed\n";
      _sp.save(std::cerr);
      assert(0);
      exit(1);
    }
  #endif
}
/*---------------------------------------------------------*/
void UmfMatrix::computeLu()
{
  init();
  #ifdef _LONG_LONG
    umfpack_dl_free_numeric (&Numeric);
    int status = umfpack_dl_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      std::ofstream file("MATRIX_NOT_OK");
      _sp.save(file);
      assert(0);
      exit(1);
    }
  #else
    umfpack_di_free_numeric (&Numeric);
    int status = umfpack_di_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control, Info);
      umfpack_di_report_status(Control, status);
      std::ofstream file("MATRIX_NOT_OK");
      _sp.save(file);
      assert(0);
      exit(1);
    }
  #endif
}
/*---------------------------------------------------------*/
void UmfMatrix::solve(armavec& x, const armavec& b) const
{
  assert( x.size() == b.size() );
  assert( x.size() == n );
  double* xb = x.begin() ;
  const double* bb = b.begin();
  #ifdef _LONG_LONG
    int status = umfpack_dl_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      std::ofstream file("MATRIX_NOT_OK");
      //      _sp.save(file, arma::arma_ascii);
      std::cerr<<"*** ERROR UmfMatrix::Solve(): umfpack_dl_solve failed\n";
      assert(0);
      exit(1);
    }
  #else
    int status = umfpack_di_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info);
    if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control, Info);
      umfpack_di_report_status(Control, status);
      std::ofstream file("MATRIX_NOT_OK");
      //      _sp.save(file, arma::arma_ascii);
      std::cerr<<"*** ERROR UmfMatrix::Solve(): umfpack_dl_solve failed\n";
      assert(0);
      exit(1);
    }
  #endif
}
