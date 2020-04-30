#ifndef __info_h
#define __info_h

class Info
{
 protected:
  const   static  char* std_txt;

  int     lastit,iputada,ierg;
  double  first_res,last_res;
  double  iresidual,itol,iglobaltol,ireduction,ireductionrate;
  double  ilastreductionrate;
  int     imaxiter,iiteration,iprintstep;

  double  compute_reduction_rate() const;

 public:
  const   char* text;
  Info(double f = 1.e-8, double r=1.e-15, int p = 10, 
	     int m = 5000, const char* txt = std_txt);

  /*   Zugriff   */

  double  residual()  const { return iresidual; }
  double& residual()        { return iresidual; }
  double  tol()       const { return itol; }
  double& tol()             { return itol; }
  double  globaltol() const { return iglobaltol; }
  double& globaltol()       { return iglobaltol; }
  double  reduction() const { return ireduction; }
  double& reduction()       { return ireduction; }
  int     maxiter()   const { return imaxiter; }
  int&    maxiter()         { return imaxiter; }
  int     iteration() const { return iiteration; }
  int&    iteration()       { return iiteration; }
  int     printstep() const { return iprintstep; }
  int&    printstep()       { return iprintstep; }
  int     putada()    const { return iputada; }
  int&    putada()          { return iputada; }
  int     erg()    const { return ierg; }
  int&    erg()          { return ierg; }
  double  reduction_rate()      const { return ireductionrate; }
  double& reduction_rate()            { return ireductionrate; }
  double  last_reduction_rate() const { return ilastreductionrate; }
  double& last_reduction_rate()       { return ilastreductionrate; }

  /*   Funktionen   */

  void    reset();
  int     check_residual(int i, double res);
};

#endif
