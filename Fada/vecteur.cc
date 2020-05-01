#include "vecteur.h"

///**************************************************/
//
//Vecteur::Vecteur(int nn, int mm)
//{
//  (*this).reinit(nn,mm);
//}
//
///**************************************************/
//
//Vecteur::Vecteur()
//{
//  pn() = 0;
//  pm() = 0;
//}
//
///**************************************************/
//
//Vecteur::Vecteur(const Vecteur& a)
//{
//  (*this).reinit(a.n(),a.m());
//  val() = a.val();
//}
//
///**************************************************/
//
//Vecteur& Vecteur::operator=(const Vecteur& a)
//{
//  val() = a.val();
//  return *this;
//}
//
///**************************************************/
//
//Vecteur& Vecteur::operator=( double a)
//{
//  val() = a;
//  return *this;
//}
//
///**************************************************/
//
//double Vecteur::operator* (const Vecteur& v) const
//{
//  double d = 0.;
//  for(int i=0;i<n()*m();i++)
//    {
//      d +=  val()(i)*v.val()(i);
//    }
//  return d;
//}
//
///**************************************************/
//
//void Vecteur::equ (double d, const Vecteur& v)
//{
//  for(int i=0;i<n()*m();i++)
//    {
//      val()(i) = d * v.val()(i);
//    }
//}
//
///**************************************************/
//
//void Vecteur::add (double d, const Vecteur& v)
//{
//  for(int i=0;i<n()*m();i++)
//    {
//      val()(i) += d * v.val()(i);
//    }
//}
//
///**************************************************/
//
//void Vecteur::sadd(double e, double d, const Vecteur& v)
//{
//  for(int i=0;i<n()*m();i++)
//    {
//      val()(i) = e * val()(i)  +  d * v.val()(i);
//    }
//}

/**************************************************/

void Vecteur::boundary(const Vecteur& v)
{
  for(int i=0;i<n();i++)
    {
      (*this)(i,0)     = v(i,0);
      (*this)(i,m()-1) = v(i,m()-1);
    }
  for(int j=0;j<m();j++)
    {
      (*this)(0,j)     = v(0,j);
      (*this)(n()-1,j) = v(n()-1,j);
    }
}

/**************************************************/

void Vecteur::boundary()
{
  for(int i=0;i<n();i++)
    {
      (*this)(i,0)     = 0.;
      (*this)(i,m()-1) = 0.;
    }			 
  for(int j=0;j<m();j++) 
    {			 
      (*this)(0,j)     = 0.;
      (*this)(n()-1,j) = 0.;
    }
}

/**************************************************/

void Vecteur::right()
{
  for(int i=0;i<n();i++)
    {
      for(int j=0;j<m();j++) 
	{			 
	  //	  (*this)(i,j)     = (1.*i+j)/(n()*m());
	  (*this)(i,j)     = 1.;///(n()*m());
	}
    }
}

///**************************************************/
//
//void Vecteur::null()
//{
//  for(int i=0;i<n()*m();i++)
//    {
//      val()(i) = 0.;
//    }
//}
//
///**************************************************/
//
//void Vecteur::output_plotmtv()
//{
//  char* name = new char[200];
//  sprintf(name,"fada.mtv");
//  FILE* fp = fopen(name,"w");
//  if(!fp) printf("Keoin output\n");
//  
//  fprintf(fp,"$ DATA=CONTOUR\n\n");
//  fprintf(fp,"#\n");
//  fprintf(fp,"# This one draws with colors \n");
//  fprintf(fp,"#\n\n");
//  fprintf(fp,"%% xmin=0.0 xmax=1.0\n");
//  fprintf(fp,"%% ymin=0.0 ymax=1.0\n");
//  fprintf(fp,"%% nx=%d ny=%d\n",n(),m());
//  fprintf(fp,"%% toplabel = \"Voici le resultat fada\"\n");
//  fprintf(fp,"%% subtitle = \"c/o roland@gaia.iwr.uni-heidelberg.de\"\n");
//  fprintf(fp,"#%% contclip=T\n");
//  fprintf(fp,"#%% meshplot=T\n");
//  fprintf(fp,"#%% interp=3\n");
//  fprintf(fp,"%% contfill \n");
//
//  int count = 0;
//  for(int i=0;i<n();i++)
//    {
//      for(int j=0;j<m();j++) 
//	{
//	  if(!((count++)%8)) fprintf(fp,"\n");
//	  fprintf(fp,"%g ",(*this)(j,i));
//	}
//    }
//
//  fprintf(fp,"\n\n$ END");
//
//  fclose(fp);
//  delete[] name;
//}
