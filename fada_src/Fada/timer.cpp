//
//  timer.cpp
//  Fada
//
//  Created by Roland Becker on 18/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <iostream>
#include  <iomanip>
#include  "timer.hpp"

/*-------------------------------------------------*/
Timer::~Timer()
{
  if(_print) print(std::cout);
}

/*-------------------------------------------------*/
void Timer::enrol(std::string name, bool sum)
{
  (*this)[name] = 0.0;
  _sum[name] = true;
}
/*-------------------------------------------------*/
void Timer::start(std::string name)
{
  if(_debug) std::cerr << "Timer: start " << name << "\n";
  _temp[name] = seconds();
}
/*-------------------------------------------------*/
void Timer::stop(std::string name)
{
  if(_debug) std::cerr << "Timer: stop " << name << "\n";
  (*this)[name] += seconds()-_temp[name];
}
/*-------------------------------------------------*/
double Timer::get(std::string name)
{
  return (*this)[name];
}
/*-------------------------------------------------*/
double Timer::total() const
{
  double d = 0.;
  for(const_iterator p = begin(); p != end(); p++)
  {
    if(_sum.find(p->first)->second)
    {
      d += p->second;
    }
  }
  return d;
}

/*-------------------------------------------------*/
void Timer::print(std::ostream& os) const
{
  double totaltime = total();
  os << "----- Total -----:  "  << std::setiosflags(std::ios::fixed)<< std::setprecision(3) << totaltime << " s\n";
  for(const_iterator p = begin(); p != end(); p++)
  {
    double singletime = p->second;
    if(_sum.find(p->first)->second)
    {
      os << std::setiosflags(std::ios::left);
      os << std::setw(15) << p->first << ": " << std::setiosflags(std::ios::fixed) << std::setprecision(3);
      os << std::resetiosflags(std::ios::left);
      os << std::setw(5) << singletime  << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(1)<<  100.0*singletime/totaltime << " \% \n";
    }
    else
    {
      os << std::setiosflags(std::ios::left);
      os << std::setw(15) << p->first << ": " << std::setiosflags(std::ios::fixed) << std::setprecision(3);
      os << std::resetiosflags(std::ios::left);
      os << std::setw(5) << singletime  << std::setw(12) << "---\n";
    }
  }
  os << "\n";
}
