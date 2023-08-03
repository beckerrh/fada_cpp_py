//
//  timer.hpp
//  Fada
//
//  Created by Roland Becker on 18/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef timer_hpp
#define timer_hpp

#include  <stdio.h>
#include  <ctime>
#include  <map>
#include  <set>

/*-------------------------------------------------*/
class Timer : public std::map<std::string, double>
{
protected:
    bool _print, _debug;
    std::map<std::string, double> _sum;
    std::map<std::string, double> _temp;
    inline double seconds(void)
    {
        static const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
        return ( (double) clock() ) * secs_per_tick;
    }
    double total() const;
    typedef std::map<std::string, double>::const_iterator const_iterator;
    typedef std::map<std::string, double>::iterator iterator;

public:
    ~Timer();
    Timer(bool print=true, bool debug=false, std::set<std::string> names={}) : std::map<std::string, double>(), _print(print), _debug(debug) 
    {
        for(auto p: names)
        {
            enrol(p);
        }
    }
    void enrol(std::string name, bool sum = true);
    void start(std::string name);
    void stop(std::string name);
    double get(std::string name);
    void print(std::ostream& os) const;
};

#endif /* timer_hpp */
