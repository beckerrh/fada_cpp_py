//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  parametermap_hpp
#define  parametermap_hpp

#include  <map>
#include  <unordered_map>
#include  <string>
#include  <any>
#include  <iostream>
#include  <vector>
#include  <sstream>

std::map<std::string,std::string> args2map(int argc, char** argv);

/*-------------------------------------------------*/
class ParameterMap
{
protected:
    std::unordered_map<std::string, std::any> _values;
    std::unordered_map<std::string, std::string> _types;
    std::unordered_map<std::string, std::vector<std::string>> _defaults_string;
    std::unordered_map<std::string, std::vector<int>> _defaults_int;
    std::unordered_map<std::string, std::vector<double>> _defaults_double;
    void check_if_already(std::string name) const
    {
        if(_values.find(name)!=_values.end())
        {
            std::cerr << "parameter '" << name << "' already given!\n";
            exit(1);
        }        
    }
    void check_if_exists(std::string name) const
    {
        if(_values.find(name)==_values.end())
        {
            std::cerr << "parameter '" << name << "' doesn't exist!\n";
            exit(1);
        }        
    }
    void _set(std::string name, std::string value);

public:
    ParameterMap() {}
    ParameterMap(const ParameterMap& pmap): _values(pmap._values), _types(pmap._types) {}

    std::vector<std::string> get_keys() const
    {
        std::vector<std::string> keys;
        for(std::unordered_map<std::string, std::any>::const_iterator it = _values.begin(); it != _values.end(); ++it) { keys.push_back(it->first);}
        return keys;
    }
    const std::string& get_type(std::string name) const {return _types.at(name);}
    bool has_key(std::string name) const {return _values.find(name)!=_values.end();}
    
    
    template<typename T>
    T const get(std::string name) const {check_if_exists(name); return std::any_cast<T>(_values.at(name));}
    
    void add(std::string name, const std::string& defaultval)
    {
        check_if_already(name);
        _values[name] = defaultval;
        _types[name] = "string";
        _defaults_string[name] = std::vector<std::string>{defaultval};
    }
    void add(std::string name, std::vector<std::string> defaults)
    {
        check_if_already(name);
        assert(defaults.size());
        _values[name] = defaults[0];
        _types[name] = "string";
        _defaults_string[name] = defaults;
    }
    // void add(std::string name, const bool& defaultval)
    // {
    //     std::cerr << "adding " << name << " -> " << defaultval << "\n";
    //     check_if_already(name);
    //     _values[name] = defaultval;
    //     _types[name] = "bool";
    // }
    void add(std::string name, int defaultval)
    {
        check_if_already(name);
        _values[name] = defaultval;
        _types[name] = "int";
        _defaults_int[name] = std::vector<int>{defaultval};
    }
    void add(std::string name, std::vector<int> defaults)
    {
        check_if_already(name);
        assert(defaults.size());
        _values[name] = defaults[0];
        _types[name] = "int";
        _defaults_int[name] = defaults;
    }
    void add(std::string name, double defaultval)
    {
        check_if_already(name);
        _values[name] = defaultval;
        _types[name] = "double";
        _defaults_double[name] = std::vector<double>{defaultval};
    }
    void add(std::string name, std::vector<double> defaults)
    {
        check_if_already(name);
        assert(defaults.size());
        _values[name] = defaults[0];
        _types[name] = "double";
        _defaults_double[name] = defaults;
    }
    // void set(std::string name, bool defaultval)
    // {
    //     _values[name] = defaultval;
    //     _types[name] = "bool";
    // }
    void set(std::string name, std::string defaultval)
    {
        _values[name] = defaultval;
        _types[name] = "string";
    }
    void set(std::string name, int defaultval)
    {
        _values[name] = defaultval;
        _types[name] = "int";
    }
    void set(std::string name, double defaultval)
    {
        _values[name] = defaultval;
        _types[name] = "double";
    }
    // void fill_from_args(int argc, char** argv);
    void fill_from_map(const std::map<std::string,std::string>& parameters);
};
std::ostream& operator<<(std::ostream& os, const ParameterMap& p);


#endif
