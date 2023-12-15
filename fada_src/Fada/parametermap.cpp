//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//
#include  <parametermap.hpp>
#include  <cstring>

std::ostream& operator<<(std::ostream& os, const ParameterMap& p)
{
    for(std::string k : p.get_keys())
    {
        std::string type = p.get_type(k);        
        os << k << " (" <<type << ") -> ";
        if(type=="string")
        {
            os << p.get<std::string>(k);
        }        
        else if(type=="int")
        {
            os << p.get<int>(k);
        }        
        else if(type=="double")
        {
            os << p.get<double>(k);
        }        
        else
        {
            os << p.get<int>(k);
        }        
        os  << "\n";
    }
    return os;
}
/*-------------------------------------------------*/
void ParameterMap::_set(std::string p, std::string value)
{
    // std::cerr << "ParameterMap " << p << " - " << get_type(p) << " - " << value << "\n";
    if(get_type(p)=="string") _values[p] = value;
    else if(get_type(p)=="int") _values[p] = stoi(value);
    else if(get_type(p)=="double") _values[p] = stod(value);
    else
    {
        std::cerr << p << " - " <<get_type(p) << " don't understand '" << value << "' for bool\n";
    //     if(value=="true" or value=="True" or value=="1") _values[p] = true;
    //     else if(value=="false" or value=="False" or value=="0") _values[p] = false;
    //     else
    //         {
    //             std::cerr << get_type(p) << " don't understand '" << value << "' for bool\n";
    //             exit(1);
    //         }
    }
}


/*-------------------------------------------------*/
std::map<std::string,std::string> args2map(int argc, char** argv)
{
    std::map<std::string,std::string> pmap;
    for(int i=1; i<argc;i++)
    {
        /* leave out first, for '-arg' */
        pmap[std::string(argv[i]).substr(1)] = std::string(argv[i+1]);
        // for(auto& p : keys)
        // {
        //     std::string pm = "-"+p;
        //     if(!std::strcmp(argv[i], pm.c_str()))
        //     {
        //         _set(p, std::string(argv[i+1]));
        //     }
        // }
        i++;
    }
    return pmap;
}
// /*-------------------------------------------------*/
// void ParameterMap::fill_from_args(int argc, char** argv)
// {
//     std::vector<std::string> keys = get_keys();
//     for(int i=1; i<argc;i++)
//     {
//         std::cerr << "i="<<i<<" argv="<<argv[i]<<" "<<argv[i+1]<<"\n";
//         for(auto& p : keys)
//         {
//             std::string pm = "-"+p;
//             if(!std::strcmp(argv[i], pm.c_str()))
//             {
//                 _set(p, std::string(argv[i+1]));
//             }
//         }
//         i++;
//     }
// }
/*-------------------------------------------------*/
void ParameterMap::fill_from_map(const std::map<std::string,std::string>& parameters)
{
    for(auto it : parameters)
    {
        _set(it.first, it.second);
    }
}

