//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "Fada/parametermap.hpp"

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    ParameterMap parameters;
    parameters.add("toto", std::vector<std::string>{"a", "b"});
    parameters.add("tata", "abra");
    parameters.add("titi", 123);
    parameters.add("tete", 123.321);
    std::cout << parameters;
    parameters.fill_from_map(args2map(argc, argv));
    std::cout << parameters;
}
