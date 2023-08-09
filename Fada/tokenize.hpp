//
//  tokenize.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef tokenize_h
#define tokenize_h

// #define _LONG_LONG

#include  <string>
#include  <vector>
#ifdef GNU
#include  <bits/stdc++.h>
#endif

static std::vector<std::string> tokenize(std::string s, std::string del) {
    std::vector<std::string> tokens;
    int end = s.find(del); 
    while (end != -1) {
        tokens.push_back(s.substr(0, end));
        s.erase(s.begin(), s.begin() + end + 1);
        end = s.find(del);
    }
    tokens.push_back(s.substr(0, end));
    return tokens;
}

#endif
