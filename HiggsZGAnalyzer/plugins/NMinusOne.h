/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _NMinusOne_H
#define _NMinusOne_H

#include <string>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

class NMinusOne {
    public:
        NMinusOne() {};
        virtual ~NMinusOne() {};
        NMinusOne(vector<void*> cuts, vector<void*> values, vector<string> ops);
        bool                TestCut(multiType cut, multiType value, string op);
        //vector<multiType>   GetOutput();

    private:
        //input parameters
        vector<void*>  _inputValues;
        vector<void*>  _cuts;
        vector<void*>  _outputValues;
        vector<string>     _ops;


};


#endif

