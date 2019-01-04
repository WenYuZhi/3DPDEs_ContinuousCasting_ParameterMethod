#ifndef STEEL_H
#define STEEL_H
#include <iostream>
#include <map> 
#include <string>  
using namespace std;
class Steel
{
    public:
	    float pho, ce, L;
	    float lamda, lamda_liquid, lamda_solid;
		float temperature_liquid, temperature_solid;
		Steel(map<string, float>);
		void print();
		map<string, float> steel_componment;
	    void physical_parameters(float);
};
#endif
