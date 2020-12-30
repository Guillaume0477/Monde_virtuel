#include "Arbre.hpp"
#include "utils.hpp"
#include <stdlib.h>
#include <stdio.h>


#ifndef SAPINCLASS
#define SAPINCLASS

class Sapin : public Arbre
{
protected:

public:
    Sapin(){
        rayon = 10;
        nombre_attrib = 3;
    };
    double get_rayon() const {
        return rayon;
    }
    double humidity(double hum){
        double value = fonction_one(0.1,1.0,hum);
        return (value);
    };
    double slope(double slo){
        double value = fonction_one(0.0,0.3,slo);
        return (value);
    };
    double stream(double str){
        double value = fonction_one(0.0,0.05,str);
        return (value);
    };

};

#endif
