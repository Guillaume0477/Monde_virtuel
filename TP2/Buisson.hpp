#include "Arbre.hpp"
#include "utils.hpp"
#include <stdlib.h>
#include <stdio.h>

#ifndef BUISSONCLASS
#define BUISSONCLASS

class Buisson : public Arbre
{
protected:

public:
    Buisson(){
        rayon = 1;
        nombre_attrib = 3;
    };
    double get_rayon() const {
        return rayon;
    }
    double humidity(double hum){
        double value = fonction_decroissante(0.0,0.2,hum);
        return (value);
    };
    double slope(double slo){
        double value = fonction_decroissante(0.0,0.3,slo);
        return (value);
    };
    double stream(double str){
        double value = fonction_decroissante(0.0,0.01,str);
        return (value);
    };

};

#endif
