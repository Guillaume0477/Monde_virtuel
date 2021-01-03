#include <stdlib.h>
#include <stdio.h>

#ifndef ARBRECLASS
#define ARBRECLASS

class Arbre
{
protected:
    int nombre_attrib;
    int rayon;
public:
    Arbre(){
        rayon = 3; 
        nombre_attrib = 3;
    };
    Arbre(const int nb, const int r){
        nombre_attrib= nb;
        rayon = r; 
    };
    virtual double get_rayon() const {
        return rayon;
    }
    virtual double humidity(double hum){
        return hum;
    };
    virtual double slope(double slo){
        return slo;
    };
    virtual double stream(double str){
        return str;
    };

};

#endif