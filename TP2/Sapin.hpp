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
        rayon = 6;
        nombre_attrib = 3;
    };
    double get_rayon() const {
        return rayon;
    }
    double humidity(double hum){
        double value = fonction_croissante_par_morceau(0.1,1.0,hum,0.8);
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
    double acces(double acc){
        double value = fonction_croissante(0.60,1.0,acc);
        return (value);
    };

};

#endif
