#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <math.h>


#ifndef UTILS
#define UTILS


float fonction_croissante_par_morceau(float lim_inf, float lim_sup,float current_num, float start);
float fonction_croissante(float, float,float);
float fonction_decroissante(float, float,float);


bool test_dist(std::pair< std::pair<int,int> , int>, std::pair< std::pair<int,int> , int>);


#endif