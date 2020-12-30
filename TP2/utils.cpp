#include "utils.hpp"


float fonction_one(float lim_inf, float lim_sup,float current_num){
    if ((current_num > lim_inf) && (current_num < lim_sup)){
        return current_num;
    }
    else {
        return 0;
    }
}

float fonction_zeros(float lim_inf, float lim_sup,float current_num){
    if ((current_num > lim_inf) && (current_num < lim_sup)){
        return current_num;
    }
    else {
        return 1;
    }
}



bool test_dist(std::pair< std::pair<int,int> , int> couple1, std::pair< std::pair<int,int> , int> couple2){
    int diff_x = (couple1.first.first - couple2.first.first);
    int diff_y = (couple1.first.second - couple2.first.second);
    float dist = std::sqrt( diff_x*diff_x + diff_y*diff_y );
    if (dist < (couple1.second + couple2.second)){
        return true;
    }
    else
    {
        return false;
    }
    

};
