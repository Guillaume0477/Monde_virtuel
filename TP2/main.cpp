#include <stdlib.h>
#include <stdio.h>
#include <iostream>


class vec2
{
protected:
    double a,b;
public :
    vec2(const double&, const double&);
    double operator[](int);
};

vec2::vec2(const double& x, const double& y){
    a = x;
    b = y;
}

double vec2::operator[](int index){
    if ((index > 2) || (index < 0)){
        std::cout << "Index out of range for vec2 type" << std::endl;
        return -1;
    } else if (index == 0) {
        return a;
    } else {
        return b;
    }
}


class Box2 
{
protected:
    vec2 a,b;

public:
    Box2(const vec2&, const vec2&);
    bool Inside(const vec2&) const;
    bool Intersect(const Box2&) const;

};

Box2::Box2(const vec2& u, const vec2& v){
    a = u;
    b = v;
}

bool Box2::Inside(const vec2& v) const {
    if ((v[0] > b[0]) || (v[0] < a[0]) || (v[1] < b[1]) || (v[1] > a[1]){
        return false;
    }
    
    return true;
}


// class Grid2 : public Box2
// {
// protected :
//     int nx, ny;
//     vec2 diagonal;
//     vec2 celldiagonal;
// public:

// }