// include "vec3.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>


#ifndef VEC2
#define VEC2

class vec2
{
protected:
    double x, y;

public:
    //vec2() : x(0.0), y(0.0) {}
    explicit vec2(double x = 0.0, double y = 0.0) : x(x), y(y) {}
    inline vec2 operator+(const vec2 &p) { return vec2(x + p.x, y + p.y); }
    inline vec2 operator-(const vec2 &p) { return vec2(x - p.x, y - p.y); }
    inline vec2 operator-() { return vec2(-x, -y); }

    inline double operator*(const vec2 &p) { return x * p.x + y * p.y; }



    inline friend vec2 operator*(double a, const vec2& p) { 
        return vec2(a * p.x, a * p.y); 
    }
    
    inline double Length() const { return sqrt(x * x + y * y); }
    inline vec2 Normalized() const
    {
        return (1.0/ Length())* (*this);
    }


    inline vec2 Scale(const vec2 &p) { return vec2(x * p.x, y * p.y); }

    inline double &operator[](int index)
    {
        if ((index > 2) || (index < 0))
        {
            std::cout << "Index out of range for vec2 type" << std::endl;
            return x;
        }
        else if (index == 0)
        {
            return x;
        }
        else
        {
            return y;
        }
    }

    inline const double operator[](int index) const
    {
        if ((index > 2) || (index < 0))
        {
            std::cout << "Index out of range for vec2 type" << std::endl;
            return x;
        }
        else if (index == 0)
        {
            return x;
        }
        else
        {
            return y;
        }
    }

    inline friend  std::ostream& operator<<(std::ostream& os, vec2 p){
        os << p.x << ' ' << p.y;
        return os;
    }

    inline vec2 round(){ return vec2(std::round(x), std::round(y));}

    inline vec2 operator<(vec2 p){ return vec2(x < p.x, y < p.y);}
    inline vec2 operator>(vec2 p){ return vec2(x > p.x, y > p.y);}

    friend class vec3;
};
#endif