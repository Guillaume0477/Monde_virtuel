#include "vec2.hpp"
#include <stdlib.h>
#include <stdio.h>


#ifndef VEC3
#define VEC3


class vec3
{
protected:
    double x, y, z;

public:
    explicit vec3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}
    vec3(const vec2 &v, double z) : x(v.x), y(v.y), z(z) {}
    inline vec3 operator+(const vec3 &p) { return vec3(x + p.x, y + p.y, z + p.z); }
    inline vec3 operator-(const vec3 &p) { return vec3(x - p.x, y - p.y, z - p.z); }
    inline double operator*(const vec3 &p) { return x * p.x + y * p.y + z * p.z; }


    //inline friend vec3 operator/(double a, const vec3 &p) { return vec3(p.x/ a ,p.y/a , p.z/a); }
    inline friend vec3 operator*(double a, const vec3 &p) { return vec3(a * p.x, a * p.y, a * p.z); }

    inline double Length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    inline vec3 Normalized() const
    {
        return (1.0/ Length())* (*this);
    }

    inline vec3 Scale(const vec3 &p) { return vec3(x * p.x, y * p.y, y * p.z); }
    inline const double operator[](int i) const
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    }

    inline double &operator[](int i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    }

    inline friend  std::ostream& operator<<(std::ostream& os, vec3 p){
        os << p.x << ' ' << p.y << ' ' << p.z;
        return os;
    }
};

#endif