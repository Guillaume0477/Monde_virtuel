#include "vec2.hpp"
#include <stdlib.h>
#include <stdio.h>


#ifndef BOX2
#define BOX2

class Box2
{
protected:
    vec2 a, b;

public:
    Box2(const vec2 &a, const vec2 &b) : a(a), b(b){}

    bool Inside(const vec2 &v) const
    {
        if ((v[0] > b[0]) || (v[0] < a[0]) || (v[1] > b[1]) || (v[1] < a[1]))
        {
            return false;
        }
        return true;
    }

    bool Intersect(const Box2 &box) const
    {
        if ((a[0] >= box.b[0]) || (a[1] >= box.b[1]) || (b[0] <= box.a[0]) || (b[1] <= box.a[1]))
        {
            return false;
        }
        else
        {
            return true;
        }
    }

public:
    static const Box2 Empty(){ return Box2(vec2(0.0, 0.0), vec2(0.0, 0.0));}

};

#endif