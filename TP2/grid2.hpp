#include "box2.hpp"
#include <stdlib.h>
#include <stdio.h>


#ifndef GRID2
#define GRID2

class Grid2 : public Box2
{
protected:
    int nx, ny;
    vec2 diagonal;
    vec2 celldiagonal;
    vec2 inversecelldiagonal;

public:
    Grid2(const Box2 &box = Box2::Empty(), int nx = 0, int ny = 0) : Box2(box), nx(nx), ny(ny)
    {
        diagonal = b - a;
        celldiagonal = diagonal.Scale(vec2(1.0 / double(nx - 1), 1.0 / double(ny - 1)));
        inversecelldiagonal = vec2(1.0 / celldiagonal[0], 1.0 / celldiagonal[1]);
    }
    bool Inside(int i, int j) const
    {
        return ((i >= 0) && (i < nx) && (j >= 0) && (j < ny));
    }
    bool Border(int i, int j) const
    {
        return ((i == 0) || (i == nx - 1) || (j == 0) || (j == ny - 1));
    }

    int Index(int i, int j) const { return i + j * nx; }

    vec2 Vertex(int i, int j) const
    {
        double u = double(i) / (nx - 1);
        double v = double(j) / (ny - 1);
        return vec2((1 - u) * a[0] + u * b[0], (1 - v) * a[1] + v * b[1]);
    }
};

#endif