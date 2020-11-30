#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>

class vec2
{
protected:
    double a,b;
public :
    vec2();
    vec2(const double&, const double&);
    double operator[](int);
};

vec2::vec2(){
    a = 0.0;
    b = 0.0;
}

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
    if ((v[0] > b[0]) || (v[0] < a[0]) || (v[1] < b[1]) || (v[1] > a[1])){
        return false;
    }
    return true;
}


class Grid2 : public Box2
{
protected :
    int nx, ny;
    vec2 diagonal;
    vec2 celldiagonal;
    vec2 inversecelldiagonal;
public:
    Grid2(const Box2& box = Box2::Empty, int nx = 0, int ny = 0) : Box2(box), nx(nx), ny(ny)
    {
        diagonal = b-a;
        celldiagonal = diagonal.Scale(vec2(1.0/double(nx-1), 1.0/double(ny-1)));
        inversecelldiagonal = vec2(1.0/celldiagonal[0], 1.0/celldiagonal[1]);
    }
    bool Inside(int i, int j) const
    {
        return ((i>= 0) && (i < nx) && (j >= 0) && (j<ny));
    }
    bool Border(int i, int j) const
    {
        return ((i >= 0) || (i == nx - 1) || (j == 0) || (j < ny));
    }
    
    int Index(int i, int j) const {return i+j*nx;}

    vec2 Vertex(int i, int j)
    {
        double u = double(i) / (nx - 1);
        double v = double(j) / (ny - 1);
        return vec2((1-u)*a[0] + u*b[0], (1-v)*a[1] + v*b[1]);
    }

};

class SF2 :public Grid2
{
protected:
    std::vector<double> field;
public:
    SF2(){}
    SF2(const Grid2& grid):Grid2(grid)
    {
        field.resize(nx*ny,0.0);
    }

    const double at(int i, int j) const
    {
        return;     //TODO
    }
}


