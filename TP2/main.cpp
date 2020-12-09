
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
//#include <QImage>

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

    friend class vec3;
};

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
};

class Box2
{
protected:
    vec2 a, b;

public:
    Box2(const vec2 &, const vec2 &);
    // Box2(const double&, const double&);
    bool Inside(const vec2 &) const;
    bool Intersect(const Box2 &) const;

public:
    static const Box2 Empty;
};

// Box2::Box2(const vec2& u, const vec2& v){
//     a = u;
//     b = v;
// }

// Box2::Box2(const double& x, const double& y){
//     a = -vec2(x/2.0, y/2.0);
//     b = vec2(x/2.0, y/2.0);
// }

bool Box2::Inside(const vec2 &v) const
{
    if ((v[0] > b[0]) || (v[0] < a[0]) || (v[1] > b[1]) || (v[1] < a[1]))
    {
        return false;
    }
    return true;
}

bool Box2::Intersect(const Box2 &box) const
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

const Box2 Box2::Empty(vec2(0.0, 0.0), vec2(0.0, 0.0));

Box2::Box2(const vec2 &a, const vec2 &b) : a(a), b(b)
{
}

class Grid2 : public Box2
{
protected:
    int nx, ny;
    vec2 diagonal;
    vec2 celldiagonal;
    vec2 inversecelldiagonal;

public:
    Grid2(const Box2 &box = Box2::Empty, int nx = 0, int ny = 0) : Box2(box), nx(nx), ny(ny)
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

class SF2 : public Grid2
{
protected:
    std::vector<double> field;

public:
    SF2() {}
    SF2(const Grid2 &grid) : Grid2(grid)
    {
        field.resize(nx * ny, 0.0);
    }

    vec2 Gradient(int, int) const;
    double Laplacian(int, int) const;

    const double at(int i, int j) const
    {
        return field[Index(i, j)];
    }

    double &at(int i, int j)
    {
        return field[Index(i, j)];
    }
};

vec2 SF2::Gradient(int i, int j) const // df/dx,df/dy ~ ( (f(x+e,y)-f(x-e,y))/2e , ... )
{
    vec2 n;

    // Gradient along x axis
    if (i == 0)
    {
        n[0] = (at(i + 1, j) - at(i, j)) * inversecelldiagonal[0];
    }
    else if (i == nx - 1)
    {
        n[0] = (at(i, j) - at(i - 1, j)) * inversecelldiagonal[0];
    }
    else
    {
        n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inversecelldiagonal[0];
    }

    // Gradient along y axis
    if (j == 0)
    {
        n[1] = (at(i, j + 1) - at(i, j)) * inversecelldiagonal[1];
    }
    else if (j == ny - 1)
    {
        n[1] = (at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
    }
    else
    {
        n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inversecelldiagonal[1];
    }
    return n;
}

double SF2::Laplacian(int i, int j) const //d2f / dx2 ~ (f(x+e)-2f(x)+f(x+e))/(e^2)
{
    double laplacian = 0.0;

    // Divergence along x axis
    if (i == 0)
    {
        laplacian += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (celldiagonal[0] * celldiagonal[0]);
    }
    else if (i == nx - 1)
    {
        laplacian += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (celldiagonal[0] * celldiagonal[0]);
    }
    else
    {
        laplacian += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (celldiagonal[0] * celldiagonal[0]);
    }

    // Divergence along y axis
    if (j == 0)
    {
        laplacian += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (celldiagonal[1] * celldiagonal[1]);
    }
    else if (j == ny - 1)
    {
        laplacian += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (celldiagonal[1] * celldiagonal[1]);
    }
    else
    {
        laplacian += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (celldiagonal[1] * celldiagonal[1]);
    }
    return laplacian;
}

class HeighField : public SF2
{
protected:
public:
    HeighField(const SF2 &s) : SF2(s) {}
    //HeighField(const QImage &image, const Box2 &, double, double) /*:Grid2(box, image.width(),image.height()){*/
    //{
        //et interpoler entre les deux z
    //}

    double Height(int i, int j) const { return at(i, j); } // Nouveau nom
    double Slope(int i, int j) const
    {
        vec2 g = Gradient(i, j);
        return sqrt(g * g);
    }

    double AverageSlope(int i, int j) const { return 0.0; }; //TODO

    vec3 Vertex(int i, int j) const { return vec3(Grid2::Vertex(i, j), Height(i, j)); }
    vec3 Normal(int i, int j) const { return vec3(-Gradient(i, j), 1.0).Normalized(); }
};

class LayeredField : public Grid2
{
protected:
    SF2 bedrock;
    SF2 sand;

public:
    LayeredField(const SF2 &bedrock) : Grid2(bedrock), bedrock(bedrock) {}
    LayeredField(const SF2 &bedrock, const SF2 &sand) : Grid2(bedrock), bedrock(bedrock), sand(sand) {}
};




int main (int argc, char *argv[]){
    return 0;
}




/*

QImage Export() const
{
    QImage image(nx,ny, QImage::Format_ARGB32);

    const vec3 lightdir = vec3(2.0, 1.0, 4.0).Normalized;

    for (int i=0, i <nx; i++)
    {
        for (int i=0, i <ny; j++)
        {
            vec3 = Normal(i,j);
            double d =n*lightdir;
            d=(1.0+d)/2.0;
            d *= d;
            image.setPixel(i,j,...)
        }
    }
}

*/

// SF2 SlopeMap() const
// {
//     SF2 n(Grid2(Box2(a, b), nx, ny));
//     for (int i = 0, i < nx; i++)
//     {
//         for (int j = 0, j < ny; j++)
//         {
//             n.at(i, j) = Slope(i, j);
//         }
//     }
// }

// SF2 Drainage() const;

/* webex ?
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; ret.setX(i * (bound_max.x() - bound_min.x())); ret.setX(j * (bound_max.y() - bound_min.y())); ret.setZ(Height(i, j)); return ret; }15:44
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; ret.setX(((nx - i) * bound_max.x() - i * bound_min.x())); ret.setX(((ny - j) * bound_max.y() - j * bound_min.y())); ret.setZ(Height(i, j)); return ret; }15:46
l'inverse sur l'interpolation pardon15:47
attendez je refais15:47
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; double u = double(i) / nx; ret.setX((u * bound_max.x() - (1 - u) * bound_min.x())); u = double(j) / ny; ret.setY((u * bound_max.x() - (1 - u) * bound_min.x())); ret.setZ(Height(i, j)); return ret; }

*/

/* from webex

ugez que ça peut complexifier la tâche c'est pas grave13:52
merci pour votre réponse :)13:52
de eric galin à Tous mes contacts:
Vec2 Gradient(int i, int j) const // df/dx,df/dy ~ ( (f(x+e,y)-f(x-e,y))/2e , ... ) { Vec2 n; // Gradient along x axis if (i == 0) { n[0] = (at(i + 1, j) - at(i, j)) * inversecelldiagonal[0]; } else if (i == nx - 1) { n[0] = (at(i, j) - at(i - 1, j)) * inversecelldiagonal[0]; } else { n[0] = (at(i + 1, j) - at(i - 1, j)) *0.5 * inversecelldiagonal[0]; } // Gradient along y axis if (j == 0) { n[1] = (at(i, j + 1) - at(i, j)) * inversecelldiagonal[1]; } else if (j == ny - 1) { n[1] = (at(i, j) - at(i, j - 1)) * inversecelldiagonal[1]; } else { n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inversecelldiagonal[1]; } return n; }14:53
de Romain Lopez-Rostain à Tous mes contacts:
ça va c'est bon !14:54
de eric galin à Tous mes contacts:
d2f / dx2 ~ (f(x+e)-2f(x)+f(x+e))/(e^2) { double laplacian = 0.0; // Divergence along x axis if (i == 0) { laplacian += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (celldiagonal[0] * celldiagonal[0]); } else if (i == nx - 1) { laplacian += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (celldiagonal[0] * celldiagonal[0]); } else { laplacian += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (celldiagonal[0] * celldiagonal[0]); } // Divergence along y axis if (j == 0) { laplacian += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (celldiagonal[1] * celldiagonal[1]); } else if (j == ny - 1) { laplacian += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (celldiagonal[1] * celldiagonal[1]); } else { laplacian += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (celldiagonal[1] * celldiagonal[1]); } return laplacian; }14:59
de Thomas à Tous mes contacts:
non ça va15:06
de Clément MARTIN à Tous mes contacts:
C'est bon pour moi15:06
de Romain Lopez-Rostain à Tous mes contacts:
Comment on a la direction de la pente si on utilise des valeurs absolues du coup ?15:13
ok15:14
de eric galin à Tous mes contacts:
double Height(int i, int j) const { return at(i, j); } // Nouveau nom double Slope(int i, int j) const { Vec2 g = Gradient(i, j); return sqrt(g * g); } // * == dot15:18
Les deux derniers dans classe HF15:18
Reprise !15:37
de tristan guillard à Tous mes contacts:
Vec3 p(int i, int j) {return Vec3(i * )}15:42
mince on peut pas écrire du code ici directement15:42
de eric galin à Tous mes contacts:
copy/paste15:42
de tristan guillard à Tous mes contacts:
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; ret.setX(i * (bound_max.x() - bound_min.x())); ret.setX(j * (bound_max.y() - bound_min.y())); ret.setZ(Height(i, j)); return ret; }15:44
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; ret.setX(((nx - i) * bound_max.x() - i * bound_min.x())); ret.setX(((ny - j) * bound_max.y() - j * bound_min.y())); ret.setZ(Height(i, j)); return ret; }15:46
l'inverse sur l'interpolation pardon15:47
attendez je refais15:47
QVector3D HeightField::p_3d(int i, int j) const { QVector3D ret; double u = double(i) / nx; ret.setX((u * bound_max.x() - (1 - u) * bound_min.x())); u = double(j) / ny; ret.setY((u * bound_max.x() - (1 - u) * bound_min.x())); ret.setZ(Height(i, j)); return ret; }15:49
de eric galin à Tous mes contacts:
Visiter le site : http://www.reliefshading.com/16:08
Encore mieux : http://www.shadedrelief.com/16:09
Normal -> Gradient16:12
Laplacien pour amplifier les détails16:13


Vec3 Vertex(int i, int j) const { return Vec3(Grid2::Vertex(i, j), Height(i, j)); } Vec3 Normal(int i, int j) const { return Vec3(-Gradient(i, j), 1.0).Normalized(); }16:14


Esquisse du shading :16:14
QImage Export() const { QImage image(nx, ny, QImage::Format_ARGB32); // Relief shading: Vec3 lightdir const Vec3 lightdir=Vec3(2.0, 1.0, 4.0).Normalized(); for (int i = 0; i < nx; i++) { for (int j = 0; j < ny; j++) { // Phong a + d ou d=lightdir*normal Vec3 n = Normal(i, j); double d = n * lightdir; // d in [-1,1] d = (1.0 + d) / 2.0; // d in [0,1] d *= d; // Altérer la couleur en fonction de l'altitude // or Height(i,j) est l'altitude image.setPixel(i, j, ....); } } return image; }16:1



*/

/* from webex le 2/12

de Casanova Delfeil à Tous mes contacts:
J'ai du mal à comprendre comment faire les variations de couleurs en fonctions des valeurs de gradient | laplacien...

de eric galin à Tous mes contacts:
C'est une fonction de transfert
Soit f une fonction dont l'image va dans [0,1]16:35
On peut faire la couleur par introduction d'une fonction de transfert (de type palette) qui génère la couleur à partir du résultat
Exemple : Col Shade(double t) { return Col(t,0,1.0-t); }16:37
Col passe du bleu au rouge (si me trompe pas) en fonction de t...
*/
