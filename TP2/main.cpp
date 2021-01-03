
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <QImage>
#include <QVector>
#include <algorithm>
#include <fstream>
#include <map>
#include <random>

#include "LayeredField.hpp"
#include "Buisson.hpp"
#include "Sapin.hpp"
#include "utils.hpp"

/******************************************
*               Classe vec2               *
******************************************/
// class vec2
// {
// protected:
//     double x, y;

// public:
//     //vec2() : x(0.0), y(0.0) {}
//     explicit vec2(double x = 0.0, double y = 0.0) : x(x), y(y) {}
//     inline vec2 operator+(const vec2 &p) { return vec2(x + p.x, y + p.y); }
//     inline vec2 operator-(const vec2 &p) { return vec2(x - p.x, y - p.y); }
//     inline vec2 operator-() { return vec2(-x, -y); }

//     inline double operator*(const vec2 &p) { return x * p.x + y * p.y; }



//     inline friend vec2 operator*(double a, const vec2& p) { 
//         return vec2(a * p.x, a * p.y); 
//     }
    
//     inline double Length() const { return sqrt(x * x + y * y); }
//     inline vec2 Normalized() const
//     {
//         return (1.0/ Length())* (*this);
//     }


//     inline vec2 Scale(const vec2 &p) { return vec2(x * p.x, y * p.y); }

//     inline double &operator[](int index)
//     {
//         if ((index > 2) || (index < 0))
//         {
//             std::cout << "Index out of range for vec2 type" << std::endl;
//             return x;
//         }
//         else if (index == 0)
//         {
//             return x;
//         }
//         else
//         {
//             return y;
//         }
//     }

//     inline const double operator[](int index) const
//     {
//         if ((index > 2) || (index < 0))
//         {
//             std::cout << "Index out of range for vec2 type" << std::endl;
//             return x;
//         }
//         else if (index == 0)
//         {
//             return x;
//         }
//         else
//         {
//             return y;
//         }
//     }

//     inline friend  std::ostream& operator<<(std::ostream& os, vec2 p){
//         os << p.x << ' ' << p.y;
//         return os;
//     }

//     inline vec2 round(){ return vec2(std::round(x), std::round(y));}

//     inline vec2 operator<(vec2 p){ return vec2(x < p.x, y < p.y);}
//     inline vec2 operator>(vec2 p){ return vec2(x > p.x, y > p.y);}

//     friend class vec3;
// };


/******************************************
*               Classe vec3               *
******************************************/
// class vec3
// {
// protected:
//     double x, y, z;

// public:
//     explicit vec3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}
//     vec3(const vec2 &v, double z) : x(v.x), y(v.y), z(z) {}
//     inline vec3 operator+(const vec3 &p) { return vec3(x + p.x, y + p.y, z + p.z); }
//     inline vec3 operator-(const vec3 &p) { return vec3(x - p.x, y - p.y, z - p.z); }
//     inline double operator*(const vec3 &p) { return x * p.x + y * p.y + z * p.z; }


//     //inline friend vec3 operator/(double a, const vec3 &p) { return vec3(p.x/ a ,p.y/a , p.z/a); }
//     inline friend vec3 operator*(double a, const vec3 &p) { return vec3(a * p.x, a * p.y, a * p.z); }

//     inline double Length() const
//     {
//         return sqrt(x * x + y * y + z * z);
//     }

//     inline vec3 Normalized() const
//     {
//         return (1.0/ Length())* (*this);
//     }

//     inline vec3 Scale(const vec3 &p) { return vec3(x * p.x, y * p.y, y * p.z); }
//     inline const double operator[](int i) const
//     {
//         if (i == 0)
//             return x;
//         else if (i == 1)
//             return y;
//         else
//             return z;
//     }

//     inline double &operator[](int i)
//     {
//         if (i == 0)
//             return x;
//         else if (i == 1)
//             return y;
//         else
//             return z;
//     }

//     inline friend  std::ostream& operator<<(std::ostream& os, vec3 p){
//         os << p.x << ' ' << p.y << ' ' << p.z;
//         return os;
//     }
// };



/******************************************
*           Classe ScalarPoint2           *
******************************************/

class VectorField2;

// class ScalarPoint2
// {
// protected:
//     QPoint p;
//     double z;
// public:
//     ScalarPoint2() {
//         p = QPoint(0.0,0.0);
//         z=0.0;
//     }
//     ScalarPoint2(const QPoint& p1, const double& z1){
//         p = p1; 
//         z= z1;
//     };
//     friend bool operator<(const ScalarPoint2& a, const ScalarPoint2& b) {
//         return a.z < b.z;
//     };
//     QPoint Point() const {return p;};
//     double Scalar() const {return z;};

// };





/******************************************
*               Classe Box2               *
******************************************/
// class Box2
// {
// protected:
//     vec2 a, b;

// public:
//     Box2(const vec2 &, const vec2 &);
//     // Box2(const double&, const double&);
//     bool Inside(const vec2 &) const;
//     bool Intersect(const Box2 &) const;

// public:
//     static const Box2 Empty;
// };

// bool Box2::Inside(const vec2 &v) const
// {
//     if ((v[0] > b[0]) || (v[0] < a[0]) || (v[1] > b[1]) || (v[1] < a[1]))
//     {
//         return false;
//     }
//     return true;
// }

// bool Box2::Intersect(const Box2 &box) const
// {
//     if ((a[0] >= box.b[0]) || (a[1] >= box.b[1]) || (b[0] <= box.a[0]) || (b[1] <= box.a[1]))
//     {
//         return false;
//     }
//     else
//     {
//         return true;
//     }
// }

// const Box2 Box2::Empty(vec2(0.0, 0.0), vec2(0.0, 0.0));

// Box2::Box2(const vec2 &a, const vec2 &b) : a(a), b(b)
// {
// }


/******************************************
*               Classe Grid2              *
******************************************/

// class Grid2 : public Box2
// {
// protected:
//     int nx, ny;
//     vec2 diagonal;
//     vec2 celldiagonal;
//     vec2 inversecelldiagonal;

// public:
//     Grid2(const Box2 &box = Box2::Empty, int nx = 0, int ny = 0) : Box2(box), nx(nx), ny(ny)
//     {
//         diagonal = b - a;
//         celldiagonal = diagonal.Scale(vec2(1.0 / double(nx - 1), 1.0 / double(ny - 1)));
//         inversecelldiagonal = vec2(1.0 / celldiagonal[0], 1.0 / celldiagonal[1]);
//     }
//     bool Inside(int i, int j) const
//     {
//         return ((i >= 0) && (i < nx) && (j >= 0) && (j < ny));
//     }
//     bool Border(int i, int j) const
//     {
//         return ((i == 0) || (i == nx - 1) || (j == 0) || (j == ny - 1));
//     }

//     int Index(int i, int j) const { return i + j * nx; }

//     // const QPoint next[8] = { QPoint (1,0) ,QPoint (1,1) ,QPoint (0,1) ,QPoint (-1,1) ,QPoint (0,-1) ,QPoint (-1,-1) ,QPoint (-1,0), QPoint (1,-1) };
//     // const double length[8] = {1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0)};

//     vec2 Vertex(int i, int j) const
//     {
//         double u = double(i) / (nx - 1);
//         double v = double(j) / (ny - 1);
//         return vec2((1 - u) * a[0] + u * b[0], (1 - v) * a[1] + v * b[1]);
//     }
// };


/******************************************
*               Classe SF2                *
******************************************/
// class SF2 : public Grid2
// {
// protected:
//     std::vector<double> field;
//     double maxVal = -1000000;
//     double minVal = 1000000;
// public:
//     SF2() {}
//     SF2(const Grid2 &grid) : Grid2(grid)
//     {
//         field.resize(nx * ny, 0.0);
//     }
//     SF2(const Grid2 &grid, const double &hz) : Grid2(grid)
//     {
//         field.resize(nx * ny, hz);
//     }

//     vec2 Gradient(int, int) const;
//     double Laplacian(int, int) const;

//     const double at(int i, int j) const
//     {
//         return field[Index(i, j)];
//     }

//     double &at(int i, int j)
//     {
//         return field[Index(i, j)];
//     }

//     const double operator[](int i) const
//     {
//         return field[i];
//     }

//     double &operator[] (int i)
//     {
//         return field[i];
//     }

//     void UpdateMinMax(){
//         for (int k = 0; k < field.size(); k++){
//             double val = field[k];
//             if (val > maxVal){
//                 maxVal = val;
//             }
//             if (val < minVal){
//                 minVal = val;
//             }
//         }
//     }

//     const double max() const {return maxVal;}
//     const double min() const {return minVal;}

//     double getConv(int, int, float, float, float);

//     SF2 GradientNorm() const ;
//     SF2 LaplacianMap();

//     void Smooth();
//     void Blur();
//     void Normalize();
//     void Clamp(float, float);
//     QVector<ScalarPoint2> GetScalarPoints() const;
// };

// vec2 SF2::Gradient(int i, int j) const // df/dx,df/dy ~ ( (f(x+e,y)-f(x-e,y))/2e , ... )
// {
//     vec2 n;

//     // Gradient along x axis
//     if (i == 0)
//     {
//         n[0] = (at(i + 1, j) - at(i, j)) * inversecelldiagonal[0];
//     }
//     else if (i == nx - 1)
//     {
//         n[0] = (at(i, j) - at(i - 1, j)) * inversecelldiagonal[0];
//     }
//     else
//     {
//         n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inversecelldiagonal[0];
//     }

//     // Gradient along y axis
//     if (j == 0)
//     {
//         n[1] = (at(i, j + 1) - at(i, j)) * inversecelldiagonal[1];
//     }
//     else if (j == ny - 1)
//     {
//         n[1] = (at(i, j) - at(i, j - 1)) * inversecelldiagonal[1];
//     }
//     else
//     {
//         n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inversecelldiagonal[1];
//     }
//     return n;
// }

// double SF2::Laplacian(int i, int j) const //d2f / dx2 ~ (f(x+e)-2f(x)+f(x+e))/(e^2)
// {
//     double laplacian = 0.0;

//     // Divergence along x axis
//     if (i == 0)
//     {
//         laplacian += (at(i, j) - 2.0 * at(i + 1, j) + at(i + 2, j)) / (celldiagonal[0] * celldiagonal[0]);
//     }
//     else if (i == nx - 1)
//     {
//         laplacian += (at(i, j) - 2.0 * at(i - 1, j) + at(i - 2, j)) / (celldiagonal[0] * celldiagonal[0]);
//     }
//     else
//     {
//         laplacian += (at(i + 1, j) - 2.0 * at(i, j) + at(i - 1, j)) / (celldiagonal[0] * celldiagonal[0]);
//     }

//     // Divergence along y axis
//     if (j == 0)
//     {
//         laplacian += (at(i, j) - 2.0 * at(i, j + 1) + at(i, j + 2)) / (celldiagonal[1] * celldiagonal[1]);
//     }
//     else if (j == ny - 1)
//     {
//         laplacian += (at(i, j) - 2.0 * at(i, j - 1) + at(i, j - 2)) / (celldiagonal[1] * celldiagonal[1]);
//     }
//     else
//     {
//         laplacian += (at(i, j + 1) - 2.0 * at(i, j) + at(i, j - 1)) / (celldiagonal[1] * celldiagonal[1]);
//     }
//     return laplacian;
// }

// double SF2::getConv(int i, int j, float mid, float side, float diag){
//     float nb = 0;
//     float value = 0;

//     //top left
//     if ((i > 0)&&(j>0)){
//         value += diag*field[Index(i-1, j-1)];
//         nb+= diag;
//     }
//     //top
//     if (i > 0){
//         value += side*field[Index(i-1, j)];
//         nb+= side;
//     }
//     //top right
//     if ((i > 0)&&(j < ny-1)){
//         value += diag*field[Index(i-1, j+1)];
//         nb+=diag;
//     }
//     //left
//     if (j>0){
//         value += side*field[Index(i, j-1)];
//         nb+=side;
//     }
//     //right
//     if (j<ny-1){
//         value += side*field[Index(i, j+1)];
//         nb+=side;
//     }
//     //bottom left
//     if ((i < nx-1) && (j>0)){
//         value += diag*field[Index(i+1, j-1)];
//         nb+=diag;
//     }
//     //bottom
//     if (i<nx-1){
//         value += side*field[Index(i+1, j)];
//         nb+=side;
//     }
//     //bottom right
//     if ((i < nx-1) && (j<ny-1)){
//         value += diag*field[Index(i+1, j+1)];
//         nb+=diag;
//     }
//     //center
//     value += mid*field[Index(i,j)];
//     nb+=mid;

//     //normalization
//     value /= nb;

//     return value;
// }

// SF2 SF2::GradientNorm() const{
//     SF2 gradNorm = SF2(Grid2(*this));

//     for (int i = 0; i < nx; i++){
//         for (int j = 0; j < ny; j++){
//             vec2 grad = Gradient(i,j);
//             gradNorm.at(i,j) = sqrt(grad*grad);
//         }
//     }

//     return gradNorm;
// }

// SF2 SF2::LaplacianMap(){
//     SF2 LaplMap = SF2(Grid2(*this));

//     for (int i = 0; i < nx; i++){
//         for (int j = 0; j < ny; j++){
//             LaplMap.at(i,j) = Laplacian(i,j);
//         }
//     }

//     return LaplMap;
// }


// void SF2::Smooth(){
//     std::vector<double> smoothed;
//     smoothed.resize(nx*ny);

//     for (int i = 0; i < nx ; i ++){
//         for (int j = 0; j < ny ; j ++){
//             int ind = Index(i,j);
            
//             smoothed[ind] = getConv(i,j,4,2,1);
//         }
//     }

//     field = smoothed;
// }

// void SF2::Blur(){
//     std::vector<double> blured;
//     blured.resize(nx*ny);

//     for (int i = 0; i < nx ; i ++){
//         for (int j = 0; j < ny ; j ++){
//             int ind = Index(i,j);
            
//             blured[ind] = getConv(i,j,1,1,1);
//         }
//     }

//     field = blured;
// }

// void SF2::Normalize(){
//     UpdateMinMax();

//     for (int k = 0; k < field.size(); k++){
//         field[k] = (field[k] - minVal)/(maxVal - minVal);
//     }

// }

// void SF2::Clamp(float mini, float maxi){
    
//     for (int k = 0; k < field.size(); k++){
//         if (field[k] < mini){
//             field[k] = mini;
//         } else if (field[k] > maxi){
//             field[k] = maxi;
//         }
//     }
// };


// QVector<ScalarPoint2> SF2::GetScalarPoints() const
// {
//     QVector<ScalarPoint2> e(nx * ny);
//     int k=0;
//     for (int i=0; i < nx ; i++){
//         for (int j=0; j < ny ; j++){
//             e[k++] = ScalarPoint2(QPoint(i,j), at(i,j));
//             //std::cout<<"atttt"<<e.at(i+nx*j).Scalar()<<std::endl;
//         }
//     }
//     return e;
// };



// float fonction_one(float lim_inf, float lim_sup,float current_num){
//     if ((current_num > lim_inf) && (current_num < lim_sup)){
//         return current_num;
//     }
//     else {
//         return 0;
//     }
// }

// float fonction_zeros(float lim_inf, float lim_sup,float current_num){
//     if ((current_num > lim_inf) && (current_num < lim_sup)){
//         return current_num;
//     }
//     else {
//         return 1;
//     }
// }




/******************************************
*              Classe Arbre               *
******************************************/
// class Arbre
// {
// protected:
//     int nombre_attrib;
//     int rayon;
// public:
//     Arbre(){
//         rayon = 3; 
//         nombre_attrib = 3;
//     };
//     Arbre(const int nb, const int r){
//         nombre_attrib= nb;
//         rayon = r; 
//     };
//     virtual double get_rayon() const {
//     }
//     virtual double humidity(double hum){
//         return hum;
//     };
//     virtual double slope(double slo){
//         return slo;
//     };
//     virtual double stream(double str){
//         return str;
//     };

// };



/******************************************
*              Classe Sapin               *
******************************************/
// class Sapin : public Arbre
// {
// protected:

// public:
//     Sapin(){
//         rayon = 10;
//         nombre_attrib = 3;
//     };
//     double get_rayon() const {
//         return rayon;
//     }
//     double humidity(double hum){
//         double value = fonction_one(0.1,1.0,hum);
//         return (value);
//     };
//     double slope(double slo){
//         double value = fonction_one(0.0,0.3,slo);
//         return (value);
//     };
//     double stream(double str){
//         double value = fonction_one(0.0,0.05,str);
//         return (value);
//     };

// };

/******************************************
*             Classe Buisson              *
******************************************/
// class Buisson : public Arbre
// {
// protected:

// public:
//     Buisson(){
//         rayon = 2;
//         nombre_attrib = 3;
//     };
//     double get_rayon() const {
//         return rayon;
//     }
//     double humidity(double hum){
//         double value = fonction_one(0.0,0.1,hum);
//         return (value);
//     };
//     double slope(double slo){
//         double value = fonction_one(0.0,0.3,slo);
//         return (value);
//     };
//     double stream(double str){
//         double value = fonction_one(0.0,0.05,str);
//         return (value);
//     };

// };





/******************************************
*           Classe HeighField2            *
******************************************/
// class HeighField : public SF2
// {
// protected:
// public:
//     HeighField(const SF2 &s) : SF2(s) {}
//     HeighField(const QImage &image, const Box2 &box, double nx, double ny) /*:Grid2(box, image.width(),image.height()){*/
//     {
//         //et interpoler entre les deux z
        
//         //Initialisation
//         Grid2 grid = Grid2(box, nx, ny);
//         *this = SF2(grid);

//         float normFact = 500.0;

//         //Remplissage des hauteurs
//         for (int x = 0; x < nx; x ++){
//             for (int y = 0; y < ny ; y++){
//                 field[Index(x,y)] = ((qRed(image.pixel(x,y)))/255.0 * normFact);
//                 // std::cout << field.back() << std::endl;
//             }
//         }

//         UpdateMinMax();
//         //std::cout << field.size() << ' ' << nx*ny << ' ' << nx << ' ' << ny << std::endl;
//     }

//     double Height(int i, int j) const { return at(i, j); } 
//     double Slope(int i, int j) const
//     {
//         vec2 g = Gradient(i, j);
//         return sqrt(g * g);
//     }

//     bool intersectRay(vec3, double&, vec3);
//     double access(int,int, int N = 20);
//     SF2 accessMap(){    
//         SF2 accMap = SF2(Grid2(*this));

//         for (int i = 0; i < nx; i++){
//             for (int j = 0; j < ny; j++){
//                 accMap.at(i,j) = access(i,j);
//             }
//         }

//         return accMap;
//     }


//     double AverageSlope(int, int) const; 

//     vec3 Vertex(int i, int j) const { return vec3(Grid2::Vertex(i, j), Height(i, j)); }
//     vec3 Normal(int i, int j) const { return vec3(-Gradient(i, j), 1.0).Normalized(); }

//     SF2 SlopeMap() const {return GradientNorm();}
//     SF2 AVGSlopeMap(){    
//         SF2 avgMap = SF2(Grid2(*this));

//         for (int i = 0; i < nx; i++){
//             for (int j = 0; j < ny; j++){
//                 avgMap.at(i,j) = AverageSlope(i,j);
//             }
//         }
//         return avgMap;
//     }
//     int CheckFlowSlope( const QPoint& p, QPoint* point, double* height, double* slope, double* nslone, int& mask) const;
//     SF2 StreamAreaStreepest() const;
//     SF2 StreamArea() const;
//     SF2 StreamPower() const;
//     SF2 WetNessIndex() const;
//     //SF2 densite_sapin(Arbre arbre) const;
//     SF2 densite_arbre(Arbre& arbre) const; //not working
//     SF2 densite_sapin() const;
//     SF2 densite_buisson() const;
//     SF2 sapin_raw_distribution() const;
//     SF2 raw_dart_throwing(Arbre& arbre) const;
//     SF2 raw_distribution(Arbre& arbre) const;
//     SF2 double_raw_distribution(Arbre& arbre1, Arbre& arbre2) const;


//     //Exportation sous format d'image ! 
//     QImage Export(SF2, bool) const;
//     QImage Shade(SF2 mapToExport) const {return Export(mapToExport, true);}

//     void ExportOBJ(char*);



// };

// QImage HeighField::Export(SF2 mapToExport, bool vis = false) const
// {
//     QImage image(nx,ny, QImage::Format_ARGB32);

//     const vec3 lightdir = vec3(2.0, 1.0, 4.0).Normalized();
//     mapToExport.Normalize();
//     for (int i=0; i <nx; i++)
//     {
//         for (int j=0; j <ny; j++)
//         {
//             int mapVal = (mapToExport.at(i,j)*255);
//             if (vis){
//                 vec3 n = Normal(i,j);
//                 double d =n*lightdir;
//                 d=(1.0+d)/2.0;
//                 d *= d;
//                 mapVal *= d;
//             }

//             image.setPixel(i,j,qRgb(mapVal, mapVal, mapVal));

//         }
//     }

//     return image;

// }

// bool HeighField::intersectRay(vec3 rayDir, double& t, vec3 origin){

//     double epsilon = 1.0;
//     vec3 ray = origin + t*rayDir;

//     while (Inside(std::round(ray[0]), std::round(ray[1])) && (ray[2] < maxVal) && (ray[2]>0)){
//         //get z value by finding the maximum heigh between the 3 closest points
//         vec2 proj = vec2(ray[0], ray[1]);
//         vec2 roundProj = proj.round();


//         //vec2 diff = (proj > roundProj) * 2 - vec2(1.0, 1.0);

//         double height = at(roundProj[0], roundProj[1]);

//         if (height > ray[2]){
//             return true;
//         } else {
//             t += epsilon;
//             ray = origin + t*rayDir;
//         }
//     }

//     return false;
// }

// double HeighField::access(int i, int j, int Nray){
//     double accessVal = 0.0;

//     vec3 origin = vec3(i,j,at(i,j));
//     for (int r = 0; r < Nray ; r++){
//         //Create ray direction randomly on hemisphere (formula (34) GI Compedium)
//         vec3 rayDir;
//         double t = 1;

//         double r1,r2; 
//         r1 = double(rand())/RAND_MAX;
//         r2 = double(rand())/RAND_MAX;
//         //Direction aléatoire dans l'espace tangent
//         vec3 rayDirTan = vec3(cos(2*M_PI*r1) * sqrt(1-r2*r2), sin(2*M_PI*r1)*sqrt(1-r2*r2), r2);

//         vec3 ta, bi, n;
//         n = Normal(i,j);

//         //Code from J-C.IEHL to transfer a vector from tangent space to world space
//         float sign= n[2] < 0 ? -1 : 1;             
//         float a= -1.0 / (sign + n[2]);
//         float d= n[0] * n[1] * a;
//         ta= vec3(1.0 + sign * n[0] * n[0] * a, sign * d, -sign * n[0]);
//         bi= vec3(d, sign + n[1] * n[1] * a, -n[1]);

//         rayDir = rayDirTan[0] * ta.Normalized() + rayDirTan[1] * bi.Normalized() + rayDirTan[2] * n.Normalized();
//         //End of the code inspired from J-C.IEHL

//         //Check the visibility
//         if (intersectRay(rayDir.Normalized(), t, origin) == false){
//             accessVal++;
//         }

//     }

//     return accessVal/double(Nray);
// }




// double HeighField::AverageSlope(int i, int j) const{
//     double avgSlope = 0.0;  
//     double slopeVal;  
//     float nb = 0;

//     float curHeight = at(i,j);

//     float diagVal = sqrt(celldiagonal[0]*celldiagonal[0] + celldiagonal[1]*celldiagonal[1]);

//     //top left
//     if ((i > 0)&&(j>0)){
//         slopeVal =  (curHeight - at(i-1,j-1)) / diagVal;
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //top
//     if (i > 0){
//         slopeVal= (curHeight - at(i-1, j)) * inversecelldiagonal[0];
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //top right
//     if ((i > 0)&&(j < ny-1)){
//         slopeVal = (curHeight - at(i-1,j+1)) / diagVal;
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //left
//     if (j>0){
//         slopeVal= (curHeight - at(i, j-1)) * inversecelldiagonal[1];
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //right
//     if (j<ny-1){
//         slopeVal= (curHeight - at(i, j+1)) * inversecelldiagonal[1];
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //bottom left
//     if ((i < nx-1) && (j>0)){
//         slopeVal = (curHeight - at(i+1,j-1)) / diagVal;
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //bottom
//     if (i<nx-1){
//         slopeVal = (curHeight - at(i+1, j)) * inversecelldiagonal[0];
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
//     //bottom right
//     if ((i < nx-1) && (j<ny-1)){
//         slopeVal = (curHeight - at(i+1,j+1)) / diagVal;
//         avgSlope += sqrt(slopeVal*slopeVal);
//         nb++;
//     }
    
//     //normalization
//     avgSlope /= nb;

//     return avgSlope;
// }

// void HeighField::ExportOBJ(char* filename){
//     std::ofstream objFile;
//     objFile.open(filename);

//     //Write vertices
//     for (int i = 0; i < nx; i ++){
//         for (int j = 0; j < ny ; j++){
//             vec3 vertex = Vertex(i,j);
//             objFile << "v " << vertex << std::endl;
//         }
//     }
    
//     objFile << std::endl;

//     //Write normals
//     for (int i = 0; i < nx; i ++){
//         for (int j = 0; j < ny ; j++){
//             vec3 normal = Normal(i,j);
//             objFile << "vn " << normal << std::endl;
//         }
//     }
    
//     objFile << std::endl;
    
//     //Write texture coords
//     for (int i = 0; i < nx; i ++){
//         for (int j = 0; j < ny ; j++){
//             vec2 texCoord = vec2(float(i)/float(nx), float(j)/float(ny));
//             objFile << "vt " << texCoord << std::endl;
//         }
//     }

//     objFile << std::endl;

//     //Write faces
//     for (int i = 0; i < nx-1; i ++){
//         for (int j = 0; j < ny-1 ; j++){
//             int ind = Index(i,j);
//             int indr = Index(i, j+1);
//             int indb = Index(i+1, j);
//             int indd = Index(i+1, j+1);
//             objFile << "f " << ind << '/' << ind << '/' << ind << ' '
//             << indr << '/' << indr << '/' << indr << ' ' 
//             << indb << '/' << indb << '/' << indb << ' ' << std::endl;

//             objFile << "f " << indr << '/' << indr << '/' << indr << ' '
//             << indb << '/' << indb << '/' << indb << ' ' 
//             << indd << '/' << indd << '/' << indd << ' ' << std::endl;


//         }
//     }    


// }


// int HeighField::CheckFlowSlope( const QPoint& p, QPoint* point, double* height, double* slope, double* nslone, int& mask) const {

//     int n = 0;
//     double zp = at(p.x(),p.y());
//     double slopesum = 0.0;

//     const QPoint next[8] = { QPoint (1,0) ,QPoint (1,1) ,QPoint (0,1) ,QPoint (-1,1) ,QPoint (0,-1) ,QPoint (-1,-1) ,QPoint (-1,0), QPoint (1,-1) };
//     const double length[8] = {1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0)};

//     //std::cout<<"zp "<<zp<<std::endl;

//     mask = 0;
//     for (int i=0; i<8; i++){
//         QPoint b = p + next[i];
//         // std::cout<<"next"<<next[i].x()<<" "<<next[i].y() <<std::endl;
//         //std::cout<<"b "<<b.x()<<" "<<b.y() <<std::endl;

//         double step = at(b.x(), b.y()) - zp;
        
//         if (b.x()>=nx || b.y()>=ny || b.x()<=0.0 || b.y()<=0.0){
//             continue;
//         }
        
//         //if(!Box2::Inside(vec2(b.x(),b.y()))){continue;};


        
//         // std::cout<<"at_b "<<at(b.x(), b.y())<<std::endl;
//         // std::cout<<"step "<<step<<std::endl;

//         if (step <0.0)
//         {
//             point[n] = b;
//             //std::cout<<"point "<<point[n].x()<<" "<<point[n].y() <<std::endl;
//             height[n] = -step;
//             slope[n] = - step/ length[i];
//             slopesum += slope[n];
//             n++;
//             mask |= 1 << i;
//         }

//     }

//     for (int k=0; k<n; k++){
//         nslone[k] = slope[k] / slopesum;
//     }
//     return n;

// }

// SF2 HeighField::StreamAreaStreepest() const{

//     SF2 stream = SF2(Grid2(Box2(a,b),nx,ny),1.0);

//     // SF2 stream = SF2(Grid2(*this));

//     // for (int i = 0; i < nx; i++){
//     //     for (int j = 0; j < ny; j++){
//     //         stream.at(i,j) = 1.0;
//     //     }
//     // }
//     // stream.at(2,2) = 0.01;

//     QVector<ScalarPoint2> QEE = GetScalarPoints();

//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
//     // for (int i = QEE.size() -1; i>=0; i--){
//     //     std::cout<<QEE.at(i).Scalar()<<std::endl;
//     //     //stream.at(q[k].x(), q[k].y()) += sp;
//     // }
//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    
//     // QVector<ScalarPoint2> QEE = QVector(ScalarPoint2(next,length))
//     std::sort(QEE.begin(), QEE.end());

//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
//     // for (int i = QEE.size() -1; i>=0; i--){
//     //     std::cout<<QEE.at(i).Scalar()<<std::endl;
//     //     //stream.at(q[k].x(), q[k].y()) += sp;
//     // }
//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;

//     for (int i = QEE.size() -1; i>=0; i--){
        
//         // QPoint q = QPoint(2,2);
//         // double z = 5;
//         // ScalarPoint2 p = ScalarPoint2(q,z);
//         // p.Point(); 
//         //QEE.at(i).Point();
//         QPoint p = QEE.at(i).Point();


        

//         QPoint q[8];
//         double h[8];
//         double s[8];
//         double sn[8];
//         int m;

//         int n = CheckFlowSlope(p,q,h,s,sn,m);

//         // std::cout<<"p "<<p.x()<<" , "<<p.y()<<std::endl;
//         // for (int j=0; j<8; j++){
//         //     std::cout<<"q "<<j<<q[j].x()<<" "<<q[j].y()<<std::endl;
//         //     std::cout<<"h "<<j<<h[j]<<std::endl;
//         //     std::cout<<"s "<<j<<s[j]<<std::endl;
//         //     std::cout<<"sn "<<j<<sn[j]<<std::endl;
//         // }

//         if (n>0)
//         {
//             double ss = s[0]; //calcul max s
//             int k=0;
//             for (int j=0; j<n;j++)
//             {
//                 if (s[j] > ss) {
//                     k=j;
//                     ss=s[j];
//                 }
//             }
//             const double sp = stream.at(p.x(),p.y());
//             stream.at(q[k].x(), q[k].y()) += sp;
//         }

//     }
//     return stream;
// }




// SF2 HeighField::StreamArea() const{

//     SF2 stream = SF2(Grid2(Box2(a,b),nx,ny),1.0);

//     // SF2 stream = SF2(Grid2(*this));

//     // for (int i = 0; i < nx; i++){
//     //     for (int j = 0; j < ny; j++){
//     //         stream.at(i,j) = 1.0;
//     //     }
//     // }
//     // stream.at(2,2) = 0.01;

//     QVector<ScalarPoint2> QEE = GetScalarPoints();

//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
//     // for (int i = QEE.size() -1; i>=0; i--){
//     //     std::cout<<QEE.at(i).Scalar()<<std::endl;
//     //     //stream.at(q[k].x(), q[k].y()) += sp;
//     // }
//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    
//     // QVector<ScalarPoint2> QEE = QVector(ScalarPoint2(next,length))
//     std::sort(QEE.begin(), QEE.end());

//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
//     // for (int i = QEE.size() -1; i>=0; i--){
//     //     std::cout<<QEE.at(i).Scalar()<<std::endl;
//     //     //stream.at(q[k].x(), q[k].y()) += sp;
//     // }
//     // std::cout<<"QEE.at(i).Scalar()"<<std::endl;

//     for (int i = QEE.size() -1; i>=0; i--){
        
//         // QPoint q = QPoint(2,2);
//         // double z = 5;
//         // ScalarPoint2 p = ScalarPoint2(q,z);
//         // p.Point(); 
//         //QEE.at(i).Point();
//         QPoint p = QEE.at(i).Point();


        

//         QPoint q[8];
//         double h[8];
//         double s[8];
//         double sn[8];
//         int m;

//         int n = CheckFlowSlope(p,q,h,s,sn,m);

//         // std::cout<<"p "<<p.x()<<" , "<<p.y()<<std::endl;
//         // for (int j=0; j<8; j++){
//         //     std::cout<<"q "<<j<<" "<<q[j].x()<<" "<<q[j].y()<<std::endl;
//         //     std::cout<<"h "<<j<<" "<<h[j]<<std::endl;
//         //     std::cout<<"s "<<j<<" "<<s[j]<<std::endl;
//         //     std::cout<<"sn "<<j<<" "<<sn[j]<<std::endl;
//         // }
//         // std::cout<<"N "<<n<<std::endl;
//         if (n>0)
//         {
//             const double sp = stream.at(p.x(),p.y());

//             for (int k=0; k<n;k++)
//             {
//                 stream.at(q[k].x(), q[k].y()) += sp * sn[k];
//                 //std::cout<<"sn "<<k<<" "<<sn[k]<<std::endl;


//             }
            

//         }

//     }
//     return stream;
// };






// SF2 HeighField::StreamPower() const{
//     SF2 stream = StreamArea();
//     SF2 slope = SlopeMap();
//     SF2 res(Grid2(*this));
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             res.at(i, j) = sqrt(stream.at(i, j)) * slope.at(i, j);
//         }
//     }
//     return res;
// }

// SF2 HeighField::WetNessIndex() const{
//     SF2 stream = StreamArea();
//     SF2 slope = SlopeMap();
//     SF2 res(Grid2(*this));
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             res.at(i, j) = log(stream.at(i, j)) * slope.at(i, j);
//         }
//     }
//     return res;
// }



// SF2 HeighField::densite_arbre(Arbre& arbre) const{


//     SF2 stream = StreamArea();
//     stream.Normalize();
//     SF2 humidity = WetNessIndex();
//     humidity.Normalize();
//     SF2 slope = SlopeMap();
//     slope.Normalize();
//     SF2 res = SF2(Grid2(Box2(a,b),nx,ny),1.0); //init 1
//     int rayon = arbre.get_rayon();

//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             for (int k=0; k<1; k++){
//             res.at(i, j) = std::min( arbre.humidity(humidity.at(i,j)) , res.at(i, j) );
//             res.at(i, j) = std::min( arbre.slope(slope.at(i,j)) , res.at(i, j) );
//             res.at(i, j) = std::min( arbre.stream(stream.at(i,j)) , res.at(i, j) );
//             }
//         }
//     }
//     return res;
// }



// bool test_dist(std::pair<int,int> couple1, std::pair<int,int> couple2, float rayon){
//     int diff_x = (couple1.first - couple2.first);
//     int diff_y = (couple1.second - couple2.second);
//     float dist = std::sqrt( diff_x*diff_x + diff_y*diff_y );
//     if (dist < 2*rayon){
//         return true;
//     }
//     else
//     {
//         return false;
//     }
    

// }

// bool test_dist(std::pair< std::pair<int,int> , int> couple1, std::pair< std::pair<int,int> , int> couple2){
//     int diff_x = (couple1.first.first - couple2.first.first);
//     int diff_y = (couple1.first.second - couple2.first.second);
//     float dist = std::sqrt( diff_x*diff_x + diff_y*diff_y );
//     if (dist < (couple1.second + couple2.second)){
//         return true;
//     }
//     else
//     {
//         return false;
//     }
    

// };



// SF2 HeighField::raw_dart_throwing(Arbre& arbre) const{

    
//     SF2 dens_arbre = densite_arbre(arbre);
//     dens_arbre.Normalize();

//     int rayon_arbre = arbre.get_rayon();
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
//     std::cout<<rayon_arbre<<std::endl;
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

//     SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

//     std::list< std::pair< std::pair<int,int> , int > >  list_arbre;

//     for (int k=0; k<100000; k++){
//         int rand_pos_x = rand()%nx;
//         int rand_pos_y = rand()%ny;


//         bool placement_not_possible = false;


//         for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
//             if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
//                 placement_not_possible = true;
//                 break;
//             }
//         }


//         if (placement_not_possible == false){
//             list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
//         }
//     }

//     for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        

//         bool test_dens_sapin = true;


//         // float rand_test = ((double) rand() / (RAND_MAX));
//         int rand_pos_x = (*it).first.first;
//         int rand_pos_y = (*it).first.second;

//         // if (rand_test <= dens_arbre.at(rand_pos_x,rand_pos_y)){
//         //     test_dens_sapin=true;
//         // }
//         // else{
//         //     continue;
//         // }

//         if (test_dens_sapin == true){
//             res.at(rand_pos_x, rand_pos_y) = 1;
//             res.at(rand_pos_x+2, rand_pos_y) = 1;
//             res.at(rand_pos_x-2, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+2) = 1;
//             res.at(rand_pos_x, rand_pos_y-2) = 1;
//             res.at(rand_pos_x+1, rand_pos_y) = 1;
//             res.at(rand_pos_x-1, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+1) = 1;
//             res.at(rand_pos_x, rand_pos_y-1) = 1;
//         }


//     }

//     return res;
// };




// SF2 HeighField::raw_distribution(Arbre& arbre) const{

    
//     SF2 dens_arbre = densite_arbre(arbre);
//     dens_arbre.Normalize();

//     int rayon_arbre = arbre.get_rayon();
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
//     std::cout<<rayon_arbre<<std::endl;
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

//     SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

//     std::list< std::pair< std::pair<int,int> , int > >  list_arbre;

//     for (int k=0; k<100000; k++){
//         int rand_pos_x = rand()%nx;
//         int rand_pos_y = rand()%ny;


//         bool placement_not_possible = false;


//         for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
//             if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
//                 placement_not_possible = true;
//                 break;
//             }
//         }


//         if (placement_not_possible == false){
//             list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
//         }
//     }

//     for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        

//         bool test_dens_sapin = true;


//         float rand_test = ((double) rand() / (RAND_MAX));
//         int rand_pos_x = (*it).first.first;
//         int rand_pos_y = (*it).first.second;

//         if (rand_test <= dens_arbre.at(rand_pos_x,rand_pos_y)){
//             test_dens_sapin=true;
//         }
//         else{
//             continue;
//         }

//         if (test_dens_sapin == true){
//             res.at(rand_pos_x, rand_pos_y) = 1;
//             res.at(rand_pos_x+2, rand_pos_y) = 1;
//             res.at(rand_pos_x-2, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+2) = 1;
//             res.at(rand_pos_x, rand_pos_y-2) = 1;
//             res.at(rand_pos_x+1, rand_pos_y) = 1;
//             res.at(rand_pos_x-1, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+1) = 1;
//             res.at(rand_pos_x, rand_pos_y-1) = 1;
//         }


//     }

//     return res;
// };


// SF2 HeighField::double_raw_distribution(Arbre& arbre1, Arbre& arbre2) const{

    
//     SF2 dens_arbre1 = densite_arbre(arbre1);
//     dens_arbre1.Normalize();

//     SF2 dens_arbre2 = densite_arbre(arbre2);
//     dens_arbre2.Normalize();

//     int rayon_arbre1 = arbre1.get_rayon();
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
//     std::cout<<rayon_arbre1<<std::endl;
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

//     int rayon_arbre2 = arbre2.get_rayon();
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
//     std::cout<<rayon_arbre2<<std::endl;
//     std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

//     SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

//     std::list< std::pair< std::pair<int,int> , int > >  list_arbre;

//     for (int k=0; k<100000; k++){
//         int rand_pos_x = rand()%nx;
//         int rand_pos_y = rand()%ny;

//         float rand_test = ((double) rand() / (RAND_MAX)); //it fix memory issue !!??

//         bool placement_not_possible = false;


//         for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
//             if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre1) , *it ) ){
//                 placement_not_possible = true;
//                 break;
//             }
//         }
//         if (placement_not_possible == false){
//             list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre1));
//         }
//     }

//     for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        
//         bool test_dens_arbre1 = false;

//         float rand_test = ((double) rand() / (RAND_MAX));
//         int rand_pos_x = (*it).first.first;
//         int rand_pos_y = (*it).first.second;


//         if (rand_test <= dens_arbre1.at(rand_pos_x,rand_pos_y)){
//             test_dens_arbre1=true;
//         }
//         else{
//             continue;
//         }


//         if (test_dens_arbre1 == true){
//             res.at(rand_pos_x, rand_pos_y) = 1;
//             res.at(rand_pos_x+2, rand_pos_y) = 1;
//             res.at(rand_pos_x-2, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+2) = 1;
//             res.at(rand_pos_x, rand_pos_y-2) = 1;
//             res.at(rand_pos_x+1, rand_pos_y) = 1;
//             res.at(rand_pos_x-1, rand_pos_y) = 1;
//             res.at(rand_pos_x, rand_pos_y+1) = 1;
//             res.at(rand_pos_x, rand_pos_y-1) = 1;
//         }

//     }



//     for (int k=0; k<100000; k++){
//         int rand_pos_x = rand()%nx;
//         int rand_pos_y = rand()%ny;

//         float rand_test = ((double) rand() / (RAND_MAX)); //it fix memory issue !!??

//         bool placement_not_possible = false;


//         for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
//             if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre2) , *it ) ){
//                 placement_not_possible = true;
//                 break;
//             }
//         }
//         if (placement_not_possible == false){
//             list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre2));
//         }
//     }

//     for (std::list< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        
//         bool test_dens_arbre2 = false;
//         int rand_pos_x = (*it).first.first;
//         int rand_pos_y = (*it).first.second;

//         float rand_test = ((double) rand() / (RAND_MAX));

//         if (rand_test <= dens_arbre1.at(rand_pos_x,rand_pos_y)){
//             test_dens_arbre2=true;
//         }
//         else{
//             continue;
//         }

//         if (test_dens_arbre2 == true){
//             res.at(rand_pos_x, rand_pos_y) = 1;
//             res.at(rand_pos_x+2, rand_pos_y+2) = 1;
//             res.at(rand_pos_x-2, rand_pos_y-2) = 1;
//             res.at(rand_pos_x-2, rand_pos_y+2) = 1;
//             res.at(rand_pos_x+2, rand_pos_y-2) = 1;
//             res.at(rand_pos_x+1, rand_pos_y+1) = 1;
//             res.at(rand_pos_x-1, rand_pos_y-1) = 1;
//             res.at(rand_pos_x-1, rand_pos_y+1) = 1;
//             res.at(rand_pos_x-1, rand_pos_y-1) = 1;
//         }

//     }
    

//     return res;
// }



/******************************************
*          Classe LayeredField            *
******************************************/

// class LayeredField : public Grid2
// {
// protected:
//     SF2 bedrock;
//     SF2 sand;

// public:
//     LayeredField(const SF2 &bedrock) : Grid2(bedrock), bedrock(bedrock) {}
//     LayeredField(const SF2 &bedrock, const SF2 &sand) : Grid2(bedrock), bedrock(bedrock), sand(sand) {}

//     LayeredField(const HeighField H, const double percentSand) : Grid2(H){
//         double sandHeight = percentSand * H.max();

//         bedrock = SF2(Grid2(H));
//         sand = SF2(Grid2(H));

//         for (int k = 0 ; k < nx*ny ; k ++){
//             sand[k] = std::min(sandHeight, H[k]);
//             bedrock[k] = H[k] - sand[k];
//         }
//         bedrock.UpdateMinMax();
//         sand.UpdateMinMax();

//     }

//     HeighField toHeighField(){
//         HeighField H = HeighField(SF2(Grid2(*this)));

//         for (int k = 0; k < nx*ny ; k++ ){
//             H[k] = bedrock[k] + sand[k];
//         }
//         H.UpdateMinMax();
//         return H;
//     }

//     void TermalErosion(int);

// };


// void LayeredField::TermalErosion(int nbEpoch){

//     double tanTheta = 2000*tan(M_PI/4.0);
//     double k = 0.001;

//     for (int n = 0; n < nbEpoch ; n++){
//         std::multimap<double, std::pair<int,int>> myMap;
//         std::multimap<double, std::pair<int,int>>::iterator it = myMap.begin();
//         for (int i = 0; i < nx ; i ++){
//             for (int j = 0; j < ny ; j++){
//                 std::pair<double, std::pair<int,int>> pairToAdd = std::pair<double, std::pair<int,int>>(bedrock.at(i,j)+sand.at(i,j), std::pair<int,int>(i,j));
//                 myMap.insert(pairToAdd);           
//                 it++;
//             }
//         }    

//         for(auto it = myMap.begin(); it != myMap.end(); it++){
//             int i,j;
//             i = it->second.first;
//             j = it->second.second;

//             //recup les voisins tq s > tan theta
//             vec2 grad = bedrock.Gradient(i,j) + sand.Gradient(i,j);
//             double s = sqrt(grad*grad);
//             // std::cout << s << std::endl;
//             if (s > tanTheta){
//                 //Donne à chaque voisin un pourcentage du sable
//                 vec2 rGrad = grad.Normalized().round();

//                 int ibis = i - rGrad[0];
//                 int jbis = j - rGrad[1];
                
//                 // std::cout <<"ij " << i << ' ' << j << std::endl;
//                 // std::cout << ibis << ' ' << jbis << std::endl;

//                 sand.at(i,j) -= std::min(sand.at(i,j), k*(s-tanTheta));

//                 if ((ibis > 0) && (jbis > 0) && (ibis < nx) && (jbis < ny)){
//                     sand.at(ibis,jbis) += std::min(sand.at(i,j), k*(s-tanTheta));
//                 }
//             }
//         }

//     }

// }

/******************************************
*           Useful fonctions              *
******************************************/


void Compute_params( HeighField hf, QString s){

    // SF2 GRAD = hf.GradientNorm();
    // SF2 LAP = hf.LaplacianMap();
    SF2 SLOPE = hf.SlopeMap();
    // SF2 AVSLOPE = hf.AVGSlopeMap();
    SF2 ACCESS = hf.accessMap();
    // SF2 AreaStreepest = hf.StreamAreaStreepest();
    // SF2 Area = hf.StreamArea();
    // SF2 Power = hf.StreamPower();
    SF2 WET = hf.WetNessIndex();

    Sapin sapin = Sapin();
    Buisson buisson = Buisson();
    Buisson buisson2 = Buisson();
    Sapin sapin2 = Sapin();
    Sapin sapin3 = Sapin();
    Buisson buisson3 = Buisson();
    Sapin sapin4 = Sapin();
    Buisson buisson4 = Buisson();
    Sapin sapin5 = Sapin();
    Buisson buisson5 = Buisson();

    std::cout<<"densdensdensdensdensdensdensdensdensdens"<<std::endl;



    std::cout<<"distridistridistridistridistridistridistri"<<std::endl;


    SF2 DISTRI_BUISSON = hf.raw_distribution(buisson2);
    SF2 DISTRI_SAPIN = hf.raw_distribution(sapin2);
    SF2 DOUBLE_DISTRI = hf.double_raw_distribution(sapin4,buisson4);

    std::cout<<"sthrowsthrowsthrowsthrowsthrowsthrowsthrow"<<std::endl;

    SF2 THROW_SAPIN = hf.raw_dart_throwing(sapin3);
    SF2 THROW_BUISSON = hf.raw_dart_throwing(buisson3);


    std::cout<<"endsthrowsthrowsthrowsthrowsthrowsthrowsthrow"<<std::endl;

    SF2 SAPIN = hf.densite_arbre(sapin);
    //SF2 SAPIN_BUG = hf.densite_arbre(sapin5);
    //SF2 BUISSON = hf.densite_arbre(buisson);
    SF2 BUISSON = hf.densite_arbre(buisson);
    //SF2 BUISSON_BUG = hf.densite_arbre(buisson5);
    // QImage hauteur_phong = hf.Shade(hf);
    // QImage hauteur = hf.Export(hf);
    // QImage gradient = hf.Export(GRAD);
    // QImage laplacian = hf.Export(LAP);
    QImage slope = hf.Export(SLOPE);
    // QImage avslope = hf.Export(AVSLOPE);
    QImage acc = hf.Export(ACCESS);
    // QImage StreamAreaStreepest = hf.Export(AreaStreepest);
    // QImage StreamArea = hf.Export(Area);
    // QImage StreamPower = hf.Export(Power);
    QImage WetNessIndex = hf.Export(WET);
    QImage densite_sapin = hf.Export(SAPIN);
    //QImage densite_sapin_bug = hf.Export(SAPIN_BUG);
    QImage raw_throw_sapin = hf.Export(THROW_SAPIN);
    QImage sapin_raw_distribution = hf.Export(DISTRI_SAPIN);

    QImage densite_buisson = hf.Export(BUISSON);
    //QImage densite_buisson_bug = hf.Export(BUISSON_BUG);
    QImage raw_throw_buisson = hf.Export(THROW_BUISSON);
    QImage buisson_raw_distribution = hf.Export(DISTRI_BUISSON);
    QImage double_raw_distribution = hf.ExportColored(DOUBLE_DISTRI,2);




    // //std::cout <<" hauteur_phong "<<s<< std::endl;
    // hauteur_phong.save("Images/hauteur_phong"+s+".png");
    // //std::cout <<" hauteur "<<s<< std::endl;
    // hauteur.save("Images/hauteur"+s+".png");
    // //std::cout <<" gradient "<<s<< std::endl;
    // gradient.save("Images/gradient"+s+".png");
    // //std::cout <<" laplacian "<<s<< std::endl;
    // laplacian.save("Images/lapla"+s+".png");
    // //std::cout <<" slope "<<s<< std::endl;
    slope.save("Images/slope"+s+".png");
    // //std::cout <<" avslope "<<s<< std::endl;
    // avslope.save("Images/avslope"+s+".png");

    acc.save("Images/access" + s + ".png");

    // StreamAreaStreepest.save("Images/StreamAreaStreepest"+s+".png");
    // StreamArea.save("Images/StreamArea"+s+".png");
    // StreamPower.save("Images/StreamPower"+s+".png");
    WetNessIndex.save("Images/WetNessIndex"+s+".png");
    densite_sapin.save("Images/densite_sapin"+s+".png");
    //densite_sapin_bug.save("Images/densite_sapin_bug"+s+".png");
    densite_buisson.save("Images/densite_buisson"+s+".png");
    //densite_buisson_bug.save("Images/densite_buisson_bug"+s+".png");

    raw_throw_sapin.save("Images/throw_sapin"+s+".png");
    raw_throw_buisson.save("Images/throw_buisson"+s+".png");

    sapin_raw_distribution.save("Images/sapin_raw_distribution"+s+".png");
    buisson_raw_distribution.save("Images/buisson_raw_distribution"+s+".png");
    double_raw_distribution.save("Images/double_raw_distribution"+s+".png");


}


/******************************************
*             Focntion main               *
******************************************/

int main (int argc, char *argv[]){

    /*******************
    //Batterie de tests
    ********************/
    // vec3 pilou = vec3(0.0, 1.0, 2.0);

    // vec3 normalized = pilou.Normalized();

    //std::cout << normalized[0] << normalized[1] << normalized[2] << std::endl;

    QImage im;
    im.load("heightmap3.jpeg");
    //im.load("heightmap3.png");
    //im.load("best.png");

    HeighField hf = HeighField(im, Box2(vec2(0,0), vec2(100,100)), im.width(), im.height(), 100.0);
    // SF2 pente = hf.SlopeMap();
    // pente.UpdateMinMax();

    // hf.Smooth();
    // std::cout << pente.max() << ' ' << pente.min() << std::endl;
    // LayeredField lf = LayeredField(hf, 0.1);

    // lf.TermalErosion(10);
    
    // HeighField h = lf.toHeighField();
    // std::cout << h.max() << ' ' << h.min() << std::endl;
    // std::cout << hf.max() << ' ' << hf.min() << std::endl;
    // std::cout << hf.at(422,422) << ' ' << hf.at(422, 423) << std::endl;
    // std::cout << hf.Gradient(422,422) << std::endl;
    // std::cout << pente.at(422, 422) << std::endl;

    //Compute_params(hf, "");


    // hf.Clamp(4, 7);
    // Compute_params(hf, "_Clamp");

    hf.Smooth();
    Compute_params(hf, "_Smooth");

    // hf.ExportOBJ("Hf.obj");

    // SF2 acc = hf.accessMap();
    // acc.UpdateMinMax();

    QImage myIm = hf.Export(hf);
    // QImage myImMap = h.Export(h);
    myIm.save("pilou.png");
    // myImMap.save("pilouTrue.png");



    // std::ofstream myFile;
    // myFile.open("test.txt");

    // myFile << vec3(0.0, 1.0, 12.0);
    

    // hf.Blur();
    // Compute_params(hf, "_Blur");


    return 0;
}


