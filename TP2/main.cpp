
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <QImage>
#include <algorithm>
#include <fstream>

/******************************************
*               Classe vec2               *
******************************************/
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

    inline friend  std::ostream& operator<<(std::ostream& os, vec2 p){
        os << p.x << ' ' << p.y;
        return os;
    }

    inline vec2 round(){ return vec2(std::round(x), std::round(y));}

    inline vec2 operator<(vec2 p){ return vec2(x < p.x, y < p.y);}
    inline vec2 operator>(vec2 p){ return vec2(x > p.x, y > p.y);}

    friend class vec3;
};


/******************************************
*               Classe vec3               *
******************************************/
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


/******************************************
*               Classe Box2               *
******************************************/
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


/******************************************
*               Classe Grid2              *
******************************************/

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


/******************************************
*               Classe SF2                *
******************************************/
class SF2 : public Grid2
{
protected:
    std::vector<double> field;
    double maxVal = -1000000;
    double minVal = 1000000;
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

    void UpdateMinMax(){
        for (int k = 0; k < field.size(); k++){
            double val = field[k];
            if (val > maxVal){
                maxVal = val;
            }
            if (val < minVal){
                minVal = val;
            }
        }
    }

    const double max() const {return maxVal;}
    const double min() const {return minVal;}

    double getConv(int, int, float, float, float);

    SF2 GradientNorm();
    SF2 LaplacianMap();

    void Smooth();
    void Blur();
    void Normalize();
    void Clamp(float, float);

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

double SF2::getConv(int i, int j, float mid, float side, float diag){
    float nb = 0;
    float value = 0;

    //top left
    if ((i > 0)&&(j>0)){
        value += diag*field[Index(i-1, j-1)];
        nb+= diag;
    }
    //top
    if (i > 0){
        value += side*field[Index(i-1, j)];
        nb+= side;
    }
    //top right
    if ((i > 0)&&(j < ny-1)){
        value += diag*field[Index(i-1, j+1)];
        nb+=diag;
    }
    //left
    if (j>0){
        value += side*field[Index(i, j-1)];
        nb+=side;
    }
    //right
    if (j<ny-1){
        value += side*field[Index(i, j+1)];
        nb+=side;
    }
    //bottom left
    if ((i < nx-1) && (j>0)){
        value += diag*field[Index(i+1, j-1)];
        nb+=diag;
    }
    //bottom
    if (i<nx-1){
        value += side*field[Index(i+1, j)];
        nb+=side;
    }
    //bottom right
    if ((i < nx-1) && (j<ny-1)){
        value += diag*field[Index(i+1, j+1)];
        nb+=diag;
    }
    //center
    value += mid*field[Index(i,j)];
    nb+=mid;

    //normalization
    value /= nb;

    return value;
}

SF2 SF2::GradientNorm(){
    SF2 gradNorm = SF2(Grid2(*this));

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            vec2 grad = Gradient(i,j);
            gradNorm.at(i,j) = sqrt(grad*grad);
        }
    }

    return gradNorm;
}

SF2 SF2::LaplacianMap(){
    SF2 LaplMap = SF2(Grid2(*this));

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            LaplMap.at(i,j) = Laplacian(i,j);
        }
    }

    return LaplMap;
}


void SF2::Smooth(){
    std::vector<double> smoothed;
    smoothed.resize(nx*ny);

    for (int i = 0; i < nx ; i ++){
        for (int j = 0; j < ny ; j ++){
            int ind = Index(i,j);
            
            smoothed[ind] = getConv(i,j,4,2,1);
        }
    }

    field = smoothed;
}

void SF2::Blur(){
    std::vector<double> blured;
    blured.resize(nx*ny);

    for (int i = 0; i < nx ; i ++){
        for (int j = 0; j < ny ; j ++){
            int ind = Index(i,j);
            
            blured[ind] = getConv(i,j,1,1,1);
        }
    }

    field = blured;
}

void SF2::Normalize(){
    UpdateMinMax();

    for (int k = 0; k < field.size(); k++){
        field[k] = (field[k] - minVal)/(maxVal - minVal);
    }

}

void SF2::Clamp(float mini, float maxi){
    
    for (int k = 0; k < field.size(); k++){
        if (field[k] < mini){
            field[k] = mini;
        } else if (field[k] > maxi){
            field[k] = maxi;
        }
    }


}

/******************************************
*           Classe HeighField2            *
******************************************/
class HeighField : public SF2
{
protected:
public:
    HeighField(const SF2 &s) : SF2(s) {}
    HeighField(const QImage &image, const Box2 &box, double nx, double ny) /*:Grid2(box, image.width(),image.height()){*/
    {
        //et interpoler entre les deux z
        
        //Initialisation
        Grid2 grid = Grid2(box, nx, ny);
        *this = SF2(grid);

        float normFact = 10.0;

        //Remplissage des hauteurs
        for (int x = 0; x < nx; x ++){
            for (int y = 0; y < ny ; y++){
                field[Index(x,y)] = ((qRed(image.pixel(x,y)))/255.0 * normFact);
                // std::cout << field.back() << std::endl;
            }
        }

        UpdateMinMax();
        //std::cout << field.size() << ' ' << nx*ny << ' ' << nx << ' ' << ny << std::endl;
    }

    double Height(int i, int j) const { return at(i, j); } 
    double Slope(int i, int j) const
    {
        vec2 g = Gradient(i, j);
        return sqrt(g * g);
    }

    bool intersectRay(vec3, double&, vec3);

    double AverageSlope(int, int) const; 

    vec3 Vertex(int i, int j) const { return vec3(Grid2::Vertex(i, j), Height(i, j)); }
    vec3 Normal(int i, int j) const { return vec3(-Gradient(i, j), 1.0).Normalized(); }

    SF2 SlopeMap(){return GradientNorm();}
    SF2 AVGSlopeMap(){    
        SF2 avgMap = SF2(Grid2(*this));

        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                avgMap.at(i,j) = AverageSlope(i,j);
            }
        }

        return avgMap;
    }

    //Exportation sous format d'image ! 
    QImage Export(SF2, bool) const;
    QImage Shade(SF2 mapToExport) const {return Export(mapToExport, true);}

    void ExportOBJ(char*);

};

QImage HeighField::Export(SF2 mapToExport, bool vis = false) const
{
    QImage image(nx,ny, QImage::Format_ARGB32);

    const vec3 lightdir = vec3(2.0, 1.0, 4.0).Normalized();
    mapToExport.Normalize();
    for (int i=0; i <nx; i++)
    {
        for (int j=0; j <ny; j++)
        {
            int mapVal = (mapToExport.at(i,j)*255);
            if (vis){
                vec3 n = Normal(i,j);
                double d =n*lightdir;
                d=(1.0+d)/2.0;
                d *= d;
                mapVal *= d;
            }

            image.setPixel(i,j,qRgb(mapVal, mapVal, mapVal));

        }
    }

    return image;

}

bool HeighField::intersectRay(vec3 rayDir, double& t, vec3 origin){

    double epsilon = 1.0;
    vec3 ray = origin + t*rayDir;

    while (Inside(ceil(ray[0]), ceil(ray[1])) && (ray[2] < maxVal)){
        //get z value by finding the maximum heigh between the 3 closest points
        vec2 proj = vec2(ray[0], ray[1]);
        vec2 roundProj = proj.round();
        //vec2 diff = (proj > roundProj) * 2 - vec2(1.0, 1.0);

        double height = at(roundProj[0], roundProj[1]);
        
        if (height > ray[2]){
            return true;
        } else {
            t += epsilon;
            ray = origin + t*rayDir;
        }
    }


    return false;
}


double HeighField::AverageSlope(int i, int j) const{
    double avgSlope = 0.0;  
    double slopeVal;  
    float nb = 0;

    float curHeight = at(i,j);

    float diagVal = sqrt(celldiagonal[0]*celldiagonal[0] + celldiagonal[1]*celldiagonal[1]);

    //top left
    if ((i > 0)&&(j>0)){
        slopeVal =  (curHeight - at(i-1,j-1)) / diagVal;
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //top
    if (i > 0){
        slopeVal= (curHeight - at(i-1, j)) * inversecelldiagonal[0];
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //top right
    if ((i > 0)&&(j < ny-1)){
        slopeVal = (curHeight - at(i-1,j+1)) / diagVal;
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //left
    if (j>0){
        slopeVal= (curHeight - at(i, j-1)) * inversecelldiagonal[1];
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //right
    if (j<ny-1){
        slopeVal= (curHeight - at(i, j+1)) * inversecelldiagonal[1];
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //bottom left
    if ((i < nx-1) && (j>0)){
        slopeVal = (curHeight - at(i+1,j-1)) / diagVal;
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //bottom
    if (i<nx-1){
        slopeVal = (curHeight - at(i+1, j)) * inversecelldiagonal[0];
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    //bottom right
    if ((i < nx-1) && (j<ny-1)){
        slopeVal = (curHeight - at(i+1,j+1)) / diagVal;
        avgSlope += sqrt(slopeVal*slopeVal);
        nb++;
    }
    
    //normalization
    avgSlope /= nb;

    return avgSlope;
}

void HeighField::ExportOBJ(char* filename){
    std::ofstream objFile;
    objFile.open(filename);

    //Write vertices
    for (int i = 0; i < nx; i ++){
        for (int j = 0; j < ny ; j++){
            vec3 vertex = Vertex(i,j);
            objFile << "v " << vertex << std::endl;
        }
    }
    
    objFile << std::endl;

    //Write normals
    for (int i = 0; i < nx; i ++){
        for (int j = 0; j < ny ; j++){
            vec3 normal = Normal(i,j);
            objFile << "vn " << normal << std::endl;
        }
    }
    
    objFile << std::endl;
    
    //Write texture coords
    for (int i = 0; i < nx; i ++){
        for (int j = 0; j < ny ; j++){
            vec2 texCoord = vec2(float(i)/float(nx), float(j)/float(ny));
            objFile << "vt " << texCoord << std::endl;
        }
    }

    objFile << std::endl;

    //Write faces
    for (int i = 0; i < nx-1; i ++){
        for (int j = 0; j < ny-1 ; j++){
            int ind = Index(i,j);
            int indr = Index(i, j+1);
            int indb = Index(i+1, j);
            int indd = Index(i+1, j+1);
            objFile << "f " << ind << '/' << ind << '/' << ind << ' '
            << indr << '/' << indr << '/' << indr << ' ' 
            << indb << '/' << indb << '/' << indb << ' ' << std::endl;

            objFile << "f " << indr << '/' << indr << '/' << indr << ' '
            << indb << '/' << indb << '/' << indb << ' ' 
            << indd << '/' << indd << '/' << indd << ' ' << std::endl;


        }
    }    


}


/******************************************
*          Classe LayeredField            *
******************************************/

class LayeredField : public Grid2
{
protected:
    SF2 bedrock;
    SF2 sand;

public:
    LayeredField(const SF2 &bedrock) : Grid2(bedrock), bedrock(bedrock) {}
    LayeredField(const SF2 &bedrock, const SF2 &sand) : Grid2(bedrock), bedrock(bedrock), sand(sand) {}
};

/******************************************
*           Useful fonctions              *
******************************************/


void Compute_params( HeighField hf, QString s){

    SF2 GRAD = hf.GradientNorm();
    SF2 LAP = hf.LaplacianMap();
    SF2 SLOPE = hf.SlopeMap();
    SF2 AVSLOPE = hf.AVGSlopeMap();

    QImage hauteur_phong = hf.Shade(hf);
    QImage hauteur = hf.Export(hf);
    QImage gradient = hf.Export(GRAD);
    QImage laplacian = hf.Export(LAP);
    QImage slope = hf.Export(SLOPE);
    QImage avslope = hf.Export(AVSLOPE);

    //std::cout <<" hauteur_phong "<<s<< std::endl;
    hauteur_phong.save("Images/hauteur_phong"+s+".png");
    //std::cout <<" hauteur "<<s<< std::endl;
    hauteur.save("Images/hauteur"+s+".png");
    //std::cout <<" gradient "<<s<< std::endl;
    gradient.save("Images/gradient"+s+".png");
    //std::cout <<" laplacian "<<s<< std::endl;
    laplacian.save("Images/lapla"+s+".png");
    //std::cout <<" slope "<<s<< std::endl;
    slope.save("Images/slope"+s+".png");
    //std::cout <<" avslope "<<s<< std::endl;
    avslope.save("Images/avslope"+s+".png");
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
    std::cout << im.width() << std::endl;
    HeighField hf = HeighField(im, Box2(vec2(0,0), vec2(1,1)), im.width(), im.height());

    Compute_params(hf, "");

    // hf.Clamp(4, 7);
    // Compute_params(hf, "_Clamp");

    // hf.Smooth();
    // Compute_params(hf, "_Smooth");

    // hf.ExportOBJ("Hf.obj");


    // QImage myIm = hf.Export(hf, true);
    // QImage myImMap = hf.Export(hf);
    // myIm.save("pilou.png");
    // myImMap.save("pilouTrue.png");

    // std::ofstream myFile;
    // myFile.open("test.txt");

    // myFile << vec3(0.0, 1.0, 12.0);
    

    // hf.Blur();
    // Compute_params(hf, "_Blur");


    return 0;
}



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
