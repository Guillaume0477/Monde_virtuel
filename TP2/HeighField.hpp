#include "SF2.hpp"
#include "Sapin.hpp"
#include "Buisson.hpp"
#include "vec3.hpp"

#include <QImage>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#ifndef HEIGHFIELD
#define HEIGHFIELD

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

        float normFact = 500.0;

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
    double access(int,int, int N = 20);
    SF2 accessMap(){    
        SF2 accMap = SF2(Grid2(*this));

        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                accMap.at(i,j) = access(i,j);
            }
        }

        return accMap;
    }


    double AverageSlope(int, int) const; 

    vec3 Vertex(int i, int j) const { return vec3(Grid2::Vertex(i, j), Height(i, j)); }
    vec3 Normal(int i, int j) const { return vec3(-Gradient(i, j), 1.0).Normalized(); }

    SF2 SlopeMap() const {return GradientNorm();}
    SF2 AVGSlopeMap(){    
        SF2 avgMap = SF2(Grid2(*this));

        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                avgMap.at(i,j) = AverageSlope(i,j);
            }
        }
        return avgMap;
    }
    int CheckFlowSlope( const QPoint& p, QPoint* point, double* height, double* slope, double* nslone, int& mask) const;
    SF2 StreamAreaStreepest() const;
    SF2 StreamArea() const;
    SF2 StreamPower() const;
    SF2 WetNessIndex() const;
    //SF2 densite_sapin(Arbre arbre) const;
    SF2 densite_arbre(Arbre& arbre) const; //not working
    SF2 densite_sapin() const;
    SF2 densite_buisson() const;
    SF2 sapin_raw_distribution() const;
    SF2 raw_dart_throwing(Arbre& arbre) const;
    SF2 raw_distribution(Arbre& arbre) const;
    SF2 double_raw_distribution(Arbre& arbre1, Arbre& arbre2) const;


    //Exportation sous format d'image ! 
    QImage Export(SF2, bool = false) const;
    QImage Shade(SF2 mapToExport) const {return Export(mapToExport, true);}

    void ExportOBJ(char*);



};

#endif