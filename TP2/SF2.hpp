#include "grid2.hpp"
#include "ScalarPoint2.hpp"
#include <QVector>
#include <stdlib.h>
#include <stdio.h>


#ifndef SF2D
#define SF2D

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
    SF2(const Grid2 &grid, const double &hz) : Grid2(grid)
    {
        field.resize(nx * ny, hz);
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

    const double operator[](int i) const
    {
        return field[i];
    }

    double &operator[] (int i)
    {
        return field[i];
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

    SF2 GradientNorm() const ;
    SF2 LaplacianMap();

    void Smooth();
    void Blur();
    void Normalize();
    void Clamp(float, float);
    QVector<ScalarPoint2> GetScalarPoints() const;
};

#endif