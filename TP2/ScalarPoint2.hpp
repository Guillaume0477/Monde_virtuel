#include <QPoint>
#include <stdlib.h>
#include <stdio.h>


#ifndef SCALARPOINT2
#define SCALARPOINT2

class ScalarPoint2
{
protected:
    QPoint p;
    double z;
public:
    ScalarPoint2() {
        p = QPoint(0.0,0.0);
        z=0.0;
    }
    ScalarPoint2(const QPoint& p1, const double& z1){
        p = p1; 
        z= z1;
    };
    friend bool operator<(const ScalarPoint2& a, const ScalarPoint2& b) {
        return a.z < b.z;
    };
    QPoint Point() const {return p;};
    double Scalar() const {return z;};

};

#endif