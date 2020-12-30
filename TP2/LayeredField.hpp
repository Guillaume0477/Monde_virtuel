#include "HeighField.hpp"
#include <map>
#include <stdlib.h>
#include <stdio.h>

#ifndef LAYEREDFIELD
#define LAYEREDFIELD

class LayeredField : public Grid2
{
protected:
    SF2 bedrock;
    SF2 sand;

public:
    LayeredField(const SF2 &bedrock) : Grid2(bedrock), bedrock(bedrock) {}
    LayeredField(const SF2 &bedrock, const SF2 &sand) : Grid2(bedrock), bedrock(bedrock), sand(sand) {}

    LayeredField(const HeighField H, const double percentSand) : Grid2(H){
        double sandHeight = percentSand * H.max();

        bedrock = SF2(Grid2(H));
        sand = SF2(Grid2(H));

        for (int k = 0 ; k < nx*ny ; k ++){
            sand[k] = std::min(sandHeight, H[k]);
            bedrock[k] = H[k] - sand[k];
        }
        bedrock.UpdateMinMax();
        sand.UpdateMinMax();

    }

    HeighField toHeighField(){
        HeighField H = HeighField(SF2(Grid2(*this)));

        for (int k = 0; k < nx*ny ; k++ ){
            H[k] = bedrock[k] + sand[k];
        }
        H.UpdateMinMax();
        return H;
    }

    void TermalErosion(int);

};

#endif