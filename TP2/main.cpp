
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
*           Useful fonctions              *
******************************************/


void Compute_params( HeighField hf, QString s){
    bool quicker = true;

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



    // // std::cout << "Enregistrement du terrain" << std::endl;
    
    // // std::cout << "Hauteurs" << std::endl;
    // // QImage hauteur = hf.Export(hf);
    // // hauteur.save("Images/hauteur"+s+".png");

    // std::cout << "Hauteurs Phong" << std::endl;
    // QImage hauteur_phong = hf.Shade(hf);
    // hauteur_phong.save("Images/hauteur_phong"+s+".png");



    // std::cout<<"Calculs pour affichage des propriétés de terrain"<<std::endl;
    
    // // std::cout<<"Gradient"<<std::endl;
    // // SF2 GRAD = hf.GradientNorm();
    // // QImage gradient = hf.Export(GRAD);
    // // gradient.save("Images/gradient"+s+".png");

    // // std::cout << "Laplacien" << std::endl;
    // // SF2 LAP = hf.LaplacianMap();
    // // QImage laplacian = hf.Export(LAP);
    // // laplacian.save("Images/lapla"+s+".png");

    std::cout << "Pente" << std::endl;
    SF2 SLOPE = hf.SlopeMap();
    QImage slope = hf.Export(SLOPE);
    slope.save("Images/slope"+s+".png");

    std::cout << "Pente Moyenne" << std::endl;
    SF2 AVSLOPE = hf.AVGSlopeMap();
    QImage avslope = hf.Export(AVSLOPE);
    avslope.save("Images/avslope"+s+".png");

    std::cout << "Accessibilite" << std::endl;
    SF2 ACCESS = hf.accessMap();
    QImage acc = hf.Export(ACCESS);
    acc.save("Images/access" + s + ".png");

    // // std::cout << "Stream Area Streepest" << std::endl;
    // // SF2 AreaStreepest = hf.StreamAreaStreepest();
    // // QImage StreamAreaStreepest = hf.Export(AreaStreepest);
    // // StreamAreaStreepest.save("Images/StreamAreaStreepest"+s+".png");

    // // std::cout << "Stream Area" << std::endl;
    // // SF2 Area = hf.StreamArea();
    // // QImage StreamArea = hf.Export(Area);
    // // StreamArea.save("Images/StreamArea"+s+".png");

    std::cout << "Stream Power" << std::endl;
    SF2 Power = hf.StreamPower();
    QImage StreamPower = hf.Export(Power);
    StreamPower.save("Images/StreamPower"+s+".png");

    std::cout << "Wetness Index" << std::endl;
    SF2 WET = hf.WetNessIndex();
    QImage WetNessIndex = hf.Export(WET);
    WetNessIndex.save("Images/WetNessIndex"+s+".png");


    std::cout<<"Calculs pour affichage des distributions"<<std::endl;
    
    std::cout<<"Buissons"<<std::endl;
    SF2 DISTRI_BUISSON = hf.raw_distribution(buisson2, quicker);
    QImage buisson_raw_distribution = hf.Export(DISTRI_BUISSON);
    buisson_raw_distribution.save("Images/buisson_raw_distribution"+s+".png");

    std::cout <<"Sapins" << std::endl;
    SF2 DISTRI_SAPIN = hf.raw_distribution(sapin2, quicker);
    QImage sapin_raw_distribution = hf.Export(DISTRI_SAPIN);
    sapin_raw_distribution.save("Images/sapin_raw_distribution"+s+".png");

    std::cout << "Combinaison sapins et buissons" << std::endl;
    SF2 DOUBLE_DISTRI = hf.double_raw_distribution_quicker(sapin4,buisson4, quicker);
    //hf.Normalize();
    QImage double_raw_distribution = hf.ExportColored(DOUBLE_DISTRI,2,true);
    double_raw_distribution.save("Images/double_raw_distribution"+s+".png");

    // Ancienne version de la double distribution
    // SF2 DOUBLE_DISTRI =  hf.double_raw_distribution(sapin4,buisson4);// hf.double_raw_distribution_quicker(sapin4,buisson4);


    // std::cout<<"Calcul pour affichage du lancer de fléchettes"<<std::endl;
    
    std::cout << "Sapins" << std::endl;
    SF2 THROW_SAPIN = hf.raw_dart_throwing(sapin3, quicker);
    QImage raw_throw_sapin = hf.Export(THROW_SAPIN);
    raw_throw_sapin.save("Images/throw_sapin"+s+".png");

    std::cout << "Buissons" << std::endl;
    SF2 THROW_BUISSON = hf.raw_dart_throwing(buisson3, quicker);
    QImage raw_throw_buisson = hf.Export(THROW_BUISSON);
    raw_throw_buisson.save("Images/throw_buisson"+s+".png");

    // std::cout<<"Calculs pour affichage des densités"<<std::endl;
    
    // std::cout << "Sapins" << std::endl;
    // SF2 SAPIN = hf.densite_arbre(sapin);
    // QImage densite_sapin = hf.Export(SAPIN);
    // densite_sapin.save("Images/densite_sapin"+s+".png");

    // // // Tests
    // // SF2 SAPIN_BUG = hf.densite_arbre(sapin5);
    // // QImage densite_sapin_bug = hf.Export(SAPIN_BUG);
    // // densite_sapin_bug.save("Images/densite_sapin_bug"+s+".png");

    // std::cout << "Buissons" << std::endl;
    // SF2 BUISSON = hf.densite_arbre(buisson);
    // QImage densite_buisson = hf.Export(BUISSON);
    // densite_buisson.save("Images/densite_buisson"+s+".png");

    // // // Tests
    // // SF2 BUISSON_BUG = hf.densite_arbre(buisson5);
    // // QImage densite_buisson_bug = hf.Export(BUISSON_BUG);
    // // densite_buisson_bug.save("Images/densite_buisson_bug"+s+".png");





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
    //im.load("ImagesToTest/im500x500.jpeg");
    //im.load("montagne.png");
    //im.load("best.png");

    HeighField hf = HeighField(im, Box2(vec2(0,0), vec2(im.width(),im.height())), im.width(), im.height(), 100.0);
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


