#include "LayeredField.hpp"

void LayeredField::TermalErosion(int nbEpoch){

    double tanTheta = 2000*tan(M_PI/4.0);
    double k = 0.001;

    for (int n = 0; n < nbEpoch ; n++){
        std::multimap<double, std::pair<int,int>> myMap;
        std::multimap<double, std::pair<int,int>>::iterator it = myMap.begin();
        for (int i = 0; i < nx ; i ++){
            for (int j = 0; j < ny ; j++){
                std::pair<double, std::pair<int,int>> pairToAdd = std::pair<double, std::pair<int,int>>(bedrock.at(i,j)+sand.at(i,j), std::pair<int,int>(i,j));
                myMap.insert(pairToAdd);           
                it++;
            }
        }    

        for(auto it = myMap.begin(); it != myMap.end(); it++){
            int i,j;
            i = it->second.first;
            j = it->second.second;

            //recup les voisins tq s > tan theta
            vec2 grad = bedrock.Gradient(i,j) + sand.Gradient(i,j);
            double s = sqrt(grad*grad);
            // std::cout << s << std::endl;
            if (s > tanTheta){
                //Donne à chaque voisin un pourcentage du sable
                vec2 rGrad = grad.Normalized().round();

                int ibis = i - rGrad[0];
                int jbis = j - rGrad[1];
                
                // std::cout <<"ij " << i << ' ' << j << std::endl;
                // std::cout << ibis << ' ' << jbis << std::endl;

                sand.at(i,j) -= std::min(sand.at(i,j), k*(s-tanTheta));

                if ((ibis > 0) && (jbis > 0) && (ibis < nx) && (jbis < ny)){
                    sand.at(ibis,jbis) += std::min(sand.at(i,j), k*(s-tanTheta));
                }
            }
        }

    }

}