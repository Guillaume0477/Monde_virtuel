#include "HeighField.hpp"


QImage HeighField::ExportColored(SF2 mapToExport, int nbColors, bool vis) const
{
    QImage image(nx,ny, QImage::Format_ARGB32);

    const vec3 lightdir = vec3(2.0, 1.0, 4.0).Normalized();
    mapToExport.Normalize();

    for (int i=0; i <nx; i++)
    {
        for (int j=0; j <ny; j++)
        {
            float value = mapToExport.at(i,j);
            int mapVal = (value*255);
            if (vis){
                vec3 n = Normal(i,j);
                double d =n*lightdir;
                d=(1.0+d)/2.0;
                d *= d;
                mapVal *= d;
            }

            if (nbColors == 1){
                image.setPixel(i,j,qRgb(mapVal, 0.0, 0.0));
            } else if (nbColors == 2){
                if (value > 2.0/(nbColors+1)){
                    image.setPixel(i,j,qRgb(124,252,0));
                } else if ((value > 0.0)&&(value<2.0/(1.0+float(nbColors)))){
                    image.setPixel(i,j,qRgb(218,242,0));
                } else {
                    image.setPixel(i,j,qRgb(0, 0, 0));
                }
            }

        }
    }

    return image;

}

QImage HeighField::Export(SF2 mapToExport, bool vis) const
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

bool HeighField::intersectRay(vec3 rayDir, double& t, vec3 origin) const{

    double epsilon = 1.0;
    vec3 ray = origin + t*rayDir;

    while (Inside(std::round(ray[0]), std::round(ray[1])) && (ray[2] < maxVal) && (ray[2]>0)){
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

double HeighField::access(int i, int j, int Nray) const {
    double accessVal = 0.0;

    vec3 origin = vec3(i,j,at(i,j));
    for (int r = 0; r < Nray ; r++){
        //Create ray direction randomly on hemisphere (formula (34) GI Compedium)
        vec3 rayDir;
        double t = 1;

        double r1,r2; 
        r1 = double(rand())/RAND_MAX;
        r2 = double(rand())/RAND_MAX;
        //Direction aléatoire dans l'espace tangent
        vec3 rayDirTan = vec3(cos(2*M_PI*r1) * sqrt(1-r2*r2), sin(2*M_PI*r1)*sqrt(1-r2*r2), r2);

        vec3 ta, bi, n;
        n = Normal(i,j);

        //Code from J-C.IEHL to transfer a vector from tangent space to world space
        float sign= n[2] < 0 ? -1 : 1;             
        float a= -1.0 / (sign + n[2]);
        float d= n[0] * n[1] * a;
        ta= vec3(1.0 + sign * n[0] * n[0] * a, sign * d, -sign * n[0]);
        bi= vec3(d, sign + n[1] * n[1] * a, -n[1]);

        rayDir = rayDirTan[0] * ta.Normalized() + rayDirTan[1] * bi.Normalized() + rayDirTan[2] * n.Normalized();
        //End of the code inspired from J-C.IEHL

        //Check the visibility
        if (intersectRay(rayDir.Normalized(), t, origin) == false){
            accessVal++;
        }

    }

    return accessVal/double(Nray);
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


int HeighField::CheckFlowSlope( const QPoint& p, QPoint* point, double* height, double* slope, double* nslone, int& mask) const {

    int n = 0;
    double zp = at(p.x(),p.y());
    double slopesum = 0.0;

    const QPoint next[8] = { QPoint (1,0) ,QPoint (1,1) ,QPoint (0,1) ,QPoint (-1,1) ,QPoint (0,-1) ,QPoint (-1,-1) ,QPoint (-1,0), QPoint (1,-1) };
    const double length[8] = {1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0),1.0 , sqrt(2.0)};

    //std::cout<<"zp "<<zp<<std::endl;

    mask = 0;
    for (int i=0; i<8; i++){
        QPoint b = p + next[i];
        // std::cout<<"next"<<next[i].x()<<" "<<next[i].y() <<std::endl;
        //std::cout<<"b "<<b.x()<<" "<<b.y() <<std::endl;

        double step = at(b.x(), b.y()) - zp;
        
        if (b.x()>=nx || b.y()>=ny || b.x()<=0.0 || b.y()<=0.0){
            continue;
        }
        
        //if(!Box2::Inside(vec2(b.x(),b.y()))){continue;};


        
        // std::cout<<"at_b "<<at(b.x(), b.y())<<std::endl;
        // std::cout<<"step "<<step<<std::endl;

        if (step <0.0)
        {
            point[n] = b;
            //std::cout<<"point "<<point[n].x()<<" "<<point[n].y() <<std::endl;
            height[n] = -step;
            slope[n] = - step/ length[i];
            slopesum += slope[n];
            n++;
            mask |= 1 << i;
        }

    }

    for (int k=0; k<n; k++){
        nslone[k] = slope[k] / slopesum;
    }
    return n;

}

SF2 HeighField::StreamAreaStreepest() const{

    SF2 stream = SF2(Grid2(Box2(a,b),nx,ny),1.0);

    // SF2 stream = SF2(Grid2(*this));

    // for (int i = 0; i < nx; i++){
    //     for (int j = 0; j < ny; j++){
    //         stream.at(i,j) = 1.0;
    //     }
    // }
    // stream.at(2,2) = 0.01;

    QVector<ScalarPoint2> QEE = GetScalarPoints();

    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    // for (int i = QEE.size() -1; i>=0; i--){
    //     std::cout<<QEE.at(i).Scalar()<<std::endl;
    //     //stream.at(q[k].x(), q[k].y()) += sp;
    // }
    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    
    // QVector<ScalarPoint2> QEE = QVector(ScalarPoint2(next,length))
    std::sort(QEE.begin(), QEE.end());

    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    // for (int i = QEE.size() -1; i>=0; i--){
    //     std::cout<<QEE.at(i).Scalar()<<std::endl;
    //     //stream.at(q[k].x(), q[k].y()) += sp;
    // }
    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;

    for (int i = QEE.size() -1; i>=0; i--){
        
        // QPoint q = QPoint(2,2);
        // double z = 5;
        // ScalarPoint2 p = ScalarPoint2(q,z);
        // p.Point(); 
        //QEE.at(i).Point();
        QPoint p = QEE.at(i).Point();


        

        QPoint q[8];
        double h[8];
        double s[8];
        double sn[8];
        int m;

        int n = CheckFlowSlope(p,q,h,s,sn,m);

        // std::cout<<"p "<<p.x()<<" , "<<p.y()<<std::endl;
        // for (int j=0; j<8; j++){
        //     std::cout<<"q "<<j<<q[j].x()<<" "<<q[j].y()<<std::endl;
        //     std::cout<<"h "<<j<<h[j]<<std::endl;
        //     std::cout<<"s "<<j<<s[j]<<std::endl;
        //     std::cout<<"sn "<<j<<sn[j]<<std::endl;
        // }

        if (n>0)
        {
            double ss = s[0]; //calcul max s
            int k=0;
            for (int j=0; j<n;j++)
            {
                if (s[j] > ss) {
                    k=j;
                    ss=s[j];
                }
            }
            const double sp = stream.at(p.x(),p.y());
            stream.at(q[k].x(), q[k].y()) += sp;
        }

    }
    return stream;
}




SF2 HeighField::StreamArea() const{

    SF2 stream = SF2(Grid2(Box2(a,b),nx,ny),1.0);

    // SF2 stream = SF2(Grid2(*this));

    // for (int i = 0; i < nx; i++){
    //     for (int j = 0; j < ny; j++){
    //         stream.at(i,j) = 1.0;
    //     }
    // }
    // stream.at(2,2) = 0.01;

    QVector<ScalarPoint2> QEE = GetScalarPoints();

    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    // for (int i = QEE.size() -1; i>=0; i--){
    //     std::cout<<QEE.at(i).Scalar()<<std::endl;
    //     //stream.at(q[k].x(), q[k].y()) += sp;
    // }
    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    
    // QVector<ScalarPoint2> QEE = QVector(ScalarPoint2(next,length))
    std::sort(QEE.begin(), QEE.end());

    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;
    // for (int i = QEE.size() -1; i>=0; i--){
    //     std::cout<<QEE.at(i).Scalar()<<std::endl;
    //     //stream.at(q[k].x(), q[k].y()) += sp;
    // }
    // std::cout<<"QEE.at(i).Scalar()"<<std::endl;

    for (int i = QEE.size() -1; i>=0; i--){
        
        // QPoint q = QPoint(2,2);
        // double z = 5;
        // ScalarPoint2 p = ScalarPoint2(q,z);
        // p.Point(); 
        //QEE.at(i).Point();
        QPoint p = QEE.at(i).Point();


        

        QPoint q[8];
        double h[8];
        double s[8];
        double sn[8];
        int m;

        int n = CheckFlowSlope(p,q,h,s,sn,m);

        // std::cout<<"p "<<p.x()<<" , "<<p.y()<<std::endl;
        // for (int j=0; j<8; j++){
        //     std::cout<<"q "<<j<<" "<<q[j].x()<<" "<<q[j].y()<<std::endl;
        //     std::cout<<"h "<<j<<" "<<h[j]<<std::endl;
        //     std::cout<<"s "<<j<<" "<<s[j]<<std::endl;
        //     std::cout<<"sn "<<j<<" "<<sn[j]<<std::endl;
        // }
        // std::cout<<"N "<<n<<std::endl;
        if (n>0)
        {
            const double sp = stream.at(p.x(),p.y());

            for (int k=0; k<n;k++)
            {
                stream.at(q[k].x(), q[k].y()) += sp * sn[k];
                //std::cout<<"sn "<<k<<" "<<sn[k]<<std::endl;


            }
            

        }

    }
    return stream;
};






SF2 HeighField::StreamPower() const{
    SF2 stream = StreamArea();
    SF2 slope = SlopeMap();
    SF2 res(Grid2(*this));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            res.at(i, j) = sqrt(stream.at(i, j)) * slope.at(i, j);
        }
    }
    return res;
}

SF2 HeighField::WetNessIndex() const{
    SF2 stream = StreamArea();
    SF2 slope = SlopeMap();
    SF2 res(Grid2(*this));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            res.at(i, j) = log(stream.at(i, j)) * slope.at(i, j);
        }
    }
    return res;
}



SF2 HeighField::densite_arbre(Arbre& arbre) const{


    SF2 stream = StreamArea();
    stream.Normalize();
    SF2 humidity = WetNessIndex();
    humidity.Normalize();
    SF2 slope = SlopeMap();
    slope.Normalize();
    SF2 access = accessMap();
    access.Normalize();
    SF2 res = SF2(Grid2(Box2(a,b),nx,ny),1.0); //init 1

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k=0; k<1; k++){
            res.at(i, j) = std::min( arbre.humidity(humidity.at(i,j)) , res.at(i, j) );
            res.at(i, j) = std::min( arbre.slope(slope.at(i,j)) , res.at(i, j) );
            res.at(i, j) = std::min( arbre.stream(stream.at(i,j)) , res.at(i, j) );
            res.at(i, j) = std::min( arbre.acces(access.at(i,j)) , res.at(i, j) );
            }
        }
    }
    return res;
}



std::vector< std::pair< std::pair<int,int> , int > > HeighField::make_dart_throwing(Arbre& arbre) const{


    int rayon_arbre = arbre.get_rayon();
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    std::cout<<rayon_arbre<<std::endl;
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    
    std::vector< std::pair< std::pair<int,int> , int > >  list_arbre;

    for (int k=0; k<100000; k++){
        int rand_pos_x = rand()%nx;
        int rand_pos_y = rand()%ny;


        bool placement_not_possible = false;


        for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
            if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
                placement_not_possible = true;
                break;
            }
        }


        if (placement_not_possible == false){
            list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
        }
    }

    return list_arbre;
};

std::vector< std::pair< std::pair<int,int> , int > > HeighField::make_dart_throwing_quicker(Arbre& arbre) const{


    int rayon_arbre = arbre.get_rayon();
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    std::cout<<rayon_arbre<<std::endl;
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    
    std::vector< std::pair< std::pair<int,int> , int > >  list_arbre;

    int val = arbre.get_rayon()*2 +2;
    int nbx = int(nx/val);

    int nby = int(ny/val);
    for (int i = 0; i < nbx; i++){
        for (int j = 0; j < nby ; j++){
            bool done = false;
            int x, y;
            while (done == false){
                done = true;
                x = int(float(rand())/RAND_MAX*val) + val*i;
                y = int(float(rand())/RAND_MAX * val) + val * j;

                //Test si collision avec les voisins
                for(int k = 0; k < 4; k++){
                    int a = i - 1 + floor(k/3);
                    int b = j - 1 + k%3;

                    if ((a<0) || (b < 0) || (a >= nbx) || (b>=nby)){
                        continue;
                    }
                    std::pair<std::pair<int,int>,int> neighbour = list_arbre[a*nby+b];

                    float dist = sqrt((neighbour.first.first - x)*(neighbour.first.first - x) + (neighbour.first.second - y)*(neighbour.first.second - y));

                    if(dist < 2*arbre.get_rayon()){
                        done = false;
                        break;
                    }

                }
            }
            // std::cout << x << ' ' << y << std::endl;

            list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(x,y),arbre.get_rayon()));

        }
    }
    std::cout << "pilou" << std::endl;

    // for (int k=0; k<100000; k++){
    //     int rand_pos_x = rand()%nx;
    //     int rand_pos_y = rand()%ny;


    //     bool placement_not_possible = false;


    //     for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
    //         if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
    //             placement_not_possible = true;
    //             break;
    //         }
    //     }


    //     if (placement_not_possible == false){
    //         list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
    //     }
    // }
    return list_arbre;
};

SF2 HeighField::raw_dart_throwing(Arbre& arbre, bool dil) const{

    
    SF2 dens_arbre = densite_arbre(arbre);
    dens_arbre.Normalize();

    // int rayon_arbre = arbre.get_rayon();
    // std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    // std::cout<<rayon_arbre<<std::endl;
    // std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

    SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

    std::vector< std::pair< std::pair<int,int> , int > >  list_arbre;

    // for (int k=0; k<100000; k++){
    //     int rand_pos_x = rand()%nx;
    //     int rand_pos_y = rand()%ny;


    //     bool placement_not_possible = false;


    //     for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
    //         if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
    //             placement_not_possible = true;
    //             break;
    //         }
    //     }


    //     if (placement_not_possible == false){
    //         list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
    //     }
    // }

    list_arbre = make_dart_throwing_quicker(arbre);

    for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        

        bool test_dens_sapin = true;


        // float rand_test = ((double) rand() / (RAND_MAX));
        int rand_pos_x = (*it).first.first;
        int rand_pos_y = (*it).first.second;

        // if (rand_test <= dens_arbre.at(rand_pos_x,rand_pos_y)){
        //     test_dens_sapin=true;
        // }
        // else{
        //     continue;
        // }

        //if (test_dens_sapin == true){
        res.at(rand_pos_x, rand_pos_y) = 1;
        // res.at(rand_pos_x+2, rand_pos_y) = 1;
        // res.at(rand_pos_x-2, rand_pos_y) = 1;
        // res.at(rand_pos_x, rand_pos_y+2) = 1;
        // res.at(rand_pos_x, rand_pos_y-2) = 1;
        // res.at(rand_pos_x+1, rand_pos_y) = 1;
        // res.at(rand_pos_x-1, rand_pos_y) = 1;
        // res.at(rand_pos_x, rand_pos_y+1) = 1;
        // res.at(rand_pos_x, rand_pos_y-1) = 1;
        //}


    }
    if(dil){
        res.UpdateMinMax();
        res.Dilate(arbre.get_rayon(), res.max());
    }
    return res;
};




SF2 HeighField::raw_distribution(Arbre& arbre, bool dil) const{

    
    SF2 dens_arbre = densite_arbre(arbre);
    dens_arbre.Normalize();

    // int rayon_arbre = arbre.get_rayon();
    // std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    // std::cout<<rayon_arbre<<std::endl;
    // std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

    SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

    std::vector< std::pair< std::pair<int,int> , int > >  list_arbre;


    std::chrono::high_resolution_clock::time_point a= std::chrono::high_resolution_clock::now();
    
    list_arbre = make_dart_throwing_quicker(arbre);

 
    std::chrono::high_resolution_clock::time_point b= std::chrono::high_resolution_clock::now();
    
    // mesurer la difference, et l'exprimer en microsecondes 
    unsigned int time= std::chrono::duration_cast<std::chrono::microseconds>(b - a).count();

    std::cout<<"Temps Dart Throwing : "<<time<<" ms"<<std::endl;
    
    // for (int k=0; k<100000; k++){
    //     int rand_pos_x = rand()%nx;
    //     int rand_pos_y = rand()%ny;


    //     bool placement_not_possible = false;


    //     for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
    //         if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre) , *it ) ){
    //             placement_not_possible = true;
    //             break;
    //         }
    //     }


    //     if (placement_not_possible == false){
    //         list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre));
    //     }
    // }

    for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        

        bool test_dens_sapin = true;


        float rand_test = ((double) rand() / (RAND_MAX));
        int rand_pos_x = (*it).first.first;
        int rand_pos_y = (*it).first.second;

        if (rand_test <= dens_arbre.at(rand_pos_x,rand_pos_y)){
            test_dens_sapin=true;
        }
        else{
            continue;
        }

        if (test_dens_sapin == true){
            res.at(rand_pos_x, rand_pos_y) = 1;
            // res.at(rand_pos_x+2, rand_pos_y) = 1;
            // res.at(rand_pos_x-2, rand_pos_y) = 1;
            // res.at(rand_pos_x, rand_pos_y+2) = 1;
            // res.at(rand_pos_x, rand_pos_y-2) = 1;
            // res.at(rand_pos_x+1, rand_pos_y) = 1;
            // res.at(rand_pos_x-1, rand_pos_y) = 1;
            // res.at(rand_pos_x, rand_pos_y+1) = 1;
            // res.at(rand_pos_x, rand_pos_y-1) = 1;
        }


    }
    if (dil){
        res.UpdateMinMax();
        res.Dilate(arbre.get_rayon(), res.max());

    }
    
    return res;
};


SF2 HeighField::double_raw_distribution(Arbre& arbre1, Arbre& arbre2) const{

    
    SF2 dens_arbre1 = densite_arbre(arbre1);
    dens_arbre1.Normalize();

    SF2 dens_arbre2 = densite_arbre(arbre2);
    dens_arbre2.Normalize();

    int rayon_arbre1 = arbre1.get_rayon();
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    std::cout<<rayon_arbre1<<std::endl;
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

    int rayon_arbre2 = arbre2.get_rayon();
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    std::cout<<rayon_arbre2<<std::endl;
    std::cout<<"rayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbrerayon_arbre"<<std::endl;
    

    SF2 res = SF2(Grid2(Box2(a,b),nx,ny),0.0); //init 0

    std::vector< std::pair< std::pair<int,int> , int > >  list_arbre;

    for (int k=0; k<100000; k++){
        int rand_pos_x = rand()%nx;
        int rand_pos_y = rand()%ny;

        float rand_test = ((double) rand() / (RAND_MAX)); //it fix memory issue !!??

        bool placement_not_possible = false;


        for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
            if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre1) , *it ) ){
                placement_not_possible = true;
                break;
            }
        }
        if (placement_not_possible == false){
            list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre1));
        }
    }

    for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        
        bool test_dens_arbre1 = false;

        float rand_test = ((double) rand() / (RAND_MAX));
        int rand_pos_x = (*it).first.first;
        int rand_pos_y = (*it).first.second;


        if (rand_test <= dens_arbre1.at(rand_pos_x,rand_pos_y)){
            test_dens_arbre1=true;
        }
        else{
            continue;
        }


        if (test_dens_arbre1 == true){
            res.at(rand_pos_x, rand_pos_y) = 1;
            // res.at(rand_pos_x+2, rand_pos_y) = 1;
            // res.at(rand_pos_x-2, rand_pos_y) = 1;
            // res.at(rand_pos_x, rand_pos_y+2) = 1;
            // res.at(rand_pos_x, rand_pos_y-2) = 1;
            // res.at(rand_pos_x+1, rand_pos_y) = 1;
            // res.at(rand_pos_x-1, rand_pos_y) = 1;
            // res.at(rand_pos_x, rand_pos_y+1) = 1;
            // res.at(rand_pos_x, rand_pos_y-1) = 1;
        }

    }



    for (int k=0; k<100000; k++){
        int rand_pos_x = rand()%nx;
        int rand_pos_y = rand()%ny;

        float rand_test = ((double) rand() / (RAND_MAX)); //it fix memory issue !!??

        bool placement_not_possible = false;


        for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
            if ( test_dist( std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre2) , *it ) ){
                placement_not_possible = true;
                break;
            }
        }
        if (placement_not_possible == false){
            list_arbre.push_back(std::pair<std::pair<int,int>,int>(std::pair<int,int>(rand_pos_x,rand_pos_y), rayon_arbre2));
        }
    }

    for (std::vector< std::pair< std::pair<int,int> , int > >::iterator it = list_arbre.begin(); it != list_arbre.end(); it++){
        
        bool test_dens_arbre2 = false;
        int rand_pos_x = (*it).first.first;
        int rand_pos_y = (*it).first.second;

        float rand_test = ((double) rand() / (RAND_MAX));

        if (rand_test <= dens_arbre2.at(rand_pos_x,rand_pos_y)){
            test_dens_arbre2=true;
        }
        else{
            continue;
        }

        if (test_dens_arbre2 == true){
            res.at(rand_pos_x, rand_pos_y) = 0.5;
            // res.at(rand_pos_x+2, rand_pos_y+2) = 1;
            // res.at(rand_pos_x-2, rand_pos_y-2) = 1;
            // res.at(rand_pos_x-2, rand_pos_y+2) = 1;
            // res.at(rand_pos_x+2, rand_pos_y-2) = 1;
            // res.at(rand_pos_x+1, rand_pos_y+1) = 1;
            // res.at(rand_pos_x-1, rand_pos_y-1) = 1;
            // res.at(rand_pos_x-1, rand_pos_y+1) = 1;
            // res.at(rand_pos_x-1, rand_pos_y-1) = 1;
        }

    }
    
    res.Dilate(arbre1.get_rayon(), 1.0);
    res.Dilate(arbre2.get_rayon(), 0.5);

    return res;
}


SF2 HeighField::double_raw_distribution_quicker(Arbre& arbre1, Arbre& arbre2) const{

    SF2 sapinDist = raw_distribution(arbre1, false);
    SF2 buissonDist = raw_distribution(arbre2, false);

    SF2 avoidArea = sapinDist;
    avoidArea.Dilate(arbre1.get_rayon()+arbre2.get_rayon(), 1.0);
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny ; j++){
            if ((buissonDist.at(i,j) == 1) && (avoidArea.at(i,j) == 0)){
                sapinDist.at(i,j) = 0.5;
            }
        }
    }

    sapinDist.Dilate(arbre1.get_rayon(), 1.0);
    sapinDist.Dilate(arbre2.get_rayon(), 0.5);


    return sapinDist;
}