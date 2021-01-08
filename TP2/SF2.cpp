#include "SF2.hpp"

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

SF2 SF2::GradientNorm() const{
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
};


QVector<ScalarPoint2> SF2::GetScalarPoints() const
{
    QVector<ScalarPoint2> e(nx * ny);
    int k=0;
    for (int i=0; i < nx ; i++){
        for (int j=0; j < ny ; j++){
            e[k++] = ScalarPoint2(QPoint(i,j), at(i,j));
        }
    }
    return e;
};

void SF2::Dilate(int rayon, double value){

    std::vector<double> fieldMem;
    fieldMem.resize(field.size());
    int ray = 2*rayon;
    if (ray%2 == 0){
        ray+=1;
    }

    for (int i = 0; i < nx ; i++){
        for (int j = 0; j < ny ; j++){
            if (at(i,j) == value){
                for (int k = 0; k < (ray*ray) ; k++){
                    int a,b;

                    a = i - floor(ray/2) + floor(k/(ray));
                    b = j - floor(ray/2) + k%(ray);
                    
                    int diff_x = i-a;
                    int diff_y = j-b;
                    float dist = std::sqrt( diff_x*diff_x + diff_y*diff_y );
  

                    if (dist<rayon){
                        if ((a >= 0) && (b >= 0) && (a < nx) && (b < ny)){
                            fieldMem[Index(a,b)] = at(i,j);
                        }
                    }
                }
            } else if (fieldMem[Index(i,j)] != value){
                fieldMem[Index(i,j)] = at(i,j);
            }
        }
    }

    field = fieldMem;

}