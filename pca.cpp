#include "pca.h"

PCA::PCA(Mat A, int k)
{
    float sum = 0;
    for(int i = 1; i <= A.getLines(); i++){
        for(int j = 1; j <= A.getColumns(); j++){
            sum = sum + A(i,j);
        }
    }
    float avg = sum/(float)(A.getLines()*A.getColumns());
    std::cout << "\navg:" << avg << "\n";
    for(int i = 1; i <= A.getLines(); i++){
        for(int j = 1; j <= A.getColumns(); j++){
            A(i,j) = A(i,j)-avg;
        }
    }
    SVDec svd = SVDec(A);
    Mat U2 = svd.getU().slice(1,svd.getU().getLines(),1,k);
    std:: cout << "\nU2:\n" << U2;
    Mat S2 = svd.getSigma().slice(1,k,1,k);
    std:: cout << "\nS2:\n" << S2;
    T = U2*S2;
}

Mat PCA::getT(){
    return T;
}
