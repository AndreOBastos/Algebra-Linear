#include "svdec.h"

SVDec::SVDec()
{

}

SVDec::SVDec(Mat A)
{
    std::vector<Eigen> eigenU = Method::similarityTransformation(A*Mat::transpose(A));
    std::vector<Eigen> eigenV = Method::similarityTransformation(Mat::transpose(A)*A);

    Mat tempU = eigenU[0].getEigenmatrix();
    for(unsigned int i = 1; i < eigenU.size(); i++){
        tempU = Mat::concatenate(tempU,eigenU[i].getEigenmatrix());
    }
    U = tempU;

    Mat tempV = eigenV[0].getEigenmatrix();
    for(unsigned int i = 1; i < eigenV.size(); i++){
        tempV = Mat::concatenate(tempV,eigenV[i].getEigenmatrix());
    }
    V = Mat::transpose(tempV);

    Mat SigmaTemp = Mat(U.getColumns(), V.getLines(), false);
    for(unsigned int i = 0; i < eigenV.size(); i++){
        SigmaTemp(i+1,i+1) = sqrt(eigenV[i].getEigenvalue());
    }
    Sigma = SigmaTemp;

    std::cout << "\nMatriz U:\n" << U << "\nMatriz Sigma:\n" << Sigma << "\nMatriz V^T:\n" << V << "\n";
}

Mat SVDec::getU() const
{
    return U;
}

void SVDec::setU(const Mat &value)
{
    U = value;
}

Mat SVDec::getSigma() const
{
    return Sigma;
}

void SVDec::setSigma(const Mat &value)
{
    Sigma = value;
}

Mat SVDec::getV() const
{
    return V;
}

void SVDec::setV(const Mat &value)
{
    V = value;
}
