#ifndef SVDEC_H
#define SVDEC_H

#include "Mat.h"
#include "method.h"

class SVDec
{
public:
    SVDec();
    SVDec(Mat A);
    Mat getU() const;
    void setU(const Mat &value);
    Mat getSigma() const;
    void setSigma(const Mat &value);
    Mat getV() const;
    void setV(const Mat &value);

    friend std::ostream& operator<<(std::ostream& os, const SVDec svd){
        os << "\nU:\n" << svd.getU() << "Sigma:\n" << svd.getSigma() << "V:\n" << svd.getV() <<"\n";

        return os;
    };

private:
    Mat U;
    Mat Sigma;
    Mat V;
};

#endif // SVDEC_H
