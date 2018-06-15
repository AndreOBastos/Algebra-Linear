#ifndef PCA_H
#define PCA_H
#include "svdec.h"

class PCA
{
public:
    PCA(Mat A, int k);
    Mat getT();
private:
    Mat T;
};

#endif // PCA_H
