#include <iostream>
#include "Mat.h"
#include "method.h"
#include "eigen.h"
#include "vector"
#include "svdec.h"
#include "pca.h"

int main()
{
    double d[] = {1,2,3,2,5,3,1,0,8};
    double f[] = {1,2,3,2,4,5,3,5,6};
    double g[] = {1,0,-1,2,1,-2,1,1,0,1,1,-1};
    double g2[] = {6,0,9,3};
    double cholesky[] = {4,12,-16,12,37,-43,-16,-43,98};
    double iterative[] = {10,-1,2,0,-1,11,-1,3,2,-1,10,-1,0,3,-1,8};
    double iterative2[] = {6,25,-11,15};
    double eigen[] = {6,2,-1,2,5,1,-1,1,4};
    double e[] = {-1,4,0};
    double t[] = {1,1,1};
    Mat m2 = Mat(3,3,d);
    Mat m3 = Mat(3,1,t);
    Mat m4 = Mat(3,3,f);
    Mat m5 = Mat(4,3,g);
    Mat m52 = Mat(4,1,g2);
    Mat c = Mat(3,3,cholesky);
    Mat i = Mat(4,4,iterative);
    Mat i2 = Mat(4,1,iterative2);
    Mat ei = Mat(3,3,eigen);

    Method::similarityTransformation(ei);
    Method::powerIterationMethod(ei,m3,0.000001);
    Method::inverseIterationMethod(ei,m3,0.000001);
    Method::shiftedIterationMethod(ei,m3,4,0.00001);
}
