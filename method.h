#ifndef METHOD_H
#define METHOD_H
#include "Mat.h"
#include "eigen.h"
#include "vector"


class Method
{
public:
    Method();
    static Mat gaussSubstitution(Mat a, Mat b);
    static Mat gaussRetroSubstitution(Mat a, Mat b);
    static Mat gaussElimination(Mat a, Mat b);
    static Mat gaussEliminationPivoting(Mat a, Mat b);
    static Mat decompositionCholesky(Mat a, Mat b);
    static Mat decompositionLU(Mat a, Mat b);
    static Mat quadraticMinimum(Mat a, Mat b);
    static Mat ortoGramSchmidt(Mat a, Mat b);
    static Mat ortoHouseHolder(Mat a, Mat b);
    static Mat ortoJacobi(Mat a, Mat b);
    static Mat jacobiIteration(Mat A, Mat b, Mat x);
    static Mat iterativeJacobi(Mat A, Mat b, float tolerance);
    static Mat gaussSeidelIteration(Mat A, Mat b, Mat x);
    static Mat iterativeGaussSeidel(Mat A, Mat b, float tolerance);
    static Mat sucessiveOverRelaxation(Mat A, Mat b, float omega, float tolerance);
    static Mat steepestDescent(Mat A, Mat b, float tolerance);
    static Eigen powerIterationMethod (Mat A, Mat x, double e);
    static Eigen inverseIterationMethod (Mat A, Mat x, double e);
    static Eigen shiftedIterationMethod (Mat A, Mat x, double mi, double e);
    static std::vector<Eigen> similarityTransformation(Mat A);
};

#endif // METHOD_H
