#ifndef EIGEN_H
#define EIGEN_H
#include "Mat.h"

class Eigen
{
public:
    Eigen();
    Eigen(Mat e, double ev);

    Mat getEigenmatrix() const;
    void setEigenmatrix(const Mat &value);
    double getEigenvalue() const;
    void setEigenvalue(double value);

    friend std::ostream& operator<<(std::ostream& os, const Eigen e){
        os << "\nEigenvector:\n" << e.getEigenmatrix() << "Eigenvalue:\n" << e.getEigenvalue() << '\n';

        return os;
    };

private:
    Mat eigenmatrix;
    double eigenvalue;
};

#endif // EIGEN_H
