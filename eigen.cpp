#include "eigen.h"

Eigen::Eigen()
{

}

Eigen::Eigen(Mat e, double ev)
{
    eigenmatrix = e;
    eigenvalue = ev;
}

Mat Eigen::getEigenmatrix() const
{
    return eigenmatrix;
}

void Eigen::setEigenmatrix(const Mat &value)
{
    eigenmatrix = value;
}

double Eigen::getEigenvalue() const
{
    return eigenvalue;
}

void Eigen::setEigenvalue(double value)
{
    eigenvalue = value;
}
