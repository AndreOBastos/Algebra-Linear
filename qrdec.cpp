#include "qrdec.h"
#include "Mat.h"

QRDec::QRDec()
{

}

QRDec QRDec::QRGramSchmidt(Mat A)
{
    Q = Mat(A.getLines(), A.getColumns(), false);
    for (int i = 1; i <= A.getColumns(); i++){
        Mat sum = Mat(A.getLines(),1,false);
        for (int k = 1; k <= i-1; k++){
            sum = sum + (Q.getColumn(k))*(A.getColumn(i)%Q.getColumn(k));
        }

        Mat v = A.getColumn(i) - sum;
        v = v/Mat::normQuadratic(v);

        for (int j = 1; j <= Q.getLines(); j++){
            Q(j,i) = v(j,1);
        }
    }

    R = Mat::transpose(Q)*A;

    return *this;
}

Mat QRDec::constructHj(Mat Aj, int j)
{
    Mat H = Mat(Aj.getLines(), Aj.getLines());
    Mat V = Mat(Aj.getLines(),1,false);
    for (int i = j; i <= Aj.getLines(); i++){ V(i,1) = Aj(i,j); }
    double normaV = Mat::normQuadratic(V);
    Mat V2 = Mat(Aj.getLines(),1,false);
    if (V(j,1) >= 0){
        V2(j,1) = -1*normaV;
    } else {
        V2(j,1) = normaV;
    }
    Mat N = V-V2;
    Mat n = N/Mat::normQuadratic(N);
    Mat r = H - n*Mat::transpose(n)*2;
    return r;
}

QRDec QRDec::QRHouseHolder(Mat A)
{
    R = A;
    Q = Mat(A.getLines(),A.getLines());
    for (int j = 1; j <= R.getColumns(); j++){
        Mat Hj = constructHj(R,j);
        R = Hj*R;
        Q = Hj*Q;
    }

    return *this;
}

Mat QRDec::constructJij(Mat Aj, int i, int j)
{
    double angle;
    if(Aj(j,j) != 0){
        angle = atan(Aj(i,j)/Aj(j,j));
    } else {
        angle = 3.14159265358979323846/2;
    }

    Mat J = Mat(Aj.getLines(),Aj.getLines());
    J(i,i) = cos(angle);
    J(i,j) = -1*sin(angle);
    J(j,i) = sin(angle);
    J(j,j) = cos(angle);

    return J;
}

QRDec QRDec::QRJacobi(Mat A)
{
    R = A;
    Q = Mat(A.getLines(), A.getLines());

        for (int j = 1; j <= A.getColumns(); j++){
            for (int i = A.getLines(); i >=j+1; i--){
                Mat J = constructJij(R,i,j);
                R = J*R;
                Q = J*Q;
            }
        }


    return *this;
}
