#ifndef QRDEC_H
#define QRDEC_H
#include "Mat.h"


class QRDec
{
public:
    QRDec();
    QRDec QRGramSchmidt(Mat A);
    static Mat constructHj(Mat Aj, int j);
    QRDec QRHouseHolder(Mat A);
    static Mat constructJij(Mat Aj, int i, int j);
    QRDec QRJacobi(Mat A);

    Mat GetQ() { return Q; }
    void SetQ(Mat val) { Q = val; }
    Mat GetR() { return R; }
    void SetR(Mat val) { R = val; }

private:
    Mat Q;
    Mat R;
};

#endif // QRDEC_H
