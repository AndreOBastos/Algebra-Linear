#ifndef LUDEC_H
#define LUDEC_H
#include "Mat.h"

class LUDec
{
    public:
        LUDec();
        LUDec(Mat dec);
        virtual ~LUDec();

        Mat GetL() { return L; }
        void SetL(Mat val) { L = val; }
        Mat GetU() { return U; }
        void SetU(Mat val) { U = val; }

    protected:

    private:
        Mat L;
        Mat U;
};

#endif // LUDEC_H
