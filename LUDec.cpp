#include "LUDec.h"
#include "Mat.h"

LUDec::LUDec()
{
    //ctor
}
LUDec::LUDec(Mat dec){
    L = Mat(dec.getLines(), dec.getColumns());
    U = Mat(dec.getLines(), dec.getColumns(), false);
    for(int j = 1; j <= dec.getColumns(); j++){
        for (int i = 1; i <= j; i++){  // Gerando U
            double sum = 0;
            for (int s = 1; s <= i-1; s++){
                sum = sum + L(i,s) * U(s,j);
            }
            U(j,i) = 0;
            U(i,j) = dec(i,j) - sum;
        }

        for (int i = j+1; i <= dec.getColumns(); i++){ //Gerando L
            int sum = 0;
            for (int s = 1; s <= j-1; s++){
                sum = sum + L(i,s) * U(s,j);
            }
            if(U(j,j) == 0){
                std::cout << "Divisao por zero detectada, parando o programa.\n";
                exit(0);
            }
            L(j,i) = 0;
            L(i,j) = (dec(i,j) - sum)/U(j,j);
        }
    }
}

LUDec::~LUDec()
{
    //dtor
}
