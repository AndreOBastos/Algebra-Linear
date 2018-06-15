#ifndef MAT_H
#define MAT_H
#include <iostream>
#include <new>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <time.h>

class Mat
{
    public:

        //Sem entrada, Mat constroi uma matriz de um valor só igual a 0.
        Mat() {
            lines = 1;
            columns = 1;
            data = new double[1];
            data[0] = 0;
        }

        //Construtor que copia os valores de uma matriz para outra.
        Mat(const Mat& other){
            lines = other.getLines();
            columns = other.getColumns();
            data = other.getData();
        }

        //Construtor que, dada linhas e colunas, retorna a identidade.
        Mat(int l,  int c){
            lines = l;
            columns = c;
            data = new double[l*c];
            int line = 1;
            int column = 1;
            for (int i = 0; i < l*c; i++){
                if(line == column){
                    data[i] = 1;
                }else{
                    data[i] = 0;
                }
                if (column<c){
                    column++;
                } else {
                    column = 1;
                    line++;
                }
            }
        }

        //Construtor onde, se random for true, retorna uma matriz com valores aleatórios, e se for falso, retorna uma matriz com todos os valores iguais a 0.
        Mat( int l,  int c, bool random){
            lines = l;
            columns = c;
            data = new double[l*c];
            if(random){
                srand(time(NULL));
                for (int i = 0; i < l*c; i++){
                    data[i] = rand()%100;
                }
            } else {
                for (int i = 0; i < l*c; i++){
                    data[i] = 0;
                }
            }
        }

        //Construtor que gera a matriz que você quiser.
        Mat( int l,  int c, double d[]){
            lines = l;
            columns = c;
            data = d;
        }

        //Funções básicas de get and set

         int getLines() const { return lines; }
        void setLines( int val) { lines = val; }
         int getColumns() const { return columns; }
        void setColumns( int val) { columns = val; }
        double* getData() const { return data; }


        //Funções básicas para acessar e modificar valores da matriz

        double getValue( int l,  int c) const{
            return data[((l-1)*columns) + (c-1)];
        }

        void setValue( int l,  int c, double val) const{
            std::memcpy(&data[((l-1)*columns) + (c-1)], &val, sizeof(double));
        }

        //Função para imprimir a matriz

        void toString(){
            for (int i = 1; i <= lines; i++){
                for (int j = 1; j <= columns; j++){
                    std::cout  << getValue(i,j) << '\t';
                }
                std::cout << '\n';
            }
            std::cout << '\n';
        }

        //Funções que retornam a matriz coluna e matriz linha que você quiser de uma matriz.

        Mat getColumn(int c){
            Mat b = Mat(this->getLines(), 1, false);
            for (int i = 1; i <= this->getLines(); i++){
                b(i,1) = this->getValue(i,c);
            }
            return b;
        }

        Mat getLine(int l){
            Mat b = Mat(1, this->getColumns(), false);
            for (int i = 1; i <= this->getColumns(); i++){
                b(1,i) = this->getValue(l,i);
            }
            return b;
        }

        //Soma de matrizes

        Mat operator+(const Mat& b){
            if((this->getLines()!=b.getLines())||(this->getColumns()!=b.getColumns())){
                std::cout << "Erro, soma indevida de matrizes.";
                exit(0);
            }
            Mat c = Mat(this->getLines(), this->getColumns(), false);
            for (int i = 1; i <= this->getLines(); i++){
                for (int j = 1; j <= this->getColumns(); j++){
                    c(i,j) = this->getValue(i,j) + b(i,j);
                }
            }
            return c;
        }

        //Subtração de matrizes

        Mat operator-(const Mat& b){
            if((this->getLines()!=b.getLines())||(this->getColumns()!=b.getColumns())){
                std::cout << "Erro, subtração indevida de matrizes.";
                exit(0);
            }
            Mat c = Mat(this->getLines(), this->getColumns(), false);
            for (int i = 1; i <= this->getLines(); i++){
                for (int j = 1; j <= this->getColumns(); j++){
                    c(i,j) = this->getValue(i,j) - b(i,j);
                }
            }
            return c;
        }

        //Produto de matriz por uma escalar

        Mat operator*(const double& b){
            Mat c = Mat(this->getLines(), this->getColumns(), false);
            for (int i = 1; i <= this->getLines(); i++){
                for (int j = 1; j <= this->getColumns(); j++){
                    c(i,j) = this->getValue(i,j)*b;
                }
            }
            return c;
        }

        //Produto entre matrizes

        Mat operator*(const Mat& b){
            if(this->getColumns()!= b.getLines()){
                std::cout << "Erro, multiplicação indevida de matrizes.";
                exit(0);
            }
            Mat c = Mat(this->getLines(), b.getColumns(), false);
            for (int i = 1; i <= this->getLines(); i++){
                for (int j = 1; j <= b.getColumns(); j++){
                    for (int k = 1; k <= this->getColumns(); k++){
                        c(i,j) = c(i,j) + this->getValue(i,k)*b(k,j);
                    }
                }
            }
            return c;
        }

        Mat operator/(const double& b){
            Mat c = Mat(this->getLines(), this->getColumns(), false);
            for (int i = 1; i <= this->getLines(); i++){
                for (int j = 1; j <= this->getColumns(); j++){
                    c(i,j) = this->getValue(i,j)/b;
                }
            }
            return c;
        }

        bool operator==(const Mat&b){
            if((this->getLines()!=b.getLines())||(this->getColumns()!=b.getColumns())){
                return false;
            }
            for(int i = 1; i <= this->getLines(); i++){
                for(int j = 1; j <= this->getColumns(); j++){
                    if(this->getValue(i,j) != b(i,j)){
                        return false;
                    }
                }
            }
            return true;
        }

        //Produto escalar entre matrizes

        double operator%(const Mat& b){
            double sum = 0;
            for (int i = 1; i <= this->getLines(); i++){
                sum = sum + (this->getValue(i,1)*b(i,1));
            }
            return sum;
        }

        Mat operator=(const Mat& b){
            if(this!= &b){

                this->setLines(b.getLines());
                this->setColumns(b.getColumns());
                this->data = b.getData();
            }
            return *this;
        }

        //Funções para acesar facilmente membros de uma matriz.

        double& operator()(int i, int j){
            return data[((i-1)*this->getColumns()) + (j-1)];
        }
        const double& operator()(int i, int j) const{
            return data[((i-1)*this->getColumns()) + (j-1)];
        }

        friend std::ostream& operator<<(std::ostream& os, const Mat a){
            for (int i = 1; i <= a.getLines(); i++){
                for (int j = 1; j <= a.getColumns(); j++){
                    os  << a.getValue(i,j) << '\t';
                }
                os << '\n';
            }
            os << '\n';

            return os;
        };

        //Funçao que retorna a transposta de uma matriz

        static Mat transpose(Mat a){
            Mat b = Mat(a.getColumns(), a.getLines(),false);
            for (int i = 1; i <=a.getColumns(); i++){
                for (int j = 1; j <=a.getLines(); j++){
                    b(i,j) = a(j,i);
                }
            }

            return b;
        }

        //Funções para a norma 1 e norma 2 de uma matriz

        static double norm(Mat a){
            double sum = 0;
            for (int i = 1; i <= a.lines; i++){
                sum = sum + a(i,1);
            }
            return sum;
        }

        static double normQuadratic(Mat a){
            double sum = 0;
            for (int i = 1; i <= a.lines; i++){
                sum = sum + (a(i,1) * a(i,1));
            }
            return sqrt(sum);
        }

        //Função que une duas matrizes lado a lado.

        static Mat concatenate(Mat a, Mat b){
            double* c_data = new double[a.getLines()*(a.getColumns()+b.getColumns())];

            for (int i = 1; i <= a.getLines(); i++){
                std::copy(a.getData() + (i-1)*(a.getColumns()),
                          a.getData() + (i-1)*(a.getColumns()) + (a.getColumns()),
                          c_data + ((i-1)*(a.getColumns()+b.getColumns())));
                std::copy(b.getData() + (i-1)*(a.getColumns()),
                          b.getData() + (i-1)*(a.getColumns()) + (a.getColumns()),
                          c_data + ((i-1)*(a.getColumns()+b.getColumns())) + (a.getColumns()));
            }

            Mat c = Mat(a.getLines(),a.getColumns()+b.getColumns(),c_data);

            return c;
        }

        //Função que retorna a inversa de uma matriz.

        static Mat inverse(Mat a){
            Mat ident = Mat(a.getLines(), a.getColumns());
            Mat A2 = Mat::concatenate(a,ident);

            //Realizando a primeira eliminação de Gauss

            for(int j=1;j<=A2.getColumns()-1;j++){
                for(int k=j+1;k<=A2.getLines();k++){
                    if(A2(j,j) == 0){
                        std::cout << "Divisao por Zero Detectada, Inversao de Gauss Jordan interrompida.";
                        Mat error = Mat(a.getLines()+1, a.getColumns()+1);
                        return(error);
                    }
                    double alpha = -1*((A2(k,j))/(A2(j,j)));
                    A2(k,j) = 0;
                    for(int l=j+1; l<=A2.getColumns();l++){
                        A2(k,l) = A2(k,l)+(alpha * A2(j,l));
                    }
                }
            }

            //Fazendo todos os membros da diagonal igual a 1

            for (int i = 1; i <= A2.getLines(); i++){
                if (A2(i,i) != 1){
                    double pivot = A2(i,i);
                    for (int j = i; j <= (A2.getColumns()/2); j++){
                            A2(i,j) = A2(i,j) / pivot;
                    }
                }
            }

            //Realizando a segunda eliminação de Gauss de baixo para cima

            for(int j = a.getColumns(); j > 1; j--){
                for(int k = j-1; k > 0; k--){
                    double alpha = -1*((A2(k,j))/(A2(j,j)));
                    A2(k,j) = 0;
                    for(int l = j+1; l <= A2.getColumns(); l++){
                        A2(k,l) = A2(k,l) + (alpha * A2(j,l));
                    }
                }
            }

            //Separando a matriz concatenada

            double* inv_data = new double[a.getLines()*a.getColumns()];
            for (int i = 1; i <= A2.getLines(); i++){
                std::copy(A2.getData() + (i-1)*(A2.getColumns()) + (a.getColumns()),
                          A2.getData() + (i-1)*(A2.getColumns()) + 2*(a.getColumns()),
                          inv_data + ((i-1)*(a.getColumns())));
            }

            Mat inverse = Mat(a.getLines(),a.getColumns(),inv_data);

            return inverse;
        }

        Mat slice(int l1, int l2, int c1, int c2){
                    Mat a = Mat(l2-l1+1, c2-c1+1, false);
                    for (int i = l1; i <= l2; i++){
                        for (int j = c1; j <= c2; j++){
                            a(i-l1+1,j-c1+1) = this->getValue(i,j);
                        }
                    }
                    return a;
                }

    protected:

    private:
         int lines;
         int columns;
        double* data;
};

#endif // MAT_H
