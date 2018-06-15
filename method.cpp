#include "method.h"
#include "LUDec.h"
#include "Mat.h"
#include "qrdec.h"
#include "eigen.h"

Method::Method()
{

}

Mat Method::gaussSubstitution(Mat a, Mat b)
{
    Mat x = Mat(a.getLines(),1,false);

    for (int i = a.getLines(); i >= 1; i--){
        double sum = 0;
        for (int j = i+1; j <= a.getColumns(); j++){
            sum = sum + a(i,j)*x(j,1);
        }
        x(i,1) = ((b(i,1) - sum)/a(i,i));
    }
    return x;
}

Mat Method::gaussRetroSubstitution(Mat a, Mat b)
{
    Mat x = Mat(a.getLines(),1,false);

    for (int i = 1; i <= a.getLines(); i++){
        double sum = 0;
        for (int j = i-1; j >= 1; j--){
            sum = sum + a(i,j)*x(j,1);
        }
        x(i,1) = ((b(i,1) - sum)/a(i,i));
    }
    return x;
}

Mat Method::gaussElimination(Mat a, Mat b)
{
    Mat a2 = Mat(a);
    Mat b2 = Mat(b);

    for(int j=1;j<=a2.getColumns()-1;j++){
        for(int k=j+1;k<=a2.getLines();k++){
            if(a2(j,j) == 0){
                std::cout << "Divisao por Zero Detectada, Eliminacao de Gauss Interrompida.";
                exit(0);
            }
            double alpha = -1*((a2(k,j))/(a2(j,j)));
            a2(k,j) = 0;
            for(int l=j+1; l<=a2.getColumns();l++){
                a2(k,l) = a2(k,l)+(alpha * a2(j,l));
            }
            b2(k,1) = b2(k,1)+(alpha * b2(j,1));
        }
    }
    std::cout << "\nMatriz A na forma escalonada: \n" << a2 << "\nMatriz b na forma escalonada: \n" << b2 << "\n";

    return gaussSubstitution(a2,b2);
}

Mat Method::gaussEliminationPivoting(Mat a, Mat b)
{
    Mat a2 = Mat(a);
    Mat b2 = Mat(b);

    for(int j=1;j<=a2.getColumns()-1;j++){
        for(int k=j+1;k<=a2.getLines();k++){
            if((a2(j,j) == 0)&&(k < a2.getLines())){
                for(int l = 1; l <=a2.getColumns(); l++){
                    float temp = a2(j,l);
                    a2(j,l) = a2(j+1,l);
                    a2(j+1,l) = temp;
                }
            }
            double alpha = -1*((a2(k,j))/(a2(j,j)));
            a2(k,j) = 0;
            for(int l=j+1; l<=a2.getColumns();l++){
                a2(k,l) = a2(k,l)+(alpha * a2(j,l));
            }
            b2(k,1) = b2(k,1)+(alpha * b2(j,1));
        }
    }
    std::cout << "\nMatriz A na forma escalonada: \n" << a2 << "\nMatriz b na forma escalonada: \n" << b2 << "\n";

    return gaussSubstitution(a2,b2);
}

Mat Method::decompositionCholesky(Mat a, Mat b)
{
    Mat A2 = Mat(a);
    Mat B2 = Mat(b);

    Mat S = Mat(A2.getLines(), A2.getColumns(), false);
    //Calculando S

    for(int i = 1; i <= S.getLines(); i++){
        for (int j = 1; j < i; j++){
            double sum = 0;
            for (int k = 1; k < j; k++){
                sum += S(i,k)*S(j,k);
            }
            if(S(j,j) == 0){
                std::cout << "\nDivisao por zero detectada, matriz impropria.\n";
                exit(0);
            } else {
                S(i,j) = ((A2(i,j) - sum))/S(j,j);
            }
        }
        double sum = 0;
        for (int k = 1; k < i; k++){
            sum += S(i,k)*S(i,k);
        }
        if((A2(i,i) - sum) < 0){
            std::cout << "\nRaiz de negativo detectada, matriz impropria\n";
            exit(0);
        } else {
            S(i,i) = sqrt(A2(i,i) - sum);
        }
    }

    std::cout << "\nMatriz S:\n" << S << "\nMatriz S transposta:\n" << Mat::transpose(S);

    //Resolvendo S*y = b
    Mat y = gaussRetroSubstitution(S,B2);

    //Resolvendo ST*x = y;
    return gaussSubstitution(Mat::transpose(S),y);
}

Mat Method::decompositionLU(Mat a, Mat b)
{
    Mat A2 = Mat(a);
    Mat B2 = Mat(b);

    //Construindo L inicial como matriz identidae
    LUDec lu_dec = LUDec(A2);

    Mat L = lu_dec.GetL();
    Mat U = lu_dec.GetU();

    std::cout << "\nMatriz L:\n" << L << "\nMatriz U:\n" << U;

    //Resolvendo L*y = b
    Mat y = gaussRetroSubstitution(L,B2);

    //Resolvendo U*x = y
    return gaussSubstitution(U,y);
}

Mat Method::quadraticMinimum(Mat a, Mat b)
{
    Mat a_t = Mat::transpose(a);

    Mat a2 = a_t*a;

    Mat b2 = a_t*b;

    std::cout << "\nMatriz (A^T)A:\n" << a2 << "\nMatriz (A^T)b:\n" << b2;

    return gaussEliminationPivoting(a2,b2);
}

Mat Method::ortoGramSchmidt(Mat a, Mat b)
{
    QRDec QR;
    QR.QRGramSchmidt(a);

    std::cout << "\nMatriz Q:\n" << QR.GetQ() << "\nMatriz R:\n" << QR.GetR();

    return quadraticMinimum(QR.GetR(), Mat::transpose(QR.GetQ())*b);
}

Mat Method::ortoHouseHolder(Mat a, Mat b)
{
    QRDec QR;
    QR.QRHouseHolder(a);

    std::cout << "\nMatriz Q:\n" << QR.GetQ() << "\nMatriz R:\n" << QR.GetR();

    return quadraticMinimum(QR.GetR(),QR.GetQ()*b);
}

Mat Method::ortoJacobi(Mat a, Mat b)
{
    QRDec QR;
    QR.QRJacobi(a);

    std::cout << "\nMatriz Q:\n" << QR.GetQ() << "\nMatriz R:\n" << QR.GetR();

    return quadraticMinimum(QR.GetR(),QR.GetQ()*b);
}

Mat Method::jacobiIteration(Mat A, Mat b, Mat x)
{
    Mat x_new = Mat(A.getLines(), 1, false);
    for(int i = 1; i <= A.getLines(); i++){
        float sum = 0;
        for(int j = 1; j <= A.getColumns(); j++){
            if(i!=j){
                sum = sum + A(i,j)*x(j,1);
            }
        }
        x_new(i,1) = (b(i,1)-sum)/A(i,i);
    }
    return x_new;
}

Mat Method::iterativeJacobi(Mat A, Mat b, float tolerance)
{
    Mat x_before, x_after;
    x_before = Mat(A.getLines(),1,false);
    std::cout << "x inicial:\n" << x_before;
    int n_iterations = 1000;
    float error;
    while (n_iterations > 0){
        x_after = jacobiIteration(A,b,x_before);
        error = Mat::normQuadratic(x_after - x_before);
        if(error < tolerance) break;
        n_iterations--;
        x_before = x_after;
    }
    if(n_iterations == 0){
        std::cout << "\n O resultado divergiu do valor correto\n" << x_after;
    }
    std::cout << "\nx apos " << 1000 - n_iterations << " iteracoes:\n" << x_after;
    return x_after;
}

Mat Method::gaussSeidelIteration(Mat A, Mat b, Mat x)
{
    Mat x_new = Mat(A.getLines(), 1, false);
    for(int i = 1; i <= A.getLines(); i++){
        float sum = 0;
        for(int j = 1; j <= A.getColumns(); j++){
            if(i > j){
                sum = sum - A(i,j)*x_new(j,1);
            }
            if(i < j){
                sum = sum - A(i,j)*x(j,1);
            }
        }
        x_new(i,1) = (sum+b(i,1))/A(i,i);
    }
    return x_new;
}

Mat Method::iterativeGaussSeidel(Mat A, Mat b, float tolerance)
{
    Mat x_before, x_after;
    x_before = Mat(A.getLines(),1,false);
    std::cout << "x inicial:\n" << x_before;
    int n_iterations = 1000;
    float error;
    while (n_iterations > 0){
        x_after = gaussSeidelIteration(A,b,x_before);
        error = Mat::normQuadratic(x_after - x_before);
        if(error < tolerance) break;
        n_iterations--;
        x_before = x_after;
    }
    if(n_iterations == 0){
        std::cout << "\n O resultado divergiu do valor correto.\n" << x_after;
    }
    std::cout << "\nx apos " << 1000 - n_iterations << " iteracoes:\n" << x_after;
    return x_after;
}

Mat Method::sucessiveOverRelaxation(Mat A, Mat b, float omega, float tolerance)
{
    Mat x_before, x_after;
    x_before = Mat(A.getLines(),1,false);
    std::cout << "x inicial:\n" << x_before;
    int n_iterations = 1000;
    float error;
    while (n_iterations > 0){
        x_after = x_before*(1-omega) + gaussSeidelIteration(A,b,x_before)*omega;
        error = Mat::normQuadratic(x_after - x_before);
        if(error < tolerance) break;
        n_iterations--;
        x_before = x_after;
    }
    if(n_iterations == 0){
        std::cout << "\n O resultado divergiu do valor correto.\n" << x_after;
    }
    std::cout << "\nx apos " << 1000 - n_iterations << " iteracoes:\n" << x_after;
    return x_after;
}

Mat Method::steepestDescent(Mat A, Mat b, float tolerance)
{
    Mat x_before, x_after,r;
    x_before = Mat(A.getLines(),1,false);
    std::cout << "x inicial:\n" << x_before;
    int n_iterations = 1000;
    float error;
    double alpha = 0;
    if(A==Mat::transpose(A)){}
    else {
        A = Mat::transpose(A)*A;
        b = Mat::transpose(A)*b;
    }
    while (n_iterations > 0){
        r = b - A*x_before;
        alpha = r%r/(r%(A*r));
        x_after = x_before + r*alpha;
        error = Mat::normQuadratic(r);
        if(error < tolerance) break;
        n_iterations--;
        x_before = x_after;
    }
    if(n_iterations == 0){
        std::cout << "\n O resultado divergiu do valor correto.\n" << x_after;
    }
    std::cout << "\nx apos " << 1000 - n_iterations << " iteracoes:\n" << x_after;
    return x_after;
}

Eigen Method::powerIterationMethod(Mat A, Mat x, double e)
{
    Mat q = x/Mat::normQuadratic(x);
    double last_eigen_value = 0;
    double new_eigen_value = 0;
    double error = 1;
    while (error > e){
        last_eigen_value = new_eigen_value;
        x = A*q;
        q = x/Mat::normQuadratic(x);
        new_eigen_value = (q%(A*q))/(q%q);
        error = (new_eigen_value - last_eigen_value);
    }

    std::cout << "\nMaior autovalor e autovetor correspondente:\n" << new_eigen_value << "\n\n" << q;

    Eigen result = Eigen(q,new_eigen_value);
    return result;
}

Eigen Method::inverseIterationMethod(Mat A, Mat x, double e)
{
    x = x/Mat::normQuadratic(x);
    double last_eigen_value = 0;
    double new_eigen_value = 0;
    double error = 1;
    LUDec LU = LUDec(A);
    Mat L = LU.GetL();
    Mat U = LU.GetU();
    Mat z = gaussRetroSubstitution(L,x);
    Mat y = gaussSubstitution(U,z);
    while (error > e){
        last_eigen_value = new_eigen_value;
        x = y/Mat::normQuadratic(y);
        z = gaussRetroSubstitution(L,x);
        y = gaussSubstitution(U,z);
        new_eigen_value = (x%y)/(x%x);
        error = (new_eigen_value - last_eigen_value);
    }

    std::cout << "\nMenor autovalor e autovetor correspondente:\n" << 1/new_eigen_value << "\n\n" << x;

    Eigen result = Eigen(x,1/new_eigen_value);
    return result;
}

Eigen Method::shiftedIterationMethod(Mat A, Mat x, double mi, double e)
{
    Mat I = Mat(A.getLines(),A.getColumns());
    Eigen result = inverseIterationMethod(A-(I*mi),x,e);
    result.setEigenvalue(result.getEigenvalue()+mi);

    std::cout << "\nAutovalor e autovetor correspondente:\n" << result.getEigenvalue() << "\n\n" << result.getEigenmatrix();

    return result;
}

std::vector<Eigen> Method::similarityTransformation(Mat A)
{
    Mat Aj = A;
    Mat X = Mat(A.getLines(),A.getColumns());

    //Houselholder Transformation
    Mat H = Mat(A.getLines(),A.getColumns());
    for (int j = 1; j <= A.getColumns()-2; j++){
        Mat Vj = Mat(A.getLines(),1,false);
        for (int k = j+1; k <= A.getLines();k++){
            Vj(k,1) = Aj(k,j);
        }
        double nVj = Mat::normQuadratic(Vj);
        Mat Vl = Mat(A.getLines(),1,false);
        if(Vj(j+1,1) > 0){
            Vl(j+1,1) = -1*nVj;
        } else {
            Vl(j+1,1) = nVj;
        }
        Mat N = Vj - Vl;
        Mat n = N/Mat::normQuadratic(N);
        Mat I = Mat(A.getLines(),A.getColumns());
        Mat Hj = I - (n*Mat::transpose(n))*2;
        Aj = Mat::transpose(Hj)*Aj*Hj;
        H = H*Hj;
    }
    X = X*H;

    bool terminado = false;
    while (!terminado){
        Mat Qt = Mat(A.getLines(),A.getColumns());
        for(int j = 1; j <= A.getColumns()-1; j++){
            Mat J = QRDec::constructJij(Aj,j+1,j);
            Aj = J*Aj;
            Qt = J*Qt;
        }
        Mat R = Aj;
        Mat Q = Mat::transpose(Qt);
        Aj = R*Q;
        X = X*Q;
        double sum = 0;
        for(int i = 1; i <= Aj.getLines(); i++){
            for(int j = 1; j <= Aj.getColumns(); j++){
                if(i>j){
                    sum = sum + fabs(Aj(i,j));
                }
            }
        }
        if(sum <= 0.00000001){
            terminado = true;
        }
    }
    std::vector<Eigen> result;
    std::cout << "\nAutovalores e autovetores:\n";
    for(int j = 1; j <= A.getColumns(); j++){
        Eigen e = Eigen(X.getColumn(j), Aj(j,j));
        std::cout << e.getEigenvalue() << "\n\n" << e.getEigenmatrix() << "\n";
        result.push_back(e);
    }
    return result;
}
