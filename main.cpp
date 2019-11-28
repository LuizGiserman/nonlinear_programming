#include <math.h>
#include <iostream>
#include <vector>


#define OK                      0
#define ERRO_DIMENSAO_MATRIZ    1
using namespace std;

/**/

double funcao (double x1, double x2)
{
    return (double) ( sin ( pow(x1, 2) - pow(x2, 2) ) * cos ( pow(x1, 2) + pow(x2, 2) ) );
}

/*retorna o valor da derivada de f em x1 no ponto (x1,x2)*/
double DerivadaX1 (double x1, double x2)
{
    double termoCosNegativo, termoCosPositivo, termoSenoNegativo, termoSenoPositivo;
    termoCosNegativo = (double) cos (pow(x1, 2) - pow(x2, 2));
    termoCosPositivo = (double) cos (pow(x1, 2) + pow(x2, 2));
    termoSenoNegativo = (double) sin (pow(x1, 2) - pow(x2, 2));
    termoSenoPositivo = (double) sin (pow(x1, 2) + pow(x2, 2));

    return (double) 2* x1 * ( (termoCosNegativo * termoCosPositivo) - (termoSenoNegativo * termoSenoPositivo));
}

/*retorna o valor da derivada de f em x2 no ponto (x1,x2)*/
double DerivadaX2 (double x1, double x2)
{
    double termoCosNegativo, termoCosPositivo, termoSenoNegativo, termoSenoPositivo;
    termoSenoNegativo = (double) sin (pow(x1, 2) - pow(x2, 2));
    termoSenoPositivo = (double) sin (pow(x1, 2) + pow(x2, 2));
    termoCosNegativo = (double) cos(pow (x1, 2) - pow(x2, 2));
    termoCosPositivo = (double) cos(pow(x1, 2) + pow (x2, 2));
    return (double) -2*x2*((termoSenoNegativo*termoSenoPositivo) + (termoCosNegativo*termoCosPositivo));
}

/*Calcula o gradiente de F no ponto (x1, x2)*/
void GradienteF (double x1, double x2, vector<double> &gradiente)
{
    gradiente.resize(2);
    gradiente[0] = DerivadaX1(x1, x2);
    gradiente[1] = DerivadaX2(x1, x2);
}
/*funcao auxiliar para o metodo da secao aurea*/
double PhiDeT (vector<double> pontox, double t, vector<double> direcao)
{
    int k = 0;
    for (auto const &dir: direcao)
    {
        pontox[k]  += dir*t;
        k++;
    }
    return funcao(pontox[0], pontox[1]);
}

double MetodoGradiente (vector<double> pontox)
{
    unsigned k = 0;
    vector<double> gradiente;
    vector<vector<double>> direcao;
    vector<double> auxiliar;
    double t = 0;
    unsigned index;

    GradienteF(pontox[0], pontox[1], gradiente);
    /*Gradiente negativo*/
    while(gradiente[0] != 0 && gradiente[1] != 0)
    {
        /*d = -gradiente transposto*/
        for (auto const &var: gradiente)
        {
            auxiliar.push_back(var*-1);
            direcao.push_back(auxiliar);
            auxiliar.clear();
        }
        /*metodo armijo ou secao aurea*/
        for (index = 0; index < (unsigned) pontox.size(); index++)
            pontox[index] += t * direcao[index][0];
        k++;
        GradienteF(pontox[0], pontox[1], gradiente);
    }

    return 2.0;

}


double BuscaArmijo (vector<double> x)
{
    double t = 1;
    /*  while (f(x + t*d) > f(x) + n*t*GradienteF(x) * d )
            t = gama * t;
    */
    return t;
}

double SecaoAurea (double ro, double epsolon, vector<double> pontox, vector<double> direcao)
{
    double a = 0;
    double s = ro;
    double b = 2*ro;
    /*obtencao do intervalo [a,b]*/
    while (PhiDeT(pontox, b, direcao) < PhiDeT(pontox, s, direcao))
    {
        a = s;
        s = b;
        b = 2*b;
    }

    /*Obtencao de t*/
    return a;
}

int MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado)
{

    unsigned linha, coluna, index;
    if (m1[0].size() != m2.size())
        return ERRO_DIMENSAO_MATRIZ;

    //mxp * pxn = mxn
    //pra linha da primeira
    //pra colunas da segunda
    resultado.resize(m1.size());
    for (linha = 0; linha < (unsigned) m1.size(); linha++)
    {
        resultado[linha].resize(m2[0].size(), 0);
        for (coluna = 0; coluna < (unsigned) m2[0].size(); coluna++)
            for (index = 0; index < (unsigned) m2.size(); index++)
                resultado[linha][coluna] += m1[linha][index] * m2[index][linha];
    }
    return OK;
}

void PrintarMatriz (vector<vector<double>> m1)
{
    for (auto const &linha: m1)
    {
        for (auto const &coluna: linha)
            cout << coluna << "\t";
        cout << endl;
    }
}

int main ()
{
    vector<vector<double>> m1;
    vector<vector<double>> m2;
    vector<vector<double>> resultado;
    vector<double> v1;
    vector<double> v2;
    v1.resize(4, 2);
    v2.resize(3, 1);

    m1.push_back(v1);
    m1.push_back(v1);

    m2.push_back(v2);
    m2.push_back(v2);
    m2.push_back(v2);
    m2.push_back(v2);


    cout << "m1:" << endl;
    PrintarMatriz(m1);
    cout << "m2:" << endl;
    PrintarMatriz(m2);
    cout << "resultado:" << endl;
    if (MultiplicarMatrizes(m1, m2, resultado) == ERRO_DIMENSAO_MATRIZ)
        cout << "Erro Multiplicando Matrizes" << endl;
    PrintarMatriz(resultado);

    return OK;
}
