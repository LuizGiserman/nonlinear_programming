#include <math.h>
#include <iostream>
#include <vector>


#define OK      0

using namespace std;

/**/

double funcao (double x1, double x2)
    return (double) ( sin ( pow(x1, 2) - pow(x2, 2) ) * cos ( pow(x1, 2) + pow(x2, 2) ) );

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
        pontox[k]  += dir*t
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

double secaoAurea (double ro, double epslon, vector<double> pontox, vector<double> direcao)
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
    

}



int main ()
{

    return OK;
}
