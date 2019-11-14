#include <math.h>
#include <iostream.h>


#define OK      0

using namespace std;

/**/

double DerivadaX1 (double x1, double x2)
{
    double termoCosNegativo, termoCosPositivo, termoSenoNegativo, termoSenoPositivo;
    termoCosNegativo = cos (pow(x1, 2) - pow(x2, 2));
    termoCosPositivo = cos (pow(x1, 2) + pow(x2, 2));
    termoSenoNegativo = sin (pow(x1, 2) - pow(x2, 2));
    termoSenoPositivo = sin (pow(x1, 2) + pow(x2, 2));

    return 2*( (termoCosNegativo * termoCosPositivo) - (termoSenoNegativo * termoSenoPositivo))
}

void DerivadaX2 (double x1, double x2)
{

}

/*Calcula o gradiente de F no ponto (x1, x2)*/
void GradienteF (double x1, double x2, vector<double> gradiente)
{
    gradiente.resize(2);
    gradiente[0] = DerivadaX1(x1, x2);
    gradiente[1] = DerivadaX2(x1, x2);
}

double MetodoGradiente (vector<double> x)
{

    vector<double> direcao;
    direcao.resize(2);
    GradienteF(x[0], x[1], direcao)
    /*Gradiente negativo*/
    for (auto& var: direcao)
        var *= -1;
    /*metodo armijo ou secao aurea*/


}





int main ()
{

    return OK;
}
