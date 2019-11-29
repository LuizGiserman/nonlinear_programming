
#ifndef _OPERACOES_MATRIZES_
#define _OPERACOES_MATRIZES_   "operacoesMatrizes.h"


#include <iostream>
#include <vector>

#define OK                      0
#define ERRO_DIMENSAO_MATRIZ    1
#define ERRO_INVERTENDO_MATRIZ  1

using namespace std;

double MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado);
/*Para quando a resposta esta no R1*/
double MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2);
/*matriz X vetor*/
double MultiplicarMatrizes (vector<vector<double>> m1, vector<double> v1, vector<vector<double> &resultado)

void PrintarMatriz (vector<vector<double>> m1);
int InverterDxD (vector<vector<double>> &m1);


double Modulo (double numero);

#endif
