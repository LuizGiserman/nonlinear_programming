
#ifndef _OPERACOES_MATRIZES_
#define _OPERACOES_MATRIZES_   "operacoesMatrizes.h"


#include <iostream>
#include <vector>
#include <cmath>

#define OK                      0
#define ERRO_DIMENSAO_MATRIZ    1
#define ERRO_INVERTENDO_MATRIZ  1
#define OPCAO_SOMA              1
#define OPCAO_SUBTRACAO         0

using namespace std;

int MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado);
/*Para quando a resposta esta no R1*/
double MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2);
/*matriz X vetor*/
int MultiplicarMatrizes (vector<vector<double>> m1, vector<double> v1, vector<vector<double>> &resultado);

/*vector X matriz(vector)*/
double MultiplicarMatrizes( vector<double> v1, vector<vector<double>> m2);

int OperarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado, bool soma);
void PrintarMatriz (vector<vector<double>> m1);
int InverterDxD (vector<vector<double>> &m1);

double ModuloVetor(vector<double> vetor);
double ModuloMatrizDxD (vector<vector<double>> &matriz);

double Modulo (double numero);

#endif
