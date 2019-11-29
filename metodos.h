#ifndef _METODOS_H
#define _METODOS_H  "metodos.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <vector>
#include "operacoesMatrizes.h"

#define BUSCA_ARMIJO            1
#define BUSCA_SECAO_AUREA       0
#define OK                      0
using namespace std;


/*Especificas da funcao*/
double funcao (double x1, double x2);
double DerivadaX1 (double x1, double x2);
double DerivadaX2 (double x1, double x2);
void GradienteF (double x1, double x2, vector<double> &gradiente);
double PhiDeT (vector<double> pontox, double t, vector<vector<double>> direcao);

void CalcularHessiana(vector<double> pontox, vector<vector<double>> &hessiana);

/*Metodo do gradiente*/
double MetodoGradiente(vector <double> &pontox, bool metodo, double ro, double epsolon);
double MetodoGradiente(vector <double> &pontox, bool metodo, double epsolon, double eta, double gama);
double MetodoGradiente (vector<double> &pontox, bool metodo, double ro, double epsolon, double eta, double gama );

/*Metodo de newton*/
double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double ro);
double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double gama, double eta);
double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double ro, double gama, double eta);



/*Buscas*/
double BuscaArmijo (vector<double> pontox, vector<vector<double>> direcao, vector<double> gradiente, double gama, double eta);
double SecaoAurea (double ro, double epsolon, vector<double> pontox, vector<vector<double>> direcao);



#endif
