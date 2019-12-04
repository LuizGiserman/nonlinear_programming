#include <iostream>
#include "operacoesMatrizes.h"
#include "metodos.h"

int main ()
{

    vector<double> pontox;
    double x1, x2;
    x1 = 3;
    x2 = 2;
    pontox.push_back(x1);
    pontox.push_back(x2);
    double index;
    vector<vector<double>> resultado;
    cout << "funcao ( " << pontox[0] << "," << pontox[1] << "): " << funcao(pontox[0], pontox[1]) << endl;

    // vector<vector<double>> dir;
    // dir.resize(2);
    // dir[0].push_back(1);
    // dir[1].push_back(1);

    // cout << "result: " << MultiplicarMatrizes(pontox, dir);
    // vector<vector<double>> m1, m2;
    // m1.resize(2);
    // m1[0].resize(2);
    // m1[1].resize(2);
    // m1[0][0] = -1;
    // m1[0][1] = 0;
    // m1[1][0] = 0;
    // m1[1][1] = -1;
    // m2.resize(2);
    // m2[0].push_back (6.37199);
    // m2[1].push_back(-5.54608);
    // PrintarMatriz(m1);
    // PrintarMatriz(m2);
    // MultiplicarMatrizes(m1, m2, m1);
    // PrintarMatriz(m1);
    // cout << "resposta: " << MetodoGradiente(pontox, BUSCA_ARMIJO, 0.0000001, 0.001, 0.7) << endl;
    // for (index = 0.01; index <= 0.1; index+=0.01)
    {
        // cout << "Index: " << index << endl;
        MetodoQuaseNewton(pontox, BUSCA_SECAO_AUREA,  0.0000001, 0.01, 0, 0);
        // MetodoQuaseNewton(pontox, BUSCA_ARMIJO, 0.0000001, 0.01, 0.25, 0.8);
        // MetodoNewton(pontox, BUSCA_ARMIJO, 0.0000001, 0.7, 0.25);
        // MetodoNewton(pontox, BUSCA_SECAO_AUREA, 0.0000001, 0.01);

        cout << "Pontox: (" << pontox[0] << ", " << pontox[1] << ")" << endl;
        cout << "funcao: " << funcao(pontox[0], pontox[1]);
        pontox[0] = x1;
        pontox[1] = x2;
        // MetodoQuaseNewton(pontox, BUSCA_ARMIJO, 0.0000001, 0.01, 0.25, 0.8);
        // // MetodoNewton(pontox, BUSCA_SECAO_AUREA, 0.0000001, index);
        // pontox[0] = x1;
        // pontox[1] = x2;

        cout << endl;
    }
    // cout << MetodoNewton(pontox, BUSCA_ARMIJO, 0.0000001, 0.7, 0.1) << endl;
    // cout << MetodoQuaseNewton(pontox, BUSCA_SECAO_AUREA, 0.00000001, 0.1, 0.7, 0.01);

    // MetodoGradiente(pontox, BUSCA_SECAO_AUREA, 0.03, 0.0000001);

    return OK;
}
