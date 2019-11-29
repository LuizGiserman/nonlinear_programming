#include <iostream>
#include "operacoesMatrizes.h"
#include "metodos.h"

int main ()
{

    vector<double> pontox;
    pontox.push_back(1.5);
    pontox.push_back(0.8);
    cout << "funcao ( " << pontox[0] << "," << pontox[1] << "): " << funcao(pontox[0], pontox[1]) << endl;
    // cout << "resposta: " << MetodoGradiente(pontox, BUSCA_ARMIJO, 0.0000001, 0.001, 0.7) << endl;
    cout << MetodoNewton(pontox, BUSCA_ARMIJO, 0.0000001, 0.7, 0.1) << endl;
    return OK;
}
