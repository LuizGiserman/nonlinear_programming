#include <iostream>
#include <vector>
#include "operacoesMatrizes.h"

double MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado)
{

    unsigned linha, coluna, index;
    bool alocarMemoria;
    if (m1[0].size() != m2.size())
        return ERRO_DIMENSAO_MATRIZ;

    //mxp * pxn = mxn
    //pra linha da primeira
    //pra colunas da segunda
    if(resultado.size() == 0)
    {
        alocarMemoria = true;
        resultado.resize(m1.size());
    }
    for (linha = 0; linha < (unsigned) m1.size(); linha++)
    {
        if(alocarMemoria)
            resultado[linha].resize(m2[0].size(), 0);
        for (coluna = 0; coluna < (unsigned) m2[0].size(); coluna++)
            for (index = 0; index < (unsigned) m2.size(); index++)
                resultado[linha][coluna] += m1[linha][index] * m2[index][linha];
    }
    if (m1.size() == 1 && m2[0].size() == 1 )
        return resultado[0][0];
    return OK;
}

double MultiplicarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2)
{
    vector<vector<double>> resultado;
    return MultiplicarMatrizes(m1, m2, resultado);
}

/*matriz X vetor*/
double MultiplicarMatrizes (vector<vector<double>> m1, vector<double> v1, vector<vector<double>> &resultado)
{
    vector<vector<double>> m2;
    m2.push_back(v1);
    return MultiplicarMatrizes(m1, m2, resultado);
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

double Modulo (double numero)
{
    if (numero < 0)
        numero *= -1;
    return numero;
}


int InverterDxD (vector<vector<double>> &m1)
{
    double aux;
    if (m1.empty())
        return ERRO_INVERTENDO_MATRIZ;
    if (m1.size() != 2 || m1[0].size() != 2)
        return ERRO_INVERTENDO_MATRIZ;

    aux = m1[0][0];
    m1[0][0] = m1[1][1];
    m1[1][1] = aux;

    m1[0][1] *= -1;
    m1[1][0] *= -1;

    return OK;
}

int OperarMatrizes (vector<vector<double>> m1, vector<vector<double>> m2, vector<vector<double>> &resultado, bool soma)
{

    unsigned linha, coluna;
    bool alocarMemoria;

    if (m1.size() != m2.size() || m1[0].size() != m2[0].size())
    {
        cout << "Erro ERRO_DIMENSAO_MATRIZ" << endl;
        return ERRO_DIMENSAO_MATRIZ;
    }

    if (resultado.size() == 0)
    {
        alocarMemoria = true;
        resultado.resize(m1.size());
    }
    for (linha = 0; linha < (unsigned) m1.size(); linha++)
    {
        if (alocarMemoria)
            m1[linha].resize(m1[0].size());
        for(coluna = 0; coluna < (unsigned) m1[0].size(); coluna++)
            if (soma)
                resultado[linha][coluna] = m1[linha][coluna] + m2[linha][coluna];
            else
                resultado[linha][coluna] = m1[linha][coluna] - m2[linha][coluna];

    }

    return OK;
}
