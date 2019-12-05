#include <cmath>
#include <iostream>
#include <vector>
#include "operacoesMatrizes.h"
#include "metodos.h"


/**/

double funcao (double x1, double x2)
{
    return ( sin ( pow(x1, 2) - pow(x2, 2) ) * cos ( pow(x1, 2) + pow(x2, 2) ) );
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
    gradiente.clear();
    gradiente.resize(2);
    gradiente[0] = DerivadaX1(x1, x2);
    gradiente[1] = DerivadaX2(x1, x2);

}
/*funcao auxiliar para o metodo da secao aurea*/
double PhiDeT (vector<double> pontox, double t, vector<vector<double>> direcao)
{
    int k = 0;
    for (auto const &dir: direcao)
    {
        /*Como direcao eh transposto, precisamos pegar o elemento 0 de cada vector*/
        pontox[k]  += dir[0]*t;
        k++;
    }
    return funcao(pontox[0], pontox[1]);
}

void CalcularHessiana(vector<double> pontox, vector<vector<double>> &hessiana)
{
    if (hessiana.size() == 0)
    {
        hessiana.resize(2);
        hessiana[0].resize(2);
        hessiana[1].resize(2);
    }

    double termoCosNegativo, termoCosPositivo, termoSenoNegativo, termoSenoPositivo;
    termoCosNegativo = cos ( pow(pontox[0], 2) - pow (pontox[1], 2));
    termoCosPositivo = cos ( pow(pontox[0], 2) + pow (pontox[1], 2));
    termoSenoNegativo = sin ( pow(pontox[0], 2) - pow (pontox[1], 2));
    termoSenoPositivo = sin ( pow(pontox[0], 2) + pow (pontox[1], 2));

    hessiana[0][0] =  (2 * termoCosNegativo * (termoCosPositivo - (4*pow(pontox[0], 2)*termoSenoPositivo))) - (2 * termoSenoNegativo *(termoSenoPositivo + (4* pow(pontox[0], 2) * termoCosPositivo)));
    hessiana[0][1] = 0.0;
    hessiana[1][0] = 0.0;
    hessiana[1][1] = -2*(termoSenoNegativo*(termoSenoPositivo + 4 * pow(pontox[1], 2) *termoCosPositivo) + termoCosNegativo*(termoCosPositivo - 4 * pow(pontox[1], 2) * termoSenoPositivo));

}

void SegundaDerivada (vector<vector<double>> &segundaDerivada, vector<double> pontox)
{
    vector<vector<double>> hessiana;
    CalcularHessiana(pontox, hessiana);

    if (segundaDerivada.size() != 0)
        segundaDerivada.clear();

    segundaDerivada.resize(2);
    segundaDerivada[0].push_back(hessiana[0][0]);
    segundaDerivada[1].push_back(hessiana[1][1]);
}

bool AtualizarMatrizH(vector<vector<double>> &matrizH, vector<double> pontox, vector<double> pontoxAntigo, vector<vector<double>> gradiente, vector<vector<double>> gradienteAntigo)
{
    /*BFGS*/
    vector<vector<double>> p, q;
    vector<vector<double>> pt, qt;
    vector<vector<double>> resultadoAuxiliar, resultadoAuxiliar2, resultadoParcial;
    double escalarAuxiliar;
    double denominador;
    unsigned index, colunas, coluna, linha;

    /**gerando p e pt**/
    p.resize(pontox.size());
    pt.resize(1);
    for (index = 0; index < (unsigned) pontox.size(); index++)
    {
        p[index].push_back(pontox[index] - pontoxAntigo[index]);
        pt[0].push_back(pontox[index] - pontoxAntigo[index]);
    }


    /**gerando q e qt**/
    q.resize(gradiente.size());
    qt.resize(gradiente[0].size());
    for(colunas = 0; colunas < (unsigned) gradiente[0].size(); colunas++)
        qt[colunas].resize(gradiente.size());

    for (index = 0; index < (unsigned) gradiente.size(); index++)
    {
        q[index].resize(gradiente[0].size());
        for(colunas = 0; colunas < (unsigned) gradiente[0].size(); colunas++)
        {
            q[index][colunas] = (gradiente[index][colunas] - gradienteAntigo[index][colunas]);
            // cout << "Q[" << index << "][" << colunas << "] = " << gradiente[index][colunas] << endl;
            qt[colunas][index] = q[index][colunas];
        }
    }

    // cout << "p" << endl;
    // PrintarMatriz(p);
    // cout << "pt" << endl;
    // PrintarMatriz(pt);
    // cout << "q" << endl;
    // PrintarMatriz(q);
    // cout << "qt" << endl;
    // PrintarMatriz(qt);

    /*1 + esse coisa*/
    MultiplicarMatrizes(qt, matrizH, resultadoAuxiliar);
    denominador = MultiplicarMatrizes(pt, q);
    // cout << "denominador: "  << denominador << endl;
    if(denominador == 0)
        return !OK;
    escalarAuxiliar = 1 + (MultiplicarMatrizes(resultadoAuxiliar, q) / denominador);

    /*Logo depois do parenteses*/
    MultiplicarMatrizes(p, pt, resultadoParcial);
    for (auto &vect: resultadoParcial)
        for (auto &var: vect)
        {
            var *= escalarAuxiliar;
            var /= denominador;
        }

    resultadoAuxiliar.clear();
    MultiplicarMatrizes(p,qt, resultadoAuxiliar);
    MultiplicarMatrizes(resultadoAuxiliar, matrizH, resultadoAuxiliar);
    MultiplicarMatrizes(matrizH, q, resultadoAuxiliar2);
    MultiplicarMatrizes(resultadoAuxiliar2, pt, resultadoAuxiliar2);
    OperarMatrizes(resultadoAuxiliar, resultadoAuxiliar2, resultadoAuxiliar, OPCAO_SOMA);

    for (auto &vect: resultadoAuxiliar)
        for(auto &var: vect)
            var /= denominador;
    OperarMatrizes(resultadoParcial, resultadoAuxiliar, resultadoParcial, OPCAO_SUBTRACAO);

    OperarMatrizes(matrizH, resultadoParcial, matrizH, OPCAO_SOMA);

    return OK;
}

double SecaoAurea (double ro, double epsolon, vector<double> pontox, vector<vector<double>> direcao)
{
    double a = 0;
    double s = ro;
    double b = 2*ro;
    double tetha1 = (3 - sqrt(5)) / 2.0;
    double tetha2 = 1 - tetha1;
    double u, v;

    /*obtencao do intervalo [a,b]*/
    while (PhiDeT(pontox, b, direcao) < PhiDeT(pontox, s, direcao))
    {
        a = s;
        s = b;
        b = 2*b;
        // cout << "(a, b): " << a << ", " << b << endl;
    }

    cout << "intervaloInicial: [" << a << ","  << b << "]" << endl;

    /*Obtencao de t*/
    u = a + tetha1 * (b-a);
    v = a + tetha2 * (b-a);

    while ( (b - a) > epsolon )
        if (PhiDeT(pontox, u, direcao) < PhiDeT(pontox, v, direcao))
        {
            b = v;
            v = u;
            u = a + tetha1 *(b-a);
        }
        else
        {
            a = u;
            u = v;
            v = a + tetha2*(b-a);
        }

    cout << "intervaloFinal: [" << a << ","  << b << "]" << endl;

    return (u+v)/2.0;
}


double BuscaArmijo (vector<double> pontox, vector<vector<double>> direcao, vector<double> gradiente, double gama, double eta)
{
    double t = 1;
    /*  while (f(x + t*d) > f(x) + n*t*GradienteF(x) * d )
            t = gama * t;
    */
    while (PhiDeT(pontox, t, direcao) >
    (funcao(pontox[0], pontox[1]) + eta*t*MultiplicarMatrizes(gradiente, direcao)))
    {
        t *= gama;
        // cout << "PhideT: " << PhiDeT(pontox, t, direcao) << " OutraParte: " << funcao(pontox[0], pontox[1]) + eta*t*MultiplicarMatrizes(gradiente, direcao) << " T: " << t << endl;
    }

    return t;
}

/*Gradiente para busca da secao aurea*/
double MetodoGradiente(vector <double> &pontox, bool metodo, double ro, double epsolon)
{
    double eta = 0;
    double gama = 0;
    return MetodoGradiente(pontox, metodo, ro, epsolon, eta, gama);
}

/*Gradiente para busca de armijo*/
double MetodoGradiente(vector <double> &pontox, bool metodo, double epsolon, double eta, double gama)
{
    double ro =0;
    return MetodoGradiente(pontox, metodo, ro, epsolon, eta, gama);
}

double MetodoGradiente (vector<double> &pontox, bool metodo, double ro, double epsolon, double eta, double gama)
{
    unsigned k = 0;
    vector<double> gradiente;
    vector<vector<double>> direcao;
    vector<double> auxiliar;
    double t = 0;
    unsigned index;
    vector<double> pontoxAntigo;

    /**so pra cond de parada**/
    for (auto const &x: pontox)
        pontoxAntigo.push_back(x*30);

    GradienteF(pontox[0], pontox[1], gradiente);
    /*Gradiente negativo*/
    // cout << "calculou o gradiente" << endl;

    while( Modulo(gradiente[0]) > epsolon && Modulo(gradiente[1]) > epsolon) //&&
    // (Modulo(ModuloVetor(pontox) - ModuloVetor(pontoxAntigo)) > epsolon))
    {
        pontoxAntigo.clear();
        for (auto const &x: pontox)
            pontoxAntigo.push_back(x);

        // cout << endl << endl << "Gradiente: (" << gradiente[0] << "," << gradiente[1] << ")" << endl;

        /*d = -gradiente transposto*/
        for (auto const &var: gradiente)
        {
            auxiliar.push_back(var*-1);
            direcao.push_back(auxiliar);
            auxiliar.clear();
        }


        /*metodo armijo ou secao aurea*/
        if (metodo == BUSCA_ARMIJO)
            t = BuscaArmijo(pontox, direcao, gradiente, gama, eta);
        else /*BUSCA SECAO AUREA*/
            t = SecaoAurea(ro, epsolon, pontox, direcao);

        for (index = 0; index < (unsigned) pontox.size(); index++)
            pontox[index] += t * direcao[index][0];
        k++;
        // cout << "pontox = " << pontox[0] << ", " << pontox[1] << endl;
        cout << "K = " << k-1 << " | pontox = (" << pontox[0] << "," << pontox[1] << ")" << endl;
        cout << "funcao dps de atualizar = " << funcao(pontox[0], pontox[1]) << endl;
        GradienteF(pontox[0], pontox[1], gradiente);
        direcao.clear();
    }
    // cout << "K = " << k << " | pontox = (" << pontox[0] << "," << pontox[1] << ")" << endl;
    // cout << "funcao dps de atualizar = " << funcao(pontox[0], pontox[1]) << endl;
    // cout << "gradiente" << endl;
    // cout << gradiente[0] << ", " << gradiente[1] << endl;
    return funcao(pontox[0], pontox[1]);

}

double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double gama, double eta)
{
    double ro = 0;
    return MetodoNewton(pontox, metodo, epsolon, ro, gama, eta);
}

double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double ro)
{
    double gama = 0;
    double eta = 0;
    return MetodoNewton(pontox, metodo, epsolon, ro, gama, eta);
}


double MetodoNewton (vector<double> &pontox, bool metodo, double epsolon, double ro, double gama, double eta)
{

    vector<double> auxGradiente;
    vector<vector<double>> gradiente;
    vector<vector<double>> hessianaInversa;
    vector<vector<double>> direcao;
    double t;
    unsigned index;
    unsigned k = 0;
    vector<double> pontoxAntigo;

    /**so pra cond de parada**/
    for (auto const &x: pontox)
        pontoxAntigo.push_back(x*30);

    GradienteF(pontox[0], pontox[1], auxGradiente);
    gradiente.resize(2);
    gradiente[0].push_back(auxGradiente[0]);
    gradiente[1].push_back(auxGradiente[1]);


    while ((Modulo(gradiente[0][0]) > epsolon) && (Modulo(gradiente[1][0]) > epsolon) &&
    (Modulo(ModuloVetor(pontox) - ModuloVetor(pontoxAntigo))> epsolon))
    {
        pontoxAntigo.clear();
        for (auto const &x: pontox)
            pontoxAntigo.push_back(x);

        CalcularHessiana(pontox, hessianaInversa);
        if(InverterDxD(hessianaInversa) == ERRO_INVERTENDO_MATRIZ)
            cout << "Erro invertendo matriz" << endl;

        /*2x2 X 2x1 = 2x1*/
        MultiplicarMatrizes(hessianaInversa, gradiente, direcao);
        for (auto &dir: direcao)
            dir[0] *= -1;

        if (metodo == BUSCA_SECAO_AUREA)
            t = SecaoAurea (ro, epsolon, pontox, direcao);
        else
            t = BuscaArmijo(pontox, direcao, auxGradiente, gama, eta);

        for (index = 0; index < (unsigned) pontox.size(); index++)
            pontox[index] += t * direcao[index][0];
        k++;
        cout << "t*d: " << t*direcao[0][0] << endl;
        cout << "K = " << k-1 << " | pontox = (" << pontox[0] << "," << pontox[1] << ")" << endl;
        cout << "funcao dps de atualizar = " << funcao(pontox[0], pontox[1]) << endl;

        hessianaInversa.clear();
        direcao.clear();
        GradienteF(pontox[0], pontox[1], auxGradiente);
        gradiente[0][0] = auxGradiente[0];
        gradiente[1][0] = auxGradiente[1];
        cout << "modulo: " << Modulo(ModuloVetor(pontox) - ModuloVetor(pontoxAntigo)) << endl;


    }

    cout << "K = " << k-1 << " | pontox = (" << pontox[0] << "," << pontox[1] << ")" << endl;
    cout << "funcao dps de atualizar = " << funcao(pontox[0], pontox[1]) << endl;

    return funcao (pontox[0], pontox[1]);

}


double MetodoQuaseNewton(vector<double> &pontox, bool metodo, double epsolon, double ro, double gama, double eta)
{
    vector<vector<double>> matrizH;
    vector<vector<double>> gradiente;
    vector<double> auxGradiente;
    vector<vector<double>> direcao;
    vector<double> pontoxAntigo;
    vector<vector<double>> gradienteAntigo, segundaDerivada;
    unsigned index, k = 0;
    double t;
    double moduloMatrizHAntiga = -60;

    /*-identidade*/
    matrizH.resize(2);
    matrizH[0].resize(2);
    matrizH[1].resize(2);
    matrizH[0][0] = -1;
    matrizH[0][1] = 0;
    matrizH[1][0] = 0;
    matrizH[1][1] = -1;
    //-1 0
    //0 -1



    GradienteF(pontox[0], pontox[1], auxGradiente);
    gradiente.resize(2);
    gradiente[0].push_back(auxGradiente[0]);
    gradiente[1].push_back(auxGradiente[1]);

    pontoxAntigo.clear();
    for (auto const &x: pontox)
        pontoxAntigo.push_back(x*30);


    while ((Modulo(gradiente[0][0]) > epsolon) && (Modulo(gradiente[1][0]) > epsolon) &&
     ((Modulo(ModuloMatrizDxD(matrizH) - moduloMatrizHAntiga)) > epsolon) && (Modulo(ModuloVetor(pontox) - ModuloVetor(pontoxAntigo)) > epsolon))
    {
        pontoxAntigo.clear();
        for (auto const &x: pontox)
            pontoxAntigo.push_back(x);

        gradienteAntigo.clear();
        gradienteAntigo.resize(2);
        gradienteAntigo[0].push_back(auxGradiente[0]);
        gradienteAntigo[1].push_back(auxGradiente[1]);

        // cout << "Matriz H" << endl;
        // PrintarMatriz(matrizH);
        // cout << "Gradiente" << endl;
        // PrintarMatriz(gradiente);

        MultiplicarMatrizes(matrizH, gradienteAntigo, direcao);

        // cout << "direcao" << endl;
        // PrintarMatriz(direcao);

        if (metodo == BUSCA_SECAO_AUREA)
            t = SecaoAurea (ro, epsolon, pontox, direcao);
        else
            t = BuscaArmijo(pontox, direcao, auxGradiente, gama, eta);

        cout << "K = "<< k << endl << "pontoX: (" << pontox[0] << "," << pontox[1] << ")" << endl;
        cout << "funcao = " << funcao(pontox[0], pontox[1]) << endl;

        for (index = 0; index < (unsigned) pontox.size(); index++)
            pontox[index] += t * direcao[index][0];

        k++;
        // cout << "T = " << t << endl;
        // cout << "Gradiente = (" << auxGradiente[0] << "," << auxGradiente[1] << ")" << endl;

        // cout << "K = " << k-1 << " | pontoxNovo = (" << pontox[0] << "," << pontox[1] << ")" << endl;

        gradiente.clear();
        gradiente.resize(2);
        GradienteF(pontox[0], pontox[1], auxGradiente);
        gradiente[0].push_back(auxGradiente[0]);
        gradiente[1].push_back(auxGradiente[1]);
        moduloMatrizHAntiga = ModuloMatrizDxD(matrizH);
        if(AtualizarMatrizH(matrizH, pontox, pontoxAntigo, gradiente, gradienteAntigo) == !OK)
        {
            cout << "erro denominador = 0" << endl;
            cout << "K = "<< k << endl << "pontoX: (" << pontox[0] << "," << pontox[1] << ")" << endl;
            cout << "funcao = " << funcao(pontox[0], pontox[1]) << endl;
            return funcao(pontox[0], pontox[1]);
        }

    }
    cout << "K = "<< k << endl << "pontoX: (" << pontox[0] << "," << pontox[1] << ")" << endl;
    cout << "funcao = " << funcao(pontox[0], pontox[1]) << endl;

    // cout << "K = " << k-1 << " | pontox = (" << pontox[0] << "," << pontox[1] << ")" << endl;
    // cout << "funcao dps de atualizar = " << funcao(pontox[0], pontox[1]) << endl;
    gradiente.clear();
    return funcao(pontox[0], pontox[1]);
}
