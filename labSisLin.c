#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main (){
    char c = 'a';
    SistLinear_t * sistLin;
    double Gauss_t, Jacobi_t, Seidel_t;
    int n_sistema = 1;
    while(c != EOF){
        if(c != ' '){
            sistLin = lerSistLinear();
            printf("Resolvendo Sistema %i\n", n_sistema);
            prnSistLinear(sistLin);

            //resolve sistema pela eliminacao de Gauss
            // eliminacaoGauss(sistLin, solution, &totalTime);
            // printf("Solução Gauss:\n");
            // prnVetor(solution, sistLin->n);
            printf("\n");

            //aloca memoria para vetor solucao
            real_t * solucao = malloc(sizeof(real_t) * sistLin->n);

            // printf("Resolvendo sistema com método de Jacobi\n");
            // gaussJacobi(sistLin, solucao, &Jacobi_t);
            // prnVetor(solucao, sistLin->n);
            // printf("\n");

            // resetVector(solucao, sistLin->n);
            printf("Resolvendo sistema com método de Gauss Seidel\n");
            gaussSeidel(sistLin, solucao, &Seidel_t);
            prnVetor(solucao, sistLin->n);
            n_sistema++;

            //Aplicar método de elimanacao de gauss com pivoteamento, Jacobi e Gauss-Seidel
            //Medir o tempo de cada solução
            //Calcular norma L2 ou norma euclidiana do residuo de cada solucao
            //Aplicar Refinamento à solução obtida caso a norma L2 do resíduo seja maior que 5.0;

            liberaSistLinear(sistLin);
        }
        //get the \n or EOF
        c = getchar();
    }
}
