#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main (){
    double Gauss_t, Jacobi_t, Seidel_t;
    double Gauss_norm, Jacobi_norm, Seidel_norm;

    SistLinear_t * sistLin;

    char c = 'a';
    int n_sistema = 1;
    while(c != EOF){
        if(c != ' '){
            sistLin = lerSistLinear();
            printf(">>>>>>>>>>Resolvendo Sistema %i<<<<<<<<<<\n", n_sistema);
            prnSistLinear(sistLin);

            real_t * solucao = malloc(sizeof(real_t) * sistLin->n);

            //resolve sistema pela eliminacao de Gauss
            printf("===========Eliminação de Gauss===========\n\n");
            eliminacaoGauss(sistLin, solucao, &Gauss_t);
            prnVetor(solucao, sistLin->n);
            printf("Tempo: %lf\n\n", Gauss_t);

            Gauss_norm = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: f:%1.9e\n\n", Gauss_norm);


            printf("===============Gauss Jacobi==============\n\n");
            gaussJacobi(sistLin, solucao, &Jacobi_t);
            prnVetor(solucao, sistLin->n);
            printf("Tempo: %lf\n\n", Jacobi_t);

            Jacobi_norm = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: f:%1.9e\n\n", Jacobi_norm);

            printf("===============Gauss Seidel===============\n\n");
            gaussSeidel(sistLin, solucao, &Seidel_t);
            prnVetor(solucao, sistLin->n);
            printf("Tempo: %lf\n\n", Seidel_t);

            Seidel_norm = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: f:%1.9e\n\n", Seidel_norm);

            n_sistema++;

            // Calcular a norma L2 (ou norma euclidiana) do resíduo de cada uma das soluções;
            // Aplicar o método de Refinamento à solução obtida em cada um dos outros métodos caso a norma L2 do resíduo seja maior que 5.0;

            liberaSistLinear(sistLin);
        }
        //get the \n or EOF
        c = getchar();
    }
}
