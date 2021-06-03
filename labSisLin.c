#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main (){
    double Gauss_t, Jacobi_t, Seidel_t, ref_t, norma;
    SistLinear_t * sistLin;

    char c = 'a';
    int n_sistema = 1, it;

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
            printf("Tempo: %lf\n", Gauss_t);

            norma = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: %1.9e\n\n", norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                it = refinamento(sistLin, solucao, &ref_t);
                prnVetor(solucao, sistLin->n);
                printf("Iterações: %i\n", it);
                printf("Tempo: %lf\n", ref_t);

                norma = normaL2Residuo(sistLin, solucao, NULL);
                printf("Norma L2 do residuo: %1.9e\n\n", norma);
            }


            printf("===============Gauss Jacobi==============\n\n");
            it = gaussJacobi(sistLin, solucao, &Jacobi_t);
            prnVetor(solucao, sistLin->n);
            printf("Iterações: %i\n", it);
            printf("Tempo: %lf\n", Jacobi_t);

            norma = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: %1.9e\n\n", norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                it = refinamento(sistLin, solucao, &ref_t);
                prnVetor(solucao, sistLin->n);
                printf("Iterações: %i\n", it);
                printf("Tempo: %lf\n", ref_t);

                norma = normaL2Residuo(sistLin, solucao, NULL);
                printf("Norma L2 do residuo: %1.9e\n\n", norma);
            }


            printf("===============Gauss Seidel===============\n\n");
            it = gaussSeidel(sistLin, solucao, &Seidel_t);
            prnVetor(solucao, sistLin->n);
            printf("Iterações: %i\n", it);
            printf("Tempo: %lf\n", Seidel_t);

            norma = normaL2Residuo(sistLin, solucao, NULL);
            printf("Norma L2 do residuo: %1.9e\n\n", norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                it = refinamento(sistLin, solucao, &ref_t);
                prnVetor(solucao, sistLin->n);
                printf("Iterações: %i\n", it);
                printf("Tempo: %lf\n", ref_t);

                norma = normaL2Residuo(sistLin, solucao, NULL);
                printf("Norma L2 do residuo: %1.9e\n\n", norma);
            }

            n_sistema++;

            liberaSistLinear(sistLin);
        }
        //get the \n or EOF
        c = getchar();
    }
}
