#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "SupportFunctions.h"

int main (){
    double Gauss_t, Jacobi_t, Seidel_t, ref_t, norma;
    SistLinear_t * sistLin, *sist_copy;

    char c = 'a';
    int n_sistema = 1, it;

    while(c != EOF){
        if(c != ' '){
            sistLin = lerSistLinear();
            printf("***** Sistema %i --> n = %i, erro %f\n", n_sistema, sistLin->n, sistLin->erro);
            real_t * solucao = malloc(sizeof(real_t) * sistLin->n);
            sist_copy = alocaSistLinear(sistLin->n);
            cpySist(sist_copy, sistLin);

            //resolve sistema pela eliminacao de Gauss
            it = eliminacaoGauss(sist_copy, solucao, &Gauss_t);
            norma = normaL2Residuo(sistLin, solucao, NULL);
            printf("===>Eliminação de Gauss: %lf\n", Gauss_t);
            printf("--> X:");
            prnVetor(solucao, sistLin->n);
            printf("--> Norma L2 do residuo: %lf\n\n", norma);

            if(norma >= 5.0){

                //aplica metodo de refinamento
                it = refinamento(sist_copy, solucao, &ref_t);

                //calcula a norma para a nova solucao encontrada
                norma = normaL2Residuo(sistLin, solucao, NULL);

                printf("===>Refinamento: %lf -> %i iteracoes\n", ref_t, it);
                prnSolucao(solucao, sistLin->n, norma);
            }


            it = gaussJacobi(sistLin, solucao, &Jacobi_t);
            norma = normaL2Residuo(sistLin, solucao, NULL);

            printf("===>Gauss Jacobi: lf --> %i iteracoes\n", Jacobi_t, it);
            prnSolucao(solucao, sistLin->n, norma);

            if(norma >= 5.0){
                //reseta sist_copy
                cpySist(sist_copy, sistLin);

                it = refinamento(sist_copy, solucao, &ref_t);
                norma = normaL2Residuo(sistLin, solucao, NULL);

                printf("===>Refinamento: %lf --> %i iteracoes\n", ref_t, it);
                prnSolucao(solucao, sistLin->n, norma);
            }

            //Resolve sistema com metodo gauss-seidel
            it = gaussSeidel(sistLin, solucao, &Seidel_t);

            //calcula norma da solucao
            norma = normaL2Residuo(sistLin, solucao, NULL);

            printf("===>Gauss Seidel: lf --> %i iteracoes \n", Seidel_t, it);
            prnSolucao(solucao, sistLin->n, norma);

            if(norma >= 5.0){
                //reseta sist_copy
                cpySist(sist_copy, sistLin);

                it = refinamento(sist_copy, solucao, &ref_t);

                //cacula a norma para a nova solucao
                norma = normaL2Residuo(sistLin, solucao, NULL);

                printf("===>Refinamento: %lf --> %i iteracoes\n", ref_t, it);
                prnSolucao(solucao, sistLin->n, it, norma);
            }

            n_sistema++;

            liberaSistLinear(sistLin);
            liberaSistLinear(sist_copy);
        }
        //get the \n or EOF
        c = getchar();
    }
}
