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
            printf(">>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<\n");
            printf(">>>>>>>>>>>>>>>Resolvendo Sistema %i<<<<<<<<<<<<<<<\n", n_sistema);
            printf(">>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<\n\n");

            // prnSistLinear(sistLin);

            real_t * solucao = malloc(sizeof(real_t) * sistLin->n);
            sist_copy = alocaSistLinear(sistLin->n);
            cpySist(sist_copy, sistLin);

            //resolve sistema pela eliminacao de Gauss
            printf("===========Eliminação de Gauss===========\n\n");
            eliminacaoGauss(sist_copy, solucao, &Gauss_t);
            norma = normaL2Residuo(sistLin, solucao, NULL);
            prnSolucao(solucao, sistLin->n, Gauss_t, it, norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                //reseta sist_copy
                // cpySist(sist_copy, sistLin);

                //aplica metodo de refinamento
                it = refinamento(sist_copy, solucao, &ref_t);

                //calcula a norma para a nova solucao encontrada
                norma = normaL2Residuo(sistLin, solucao, NULL);

                prnSolucao(solucao, sistLin->n, ref_t, it, norma);
            }

            printf("===============Gauss Jacobi==============\n\n");

            it = gaussJacobi(sistLin, solucao, &Jacobi_t);
            norma = normaL2Residuo(sistLin, solucao, NULL);

            prnSolucao(solucao, sistLin->n, Jacobi_t, it, norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                //reseta sist_copy
                cpySist(sist_copy, sistLin);

                it = refinamento(sist_copy, solucao, &ref_t);
                norma = normaL2Residuo(sistLin, solucao, NULL);

                prnSolucao(solucao, sistLin->n, ref_t, it, norma);
            }


            printf("===============Gauss Seidel===============\n\n");

            it = gaussSeidel(sistLin, solucao, &Seidel_t);

            //calcula norma da solucao
            norma = normaL2Residuo(sistLin, solucao, NULL);

            prnSolucao(solucao, sistLin->n, Seidel_t, it, norma);

            if(norma >= 5.0){
                printf("===============Refinamento==============\n\n");
                //reseta sist_copy
                cpySist(sist_copy, sistLin);

                it = refinamento(sist_copy, solucao, &ref_t);

                //cacula a norma para a nova solucao
                norma = normaL2Residuo(sistLin, solucao, NULL);

                prnSolucao(solucao, sistLin->n, ref_t, it, norma);
            }

            n_sistema++;

            liberaSistLinear(sistLin);
            liberaSistLinear(sist_copy);
        }
        //get the \n or EOF
        c = getchar();
    }
}
