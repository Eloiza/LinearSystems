#include "SupportFunctions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int sistemaIndeterminado(SistLinear_t * SL){
    //uma linha igual a zero
    int termosZero = 0;
    for(int i=0; i< SL->n; i++){
        termosZero = 0;
        for(int j=0; j< SL->n; j++){
            if(SL->A[i][j] == 0){
                termosZero++;
            }
        }

        if(termosZero >= SL->n){
            return 1;
        }
    }

    return 0;
}

int sistemaImpossivel(SistLinear_t * SL){
    //linha todos os coeficientes = 0, vetor b != 0
    int termosZero = 0;
    for(int i=0; i< SL->n; i++){
        for(int j=0; j< SL->n; j++){
            if(SL->A[i][j] == 0){
                termosZero++;
            }
        }
    }

    //checa se toda a matriz é igual a 0
    if(termosZero >= SL->n*SL->n){
        for(int i=0; i< SL->n; i++){
            if(SL->b[i] != 0){
                return 1;
            }
        }
    }

    return 0;
}

void prnSolucao(real_t * solucao, int n, double norma){
    printf("--> X:");
    prnVetor(solucao, sistLin->n);
    printf("--> Norma L2 do residuo: %lf\n\n", norma);
}

real_t * sumVector(real_t * a, real_t * b, unsigned int size){
    real_t * sum = malloc(sizeof(real_t)*size);
    for(int i=0; i<size; i++){
        sum[i] = a[i] + b[i];
    }

    return sum;
}
/*Copia os elementos de source para o vetor destiny.
Os vetores precisam ter o mesmo tamanho*/
void cpyVector(real_t * destiny, real_t * source, int size){
    for(int i=0; i<size; i++){
        destiny[i] = source[i];
    }
}

/*Copia os elementos de source para a matriz destiny
As matrizes precisam ser quadradas e ter mesmo tamanho*/
void cpyMatrix(real_t ** destiny, real_t ** source, int size){
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            destiny[i][j] = source[i][j];
        }
    }
}

/*Copia os atributos do sistema source para o sistema destiny.
As matrizes e vetores precisam ter o mesmo tamanho*/
void cpySist(SistLinear_t * destiny, SistLinear_t * source){
    destiny->n = source->n;
    destiny->erro = source->erro;
    cpyVector(destiny->b, source->b, destiny->n);
    cpyMatrix(destiny->A, source->A, destiny->n);
}

/*Troca duas linhas de um sistema linear. Essa operação inclui
trocas na matriz A e no vetor b*/
void swapLine(SistLinear_t * SL, int i, int j){
    real_t * aux_m;
    real_t aux_v;

    //swap lines in matrix
    aux_m = SL->A[i];
    SL->A[i] = SL->A[j];
    SL->A[j] = aux_m;

    //swap lines in vector b
    aux_v = SL->b[i];
    SL->b[i] = SL->b[j];
    SL->b[j] = aux_v;
}

/*Recebe uma matriz triangular superior e aliza a retro substituição
dos termos.*/
int retroSubst(SistLinear_t *SL, real_t * x, double * tTotal){
    for(int i= SL->n -1; i >=0; i--){
        x[i] = SL->b[i];
        for(int j = i+1; j< SL->n; j++)
            //overflow, erros de aproximacao devido subtracao
            x[i] -= SL->A[i][j] * x[j];

        //divisão por 0
        x[i] /= SL->A[i][i];
    }

    return 0;
};

/*Encontra maior valor em dada coluna da matriz A do sistema SL*/
unsigned int findMAX(SistLinear_t *SL, unsigned int lin, unsigned int col){
    real_t max_value = SL->A[lin][col];
    unsigned int max_index = lin;

    for(int i=lin; i< SL->n; i++){
        if(fabs(SL->A[i][col]) > fabs(max_value)){
            max_value = SL->A[i][col];
            max_index = i;
        }
    }

    return max_index;
}

/*Calula o erro absoluto entre dois vetores. Essa função é utilizada
nos métodos Gauss Jacobi e Gauss Seidel. Calcula max(|xi - xi+1|)*/
real_t calculateError(real_t * a, real_t * b, int length){
    real_t * result = malloc(sizeof(real_t)*length);

    for(int i=0; i< length; i++){
        result[i] = fabs(a[i] - b[i]);
    }

    real_t maior = 0;
    //encontra o maior numero no vetor
    for(int j=0; j< length; j++){
        if(result[j] > maior){
            maior = result[j];
        }
    }

    return maior;
}
