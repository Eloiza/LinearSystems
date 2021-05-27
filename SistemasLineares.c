#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{

};


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTriangulariza tempo gasto na triangularização
  \param tRetroSubst tempo gasto na retrosubstituição

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{


};

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{


};

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo gasto pelo método
  \param tIteração tempo gasto em cada iteração

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{


};


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{


};

/*!
  \brief Alocaçao de memória

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n){
    SistLinear_t * sistLin = malloc(sizeof(SistLinear_t));
    if(!sistLin){
        printf("Error to alocate memory for sistLinear_t\n");
    }
    sistLin->n = n;
    sistLin->erro = 0;
    sistLin->b = malloc(sizeof(real_t)*n);

    sistLin->A = malloc(sizeof(real_t*)*n);
    for(int i=0; i<n; i++){
        sistLin->A[i] = malloc(sizeof(real_t)*n);
    }

    return sistLin;
};

/*!
  \brief Liberaçao de memória

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{


};

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{

};


// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL){
    int i;
    printf("n: %d \nerro: %lf \n", SL->n, SL->erro);

    // print the independet terms vector
    printf("b: ");
    prnVetor(SL->b, SL->n);

    //print the coeficients matrix
    printf("A:\n");
    for(i=0; i<SL->n; i++){
        prnVetor(SL->A[i], SL->n);
    }
};

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n){
    int i = 0;
    for(i=0; i<n; i++){
        //case last position
        if(i == n - 1){
            printf("%lf ", v[i]);

        }else{
            printf("%lf, ", v[i]);
        }
    }
    printf("\n");
};
