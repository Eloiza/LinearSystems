#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"

#define BUFFER_SIZE 300


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

//Realiza retro substituição.
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

unsigned int findMAX(SistLinear_t *SL, unsigned int col){
    real_t max_value = SL->A[0][col];
    unsigned int max_index = 0;

    for(int i=0; i< SL->n; i++){
        if(fabs(SL->A[i][col]) > fabs(max_value)){
            max_value = SL->A[i][col];
            max_index = i;
        }
    }

    return max_index;
}

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
  \param tTotal tempo gasto na triangularização
  \param tRetroSubst tempo gasto na retrosubstituição

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal){
    unsigned int iPivo; //stores the pivo index
    double m;

    for(int i=0; i<SL->n; i++){
        iPivo = findMAX(SL, i);

        if(i != iPivo){
            swapLine(SL, i, iPivo);
        }

        for(int k= i+1; k< SL->n; k++){
            m = SL->A[k][i] / SL->A[i][i];

            SL->A[k][i] = 0;
            for(int j= i+1; j<SL->n; j++){
                SL->A[k][j] = SL->A[k][j] - (SL->A[i][j] * m);
            }

            SL->b[k] = SL->b[k] - (SL->b[i] * m);
        }
    }

    double tRetro = 0;
    retroSubst(SL, x, &tRetro);

    return 0;
};

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
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal){
    unsigned int it_count = 0;

    //allocate memory and set to 0
    real_t * x_old = calloc(SL->n, sizeof(real_t));
    real_t * x_new = calloc(SL->n, sizeof(real_t));

    double soma;
    real_t it_error;
    for(it_count=0; it_count< MAXIT; it_count++){

        printf("X_old: ");
        prnVetor(x_old, SL->n);

        printf("X_new: ");
        prnVetor(x_new, SL->n);


        for(int i=0; i< SL->n; i++){
            soma = 0;
            for(int j=0; j<SL->n; j++){
                if(j != i){
                    soma += (SL->A[i][j] * x_old[j]) / SL->A[i][i];
                }
                x_new[i] = (SL->b[i] / SL->A[i][i]) - soma;
            }
        }

        printf("X_old: ");
        prnVetor(x_old, SL->n);

        printf("X_new: ");
        prnVetor(x_new, SL->n);

        it_error = calculateError(x_new, x_old, SL->n);
        printf("IT_ERROR: %f\n", it_error);
        if(it_error <= SL->erro){
            break;
        }

        //x_old receives x_new values
        memcpy(x_old, x_new, SL->n);

        printf("Iterações para convergencia: %i\n", it_count);
    }


    x = x_new;
    return 0;
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

    //case of error in memory allocation
    if(!sistLin->b){
        fprintf(stderr, "%s", "Error to allocate memory for vector b!\n");
        exit(-1);
    }

    sistLin->A = malloc(sizeof(real_t*)*n);

    if(!sistLin->A){
        fprintf(stderr, "%s", "Error to allocate memory for vector A!\n");
        exit(-1);
    }

    for(int i=0; i<n; i++){
        sistLin->A[i] = malloc(sizeof(real_t)*n);

        if(!sistLin->A[i]){
            fprintf(stderr, "%s", "Error to allocate memory for vector A!\n");
            exit(-1);
        }

    }

    return sistLin;
};

/*!
  \brief Liberaçao de memória

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL){
    int i=0;

    //free the matrix lines
    for(i=0; i<SL->n; i++){
        free(SL->A[i]);
    }

    //free the left pointers
    free(SL->A);
    free(SL->b);

    //set to 0 the other atributes
    SL->erro = 0;
    SL->n = 0;

    //set to NULL all pointers
    SL->A = NULL;
    SL->b = NULL;

    //free all the structure
    free(SL);
};

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear (){
    int n;
    real_t error;

    scanf("%i", &n);
    scanf("%f", &error);
    //
    // printf("n: %d\n", n);
    // printf("error: %f\n", error);

    //get the \n
    getchar();

    SistLinear_t * sistLin = alocaSistLinear(n);
    sistLin->erro = error;

    char buffer[BUFFER_SIZE];
    char * token;
    int index = 0;

    for(int i=0; i< n + 1; i++){
        index = 0;
        //read a line and store in buffer
        if(fgets(buffer, BUFFER_SIZE, stdin) == NULL){
            fprintf(stderr, "%s", "Error to read a line!\n");
            exit(-1);

        }else{
            //get one digit from the line
            token = strtok(buffer, " ");
            while(token != NULL){
                if(i == n){
                    sistLin->b[index] = atof(token);
                    // printf("sistLin->b[%i] = %f\n", index, atof(token));
                }
                else{
                    sistLin->A[i][index] = atof(token);
                    // printf("sistLin->A[%i][%i] = %f\n", i, index, atof(token));
                }

                index++;
                token = strtok(NULL, " ");
            }
        }//end else
    }//end for

    return sistLin;
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

    printf("\n");
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
