#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"
#include "SupportFunctions.h"

#define BUFFER_SIZE 300


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

    // SistLinear_t * SL_copy = alocaSistLinear(SL->n);
    // cpySist(SL_copy, SL);

    double start_t = timestamp();
    for(int i=0; i<SL->n; i++){
        iPivo = findMAX(SL, i, i);

        //case pivo = 0 go to next iteration
        if(!SL->A[iPivo][i]){
            continue;

        }else if(i != iPivo){
            swapLine(SL, i, iPivo);
        }

        for(int k= i+1; k< SL->n; k++){
            if(!SL->A[i][i]){
                fprintf(stderr, "%s", "eliminacaoGauss - Division by 0\n");
                return 1;
            }

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

    *tTotal = timestamp() - start_t;

    return 0;
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
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal){
    unsigned int it_count = 0;

    //allocate memory and set to 0
    real_t * x_old = calloc(SL->n, sizeof(real_t));
    real_t * x_new = calloc(SL->n, sizeof(real_t));

    double soma;
    real_t it_error;
    double start_t = timestamp();
    for(it_count=0; it_count< MAXIT; it_count++){

        for(int i=0; i< SL->n; i++){
            soma = 0;
            for(int j=0; j<SL->n; j++){
                if(j != i){
                    soma += (SL->A[i][j] * x_old[j]) / SL->A[i][i];
                }
                x_new[i] = (SL->b[i] / SL->A[i][i]) - soma;
            }
        }

        it_error = calculateError(x_new, x_old, SL->n);
        if(it_error <= SL->erro){
            break;
        }

        //x_old receives x_new values
        cpyVector(x_old, x_new, SL->n);
    }

    *tTotal = timestamp() - start_t;
    cpyVector(x, x_new, SL->n);

    return it_count;
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
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal){
    unsigned int it_count;

    real_t * x_old = calloc(SL->n, sizeof(real_t));
    real_t * x_new = calloc(SL->n, sizeof(real_t));

    real_t soma, it_error;
    double start_t = timestamp();
    for(it_count=0; it_count < MAXIT; it_count++){
        for(int i=0; i< SL->n; i++){
            soma = 0;
            for(int j=0; j< i; j++){
                soma += SL->A[i][j] * x_new[j];
            }
            for(int j= i+1; j< SL->n; j++){
                soma += SL->A[i][j] * x_old[j];
            }

            x_new[i] = (SL->b[i] - soma) / SL->A[i][i];
        }

        it_error = calculateError(x_new, x_old, SL->n);
        if(it_error <= SL->erro){
            break;
        }

        cpyVector(x_old, x_new, SL->n);
    }

    *tTotal = timestamp() - start_t;

    cpyVector(x, x_new, SL->n);

    return it_count;
};

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res){

    //se o vetor res não for fornecido, aloca memoria
    if(!res){
        res = malloc(sizeof(real_t)*SL->n);
    }

    real_t soma = 0;

    //obtem o vetor residuo para a solução
    for(int i=0; i<SL->n; i++){
        soma = 0;
        for(int j=0; j<SL->n; j++){
            soma += SL->A[i][j] * x[j];
            // printf("soma(%f) = A[%i][%i](%f) * x[%i](%f) \n", soma, i, j, SL->A[i][j], i, x[i]);
        }
        res[i] = SL->b[i] - soma;
        // printf("res[%i] = b[%i](%f) - soma(%f)\n", i, i, SL->b[i], soma);
        // printf("res[%i] = %f\n", i, res[i]);
    }

    //calcula a norma do vetor residuo
    real_t norma = 0;
    for(int i=0; i< SL->n; i++){
        norma += pow(res[i], 2);
    }

    norma = sqrt(norma);

    return norma;
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
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal){
    real_t * w = calloc(SL->n, sizeof(real_t));
    real_t * r = calloc(SL->n, sizeof(real_t));

    real_t * x_old = calloc(SL->n, sizeof(real_t));
    real_t * x_new = calloc(SL->n, sizeof(real_t));

    cpyVector(x_new, x, SL->n);

    int it_count = 0;
    double seidel_t, start_t = timestamp();
    real_t norma, it_erro;
    for(it_count=0; it_count<MAXIT; it_count++){

        //calcular residuo e testa primeira condicao de parada
        norma = normaL2Residuo(SL, x_new, r);
        // printf("norma (%lf) <= erro(%lf)\n", norma, SL->erro);
        if(norma <= SL->erro){
            break;
        }

        //obter w resolvendo AW = r
        SL->b = r;
        // prnVetor(r, SL->n);
        // prnSistLinear(SL);
        eliminacaoGauss(SL, w, &seidel_t);

        //obter nova solucao x(i) + w
        x_new = sumVector(x_old, w, SL->n);

        //testar segundo criterio de parada
        it_erro = calculateError(x_old, x_new, SL->n);
        // printf("it_erro (%lf) <= erro(%lf) \n",it_erro, SL->erro);
        if(calculateError(x_old, x_new, SL->n) <= SL->erro){
            break;
        }

        //x_old = x_new
        cpyVector(x_old, x_new, SL->n);
    }

    *tTotal = timestamp() - start_t;

    cpyVector(x, x_new, SL->n);

    return it_count;
};

/*!
  \brief Alocaçao de memória

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n){
    SistLinear_t * sistLin = malloc(sizeof(SistLinear_t));
    if(!sistLin){
        fprintf(stderr, "%s", "Error to alocate memory for sistLinear_t\n");
        return NULL;
    }
    sistLin->n = n;
    sistLin->erro = 0;
    sistLin->b = malloc(sizeof(real_t)*n);

    //case of error in memory allocation
    if(!sistLin->b){
        fprintf(stderr, "%s", "Error to allocate memory for vector b!\n");
        return NULL;
    }

    sistLin->A = malloc(sizeof(real_t*)*n);

    if(!sistLin->A){
        fprintf(stderr, "%s", "Error to allocate memory for vector A!\n");
        return NULL;
    }

    for(int i=0; i<n; i++){
        sistLin->A[i] = malloc(sizeof(real_t)*n);

        if(!sistLin->A[i]){
            fprintf(stderr, "%s", "Error to allocate memory for vector A!\n");
            return NULL;
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
    printf("n: %d \nerro: %1.9e \n", SL->n, SL->erro);

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
            printf("%1.9e ", v[i]);

        }else{
            printf("%1.9e, ", v[i]);
        }
    }
    printf("\n");
};
