#include "SistemasLineares.h"

/*Escreve na tela a soluação encontrada e informações sobre a execucao
dos métodos*/
void prnSolucao(real_t * solucao, int n, double tempo, int iteracoes, double norma);

/*Soma dois vetores de mesmo tamanho. O resultado é guardado em
um terceiro vetor de mesmo tipo*/
real_t * sumVector(real_t * a, real_t * b, unsigned int size);

/*Copia os elementos de source para o vetor destiny.
Os vetores precisam ter o mesmo tamanho*/
void cpyVector(real_t * destiny, real_t * source, int size);

/*Copia os elementos de source para a matriz destiny
As matrizes precisam ser quadradas e ter mesmo tamanho*/
void cpyMatrix(real_t ** destiny, real_t ** source, int size);

/*Copia os atributos do sistema source para o sistema destiny.
As matrizes e vetores precisam ter o mesmo tamanho*/
void cpySist(SistLinear_t * destiny, SistLinear_t * source);

/*Troca duas linhas de um sistema linear. Essa operação inclui
trocas na matriz A e no vetor b*/
void swapLine(SistLinear_t * SL, int i, int j);

/*Recebe uma matriz triangular superior e aliza a retro substituição
dos termos.*/
int retroSubst(SistLinear_t *SL, real_t * x, double * tTotal);

/*Encontra maior valor em dada coluna da matriz A do sistema SL*/
unsigned int findMAX(SistLinear_t *SL, unsigned int lin, unsigned int col);

/*Calula o erro absoluto entre dois vetores. Essa função é utilizada
nos métodos Gauss Jacobi e Gauss Seidel*/
real_t calculateError(real_t * a, real_t * b, int length);
