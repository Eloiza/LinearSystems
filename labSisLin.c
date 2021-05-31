#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main (){
    //allcoate memory for a  vector of linear systems
    // SistLinear_t ** sistLin_v = (SistLinear_t **) malloc(sizeof(SistLinear_t));

    char c = 'a';
    SistLinear_t * sistLin;
    real_t * solution = malloc(sizeof(int)*3);
    double totalTime;
    while(c != EOF){
        if(c != ' '){
            sistLin = lerSistLinear();
            prnSistLinear(sistLin);
            eliminacaoGauss(sistLin, solution, &totalTime);
            printf("SOLUCAO:\n");
            prnVetor(solution, sistLin->n);

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
