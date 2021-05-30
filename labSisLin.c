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
    while(c != EOF){
        if(c != ' '){
            sistLin = lerSistLinear();
            prnSistLinear(sistLin);
            liberaSistLinear(sistLin);            
        }
        //get the \n or EOF
        c = getchar();
    }
}
