#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    SistLinear_t * sistLin = lerSistLinear();
    prnSistLinear(sistLin);
    liberaSistLinear(sistLin);
}
