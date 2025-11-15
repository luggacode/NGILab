#include <stdio.h>
#include </opt/anaconda3/envs/Brian2_NGILab/include/gsl/gsl_const_mksa.h>

int main(void) {
    const double R = GSL_CONST_MKSA_MOLAR_GAS;
    const double F = GSL_CONST_MKSA_FARADAY;

    printf("Gas constant R = %.15e J/(mol*K)\n", R);
    printf("Faraday constant F = %.15e C/mol\n", F);

    return 0;
}
