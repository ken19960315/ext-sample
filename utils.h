#include "base/abc/abc.h"

// optimize subroutine (dc2; dfraig * 5)
Abc_Ntk_t * Ntk_Optimize(Abc_Ntk_t * pNtk);

// generate single stuck-at fault
Abc_Ntk_t * Ntk_StuckGen(Abc_Ntk_t * pNtk);
