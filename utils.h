#include "base/abc/abc.h"
#include <vector>

// optimize subroutine (dc2; dfraig * 5)
Abc_Ntk_t * Ntk_Optimize(Abc_Ntk_t * pNtk);

// generate single stuck-at fault
Abc_Ntk_t * Ntk_StuckGen(Abc_Ntk_t * pNtk);

// enumerate a certain number of minterms
std::vector<std::vector<bool>> Ntk_Minterm(Abc_Ntk_t * pNtk, int nSample);
