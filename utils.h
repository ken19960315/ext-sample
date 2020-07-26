#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include <vector>
#include <string>

// optimize subroutine (dc2; dfraig * 5)
Abc_Ntk_t * Ntk_Optimize(Abc_Ntk_t * pNtk);

// generate single stuck-at fault
Abc_Ntk_t * Ntk_StuckGen(Abc_Ntk_t * pNtk);

// split the string by delim
void split(std::vector<std::string> &vs, std::string str, char delim=' ');
