#include "base/abc/abc.h"
#include <vector>
#include <string>

// optimize subroutine (dc2; dfraig * 5)
Abc_Ntk_t * Ntk_Optimize(Abc_Ntk_t * pNtk);

// generate single stuck-at fault
Abc_Ntk_t * Ntk_StuckGen(Abc_Ntk_t * pNtk);

// convert to BDD and enumerate a certain number of minterms
std::vector<int*> Ntk_Minterm(Abc_Ntk_t * pNtk, int num);

// generate samples
int genSample(Abc_Ntk_t * pNtk, int hashBits, int loThresh, int hiThresh, std::vector<int*> &vSample);

// split the string by delim
void split(std::vector<std::string> &vs, std::string str, char delim=' ');
