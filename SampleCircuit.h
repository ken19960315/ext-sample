#include <iostream>
#include <vector>

#include "base/main/main.h"
#include "base/main/mainInt.h"

using namespace std;

class SampleCircuit{

public:
    SampleCircuit() {f_init=false; f_gen=false; srand((unsigned)time(NULL));}
    SampleCircuit(int nPI, int nPO);

    void setIOnum(int nPI, int nPO);
    void setRndSeed(unsigned seed) {srand(seed);}
    Abc_Ntk_t* genCircuit();
    Abc_Ntk_t* genCircuit(Abc_Ntk_t* pNtk);
    Abc_Ntk_t* genCircuit2(Abc_Ntk_t* pNtk);
    Abc_Ntk_t* connect(Abc_Ntk_t* pNtk, char* pName="connected");

    friend ostream& operator<<(ostream& out, const SampleCircuit& obj);
    vector< vector<int> > getXOR() {return XOR;}
    Abc_Ntk_t* getCircuit() {return pAig;}

private:
    bool f_init;
    bool f_gen;
    int nPI;
    int nPO;
    unsigned seed;
    vector< vector<int> > XOR;
    Abc_Ntk_t* pAig;

    void rndXORGen(int nPI, int nPO, vector< vector<bool> > &Mat, vector<int> &pivot);
    bool checkAvailable(int nPI, int nPO, vector< vector<bool> > &Mat);
};
