#include <iostream>
#include <vector>

#include "base/main/main.h"
#include "base/main/mainInt.h"

using namespace std;

class SampleCircuit{

public:
    ~SampleCircuit();
    SampleCircuit();
    SampleCircuit(int nPI, int nPO);

    void setIOnum(int nPI, int nPO);
    void setRndSeed(unsigned seed) {srand(seed);}
    Abc_Ntk_t* genCircuit(bool fVerbose=false);
    Abc_Ntk_t* genCircuit(Abc_Ntk_t* pNtk, bool fVerbose=false);
    //Abc_Ntk_t* genCircuit2(Abc_Ntk_t* pNtk); //deprecated
    Abc_Ntk_t* connect(Abc_Ntk_t* pCkt, Abc_Ntk_t* pNtk, char* pName="sampled");

    friend ostream& operator<<(ostream& out, const SampleCircuit& obj);
    vector< vector<int> > getXOR() {return XOR;}

private:
    bool fInit;
    bool fGen;
    int nPI;
    int nPO;
    unsigned seed;
    vector< vector<int> > XOR;

    void rndXORGen(int nPI, int nPO, vector< vector<bool> > &Mat, vector<int> &pivot);
    bool checkAvailable(int nPI, int nPO, vector< vector<bool> > &Mat);
};
