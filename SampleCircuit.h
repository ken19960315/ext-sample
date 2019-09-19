#include <iostream>
#include <vector>

#include "base/main/main.h"
#include "base/main/mainInt.h"

using namespace std;

class SampleCircuit{

public:
    SampleCircuit() {f_init=false; f_gen=false;}
    SampleCircuit(int nPI, int nPO);
    ~SampleCircuit();

    void setIOnum(int nPI, int nPO);
    void setRndSeed(int seed) {srand(seed);}
    Abc_Ntk_t* genCircuit(char* pName="sample", bool f_verbose=false);
    Abc_Ntk_t* connect(Abc_Ntk_t* pNtk, char* pName="connected");

    friend ostream& operator<<(ostream& out, const SampleCircuit& obj);
    vector<int> getPivot() {return pivot;}
    vector< vector<bool> > getMat() {return Mat;}
    Abc_Ntk_t* getCircuit() {return pAig;}

private:
    bool f_init;
    bool f_gen;
    int nPI;
    int nPO;
    vector<int> pivot;
    vector< vector<bool> > Mat;
    Abc_Ntk_t* pAig;

    void rndXORGen();
    void GaussianElim();
    bool checkAvailable();
};
