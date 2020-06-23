#include "ext-sample/utils.h"
#include <iostream>
#include <cassert>
#include <random>
#include <ctime>
#include <algorithm>


#include "ext-sample/SampleCircuit.h"

using namespace std;

extern "C"{
Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );
Abc_Ntk_t * Abc_NtkDarSeqSweep2( Abc_Ntk_t * pNtk, Ssw_Pars_t * pPars );
Abc_Ntk_t * Abc_NtkDarFraig( Abc_Ntk_t * pNtk, int nConfLimit, int fDoSparse, int fProve, int fTransfer, int fSpeculate, int fChoicing, int fVerbose );
void Io_WriteAiger( Abc_Ntk_t * pNtk, char * pFileName, int fWriteSymbols, int fCompact, int fUnique );
}

Abc_Ntk_t * Ntk_Optimize(Abc_Ntk_t* pNtk)
{
    Ssw_Pars_t Pars, * pPars = &Pars;
    Ssw_ManSetDefaultParams( pPars );
    Abc_Ntk_t * pNtkOpt = Abc_NtkDup( pNtk );
    strcat(pNtkOpt->pName, "Opt");
    for (int j = 0; j < 5; j++)
    {
        if ( !Abc_NtkIsComb(pNtkOpt) )
            pNtkOpt = Abc_NtkDarSeqSweep2( pNtkOpt, pPars );
        pNtkOpt = Abc_NtkDC2( pNtkOpt, 0, 0, 1, 0, 0 );
        pNtkOpt = Abc_NtkDarFraig( pNtkOpt, 100, 1, 0, 0, 0, 0, 0 );
    }    

    return pNtkOpt;
}

Abc_Ntk_t * Ntk_StuckGen(Abc_Ntk_t * pNtk)
{
    int i, nAnd=0, picked;
    Abc_Obj_t * pObj;
    Abc_Ntk_t * pNtkRes;
    vector<Abc_Obj_t*> vAnd;

    assert(Abc_NtkIsStrash(pNtk) && Abc_NtkHasAig(pNtk));
    pNtkRes = Abc_NtkDup(pNtk);

    // const node
    Abc_Obj_t * pAigOne  = Abc_AigConst1(pNtkRes);
    Abc_Obj_t * pAigZero = Abc_ObjNot(pAigOne);
    
    // collect all nodes
    Abc_AigForEachAnd( pNtkRes, pObj, i )
    {
        nAnd++;
        vAnd.push_back(pObj);
    }

    // randomly choose a node
    random_device rd;
    default_random_engine gen = default_random_engine(rd());
    uniform_int_distribution<int> dis(0,nAnd-1);
    picked = dis(gen);

    // stuck to 0 or 1
    srand((unsigned)time(NULL));
    if (rand()%2 == 0)
        Abc_AigReplace( (Abc_Aig_t*)pNtkRes->pManFunc, vAnd[picked], pAigZero, 1 );
    else
        Abc_AigReplace( (Abc_Aig_t*)pNtkRes->pManFunc, vAnd[picked], pAigOne, 1 );

    return pNtkRes;
}

// convert a single-output AIG to BDD and enumerate the minterms
vector<int*> Ntk_Minterm(Abc_Ntk_t * pNtk, int num)
{
    int nPI, nMinterm;
    int count;
    Abc_Obj_t * pObj;
    DdManager * dd;
    DdNode * df;
    DdNode ** nodes;
    vector<int*> vMinterm;

    //assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
    assert(Abc_NtkPoNum(pNtk)==1);
    if (!Abc_NtkIsStrash(pNtk))
        pNtk = Abc_NtkStrash(pNtk, 0, 1, 0);
    nPI = Abc_NtkPiNum(pNtk);

    // convert to BDD
    Abc_NtkBuildGlobalBdds(pNtk, ABC_INFINITY, 1, 1, 0, 0);
    pObj = Abc_NtkPo(pNtk, 0);
    df = (DdNode *)Abc_ObjGlobalBdd(pObj);
    dd = (DdManager *)Abc_NtkGlobalBddMan( pNtk );

    // get witness
    nMinterm = Cudd_CountMinterm(dd, df, nPI);
	nodes = new DdNode*[nPI];
	for (int i = 0; i < nPI; i++)
		nodes[i] = Cudd_ReadVars(dd, i);
    count = 0;
	while (df != Cudd_ReadLogicZero(dd) && count < num)
    {
		DdNode * mintermNode = Cudd_bddPickOneMinterm(dd, df, nodes, nPI);
    	Cudd_Ref(mintermNode);
		int* minterm = new int[nPI];
        for (int i = 0; i < nPI; i++)
            minterm[i] = Cudd_bddLeq(dd, mintermNode, nodes[i]);
        vMinterm.push_back(minterm);
		DdNode *df2 = Cudd_bddAnd(dd, df, Cudd_Not(mintermNode));
    	Cudd_Ref(df2);
    	Cudd_RecursiveDeref(dd, df);
    	df = df2;	
		count++;
    }

    delete [] nodes;
    return vMinterm;
}

/*void genSample(Abc_Ntk_t * pNtk, int hashBits, int loThresh, int hiThresh, Abc_Ntk_t * &pCkt, std::vector<int*> &vSample)
{
    SampleCircuit sc;
	vector<int> V{hashBits, hashBits-1, hashBits-2};
	vector<int*> vMinterm;
    Abc_Ntk_t * pNtkRes;
    int nPI = Abc_NtkPiNum(pNtk);

    random_shuffle(V.begin(), V.end());
    vSample.clear();
	for (int i = 0; i < V.size(); i++)
	{
		sc.setIOnum(nPI-V[i], nPI);
        pCkt = sc.genCircuit();
		pNtkRes = sc.connect(pNtk);
    	vMinterm = Ntk_Minterm(pNtkRes, hiThresh);
		int nSize = vMinterm.size();
		if (nSize >= loThresh && nSize < hiThresh)
		{
			for (int j = 0; j < loThresh; j++)
			{
            	int * sample = Abc_NtkVerifySimulatePattern( pCkt, vMinterm[j] );
			    vSample.push_back(sample);
            }
            Abc_NtkDelete(pNtkRes);
            return;
		}
		Abc_NtkDelete(pCkt);
		Abc_NtkDelete(pNtkRes);
	}
    
    return;
}*/

void split(vector<string> &vs, string str, char delim)
{
    size_t pos;
    string s;

    vs.clear();
    while ( (pos=str.find(delim)) != string::npos)
    {
        if (pos!=0)
        {
            s = str.substr(0, pos);
            vs.push_back(s);
        }
        str = str.substr(pos+1);
    }
    if (str.length()!=0)
        vs.push_back(str);

    return;
}
