#include "ext-sample/utils.h"
#include <iostream>
#include <cassert>
#include <random>
#include <ctime>

#include "base/main/main.h"
#include "base/main/mainInt.h"

using namespace std;

extern "C"{
Abc_Ntk_t * Abc_NtkDC2( Abc_Ntk_t * pNtk, int fBalance, int fUpdateLevel, int fFanout, int fPower, int fVerbose );
Abc_Ntk_t * Abc_NtkDarSeqSweep2( Abc_Ntk_t * pNtk, Ssw_Pars_t * pPars );
Abc_Ntk_t * Abc_NtkDarFraig( Abc_Ntk_t * pNtk, int nConfLimit, int fDoSparse, int fProve, int fTransfer, int fSpeculate, int fChoicing, int fVerbose );
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
        Abc_AigReplace( (Abc_Aig_t*)pNtkRes->pManFunc, vAnd[picked], pAigZero, 0 );
    else
        Abc_AigReplace( (Abc_Aig_t*)pNtkRes->pManFunc, vAnd[picked], pAigOne, 0 );

    return pNtkRes;
}

// convert a single-output AIG to BDD and enumerate the minterms
vector<vector<bool>> Ntk_Minterm(Abc_Ntk_t * pNtk, int nSample)
{
    int nPI, nMinterm;
    int count;
    Abc_Obj_t * pObj;
    DdNode * df;
    DdNode ** nodes;
	DdManager * dd;
    vector<vector<bool>> vMinterm;

    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
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
	while (df != Cudd_ReadLogicZero(dd) && count < nSample)
    {
		DdNode * mintermNode = Cudd_bddPickOneMinterm(dd, df, nodes, nPI);
    	Cudd_Ref(mintermNode);
		vector<bool> minterm;
        for (int i = 0; i < nPI; i++)
            minterm.push_back(Cudd_bddLeq(dd, mintermNode, nodes[i]));
        vMinterm.push_back(minterm);
		DdNode *df2 = Cudd_bddAnd(dd, df, Cudd_Not(mintermNode));
    	Cudd_Ref(df2);
    	Cudd_RecursiveDeref(dd, df);
    	df = df2;	
		count++;
    }

    Abc_NtkFreeGlobalBdds( pNtk, 1 );

    return vMinterm;
}
