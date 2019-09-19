#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <queue>
#include <cassert>

#include <ext-sample/SampleCircuit.h>

using namespace std;

#define SWAP(a, b)  do {a ^= b; b ^= a; a ^= b; } while(0)

SampleCircuit::SampleCircuit(int nPI, int nPO)
{
    assert(nPI>0 && nPO>0 && nPI<nPO);
    this->nPI = nPI;
    this->nPO = nPO;

    Mat.resize(nPO-nPI);
    for (int i = 0; i < nPO-nPI; i++)
        Mat[i].resize(nPO+1, 0);

    f_init = true;
    f_gen = false;
    srand((unsigned)time(NULL));
}

SampleCircuit::~SampleCircuit()
{
    Abc_NtkDelete( pAig );
    pivot.clear();
    Mat.clear();
}

void SampleCircuit::setIOnum(int nPI, int nPO)
{
    assert(nPI>0 && nPO>0 && nPI<nPO);
    this->nPI = nPI;
    this->nPO = nPO;

    Mat.resize(nPO-nPI);
    for (int i = 0; i < nPO-nPI; i++)
        Mat[i].resize(nPO+1, 0);

    f_init = true;
}

Abc_Ntk_t* SampleCircuit::genCircuit(char* pName, bool f_verbose)
{
    Abc_Obj_t * pObj;
    
    if ( pAig )
        Abc_NtkDelete( pAig );
    pAig = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 ); 
    pAig->pName = Extra_UtilStrsav( pName );
    
    // create PI/PO 
    for ( int i = 0; i < nPI; i++ )
        pObj = Abc_NtkCreatePi( pAig );
    for ( int i = 0; i < nPO; i++ )
        pObj = Abc_NtkCreatePo( pAig );

    // const node
    Abc_Obj_t * pAigOne  = Abc_AigConst1(pAig);
    Abc_Obj_t * pAigZero = Abc_ObjNot(pAigOne);

    // generate XOR constraints
    assert(f_init);
    do{
        rndXORGen();
        GaussianElim();
    } while(!checkAvailable());    
    if (f_verbose)
    {
        cout << "After Gaussian Elimination, XOR constraints =\n";
        cout << *this;
    }
    
    // create mapping 
    vector<int>::iterator iter = pivot.begin();
    queue<int> restPO;
    unordered_map<int,int> PImap;
    int mPO, mPI=0; 
    Abc_Obj_t * pObjF;
    for (mPO = 0; mPO < nPO; mPO++)
    {
        if (iter == pivot.end() || mPO != *iter)
        {
            assert(mPI < nPI);
            PImap[mPO] = mPI;
            // connect to PO
            pObjF = Abc_AigAnd((Abc_Aig_t*)pAig->pManFunc, Abc_NtkPi(pAig, mPI++), pAigOne);
            Abc_ObjAddFanin(Abc_NtkPo(pAig, mPO), pObjF);
        }
        else
        {
            restPO.push(mPO);
            iter++;
        }
    }
    
    // create internal nodes 
    iter = pivot.begin();
    for (int i = 0; i < nPO-nPI; i++)
    {
        if (Mat[i][nPO])
            pObjF = pAigOne;
        else
            pObjF = pAigZero;
        for (int j = *iter+1; j < nPO; j++)
        {
            if (Mat[i][j])
                pObjF = Abc_AigXor((Abc_Aig_t*)pAig->pManFunc, pObjF, Abc_NtkPi(pAig,PImap[j]));
        }
        
        // connect to PO
        Abc_ObjAddFanin(Abc_NtkPo(pAig, restPO.front()), pObjF);
        restPO.pop();

        iter++;
    }

    //
    Abc_AigCleanup( (Abc_Aig_t*)pAig->pManFunc );
    Abc_NtkAddDummyPiNames( pAig );
    Abc_NtkAddDummyPoNames( pAig );
    Abc_NtkAddDummyBoxNames( pAig );

    // return constructed AIG
    if ( !Abc_NtkCheck( pAig ) )
    {
        printf( "The AIG construction has failed.\n" );
        Abc_NtkDelete( pAig );
        return NULL;
    }
    f_gen = true;
    return pAig;
}

Abc_Ntk_t* SampleCircuit::connect(Abc_Ntk_t* pNtk, char* pName)
{
    assert(pNtk && Abc_NtkHasAig(pNtk) && Abc_NtkIsStrash(pNtk));
    assert(f_gen);

    int i;
    Abc_Ntk_t * pAigNew; 
    Abc_Obj_t * pObj; 
    Abc_Obj_t * pChild0; 
    Abc_Obj_t * pChild1; 
    unordered_map<int, int> IDmap;
    
    // create new network
    pAigNew = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 ); 
    pAigNew->pName = Extra_UtilStrsav( pName );

    // map constant nodes
    Abc_AigConst1(pAig)->pCopy = Abc_AigConst1(pAigNew);
    Abc_AigConst1(pNtk)->pCopy = Abc_AigConst1(pAigNew);
    // clone CIs/CIs/boxes
    Abc_NtkForEachPi( pAig, pObj, i )
        Abc_NtkDupObj( pAigNew, pObj, 1 );
    Abc_NtkForEachPo( pNtk, pObj, i )
        Abc_NtkDupObj( pAigNew, pObj, 1 );
    Abc_NtkForEachBox( pNtk, pObj, i )
        Abc_NtkDupBox( pAigNew, pObj, 1 );
    //
    Abc_NtkForEachPi( pNtk, pObj, i )
        IDmap[pObj->Id] = i;
    // copy the AND gates
    Abc_AigForEachAnd( pAig, pObj, i )
        pObj->pCopy = Abc_AigAnd( (Abc_Aig_t *)pAigNew->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj) );
    Abc_AigForEachAnd( pNtk, pObj, i )
    {
        if (Abc_ObjIsPi(Abc_ObjFanin0(pObj)))
            pChild0 = Abc_ObjNotCond( Abc_ObjFanin0(Abc_NtkPo(pAig,IDmap[Abc_ObjFanin0(pObj)->Id]))->pCopy, Abc_ObjFaninC0(pObj) );
        else
            pChild0 = Abc_ObjChild0Copy(pObj);
        if (Abc_ObjIsPi(Abc_ObjFanin1(pObj)))
            pChild1 = Abc_ObjNotCond( Abc_ObjFanin0(Abc_NtkPo(pAig,IDmap[Abc_ObjFanin1(pObj)->Id]))->pCopy, Abc_ObjFaninC1(pObj) );
        else
            pChild1 = Abc_ObjChild1Copy(pObj);
        pObj->pCopy = Abc_AigAnd( (Abc_Aig_t *)pAigNew->pManFunc, pChild0, pChild1 );
    }
    // relink CO
    Abc_NtkForEachCo( pNtk, pObj, i )
    {
        if (Abc_ObjIsPi(Abc_ObjFanin0(pObj)))
            pChild0 = Abc_ObjNotCond( Abc_ObjFanin0(Abc_NtkPo(pAig,IDmap[Abc_ObjFanin0(pObj)->Id]))->pCopy, Abc_ObjFaninC0(pObj) );
        //else if (Abc_AigNodeIsConst(Abc_ObjFanin0(pObj)))
        //    pChild0 = Abc_ObjNotCond( Abc_AigConst1(pAigNew), Abc_ObjFaninC0(pObj) );
        else
            pChild0 = Abc_ObjChild0Copy(pObj);
        Abc_ObjAddFanin( pObj->pCopy, pChild0 );
    }

    //
    Abc_AigCleanup( (Abc_Aig_t*)pAigNew->pManFunc );
    // return constructed AIG
    if ( !Abc_NtkCheck( pAigNew ) )
    {
        printf( "The AIG construction has failed.\n" );
        Abc_NtkDelete( pAigNew );
        return NULL;
    }
    return pAigNew;
}

void SampleCircuit::rndXORGen()
{
    int nXOR, pos;

    // for shuffle
    vector<int> seq(nPO);
    for (int i = 0; i < nPO; i++)
        seq[i] = i;
    
    // generate XOR constraints
    for (int i = 0; i < nPO-nPI; i++)
    {
        nXOR = nPO * (rand()/(RAND_MAX+1.0)) + 1;
        for (int j = 0; j < nPO; j++)
        {
            pos = nPO * (rand()/(RAND_MAX+1.0));
            if (j != pos)
                SWAP(seq[j], seq[pos]);
        }
        for (int j = 0; j < nXOR; j++)
            Mat[i][seq[j]] = 1;

        Mat[i][nPO] = rand()%2;
    }

    return;
}

void SampleCircuit::GaussianElim()
{
    pivot.clear();
    int row = 0;
    for (int i = 0; i < nPO && row < nPO-nPI; i++)
    {
        // find pivot
        if (!Mat[row][i])
        {
            for (int j = row+1; j < nPO-nPI; j++)
            {
                // swap
                if (Mat[j][i])
                {
                    Mat[row].swap(Mat[j]);
                    break;
                }
            }
            if (!Mat[row][i])
                continue;
        }

        // reduction
        for (int j = 0; j < nPO-nPI; j++)
        {
            if (Mat[j][i] && j != row)
            {
                for (int k = i; k <= nPO; k++)
                    Mat[j][k] = Mat[j][k]^Mat[row][k];
            }
        }
        row += 1;
        pivot.push_back(i);
    }

    return;
}
    
ostream& operator<<(ostream& out, const SampleCircuit& obj)
{
    for (int i = 0; i < obj.nPO-obj.nPI; i++)
    {
        out << "|";
        for (int j = 0; j < obj.nPO; j++)
            out << " " << obj.Mat[i][j];
        out << " | " << obj.Mat[i][obj.nPO];
        out << " |\n";
    }

    return out;
}
    
bool SampleCircuit::checkAvailable()
{
    // check last row
    for (int i = 0; i < nPO; i++)
        if (Mat[nPO-nPI-1][i])
            return true;
    return false;
}
