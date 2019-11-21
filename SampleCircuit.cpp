#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <random>
#include <cassert>
#include <chrono>

#include <ext-sample/SampleCircuit.h>

using namespace std;

#define SWAP(a, b)  do {a ^= b; b ^= a; a ^= b; } while(0)

SampleCircuit::SampleCircuit(int nPI, int nPO)
{
    assert(nPI>0 && nPO>0 && nPI<nPO);
    this->nPI = nPI;
    this->nPO = nPO;

    f_init = true;
    f_gen = false;
    seed = chrono::system_clock::now().time_since_epoch().count();
    srand(seed);
}

void SampleCircuit::setIOnum(int nPI, int nPO)
{
    assert(nPI>0 && nPO>0 && nPI<nPO);
    this->nPI = nPI;
    this->nPO = nPO;
    
    f_init = true;
}

Abc_Ntk_t* SampleCircuit::genCircuit()
{
    vector<int> pivot;
    vector< vector<bool> > Mat;
    char* pName = "sample";
    Abc_Obj_t * pObj;
    
    // allocate memory
    XOR.resize(nPO);
    Mat.resize(nPO-nPI);

    // initialize
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
        rndXORGen(nPI, nPO, Mat, pivot);
    } while(!checkAvailable(nPI, nPO, Mat));    
    /*cout << "Matrix =\n";
    for (int i = 0; i < nPO-nPI; i++)
    {
        cout << "|";
        for (int j = 0; j < nPO; j++)
            cout << " " << Mat[i][j];
        cout << " | " << Mat[i][nPO];
        cout << " |\n";
    }*/
    
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
            vector<int> XOR_vec(2);
            XOR_vec[0] = 0;
            XOR_vec[1] = mPI;
            XOR[mPO] = XOR_vec;
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
        vector<int> XOR_vec;
        XOR_vec.push_back(Mat[i][nPO]);
        for (int j = *iter+1; j < nPO; j++)
        {
            if (Mat[i][j])
            {
                pObjF = Abc_AigXor((Abc_Aig_t*)pAig->pManFunc, pObjF, Abc_NtkPi(pAig,PImap[j]));
                XOR_vec.push_back(PImap[j]);
            }
        }
        XOR[restPO.front()] = XOR_vec;

        // connect to PO
        Abc_ObjAddFanin(Abc_NtkPo(pAig, restPO.front()), pObjF);
        restPO.pop();

        iter++;
    }

    pivot.clear();
    Mat.clear();
    
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
    return Abc_NtkDup( pAig );
}

bool cmp(vector<int> i, vector<int> j) { return (i.size()>j.size()); }
Abc_Ntk_t* SampleCircuit::genCircuit(Abc_Ntk_t* pNtk)
{
    int i, j;
    vector<int> pivot;
    vector< vector<bool> > Mat;
    char* pName = "sample";
    Abc_Obj_t * pCo, * pObj, * pObjF;
    Vec_Ptr_t * vSupp;
    vector< vector<int> > sup_vec; 
    vector<int> vSuppRef; 
    unordered_map<string, int> sup_map;

    // get support information
    sup_vec.resize(Abc_NtkPoNum(pNtk));
    vSuppRef.resize(Abc_NtkPiNum(pNtk), 0);
    Abc_NtkForEachPi( pNtk, pObj, i )
        sup_map[Abc_ObjName(pObj)] = i;
    Abc_NtkForEachCo( pNtk, pCo, i )
    {
        //cout << Abc_ObjName(pCo) << "\n";
        vSupp  = Abc_NtkNodeSupport( pNtk, &pCo, 1 );
        Vec_PtrForEachEntry( Abc_Obj_t *, vSupp, pObj, j )
        {
            sup_vec[i].push_back(sup_map[Abc_ObjName(pObj)]);
            vSuppRef[sup_map[Abc_ObjName(pObj)]]++;
        }
        //for (vector<int>::iterator iter = sup_vec[i].begin(); iter != sup_vec[i].end(); iter++)
        //    cout << *iter << "(" << Abc_ObjName(Abc_NtkPi(pNtk, *iter)) << ") ";
        //cout << "\n";
        Vec_PtrFree( vSupp );
    }

    // sort by support reference
    for (i = 0; i < sup_vec.size(); i++)
        sort(sup_vec[i].begin(), sup_vec[i].end(), [&vSuppRef](size_t i1, size_t i2) {return vSuppRef[i1] < vSuppRef[i2];});
    // sort by support number
    sort(sup_vec.begin(), sup_vec.end(), cmp);

    // initialize
    pAig = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 ); 
    pAig->pName = Extra_UtilStrsav( pName );
    
    // create PI/PO 
    for ( i = 0; i < nPI; i++ )
        pObj = Abc_NtkCreatePi( pAig );
    for ( i = 0; i < nPO; i++ )
        pObj = Abc_NtkCreatePo( pAig );

    // const node
    Abc_Obj_t * pAigOne  = Abc_AigConst1(pAig);
    Abc_Obj_t * pAigZero = Abc_ObjNot(pAigOne);

    // generate sampling circuit
    assert(f_init);
    vector<bool> vCheck(nPO, false); 
    XOR.resize(nPO);
    while (!sup_vec.empty())
    {
        // sampling the cone with the most supports
        vector<int> supp = sup_vec.front();
        if (supp.size() <= 0) break;
        int nPOc = supp.size();
        sup_vec.erase(sup_vec.begin());
        
        // deal the case nPOc <= nPI
        if (nPOc <= nPI)
        {
            // shuffle PI variable
            vector<int> var(nPI);
            for (i = 0; i < nPI; i++)
                var[i] = i;
            shuffle(var.begin(), var.end(), default_random_engine(seed));
            for (i = 0; i < nPOc; i++)
            {
                vector<int> XOR_vec(2);
                XOR_vec[0] = 0;
                XOR_vec[1] = var[i];
                XOR[supp[i]] = XOR_vec;
                pObjF = Abc_AigAnd((Abc_Aig_t*)pAig->pManFunc, Abc_NtkPi(pAig, var[i]), pAigOne);
                Abc_ObjAddFanin(Abc_NtkPo(pAig, supp[i]), pObjF);
            }
        }
        else
        {
            // allocate memory
            Mat.resize(nPOc-nPI);
            
            // generate XOR constraints
            do{
                rndXORGen(nPI, nPOc, Mat, pivot);
            } while(!checkAvailable(nPI, nPOc, Mat));    
            /*cout << "Matrix =\n";
            for (i = 0; i < nPO-nPI; i++)
            {
                cout << "|";
                for (j = 0; j < nPO; j++)
                    cout << " " << Mat[i][j];
                cout << " | " << Mat[i][nPO];
                cout << " |\n";
            }*/
            
            // create mapping 
            vector<int>::iterator iter = pivot.begin();
            queue<int> restPO;
            unordered_map<int,int> PImap;
            int mPO, mPI=0; 
            for (mPO = 0; mPO < nPOc; mPO++)
            {
                if (iter == pivot.end() || mPO != *iter)
                {
                    assert(mPI < nPI);
                    PImap[mPO] = mPI;
                    // connect to PO
                    vector<int> XOR_vec(2);
                    XOR_vec[0] = 0;
                    XOR_vec[1] = mPI;
                    XOR[supp[mPO]] = XOR_vec;
                    pObjF = Abc_AigAnd((Abc_Aig_t*)pAig->pManFunc, Abc_NtkPi(pAig, mPI++), pAigOne);
                    Abc_ObjAddFanin(Abc_NtkPo(pAig, supp[mPO]), pObjF);
                }
                else
                {
                    restPO.push(supp[mPO]);
                    iter++;
                }
            }
            
            // create internal nodes 
            iter = pivot.begin();
            for (i = 0; i < nPOc-nPI; i++)
            {
                if (Mat[i][nPOc])
                    pObjF = pAigOne;
                else
                    pObjF = pAigZero;
                vector<int> XOR_vec;
                XOR_vec.push_back(Mat[i][nPOc]);
                for (j = *iter+1; j < nPOc; j++)
                {
                    if (Mat[i][j])
                    {
                        pObjF = Abc_AigXor((Abc_Aig_t*)pAig->pManFunc, pObjF, Abc_NtkPi(pAig,PImap[j]));
                        XOR_vec.push_back(PImap[j]);
                    }
                }
                XOR[restPO.front()] = XOR_vec;

                // connect to PO
                Abc_ObjAddFanin(Abc_NtkPo(pAig, restPO.front()), pObjF);
                restPO.pop();

                iter++;
            }

            pivot.clear();
            Mat.clear();
        }
        
        // update check vector
        for (i = 0; i < supp.size(); i++)
            vCheck[supp[i]] = true;

        // remove the connected variables
        for (i = 0; i < sup_vec.size(); i++)
        {
            int erase = 0;
            for (j = 0; j < sup_vec[i].size(); j++)
            {
                if (vCheck[sup_vec[i][j]])
                {
                    if (erase != j)
                        SWAP(sup_vec[i][j], sup_vec[i][erase]);
                    erase++;
                }
            }
            sup_vec[i].erase(sup_vec[i].begin(), sup_vec[i].begin()+erase);
        }
        // sort the rest cone
        sort(sup_vec.begin(), sup_vec.end(), cmp);
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
    return Abc_NtkDup( pAig );
}

vector<size_t> argsort(const vector<int> &v) 
{
    //initialize original index locations
    vector<size_t> idx(v.size());
    for (size_t i = 0; i < idx.size(); ++i) 
        idx[i] = i;

    //sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
    [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}
Abc_Ntk_t* SampleCircuit::genCircuit2(Abc_Ntk_t* pNtk)
{
    int i, j;
    vector<int> pivot;
    vector< vector<bool> > Mat;
    char* pName = "sample";
    Abc_Obj_t * pCo, * pObj, * pObjF;
    Vec_Ptr_t * vSupp;
    vector<int> vSuppRef; 
    unordered_map<string, int> sup_map;

    // get support information
    vSuppRef.resize(Abc_NtkPiNum(pNtk), 0);
    Abc_NtkForEachPi( pNtk, pObj, i )
        sup_map[Abc_ObjName(pObj)] = i;
    Abc_NtkForEachCo( pNtk, pCo, i )
    {
        //cout << Abc_ObjName(pCo) << "\n";
        vSupp = Abc_NtkNodeSupport( pNtk, &pCo, 1 );
        Vec_PtrForEachEntry( Abc_Obj_t *, vSupp, pObj, j )
            vSuppRef[sup_map[Abc_ObjName(pObj)]]++;
        //for (vector<int>::iterator iter = sup_vec[i].begin(); iter != sup_vec[i].end(); iter++)
        //    cout << *iter << "(" << Abc_ObjName(Abc_NtkPi(pNtk, *iter)) << ") ";
        //cout << "\n";
        Vec_PtrFree( vSupp );
    }

    // initialize
    pAig = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 ); 
    pAig->pName = Extra_UtilStrsav( pName );
    
    // create PI/PO 
    for ( i = 0; i < nPI; i++ )
        pObj = Abc_NtkCreatePi( pAig );
    for ( i = 0; i < nPO; i++ )
        pObj = Abc_NtkCreatePo( pAig );

    // const node
    Abc_Obj_t * pAigOne  = Abc_AigConst1(pAig);
    Abc_Obj_t * pAigZero = Abc_ObjNot(pAigOne);

    // generate sampling matrix
    assert(f_init);
    Mat.resize(nPO-nPI);
    do{
        rndXORGen(nPI, nPO, Mat, pivot);
    } while(!checkAvailable(nPI, nPO, Mat));    
    
    // create XOR constraints 
    vector<int>::iterator iter = pivot.begin();
    queue<int> restPO;
    unordered_map<int,int> PImap;
    int mPO, mPI=0; 
    XOR.resize(nPO);
    for (mPO = 0; mPO < nPO; mPO++)
    {
        if (iter == pivot.end() || mPO != *iter)
        {
            assert(mPI < nPI);
            PImap[mPO] = mPI;
            
            vector<int> vXOR(2);
            vXOR[0] = 0;
            vXOR[1] = mPI;
            XOR[mPO] = vXOR;
            mPI++;
        }
        else
        {
            restPO.push(mPO);
            iter++;
        }
    }
    iter = pivot.begin();
    for (i = 0; i < nPO-nPI; i++)
    {
        vector<int> vXOR;
        vXOR.push_back(Mat[i][nPO]);
        for (j = *iter+1; j < nPO; j++)
        {
            if (Mat[i][j])
                vXOR.push_back(PImap[j]);
        }
        XOR[restPO.front()] = vXOR;
        restPO.pop();

        iter++;
    }
    pivot.clear();
    Mat.clear();

    // sort XOR
    sort(XOR.begin(), XOR.end(), cmp);

    // connect from the largest XOR constraint
    vector<size_t> arg = argsort(vSuppRef);
    /*for (i = 0; i < arg.size(); i++)
        cout << vSuppRef[i] << " ";
    cout << "\n";
    for (i = 0; i < arg.size(); i++)
        cout << arg[i] << " ";
    cout << "\n";*/
    assert(arg.size() == XOR.size());
    for (i = 0; i < XOR.size(); i++)
    {
        vector<int> vXOR = XOR[i];
        
        // create AIG internel nodes
        if (vXOR[0])
            pObjF = pAigOne;
        else
            pObjF = pAigZero;
        for (j = 0; j < vXOR.size(); j++)
            pObjF = Abc_AigXor((Abc_Aig_t*)pAig->pManFunc, pObjF, Abc_NtkPi(pAig,vXOR[j]));
        Abc_ObjAddFanin(Abc_NtkPo(pAig, arg[i]), pObjF);
    }
    // swap XOR to correct position
    vector< vector<int> > tmp(XOR.size());
    for (i = 0; i < XOR.size(); i++)
        tmp[arg[i]] = XOR[i];
    XOR = tmp;
    tmp.clear();

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
    return Abc_NtkDup( pAig );
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
    // clone CIs/COs/boxes
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

void SampleCircuit::rndXORGen(int nPI, int nPO, vector< vector<bool> > &Mat, vector<int> &pivot)
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
        Mat[i].resize(nPO+1, 0);
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

    // Gaussian Elimination
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
    out << "XOR constraints:\n";
    for (int i = 0; i < obj.nPO; i++)
    {
        vector<int> vec = obj.XOR[i];
        out << "\tPO" << i << " = ";
        for (vector<int>::iterator iter = vec.begin(); iter != vec.end(); iter++)
        {
            if (iter == vec.begin())
            {
                if (vec.size() == 1) out << *iter;
                else if (*iter == 1) out << "1+";
                continue;
            }
            out << "PI" << *iter;
            if (iter+1 != vec.end()) out << "+";
        }
        out << "\n";
    }

    return out;
}
    
bool SampleCircuit::checkAvailable(int nPI, int nPO, vector< vector<bool> > &Mat)
{
    // check last row
    for (int i = 0; i < nPO; i++)
        if (Mat[nPO-nPI-1][i])
            return true;
    return false;
}
