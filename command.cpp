#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "proof/fraig/fraig.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <map>
#include <unordered_map>
#include <set>
#include <bitset>
#include <algorithm>
#include <numeric>

#include "ext-sample/SampleCircuit.h"
#include "ext-sample/chisqr.h"
#include "ext-sample/utils.h"


extern "C" void Abc_NtkPrintStrSupports( Abc_Ntk_t * pNtk, int fMatrix );

namespace
{

// ABC command: print/write network infomation
int Info_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    char * filename;
    bool fIO = false;
    bool fSize = false;
    bool fRedirect = false;
    fstream file;
    Abc_Ntk_t * pNtk;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "isrh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            fIO = true;
            break;
        case 's':
            fSize = true;
            break;
        case 'r':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-r\" should be followed by an string.\n" );
                goto usage;
            }
            fRedirect = true;
            filename = argv[globalUtilOptind];
            assert(filename != NULL);
            globalUtilOptind++;
            break;
        case 'h':
            goto usage;
        default:
            goto usage;
        }
    }

    // get the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));

    // print the information
    if (fRedirect)
    {
        file.open(filename, ios::out|ios::trunc);
        assert(file.is_open());
        if (fIO)
            file << Abc_NtkPiNum(pNtk) << " " << Abc_NtkPoNum(pNtk) << "\n";
        if (fSize)
        {
            if (Abc_NtkHasAig(pNtk))
                file << Abc_NtkNodeNum(pNtk) << "\n";
            else if (Abc_NtkHasBdd(pNtk))
                file << Abc_NtkGetBddNodeNum(pNtk) << "\n";
        }
        file.close();
    }
    else
    {
        if (fIO)
            cout << Abc_NtkPiNum(pNtk) << " " << Abc_NtkPoNum(pNtk) << "\n";
        if (fSize)
        {
            if (Abc_NtkHasAig(pNtk))
                cout << Abc_NtkNodeNum(pNtk) << "\n";
            else if (Abc_NtkHasBdd(pNtk))
                cout << Abc_NtkGetBddNodeNum(pNtk) << "\n";
        }
    }

    return 0;

usage:
    Abc_Print( -2, "usage: info [-r <file>] [-ish]\n" );
    Abc_Print( -2, "\t        Print out the information of the current network\n" );
    Abc_Print( -2, "\t-r <file> : redirect the result to the given file\n");
    Abc_Print( -2, "\t-i        : print the IO num\n");
    Abc_Print( -2, "\t-s        : print the size\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: print/write network infomation
int ReadCNF_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
	char c;
	char * filename;
	char * pCircuitName = "transformedCNF";
	string line;
	int nVar, var;
	bool fVarInit = false;
	fstream fin;
	vector<string> vSplit;
	Abc_Ntk_t * pAig;
	Abc_Obj_t * pAigOne, * pAigZero;
	Abc_Obj_t * pObjA, * pObjClause, * pObjF;

	// parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "ih" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an string.\n" );
                goto usage;
            }
            filename = argv[globalUtilOptind];
            assert(filename != NULL);
            globalUtilOptind++;
            break;
        case 'h':
            goto usage;
        default:
            goto usage;
        }
    }

	// initialize AIG
    pAig = Abc_NtkAlloc( ABC_NTK_STRASH, ABC_FUNC_AIG, 1 ); 
    pAig->pName = Extra_UtilStrsav( pCircuitName );

    // const node
    pAigOne  = Abc_AigConst1(pAig);
    pAigZero = Abc_ObjNot(pAigOne);
	pObjF = pAigOne;
	
	// parse CNF & construct AIG
    fin.open(filename, ios::in);
    assert(fin.is_open());
	while (getline(fin, line))
	{
		if (line[0] == 'c')
			continue;
		else if (line[0] == 'p')
		{
			assert(!fVarInit);
			split(vSplit, line);
			nVar = stoi(vSplit[2]); 
    		// create PI/PO 
    		for ( int i = 0; i < nVar; i++ )
    		    pObjA = Abc_NtkCreatePi( pAig );
    		pObjA = Abc_NtkCreatePo( pAig );
			fVarInit = true;
		}
		else
		{
			assert(fVarInit);
			split(vSplit, line);
			if (vSplit.size() == 0) break;
			pObjClause = pAigZero;
			for (int i = 0; i < vSplit.size()-1; i++)
			{
				var = stoi(vSplit[i]);
				pObjA = Abc_NtkPi(pAig, abs(var)-1);
				if (var < 0)
					pObjA = Abc_ObjNot(pObjA);
				pObjClause = Abc_AigOr( (Abc_Aig_t*)pAig->pManFunc, pObjClause, pObjA );
			}
			pObjF = Abc_AigAnd( (Abc_Aig_t*)pAig->pManFunc, pObjF, pObjClause );
		}
	}

	// connect to PO
	Abc_ObjAddFanin(Abc_NtkPo(pAig, 0), pObjF);

	// remove dangling nodes
    Abc_AigCleanup( (Abc_Aig_t*)pAig->pManFunc );
    
	// add names
	Abc_NtkAddDummyPiNames( pAig );
	Abc_NtkAddDummyPoNames( pAig );

	// check network
    if ( !Abc_NtkCheck( pAig ) )
    {
        printf( "The AIG construction has failed.\n" );
        Abc_NtkDelete( pAig );
        return 0;
    }
	
    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pAig);

    return 0;

usage:
    Abc_Print( -2, "usage: read_cnf [-i <cnf>] [-h]\n" );
    Abc_Print( -2, "\t        Read DIMACS CNF file and construct corresponding AIG network\n" );
    Abc_Print( -2, "\t-i <cnf> : CNF file input\n");
    Abc_Print( -2, "\t-h       : print the command usage\n");
    return 0;
}

// ABC command: Generate a sampling circuit 
int SampleCkt_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 0;
    int nPO = 0;
    bool fCorrelation = false;
    bool fVerbose = false;
    Abc_Ntk_t * pNtk, * pAig;
    SampleCircuit sc;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "iocvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nPI = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPI <= 0 )
                goto usage;
            break;
        case 'o':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-o\" should be followed by an integer.\n" );
                goto usage;
            }
            nPO = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPO <= 0 )
                goto usage;
            break;
        case 'c':
            fCorrelation = true;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'h':
            goto usage;
        default:
            goto usage;
        }
    }

    // read the current network if -c flag is activate
    if (fCorrelation)
    {
        pNtk = Abc_FrameReadNtk(pAbc);
        assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
        nPO = Abc_NtkPiNum(pNtk);
    }

    // check if the number of PO/PI is available
    assert(nPI > 0 && nPO > 0);
    assert(nPI < nPO);

    // generate sampling circuit
    sc.setIOnum(nPI, nPO);
    if (fCorrelation)
        pAig = sc.genCircuit(pNtk);
    else
        pAig = sc.genCircuit();
    if (fVerbose)
    {
        Abc_Print( 2, "Generate sampling circuit w/ nPI = %d, nPO = %d\n", nPI, nPO );
        cout << sc;
    }

    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pAig);

    return 0;

usage:
    Abc_Print( -2, "usage: sampleCkt [-i <num>] [-o <num> / -c] [-vh]\n" );
    Abc_Print( -2, "\t        Generate a sampling circuit with given PI and PO number\n" );
    Abc_Print( -2, "\t-i <num>  : set the number of PI\n");
    Abc_Print( -2, "\t-o <num>  : set the number of PO\n");
    Abc_Print( -2, "\t-c        : consider supports information of the current network\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: Generate a sampling circuit 
int SampleGen_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI, nPO;
    int nSample = 0;
    bool fRedirect = false;
    bool fVerbose = false;
    int * pPattern;
    char* filename;
    fstream file;
    vector<int> vNum;
    Abc_Ntk_t * pNtk;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "srvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 's':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nSample = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nSample <= 0 )
                goto usage;
            break;
        case 'r':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-r\" should be followed by an string.\n" );
                goto usage;
            }
            fRedirect = true;
            filename = argv[globalUtilOptind];
            globalUtilOptind++;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'h':
            goto usage;
        default:
            goto usage;
        }
    }

    // read the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
    nPI = Abc_NtkPiNum(pNtk);
    nPO = Abc_NtkPoNum(pNtk);
    assert(nSample <= pow(2,nPI));

    // generate samples
    if (fRedirect)    
    {
        file.open(filename, ios::out|ios::trunc);
        assert(file.is_open());
    }
    for (int i = 0; i < int(pow(2,nPI)); i++)
        vNum.push_back(i);
    random_shuffle(vNum.begin(), vNum.end());
    for (int i = 0; i < nSample; i++)
    {
        bitset<20> bs(vNum[i]);
        pPattern = new int[nPI];
        for (int j = 0; j < nPI; j++)
            pPattern[j] = bs[j];
        int * pValues = Abc_NtkVerifySimulatePattern( pNtk, pPattern );
        if (fRedirect)
        {
            for (int j = 0; j < nPO; j++)
            {
                if (pValues[j] == 0)
                    file << "-";
                file << j+1 << " ";
            }
            file << "\n";
        }
        else
        {
            for (int j = 0; j < nPO; j++)
            {
                if (pValues[j] == 0)
                    cout << "-";
                cout << j+1 << " ";
            }
            cout << "\n";
        }
    }
    
    file.close();
    return 0;

usage:
    Abc_Print( -2, "usage: sampleGen [-s <num>] [-r <file>] [-vh]\n" );
    Abc_Print( -2, "\t        Using created sampling circuit to generate samples\n" );
    Abc_Print( -2, "\t-s <num>  : set the number of samples\n");
    Abc_Print( -2, "\t-r <file> : redirect the result to the given file\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: Generate a sampling circuit and connect it to the current network
int SampleCnt_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 0;
    int nPO;
    bool fVerbose = false;
    bool fCorrelation = false;
    Abc_Ntk_t * pNtk;
    Abc_Ntk_t * pAig, * pAigNew; 
    SampleCircuit sc;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "icvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nPI = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPI <= 0 )
                goto usage;
            break;
        case 'c':
            fCorrelation = true;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'h':
        default:
            goto usage;
        }
    }

    // get the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
    if (!Abc_NtkIsStrash(pNtk))
        pNtk = Abc_NtkStrash(pNtk, 0, 1, 0);
    if (fVerbose)
    {
        Abc_Print( -2, "Network Information\n");
        Abc_Print( -2, "Pi Num: %d\n", Abc_NtkPiNum(pNtk));
        Abc_Print( -2, "Po Num: %d\n", Abc_NtkPoNum(pNtk));
    }
    nPO = Abc_NtkPiNum(pNtk);

    // check if the number of PO/PI is available
    assert(nPI > 0 && nPO > 0);
    assert(nPI < nPO);

    // generate sample circuit & connect
    sc.setIOnum(nPI, nPO);
    if (fCorrelation)
        pAig = sc.genCircuit(pNtk);
    else
        pAig = sc.genCircuit();
    pAigNew = sc.connect(pAig, pNtk);
    if (fVerbose)
    {
        Abc_Print( 2, "Generate sample circuit w/ nPI = %d, nPO = %d\n", nPI, nPO );
        cout << sc;
        /*Abc_Print( 2, "Sampling Circuit:\n" );
        Abc_NtkPrintStrSupports( pAig, 0 );
        Abc_Print( 2, "Connect to current network\n");
        Abc_Print( 2, "Original Network:\n" );
        Abc_NtkPrintStrSupports( pNtk, 0 );
        Abc_Print( 2, "Sampling Network:\n" );
        Abc_NtkPrintStrSupports( pAigNew, 0 );*/
    }

    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pAigNew);

    return 0;

usage:
    Abc_Print( -2, "usage: sampleCnt [-i <num>] [-cvh]\n" );
    Abc_Print( -2, "\t        Generate a sampling circuit with given PI number and connect it to the current network\n" );
    Abc_Print( -2, "\t-i <num> : set the number of PI\n");
    Abc_Print( -2, "\t-c       : consider supports information of the current network\n");
    Abc_Print( -2, "\t-v       : verbose\n");
    Abc_Print( -2, "\t-h       : print the command usage\n");
    return 0;
}

// ABC command: Draw samples from the witness of the given formula
int SampleWit_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nSample;
    int nPI, nPO, nGen = 0, nAttempt = 0, nCircuit = 0, nSize;
	bool fVerbose = false, fRedirect = false;
   	char * filename;
   	char circuitname[10];
	fstream file;
    Abc_Ntk_t * pNtk, * pCkt, * pNtkRes;
	int * pIn1, * pIn2, * pOut;
    vector<int*> vMinterm;
	vector<int*> vResult;
    SampleCircuit sc;
    clock_t t1, t2;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "isrvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nPI = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPI <= 0 )
                goto usage;
            break;
        case 's':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-s\" should be followed by an integer.\n" );
                goto usage;
            }
            nSample = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nSample <= 0 )
                goto usage;
            break;
        case 'r':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-r\" should be followed by an string.\n" );
                goto usage;
            }
            fRedirect = true;
            filename = argv[globalUtilOptind];
            globalUtilOptind++;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'h':
        default:
            goto usage;
        }
    }

    // get the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
    if (!Abc_NtkIsStrash(pNtk))
        pNtk = Abc_NtkStrash(pNtk, 0, 1, 0);
    nPO = Abc_NtkPiNum(pNtk);
	
	// generate samples
    t1 = clock();
	sc.setIOnum(nPI, nPO);
    pIn1 = new int[nPI];
	while (nGen < nSample)
	{
        // symbolic sampling
        pCkt = sc.genCircuit();
        pNtkRes = sc.connect(pCkt, pNtk);        

        // count minterm
        for (int i = 0; i < pow(2,nPI); i++)
        {
            for (int j = 0; j < nPI; j++)
                pIn1[j] = (i & (1 << j)) ? 1 : 0;
            pOut = Abc_NtkVerifySimulatePattern(pNtkRes, pIn1);
            if (pOut[0] == 1)
            {
                pIn2 = Abc_NtkVerifySimulatePattern(pCkt, pIn1);
                vMinterm.push_back(pIn2);
            }
            delete pOut;
        }
       
        // draw sample 
	    nSize = vMinterm.size();
	    if (nSize > 0)
	    {
            random_shuffle(vMinterm.begin(), vMinterm.end());
	    	vResult.push_back(vMinterm[0]);
            vMinterm.clear();
            nGen++;
            /*// dump the circuit
            sprintf(circuitname, "wit%d.aig", nCircuit);
            assert(vNtk[t] != NULL);
            Io_WriteAiger( pCkt, circuitname, 0, 0, 0 ); //io name,compact,conanical
		    nCircuit++;*/
	    }
	    Abc_NtkDelete(pCkt);
	    Abc_NtkDelete(pNtkRes);
        nAttempt++;
	}
    delete pIn1;
    t2 = clock();
    if (fVerbose)
        Abc_Print( 2, "Generate %d samples with p=%f%% in %f secs.\n", nSample, ((float)nSample/nAttempt)*100, (double)(t2-t1)/CLOCKS_PER_SEC );

	// dump the result
	assert(vResult.size() == nSample);
	if (fRedirect)
	{
        file.open(filename, ios::out|ios::trunc);
        assert(file.is_open());
        file << (float)nSample/nAttempt << "\n";
		for (int i = 0; i < nSample; i++)
		{
			for (int j = 0; j < nPO; j++)
			{
				if (vResult[i][j] == 0)
					file << "-";
				file << j+1 << " ";
			}
			file << "\n";
		}
        file.close();
	}
	/*else
	{
		for (int i = 0; i < nSample; i++)
		{
			for (int j = 0; j < nPO; j++)
			{
				if (vResult[i][j] == 0)
					cout << "-";
				cout << j+1 << " ";
			}
			cout << "\n";
		}
	}*/

    vResult.clear();
	return 0;

usage:
    Abc_Print( -2, "usage: sampleWit [-s <num>] [-vh]\n" );
    Abc_Print( -2, "\t        Generate witnesses of the current network and dump the sampling circuits named by \"wit<num>.aig\"\n" );
    Abc_Print( -2, "\t-i <num>  : set the number of PI\n");
    Abc_Print( -2, "\t-s <num>  : set the number of samples\n");
    Abc_Print( -2, "\t-r <file> : redirect the result to the given file\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: Pearson Chi Square Test for sampling circuit
int SampleChiTest_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int i, j, pos;
    int nPI = 0, nPO = 0;
    int nGen = 0, nSample = 0;
    double nExpect;
    bool fCorrelation = false;
    bool fRedirect = false;
    bool fVerbose = false;
    int * pInput;
    string strPattern;
    fstream file;
    char* filename;
    Abc_Obj_t * pPi;
    Abc_Ntk_t * pNtk, * pSampleCkt;
    SampleCircuit sc;
    unordered_map<string, int> PImap;
    map<string, unsigned char> mObserve;
    map<string, unsigned char>::iterator mIter;
    double valCrit, valP;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "ioecrvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nPI = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPI <= 0 )
                goto usage;
            break;
        case 'o':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-o\" should be followed by an integer.\n" );
                goto usage;
            }
            nPO = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPO <= 0 )
                goto usage;
            break;
        case 'e':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-e\" should be followed by an integer.\n" );
                goto usage;
            }
            nSample = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nSample <= 0 )
                goto usage;
            break;
        case 'r':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-r\" should be followed by an string.\n" );
                goto usage;
            }
            fRedirect = true;
            filename = argv[globalUtilOptind];
            globalUtilOptind++;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'c':
            fCorrelation = true;
            break;
        case 'h':
        default:
            goto usage;
        }
    }

    // get the current network
    if (fCorrelation)
    {
        pNtk = Abc_FrameReadNtk(pAbc);
        assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
        assert(Abc_NtkHasAig(pNtk) && Abc_NtkIsStrash(pNtk));
        nPO = Abc_NtkPiNum(pNtk);
    
        // construct mapping between ID and name for current network
        Abc_NtkForEachPi( pNtk, pPi, i )
        {
            PImap[Abc_ObjName(pPi)] = i;
            if (fVerbose)
                cout << "Map " << Abc_ObjName(pPi) << " to " << PImap[Abc_ObjName(pPi)] << "\n";
        }
        if (fVerbose)
            cout << "\n";
    }

    // check if the number of PO/PI is available
    assert(nPI > 1 && nPI < nPO);    

    // set to default value
    if (nSample == 0)
        nSample = (int)pow(2,nPO);
    if (fVerbose)
        Abc_Print( -2, "\nConduct Pearson's chi-squared test w/ nPI = %d, nPO = %d, and nSample = %d.\n", nPI, nPO, nSample );
    pInput = new int[nPI];
    while (nGen < nSample)
    {
        if (fVerbose)
            Abc_Print( -2, "\nSampling...(%d/%d)\n", nGen, nSample );
        
        // generate sampling circuit
        sc.setIOnum(nPI, nPO);
        if (fCorrelation)
            pSampleCkt = sc.genCircuit(pNtk);
        else
            pSampleCkt = sc.genCircuit();
        
        // count patterns
        strPattern.clear();
        for (i = 0; i < pow(2,nPI); i++)
        {
            for (j = 0; j < nPI; j++)
                pInput[j] = (i & (1 << j)) ? 1 : 0;
            int * pValues = Abc_NtkVerifySimulatePattern(pSampleCkt, pInput);
            for (j = 0; j < nPO; j++)
                strPattern.push_back((char)(pValues[j]+48));
            if ((mIter = mObserve.find(strPattern)) != mObserve.end())
                mIter->second++;
            else
                mObserve[strPattern] = 1;
        }

        nGen += pow(2,nPI);
    }

    // compute chi-square value
    valCrit = 0;
    nExpect = nSample/pow(2,nPO);
    for (mIter = mObserve.begin(); mIter != mObserve.end(); mIter++)
    {
        double XSqr = mIter->second - nExpect;
        valCrit += ((XSqr * XSqr) / nExpect);
    }
    valCrit += nExpect*(pow(2,nPO)-mObserve.size());    

    // compute p-value
    valP = chisqr(pow(2,nPO)-1, valCrit);

    // print
    if (fRedirect)
    {
        file.open(filename, ios::out|ios::trunc);
        assert(file.is_open());
        file << "DoF = " << pow(2,nPO)-1 << "\n";
        file << "Chi-value = " << valCrit << "\n";
        file << "P-value = " << valP << "\n";
        file.close();
    }
    else
    {
        cout << "DoF = " << pow(2,nPO)-1 << "\n";
        cout << "Chi-value = " << valCrit << "\n";
        cout << "P-value = " << valP << "\n";
    }

    return 0;

usage:
    Abc_Print( -2, "usage: sampleChiTest [-i <num>] [-e <num>] [-o <num>/-c] [-r <file>] [-vh]\n" );
    Abc_Print( -2, "\t        Apply Pearson Chi-Square test on sampling circuit\n" );
    Abc_Print( -2, "\t-i <num>  : set the number of PI\n");
    Abc_Print( -2, "\t-o <num>  : set the number of PO\n");
    Abc_Print( -2, "\t-e <num>  : set the number of samples\n");
    Abc_Print( -2, "\t-r <file> : redirect the result to the given file\n");
    Abc_Print( -2, "\t-c       : consider supports information of the current network\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: Bias sampling test on single stuck-at fault circuit
int SampleStuckTest_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 10, nPO;
    int nGen, nSample = 1024, nTest = 1;
    int RetValue;
    int nFixed, nCount1, nCount2, nDup1, nDup2;
    int * pValues, * pModel, * pIn, * pIn2;
    bool fDump = false;
    bool fRedirect = false;
    bool fVerbose = false;
    char * filename;
    Abc_Ntk_t * pNtk, * pNtkRes;
    Abc_Ntk_t * pMiter;
    Abc_Ntk_t * pCkt, * pCnt; 
    Prove_Params_t Params, * pParams = &Params;
    SampleCircuit sc;
    set<string> setPattern;
    vector<int> vFixed, vNFixed;
    vector<int> vCount1, vCount2;
    vector<int> vDup1, vDup2;
    clock_t t1, t2;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "isndrvh" ) ) != EOF )
    {
        switch ( c )
        {
        case 'i':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-i\" should be followed by an integer.\n" );
                goto usage;
            }
            nPI = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nPI <= 0 )
                goto usage;
            break;
        case 's':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-s\" should be followed by an integer.\n" );
                goto usage;
            }
            nSample = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nSample <= 0 )
                goto usage;
            break;
        case 'n':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-n\" should be followed by an integer.\n" );
                goto usage;
            }
            nTest = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( nTest <= 0 )
                goto usage;
            break;
        case 'r':
            if ( globalUtilOptind >= argc )
            {
                Abc_Print( -1, "Command line switch \"-r\" should be followed by an string.\n" );
                goto usage;
            }
            fRedirect = true;
            filename = argv[globalUtilOptind];
            globalUtilOptind++;
            break;
        case 'd':
            fDump = true;
            break;
        case 'v':
            fVerbose = true;
            break;
        case 'h':
        default:
            goto usage;
        }
    }

    // get the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    assert(pNtk != NULL && Abc_NtkIsComb(pNtk));
    if (!Abc_NtkIsStrash(pNtk))
        pNtk = Abc_NtkStrash(pNtk, 0, 1, 0);
    nPO = Abc_NtkPiNum(pNtk);

    // check if the number of PO/PI is available
    assert(nPI > 0 && nPO > 0);
    assert(nPI < nPO);

    // conduct experiment
    pIn = new int[nPI];
    pIn2 = new int[nPO];
    t1 = clock();
    for (int i = 0; i < nTest; i++)
    {
        if (fVerbose)
            Abc_Print( 2, "Test...(%d/%d)\n", i+1, nTest );

        // generate single stuck-at fault
        pNtkRes = Ntk_StuckGen(pNtk);

        // build the miter
        pMiter = Abc_NtkMiter( pNtk, pNtkRes, 1, 0, 0, 0 );
        assert(pMiter && Abc_NtkMiterIsConstant(pMiter) != 1);

        // get the counter-example by SAT solver
        Prove_ParamsSetDefault( pParams );
        pParams->nItersMax = 5;
        RetValue = Abc_NtkIvyProve( &pMiter, pParams );
        assert(RetValue == 0); // 0:NEQ, 1:EQ, -1:Unknown
        pModel = pMiter->pModel;
        assert(pModel);

        // determine variables that make miter alter after flipped
        vFixed.clear();
        for (int j = 0; j < nPO; j++)
        {
            pModel[j] ^= 1;
            pValues = Abc_NtkVerifySimulatePattern( pMiter, pModel );
            if (pValues[0] == 0)
                vFixed.push_back(j);
            pModel[j] ^= 1;
            delete [] pValues;
        }

        // original sim
        nGen = 0;
        nCount1 = 0;
        nDup1 = 0;
        setPattern.clear();
        while (nGen < nSample)
        {
            sc.setIOnum(nPI, nPO);
            pCkt = sc.genCircuit(pNtk);
            pCnt = sc.connect(pCkt, pMiter);
            for (int j = 0; j < pow(2,nPI); j++)
            {
                // generate pattern
                for (int k = 0; k < nPI; k++)
                    pIn[k] = (j & (1 << k)) ? 1 : 0;
                
                // sim
                pValues = Abc_NtkVerifySimulatePattern(pCnt, pIn);
                
                // count satisfying patterns
                if (pValues[0] == 1)
                    nCount1++;

                // count duplicated patterns
                delete [] pValues;
                pValues = Abc_NtkVerifySimulatePattern(pCkt, pIn);
                string s;
                for (int k = 0; k < nPO; k++)
                    s.push_back((char)(pValues[k]+48));
                if (setPattern.count(s))
                    nDup1++;
                else
                    setPattern.insert(s);

                delete [] pValues; 
            }

            nGen += (int)pow(2,nPI);
        }
        
        // biased sim
        random_shuffle(vFixed.begin(), vFixed.end());
        nFixed = min((int)vFixed.size(), (int)(nPO/10));
        nGen = 0;
        nCount2 = 0;
        nDup2 = 0;
        setPattern.clear();
        while (nGen < nSample)
        {
            sc.setIOnum(nPI, nPO-nFixed);
            pCkt = sc.genCircuit();
            for (int j = 0; j < pow(2,nPI); j++)
            {
                // generate pattern
                for (int k = 0; k < nPI; k++)
                    pIn[k] = (j & (1 << k)) ? 1 : 0;
                pValues = Abc_NtkVerifySimulatePattern(pCkt, pIn);
                for (int k = 0; k < nPO; k++)
                    pIn2[k] = pModel[k];
                int m = 0, n = 0;
                for (int k = 0; k < nPO; k++)
                {
                    if (k == vFixed[m])
                        m++; 
                    else
                        pIn2[k] = pValues[n++];
                }
                
                // sim
                delete [] pValues;
                pValues = Abc_NtkVerifySimulatePattern(pMiter, pIn2);

                // count satisfying patterns
                if (pValues[0] == 1)
                    nCount2++;
   
                // count duplicated patterns
                string s;
                for (int k = 0; k < nPO; k++)
                    s.push_back((char)(pIn2[k]+48));
                if (setPattern.count(s))
                    nDup2++;
                else
                    setPattern.insert(s);
                
                delete [] pValues; 
            }

            nGen += (int)pow(2,nPI);
        }

        // save result
        vNFixed.push_back(nFixed);
        vCount1.push_back(nCount1);
        vCount2.push_back(nCount2);
        vDup1.push_back(nDup1);
        vDup2.push_back(nDup2);

        ABC_FREE( pMiter->pModel );
        Abc_NtkDelete( pMiter );
        Abc_NtkDelete( pNtkRes );
    }
    t2 = clock();
    if (fVerbose)
    {
        Abc_Print( 2, "# of tests: %d\n", nTest);
        Abc_Print( 2, "# of samples per test: %d\n", nSample);
        Abc_Print( 2, "Avg. time in secs: %.1f\n", ((double)(t2-t1)/CLOCKS_PER_SEC)/nTest);
        Abc_Print( 2, "Avg.%% Hit in origin: %.1f\n", 100*accumulate(vCount1.begin(), vCount1.end(), 0.0)/nTest/nSample);
        Abc_Print( 2, "Avg.%% Dup in origin: %f\n", 100*accumulate(vDup1.begin(), vDup1.end(), 0.0)/nTest/nSample);
        Abc_Print( 2, "Avg. Fixed variables: %.1f\n", accumulate(vNFixed.begin(), vNFixed.end(), 0.0)/nTest);
        Abc_Print( 2, "Avg.%% Hit in biased: %.1f\n", 100*accumulate(vCount2.begin(), vCount2.end(), 0.0)/nTest/nSample);
        Abc_Print( 2, "Avg.%% Dup in biased: %f\n", 100*accumulate(vDup2.begin(), vDup2.end(), 0.0)/nTest/nSample);
    }

    delete [] pIn;
    delete [] pIn2;

    return 0;

usage:
    Abc_Print( -2, "usage: sampleStuckTest [-i <num>] [-s <num>] [-n <num>] [-r <file>] [-dvh]\n" );
    Abc_Print( -2, "\t        Generate single stuck-at fault on current network and test bias sampling\n" );
    Abc_Print( -2, "\t-i <num>  : sets the number of fanins[default:10]\n");
    Abc_Print( -2, "\t-s <num>  : sets the number of samples[default:1024]\n");
    Abc_Print( -2, "\t-n <num>  : sets the number of times[default:1]\n");
    Abc_Print( -2, "\t-d        : toggles whether to dump the aig circuits[default: False]\n");
    Abc_Print( -2, "\t-r <file> : redirects the result to the file\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// called during ABC startup
void init(Abc_Frame_t* pAbc)
{
    Cmd_CommandAdd( pAbc, "Sample", "info", Info_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "read_cnf", ReadCNF_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleCkt", SampleCkt_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleGen", SampleGen_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleCnt", SampleCnt_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleWit", SampleWit_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleChiTest", SampleChiTest_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleStuckTest", SampleStuckTest_Command, 0);
}

// called during ABC termination
void destroy(Abc_Frame_t* pAbc)
{
}

// this object should not be modified after the call to Abc_FrameAddInitializer
Abc_FrameInitializer_t frame_initializer = { init, destroy };

// register the initializer a constructor of a global object
// called before main (and ABC startup)
struct registrar
{
    registrar() 
    {
        Abc_FrameAddInitializer(&frame_initializer);
    }
} sample_registrar;

} // unnamed namespace
