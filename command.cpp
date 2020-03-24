#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "proof/fraig/fraig.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <unordered_map>
#include <bitset>
#include <algorithm>

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
        return NULL;
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
    Abc_Print( -2, "\t        Generate a sample circuit with given PI and PO number\n" );
    Abc_Print( -2, "\t-i <num>  : set the number of PI\n");
    Abc_Print( -2, "\t-o <num>  : set the number of PO\n");
    Abc_Print( -2, "\t-c        : consider supports correlation\n");
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
    file.open(filename, ios::out|ios::trunc);
    assert(file.is_open());
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
    pAigNew = sc.connect(pNtk);
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
    Abc_Print( -2, "\t        Generate a sample circuit with given PI number and connect it to current AIG network\n" );
    Abc_Print( -2, "\t-i <num> : set the number of PI\n");
    Abc_Print( -2, "\t-c       : consider supports correlation\n");
    Abc_Print( -2, "\t-v       : verbose\n");
    Abc_Print( -2, "\t-h       : print the command usage\n");
    return 0;
}

// ABC command: Draw samples from the witness of the given formula
int SampleWit_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI, nSample = 0, nGen = 0, nCircuit = 0, nSize;
	int hashBits = 0, loThresh, hiThresh;
	double kappa = 0.638;
    double pivot;
	bool fVerbose = false, fRedirect = false;
   	char * filename;
   	char circuitname[10];
	fstream file;
    Abc_Ntk_t * pNtk, * pCkt, * pNtkRes;
    Abc_Obj_t * pObj;
	vector<int*> vMinterm;
	vector<int*> vSample;
	vector<int*> vResult;
    SampleCircuit sc;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "srvh" ) ) != EOF )
    {
        switch ( c )
        {
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
    nPI = Abc_NtkPiNum(pNtk);
	
	// estimate parameters
	pivot = ceil(4.03 * (1 + 1/kappa) * (1 + 1/kappa));
	hiThresh = (int)ceil(1 + sqrt(2)*(1+kappa)*pivot);	
	loThresh = (int)floor(pivot / (sqrt(2)*(1+kappa)));	
    for (int i = nPI-1; i > 0; i--)
	{
		sc.setIOnum(i, nPI);
        pCkt = sc.genCircuit();
		pNtkRes = sc.connect(pNtk);
    	vMinterm = Ntk_Minterm(pNtkRes, 61);
		Abc_NtkDelete(pCkt);
		Abc_NtkDelete(pNtkRes);
		nSize = vMinterm.size();
		if (nSize >= 1 && nSize <= 60)
		{
			hashBits = (int)round(log2(nSize)+(nPI-i)+log2(1.8)-log2(pivot));
			break;
		}
	}
	if (hashBits <= 2)
	{
		Abc_Print( 2, "Too few solutions\n");
		return 0;
	}

	if (fVerbose)
	{
		Abc_Print( 2, "hashBits: %d\n", hashBits );
		Abc_Print( 2, "loThresh: %d, hiThresh: %d\n", loThresh, hiThresh );
	}
	// generate samples
	while (nGen < nSample)
	{
		if ( (pNtkRes=genSample(pNtk, hashBits, loThresh, hiThresh, vSample)) != NULL )
		{		
			assert(vSample.size() == loThresh);
			for (int i = 0; i < vSample.size() && nGen < nSample; i++, nGen++)
				vResult.push_back(vSample[i]);
            // dump the circuit
            sprintf(circuitname, "wit%d.aig", nCircuit);
            Io_WriteAiger( pNtkRes, circuitname, 0, 0, 0 ); //io name,compact,conanical
		    nCircuit++;
        }
	}

	// dump the result
	assert(vResult.size() == nSample);
	if (fRedirect)
	{
        file.open(filename, ios::out|ios::trunc);
        assert(file.is_open());
		for (int i = 0; i < nSample; i++)
		{
			for (int j = 0; j < nPI; j++)
			{
				if (vResult[i][j] == 0)
					file << "-";
				file << j+1 << " ";
			}
			file << "\n";
		}
	}
	/*else
	{
		for (int i = 0; i < nSample; i++)
		{
			for (int j = 0; j < nPI; j++)
			{
				if (vResult[i][j] == 0)
					cout << "-";
				cout << j+1 << " ";
			}
			cout << "\n";
		}
	}*/

	return 0;

usage:
    Abc_Print( -2, "usage: sampleWit [-s <num>] [-vh]\n" );
    Abc_Print( -2, "\t        Generate witnesses of the current network\n" );
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
    int nPI = 0;
    int nPO = 0;
    int times = 0;
    bool fCorrelation = false;
    bool fRedirect = false;
    bool fVerbose = false;
    fstream file;
    char* filename;
    Abc_Obj_t * pPi;
    Abc_Ntk_t * pNtk;
    SampleCircuit sc;
    unordered_map<string, int> PImap;
    int expect = 5;
    vector<short> vObserve;
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
            expect = atoi(argv[globalUtilOptind]);
            globalUtilOptind++;
            if ( expect <= 0 )
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
        if (pNtk == NULL)
        {
            Abc_Print( -1, "Cannot get the current network.\n");
            return 0;
        }
        // check whether the network is combinational
        if (!Abc_NtkIsComb(pNtk))
        {
            Abc_Print( -1, "Network is not combinational.\n" );
            return 0;
        }
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
    if (nPO == 1)
    {
        Abc_Print( -1, "Circuit has only one PI.\n");
        return 0;
    }
    if (nPI <= 0)
    {
        Abc_Print( -1, "The number of PI should be larger than 0.\n" );
        return 0;
    }
    if (nPI >= nPO)
    {
        Abc_Print( -1, "The number of PI should be less than the number of PO.\n" );
        return 0;
    }
    if (nPO > 20)
    {
        Abc_Print( -1, "The number of PO should be less than 20 due to memory concern.\n" );
        return 0;
    }

    // generate sampling circuit
    times = expect*pow(2,nPO) / pow(2,nPI);
    if (fVerbose)
        Abc_Print( -2, "\nGenerate sample circuit w/ nPI = %d, nPO = %d, and %d times.\n", nPI, nPO, times );
    sc.setRndSeed((unsigned)time(NULL));
    vObserve.resize(pow(2,nPO));
    for (i = 0; i < pow(2,nPO); i++)
        vObserve[i] = 0;
    for (i = 0; i < times; i++)
    {
        if (fVerbose)
            Abc_Print( -2, "\nSampling...(%d/%d)\n", i+1, times );
            
        sc.setIOnum(nPI, nPO);
        sc.genCircuit();
        if (fVerbose)
            cout << sc;

        // count patterns
        vector< vector<int> > XOR = sc.getXOR();
        for (j = 0; j < pow(2,nPI); j++)
        {
            string in = bitset<20>(j).to_string();
            string out(nPO, '-');
            for (int m = 0; m < nPO; m++)
            {
                int output = XOR[m][0];
                for (int n = 1; n < XOR[m].size(); n++)
                {
                    if (in[in.size()-XOR[m][n]-1] == '1')
                        output ^= 1;
                }
                pos = fCorrelation? PImap[Abc_ObjName(Abc_NtkPi(pNtk, m))] : m;
                out[pos] = (char)(output+48);
            }
            vObserve[bitset<20>(out).to_ulong()]++;
        }
    }

    // compute chi-square value
    valCrit = 0;
    for (i = 0; i < pow(2,nPO); i++)
    {
        double XSqr = vObserve[i] - expect;
        valCrit += ((XSqr * XSqr) / expect);
    }
    
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
    Abc_Print( -2, "usage: sampleChiTest [-i <num>] [-o <num>] [-r <file>] [-cvh]\n" );
    Abc_Print( -2, "\t        Apply Pearson Chi-Square test on sampling circuit\n" );
    Abc_Print( -2, "\t-i <num>  : set the number of PI\n");
    Abc_Print( -2, "\t-o <num>  : set the number of PO\n");
    Abc_Print( -2, "\t-e <num>  : set the expected number of each pattern[default:5]\n");
    Abc_Print( -2, "\t-r <file> : redirect the result to the given file\n");
    Abc_Print( -2, "\t-c        : consider network supports correlation\n");
    Abc_Print( -2, "\t-v        : verbose\n");
    Abc_Print( -2, "\t-h        : print the command usage\n");
    return 0;
}

// ABC command: Bias sampling test on single stuck-at fault circuit
int SampleStuckTest_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 10, nPO;
    int nSample = 1024, nTest = 1;
    int RetValue;
    int * pModel;
    bool fDump = false;
    bool fRedirect = false;
    bool fVerbose = false;
    char * filename;
    Abc_Ntk_t * pNtk, * pNtkRes;
    Abc_Ntk_t * pMiter;
    Abc_Ntk_t * pAig, * pAigNew; 
    Prove_Params_t Params, * pParams = &Params;
    SampleCircuit sc;
    vector<int> vFixed;

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

    pNtkRes = Ntk_StuckGen(pNtk);
    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pNtkRes);
    
    /*for (int i = 0; i < nTest; i++)
    {
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

        // determine variables that make miter remained 0 after flipped
        for (int j = 0; j < nPO; j++)
        {
            pModel[j] ^= 1;
            int * pValues = Abc_NtkVerifySimulatePattern( pMiter, pModel );
            if (pValues[0] == 1)
                vFixed.push_back(j);
            pModel[j] ^= 1;
        }

        // generate sample circuit & connect
        
        // simulate

        ABC_FREE( pMiter->pModel );
        Abc_NtkDelete( pMiter );
    }*/
    /*sc.setIOnum(nPI, nPO);
    pAig = sc.genCircuit(pNtk);
    pAigNew = sc.connect(pNtk);
    if (fVerbose)
    {
        Abc_Print( 2, "Generate sample circuit w/ nPI = %d, nPO = %d\n", nPI, nPO );
        cout << sc;
    }*/

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
