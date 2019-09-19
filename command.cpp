#include "base/main/main.h"
#include "base/main/mainInt.h"

#include <iostream>
#include <ctime>

#include "ext-sample/SampleCircuit.h"

namespace
{

// ABC command: Generate a sample circuit 
int SampleGen_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 0;
    int nPO = 0;
    Abc_Ntk_t *pAig;
    SampleCircuit sc;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "ioh" ) ) != EOF )
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
        case 'h':
            goto usage;
        default:
            goto usage;
        }
    }

    // check if the number of PO/PI is available
    if (nPI <= 0 || nPO <= 0)
    {
        Abc_Print( -1, "The number of PI/PO should be larger than 0.\n" );
        goto usage;
    }
    if (nPI >= nPO)
    {
        Abc_Print( -1, "The number of PI should be less than the number of PO.\n" );
        goto usage;
    }

    // generate sample circuit
    Abc_Print( -2, "Generate sample circuit w/ nPI = %d, nPO = %d\n", nPI, nPO );
    sc.setIOnum(nPI, nPO);
    sc.setRndSeed((unsigned)time(NULL));
    pAig = sc.genCircuit();

    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pAig);

    return 0;

usage:
    Abc_Print( -2, "usage: sampleGen [-ioh]\n" );
    Abc_Print( -2, "\t        Generate a sample circuit with given PI and PO number\n" );
    Abc_Print( -2, "\t-i <num> : sets the number of PI\n");
    Abc_Print( -2, "\t-o <num> : sets the number of PO\n");
    Abc_Print( -2, "\t-h       : print the command usage\n");
    return 0;
}

// ABC command: Generate a sample circuit and connect to current Aig
int SampleCnt_Command( Abc_Frame_t * pAbc, int argc, char ** argv )
{
    char c;
    int nPI = 0;
    int nPO;
    bool f_verbose = false;
    bool f_auto = false;
    Abc_Ntk_t * pNtk;
    Abc_Ntk_t * pAigNew; 
    SampleCircuit sc;

    // parse arguments
    Extra_UtilGetoptReset();
    while ( ( c = Extra_UtilGetopt( argc, argv, "iavh" ) ) != EOF )
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
        case 'v':
            f_verbose = true;
            break;
        case 'a':
            f_auto = true;
            break;
        case 'h':
        default:
            goto usage;
        }
    }

    // get the current network
    pNtk = Abc_FrameReadNtk(pAbc);
    if (pNtk == NULL)
    {
        Abc_Print( -1, "Cannot get the current network.\n");
        goto usage;
    }
    assert(Abc_NtkHasAig(pNtk) && Abc_NtkIsStrash(pNtk));
    if (f_verbose)
    {
        Abc_Print( -2, "Network Information\n");
        Abc_Print( -2, "Pi Num: %d\n", Abc_NtkPiNum(pNtk));
        Abc_Print( -2, "Ci Num: %d\n", Abc_NtkCiNum(pNtk));
        Abc_Print( -2, "Box Num: %d\n", Abc_NtkBoxNum(pNtk));
        Abc_Print( -2, "Latch Num: %d\n", Abc_NtkLatchNum(pNtk));
    }
    nPO = Abc_NtkPiNum(pNtk);
    if (f_auto) 
        nPI = (nPO+1)/2;

    // check if the number of PO/PI is available
    if (nPI <= 0)
    {
        Abc_Print( -1, "The number of PI should be larger than 0.\n" );
        goto usage;
    }
    if (nPI >= nPO)
    {
        Abc_Print( -1, "The number of PI should be less than the number of PO.\n" );
        goto usage;
    }

    // generate sample circuit & connect
    if (f_verbose)
        Abc_Print( -2, "Generate sample circuit w/ nPI = %d, nPO = %d\n", nPI, nPO );
    sc.setIOnum(nPI, nPO);
    sc.setRndSeed((unsigned)time(NULL));
    sc.genCircuit();
    if (f_verbose)
        Abc_Print( -2, "Connect to current network\n");
    pAigNew = sc.connect(pNtk);

    // replace the current network
    Abc_FrameReplaceCurrentNetwork(pAbc, pAigNew);

    return 0;

usage:
    Abc_Print( -2, "usage: sampleCnt [-iavh]\n" );
    Abc_Print( -2, "\t        Generate a sample circuit with given PI number and connect it to current AIG network\n" );
    Abc_Print( -2, "\t-i <num> : sets the number of PI\n");
    Abc_Print( -2, "\t-a       : sets the number of PI equal to the half number of PO of current network\n");
    Abc_Print( -2, "\t-v       : verbose\n");
    Abc_Print( -2, "\t-h       : print the command usage\n");
    return 0;
}

// called during ABC startup
void init(Abc_Frame_t* pAbc)
{
    Cmd_CommandAdd( pAbc, "Sample", "sampleGen", SampleGen_Command, 0);
    Cmd_CommandAdd( pAbc, "Sample", "sampleCnt", SampleCnt_Command, 0);
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
