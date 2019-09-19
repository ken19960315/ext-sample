# ext-sample
A sampling circuit generator package embedded in [abc](https://github.com/berkeley-abc/abc).

# Compile
1. Put the entire "ext-sample" directory to the "src" directory under abc.
2. type "make" command at the root directory of abc

# Documentation
- command.cpp  
  Implement two commands for abc, and there is some usage of SampleCircuit class.
  1. **sampleGen \[-i \<num>] \[-o \<num>]**  
    Generate a sample circuit with specific PI and PO number in AIG format.
  2. **sampleCnt \[-i \<num>]**  
    Generate a sample circuit with specific PI number, and connect it to the current network. Notice that the current network should be strash.

- SampleCircuit.h/cpp  
  C++ interface  
  1. **Constructor**  
    can initialize the number of PI/PO here
  2. **void setIOnum(int nPI, int nPO)**  
    set the number of PI/PO
  3. **void setRndSeed(int seed)**  
    set random seed
  4. **Abc_Ntk_t\* genCircuit()**  
    generate sampling circuit with given PI/PO number and return the network(also stored in "SampleCircuit::pAig" variable).
  5. **Abc_Ntk_t\* connect(Abc_Ntk_t\* pNtk)**  
    connect to current network and return the new network
