# ext-sample
A sampling circuit generator package embedded in [abc](https://github.com/berkeley-abc/abc).
## Features
- Generate sampling circuits to perform uniform sampling on given Boolean space.
- Connect generated sampling circuits to a circuit to cast the computation from original domain to sampling domain.

# Compile
1. Put the entire "ext-sample" directory to the "src" directory under abc.
2. type "make" command at the root directory of abc

# Documentation
- **Commands Usage**  
  \[-h] option shows the detailed usage of each command
  - **sampleCkt**  
    - Generate a sampling circuit with given PI and PO number.
  - **sampleCnt**  
    - Generate a sampling circuit with given PI number and connect it to the current network.  
    - \[-c] option make the circuit size smaller after connection in our experiments.

- **Interface**    
  - SamplingCircuit.h
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
