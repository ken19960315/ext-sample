# ext-sample
A sampling circuit generator package embedded in [abc](https://github.com/berkeley-abc/abc).
## Features
- Generate sampling circuits to perform uniform sampling on given Boolean space.
- Capabale to draws samples under certain constraints presented by DIMACS CNF or BLIF format.
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
  - **sampleWit**  
    - Generate witnesses of the current network and dump the sampling circuits named by "wit\<num>.aig".

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
      
# Example
## Generate sampling circuits under constraints    
  Suppose there are 10 variables from v1 to v10, and we want to draw samples that meet the constraint (v1+v5).  
  Then we can use either DIMACS CNF or BLIF to represent the constraints like the following:
  - DIMACS CNF
    ```
    c ind 1 2 3 4 5 6 7 8 9 10 0
    p cnf 10 1
    1 5 0
    ```
  - BLIF
    ```
    .model test
    .inputs v1 v2 v3 v4 v5 v6 v7 v8 v9 v10
    .outputs out
    .names v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 out
    1--------- 1
    ----1----- 1
    .end
    ```
  After, we use the implemented commands in abc to perform sampling
  ```
  ./abc
  # abc commandline
  > read_cnf -i <cnf> or read <blif> # read constraints as a circuit
  > sampleWit -s 100                 # generate sampling circuits until contain at least 100 samples 
  ```
  If -v option is added after command "sampleWit", the program will print out the value of "loThresh", and each generated circuit contains "loThresh" number of samples. 
