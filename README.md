# ext-sample
A sampling circuit generator package embedded in [abc](https://github.com/berkeley-abc/abc).

# Compile
1. Put the entire "ext-sample" directory to the "src" directory under abc.
2. type "make" command at the root directory of abc

# Documentation
1. command.cpp  
  Implement two commands for abc
  - sampleGen \[-i \<num>] \[-o \<num>]  
    Generate a sample circuit with specific PI and PO number in AIG format.
  - sampleCnt \[-i \<num>]  
    Generate a sample circuit with specific PI number, and connect it to the current network. Notice that the current network should be strash.

2. SampleCircuit.h/cpp
