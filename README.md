# Parallel implementation of AD3


## Objective
We had a double objective.
**First**: Make an efficient implementation of AD3 and study the possibilites of parallelisation
**Second**: Apply AD3 to a new domain of application: Combinatorial Auctions.

## Implementation
Efficient implementation and parallelisation. We made a new implementation trying to exploit as much as possible the CPU resources. we could summarise our approach as follows:
- A reorganisation of the data memory layout, to improve the use of the memory system
- A reogranisation of the execution flow, so we can take profit of vectorisation (Which is a capability of modern processors. They are able to perform some operations in one single instruction given some circusntamces)
- A parallelisation of the code using openMP libraries. We studied how and where parallelise the code.
- The execution is run in parallel inside the same computer. All computational resourcres access to the same memory space. This allow to have high performance since there are no network bottlenecks and today computers can store huge amounts of data (a desktop computer could have dozens of GB)

## Domain of application

We have studied Combinatorial Auctions, where CPLEX/GUROBI is the de-facto standard. We tested PAR-AD3 extensively with distributions generated by CATS, a standard Combianatorial Auction generator. Then, we converted that output to a factor graph format. Note that we can encode that problems using only the At Most One factor.
Here you can find some examples of combinatorial auctions generated by [CATS](http://www.cs.ubc.ca/~kevinlb/CATS/) and then converted to Factor Graph. I have selected this ones from a pool of experiments, randomly picking one representative of  different distributions we are considering. (In data dir).



| FileName | Size | Size gzipped |  
|--------|--------|--------------|
| paths-9000-30000-00.normal.fg	| 2.6M |	830K	| 
| regions-npv-9000-40000-00.normal.fg |	3.7M |	1.0M	|
| arbitrary-npv-9000-40000-00.normal.fg	| 3.7M	| 1.6M	|

Considerations about these distributions are described at [this paper](http://www.cs.ubc.ca/~kevinlb/pub.php?u=EmpiricalHardness.pdf). General overview at section 3.2.

##Running the Code
I was cleaning up the code and making an ansi c version with compiles fine with gcc. You can downlowad it and compile it. You will recognise your code inside my implementation.
 
You can download the source code [here](https://raw.githubusercontent.com/parad3-wdp/ca/master/src/parAD3.ansi.c) (one single file).
 
For compiling:

`gcc -O3 -fopenmp parAD3.ansi.c -lm -Wno-unused-result -o parAD3`

For executing:

`./parAD3 --file_graphs=[IN] --file_posteriors=[OUT] (--max_iterations=[NUM]  --eta=[NUM]  --residual_threshold=[NUM] --threads=[NUM])`

which is a simplifaction of your AD3, but adding the threads parameter. For maximium performance match threads with the actual number of cores of your system (you can obtain it issuing nproc in a linux shell)
