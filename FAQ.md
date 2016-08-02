
## Frequently Asked Questions
PIPS
====
**Q:  How is PIPS different than general purpose solvers like CPLEX, Gurobi, IPOPT, etc.?**

A. PIPS can understand certain structured problems and use specific algorithms to solve them.
This feature isn't available on standard solvers.
PIPS is also a parallel solver, meaning that it can distribute work across nodes (other computers or CPU cores).

**Q: What if I have a weak notebook with a single or dual core CPU? Is PIPS still useful for me?**
s
A. On short, yes, but you will not get the parallelism benefits, and it will most likely be slower than general purpose solvers.

**Q: Can I use PIPS to solve standard (unstructured) problems?**

A. Yes, but it is not a good idea in general. Structure-exploiting solvers are generally slover for unstructured problems due to overhead and other algorithmic considerations.
The algorithms implemented in PIPS are designed to exploit the structure and perform (most of) the computations in parallel.

**Q: Does PIPS support multistage stochastic programming or only 2-stage problems?**

A. PIPS only supports 2-stage stochastic problems, however, be aware that multistage problems can be expressed as two-stage in most the cases. 

StructJuMP
====
**Q: What is StructJuMP?**

A. [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl) is a modeling environment that extends [JuMP](https://github.com/JuliaOpt/JuMP.jl) to easily express the structure of problems and interface efficiently with 
structure-exploiting solvers such as PIPS. It also works in parallel (distributed memory), and thus allows the specification of much  
larger (structured) problems than JuMP can handle.

**Q: What do you mean by structured problem?**

A. Problems that have block-angular structure, such as stochastic problems with recourse.

**Q: Does StructJuMP accept stochastic problem notation directly?**

A. No. You must translate your stochastic problem into the deterministic, scenario-based equivalent model.

**Q: But I can already solve deterministic equivalent problems with standard tools. Why would I use these instead?**

A. StructJuMP marks the structure of the problem and ensures streamlined interface with structure-exploiting solvers working in 
parallel, such as PIPS. As such, StructJuMP can handle much larger (structured) problems than any of the general purpose modeling environments.

**Q: Does StructJuMP support multistage stochastic programming or only 2-stage problems?**

A. StructJuMP supports multilevel structure, however, in the absence of solvers that solve multilevel problems, it has limited applicability.

StructJuMP+PIPS
====
**Q: Is StructJuMP+PIPS always the fastest method for solving structured problems?**

A. No. For small problems general purpose solutions CPLEX/Gurobi/Ipopt/etc+JuMP/AMPL/etc can be and usually are faster.
PIPS pays off as the problem size increases.

**Q: Increases by how much? Is there any rule of thumb?**

A. This is in general problem dependent. However, whenever your problem does not fit in the memory because 
it has too many scenarios, StructJuMP+PIPS can still handle it on a memory distributed computer.

**Q: Can you give me a quick recap on when to use StructJuMP+PIPS and when to use general purpose AMLs and solvers?**

A. You should use StrucJuMP+PIPS if:
* your problem is structured (for example, a stochastic problem)
* it is sufficiently large (for example, has many scenarios)
* you have a highly parallel computing hardware

Acknowledgments
====
Cosmin Petra would like to thank Oswaldo Melo and Alex Dowling for their invaluable input and help in preparing this document.
