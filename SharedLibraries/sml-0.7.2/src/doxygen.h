/* (c) 2008,2009 Jonathan Hogg and Andreas Grothey, University of Edinburgh
 *
 * This file is part of SML.
 *
 * SML is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, using version 3 of the License.
 *
 * SML is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */
/** 
\mainpage Structured Modelling Language for OOPS 

This is the doxygen documentation for the \ref language "Structured Modelling Language (SML)"
Interface to OOPS.

\section sec_introduction Introduction


\ref language "SML" is a structured modelling extension to AMPL, designed to act as a
preprocessor to AMPL. The driver performs the following tasks:
- read an SML input file,
- analyse the intended problem/matrix structure,
- create a separate (plain) AMPL file for each node in the matrix tree
- process the AMPL files through AMPL (therefore creating *.nl files),
- fill the internal data structures with information about the resulting
  matrix tree.

The OOPS driver backend will create OOPS structured matrices from the
objects, and set up the appropriate call back functions that fill in
the sparse leaf nodes of the matrix tree by processing the appropriate
*.nl file.

A different backend (for example for a decomposition solver) could be
provided. The backend is independent of AMPL (just the amplsolver
library needs to be called). The frontend is dependent on AMPL.


\section sec_installation Installation and usage

For information on how to install and use SML, please refer to the user guide.

\section sec_interface Solver Interface

SML can be interfaced with any structure exploiting solver as the backend. 

A detailed description of the \ref page_interface "Solver Interface"
can be found here.

\section Internals

The processing of an SML file is done in the following steps:
- \ref frontend "Frontend work" (reading the file into the internal data
  structures)
- \ref scriptfile "Writing out" the script file and submodel AMPL files
- Creating the \ref expanded "Expanded Model"
- \ref backend "Backend work" (creating OOPS data structures from the internal
  SML objects)

\section sec_problems Known problems

Read the \ref problems "Known Problems" page to see the current limitations
in SML.

*/

/**
\page frontend Frontend

ampl.l is the LEX/FLEX input file. It mainly consists of a list of tokens
that should be recognised in SML.

ampl.ypp is the YACC/BISON grammar file. It specifies the grammar of
SML in Backus-Naur form (BNF). Large parts of it follow the
specification of the AMPL grammar in the appendix of the AMPL book. In
YACC/BISON every "rule" returns a "value" that is computed from the
values of its components. In SML most rules return a pointer to a
SyntaxNode. A SyntaxNode represents an AMPL/SML operator, and an AMPL/SML
expression is thus represented as a tree of SyntaxNodes. 

The grammar processor recognizes the start and end of a block/submodel
and creates an AmplModel object for each of these. It classifies all
other lines into set/parameter/variable/constraint/objective/submodel
declarations - which are stored in a ModelComp object - and attaches
them to the appropriate (current) model.

\ref stochmodel "Stochastic Programming models" defined by the "block
stochastic" commands are treated specially. They will be read into a StochModel
object.

Every declaration is further divided into a name, indexing and
attribute (body) section. These are attached to the appropriate
ModelComp object. The indexing and attribute expressions are
represented as a tree of SyntaxNodes. The indexing expression is actually
represented by the subobject SyntaxNodeIx, which has extra fields and
methods that allow a quick recovery of the defined dummy variables and
their dimension.

The indexing and attributes trees are postprocessed: they are
scanned for name references. Every reference is compared with the
names of model components that are currently in scope (this includes
dummy variables that are used in the indexing expression), and if
found the reference in the SyntaxNode tree is replaced by a pointer to the
appropriate ModelComp (done by find_var_ref_in_context()).
\bug Currently no hashing is done in this search.

Internally all names are represented by a SyntaxNode with
SyntaxNode.opCode==ID. These are replaced by an object of subclass
SyntaxNodeIDREF with SyntaxNodeIDREF.opCode==IDREF, which carry a pointer to
the appropriate ModelComponent.

The output of the frontend is a tree of AmplModel objects (each
representing a node of the model tree), consisting of ModelComp
objects (each representing an AMPL declaration). These in turn consist
of several SyntaxNode trees representing the indexing and attribute (body)
section of the declaration

\note Also need to describe the data file parser and its classes
*/

/**
\page stochmodel Processing of Stochastic Programming blocks

A stochastic block definition is read in as a normal block definition (just
that it creates a StochModel object rather than an AmplModel object).
The difference is that StochModel carries information about the stochastic
parameters of the block (STAGES, NODES, PROBABILITY, PARENT).

The components of a StochModel are StochModelComp objects (rather than
ModelComp objects). The difference here is that a StochModelComp
carries information about which stages this component belongs to and
if a component is "deterministic" (i.e. there is only one copy of it
per stage, instead of one copy per node).

After the reading of the stochastic block is complete,
the StochModel object is translated into a chain of
AmplModel objects.  This is done by StochModel::expandToFlatModel().

\section expFlat Expansion to flat model tree

Expansion works in two passes. In the first pass the chain of
AmplModel objects is built (whose components are still StochModelComp
objects) In the second pass the StochModelComp objects are translated
to ModelComp objects (that is references to StochModelComp objects
are resolved to references to ModelComp objects and special
stochastic programming expressions such as Exp, node, stage are
translated).

The two passes are necessary since model components can refer to other
components up or down the tree. Hence the re-resolving of references
from StochModelComp to ModelComp can only be done once the AmplModel
tree is complete.

In detail the steps in the expansion procedure are as follows

PASS 1:
 - Expand the STAGES set for the StochModel (StochModel::expandStages())
   and for all its components (StochModel::expandStagesOfComp())
   This is done by setting up an ampl script that is processed by ampl and
   whose output is read in again
 - Create an AmplModel for each element in STAGES. Add a clone of all
   StochModelComp object that should be included in this stage to the model.
   Also add a ModelComp for the next stage down and its indexing expression
   to the AmplModel. The indexing
   expression is of the form
   \code
     set indS0 :={this_nd in NODES:A[this_nd]=="root"};
     block S1{ix0 in indS0}
     ...
       set indS1 :={this_nd in NODES:A[this_nd]==ix0};
       block S2{ix1 in indS1}
   \endcode
   Here indS0 is a new ModelComp that is added to the AmplModel before the
   ModelComp representing the subblock of the next stage.

PASS2:
  - StochModel::_transcribeComponents(): recursively call
    StochModelComp::transcribeToModelComp() for all components of all AmplModels
    within the chain.
    StochModelComp::transcribeToModelComp() will
    - create a deep copy of the StochModelComp
    - find all IDREF nodes in dependency and resolve them w.r.t  AmplModel chain
      (also resolving 'ancestor' references)
    - find all STAGE/NODE nodes and translate them
    - find all EXP expressions and translate them
      (EXP constraints are moved up to the correct level in the AmplModel chain
       Actually they are queued to be moved later to not mess up the recursion)
  - AmplModel::applyChanges
     Apply queued moves of model components (originating from Exp constraints)
  - Finally the dependency lists ModelComp::dependencies are rebuilt 
    (root->reassignDependencies)

*/

/**
\page scriptfile Writing the scriptfile and submodel AMPL files

\section submodel Writing the submodel AMPL files

This is done in backend.cpp by process_model() and write_ampl_for_submodel() called from it.

For every AmplModel object a separate *.mod file is generated. The
name of the *.mod file is obtained by concatenating all the block
names on the path from the root model down to this model. Every *.mod
file contains the subset of all given ModelComp declarations needed
for this block. They are obtained in the following manner:

- All set and parameter declarations in nodes *above* and including
  the current one are included. 

- All variable, constraint and objective declarations in the current
  node are included.

- All entities referenced by any included model component are also
included (there is a ModelComp::dependencies list that helps this
tasks. Also there is the possibility to recursively 'tag' all
dependencies of a given ModelComp object:
ModelComp::tagDependencies(), and to get a list of all tagged
ModelComp objects) 

- All included components are printed out in the order in which they appeared in the original AMPL/SML file. 

- Names of components and parameter list have to be changed. They are
translated into "global" names. Basically every block adds the name of
the block to the beginning of the name of every component defined
within it, and the dummy variable in the indexing expression for this
block becomes the first index of every component. For example a
variable Flow{i in ARCS} defined inside a block MCNF{j in COMM}
becomes MCNF_Flow{j in COMM, i in ARCS}. All references to Flow would
be changed from Flow[i] to Flow[j,i]

The writing out of one component of the *.mod file is done by
modified_write() Translating of local names into global names is done
by a static variable SyntaxNode::use_global_names. SyntaxNode has a
SyntaxNode::print() method that prints out the AMPL expression represented
by the tree. If SyntaxNode::use_global_names is set then all names of
IDREF nodes are replaced by their global name with the proper indexing
expression. 

A stack of applicable indexing expressions is kept in the global
variable l_addIndex/n_addIndex in backend.cpp 

The resulting *.mod file (and therefore the corresponding *.nl file
contains exactly the constraints defined in the block, but may have
more variables. Therefore in order to get the sparse matrix information 
for an OOPS Matrix tree node, some matching of variable names has to be done.

\section writescript Writing the scriptfile "script.scr"

Is done in process_model() in backend.cpp. It uses much the same logic
as the printing out of the *.mod files.

*/

/**
\page expanded Creating Expanded Model

The SML file describes the problem in "flat" form: that is only one
instance of each type of block is present (and therefore only one
instance of an AmplModel object). The indexing expression that
carries information about the number of instances of a block is just
passed through but not analysed.

In the ExpandedModel tree however there is one instance of each block
for each member of the indexing set.
Since the SML interpreter does not attempt
to understand AMPL enough to generate the indexing set member list
itself, the ampl-interpreter is used for this purpose. 

\section sec_naming_nodes Naming of Nodes

While nodes in the "flat" model are named by concatenating the block
names from root to the current model, in the expanded model nodes are
named by also concatenating the instance of each block (this is the
value of the dummy variable in the indexing expression in this branch). 
That is a model tree
\code
block MCNF{i in ARCS}:

 block RouteComm{j in COMM}:

 end block

end block
\endcode

results in flat tree nodes (and the AMPL submodel *.mod files) being
called root, root_MCNF and root_MCNF_RouteComm and expanded tree nodes
being called root_MCNF_A1, root_MCNF_A2, ... and
root_MCNF_RouteComm_A1_C1, root_MCNF_RouteComm_A1_C2, etc.

\section genExpModel Generation of expanded model instance lists/cardinalities

Lines 

\code 
print card(<indexing-set-name>) > "<name-of-node>.crd"
display <indexing-set-name> > "<name-of-node>.set" 
\endcode 

are included in the script file for every block definition within a
given (sub)model to print out the cardinality and the members for each
indexing set. <name-of-node> is here the name of the node in the
expanded tree (including instance names). So the cardinality (and list
of members) of COMM in the example above could be different depending
on the value of the first level dummy variable i. These different
values for the second level indexing set COMM would be reflected in
the files root_MCNF_RouteComm_A1.crd/set,
root_MCNF_RouteComm_A2.crd/set, giving the cardinality and instance
names of the root_MCNF_RouteComm submodel in the root_MCNF_A1/A2
branches of the expanded model

\section genExpModelTree Generation of the expanded model tree

The expanded model is represented as a tree of ExpandedTree
objects. The tree is generated by calling the constructor
ExpandedModel::ExpandedModel(AmplModel* ) which in turn calls
AmplModel::createExpandedModel(string, string).

*/

/**
\page backend OOPS Backend

This is done by SML_OOPS_driver() in sml-oops.cpp. A node in the OOPS
Matrix tree is represented by an OOPSBlock object. The OOPS matrices A
and Q are created in generateSML(ExpandedModel*). 

A node in the OOPS Matrix tree for A is basically represented by two
ExpandedModel objects (with their corresponding NlFile objects):

-# the first specifies the row node, i.e. the NlFile from which the constraints
should be taken,
-# the other specifies the col node, i.e. which variables should be taken.

These are the two parameters in the OOPSBlock constructor.

The specification of the variables is a bit tricky. The ExpandedModel
node specifying the column doesn't directly carry a list of relevant
expanded variable (names). All it contains is a list of variable
definitions (from the AmplModel object that spawned it) and an NlFile
object (that however contains more variable definitions than
needed). What is done is to get a list of applicable variables for
each ExpandedModel object by getting the list of variables defined in
its NlFile and comparing this against the variable declarations in the
originating AmplModel object. After obtaining this list of applicable variables
this is compared with the variables declared in the row ExpandedModel. 
\bug This part of the implementation needs a lot of comparing
strings against strings in nested loops and is quite inefficient.

This is all done in the constructor
OOPSBlock::OOPSBlock(ExpandedModel*, ExpandedModel*).

Once this is done, the setup of the OOPS Matrices is straightforward.
Each node in the OOPS matrix tree is generated with a pointer to the
corresponding OOPSBlock object being passed as the identifier (which
is subsequently used in the CallBack function)

All calls to amplsolver (the AMPL nl-file reader library) are done in
the class NlFile. This is mainly because the main amplsolver include
file "asl.h" defines many global variables, some of which clash with
C++ keywords (e.g. list). This way the use of these keywords only
needs to be avoided in NlFile.cpp

*/


/**
\page language The Structured Modelling Language (SML)

\section block Blocks (Submodels)

The main extension of SML over AMPL is the introduction of the keyword
'block':
\code
 set COMMODITIES;
 block Net{k in COMMODITIES}: 
    var Flow{ARCS} >=0;
    ...
 end block;
\endcode
A block groups a set of model declarations (set, var, param, subject
to, minimize/maximize) together. These can be repeated over an
indexing set (as any AMPL entity can). A block is a natural
representation of a subproblem in AMPL.

Blocks can be nested. Several blocks can be defined on the
same level, thus creating a tree of blocks.

All entities within the block are local variables to this block. They
are all repeated over the indexing set indicated in the block
command. The above piece of code actually defines a variable Net_Flow:
\code
 var Net_Flow{k in COMMODITIES, ARCS} >=0;
\endcode

Blocks can also be defined with the alternative syntax
\code
 set COMMODITIES;
 block Net{k in COMMODITIES}:{
    var Flow{ARCS} >=0;
    ...
 }
\endcode

\subsection Scoping 

From within the block all model components that were defined in the
block and its ancestor blocks can be used in definitions. Model
components defined in a sublock can be used by need to be accessed
"through" the name of the subblock. That is from outside the block
variable Flow can be referred to as
\code
 Net[k].Flow[j];
\endcode

Model components defined in sibling blocks (i.e. blocks defined on the
same level) or their child blocks cannot be used.

\section sml_sp Stochastic Programming

A stochastic programming block can be defined as
\code
 block alm stochastic using (NODES, PARENTS, PROBS, STAGES): {
   ...
 }
\endcode

Here 'alm' is the name of the block, 'STAGES' and 'NODES' are sets of
stages and nodes respectively, 'PARENTS{NODES}, PROBS{NODES}' are
parameter arrays indexed over the set NODES that give the parent
node and the conditional probability of a given node respectively.
The parent of the root node must be set to the string "null".
The PARENTS array must imply the same number of stages as are given
in the STAGES set.

@note PARENT is a (symbolic) parameter indexed over the set NODES that 
gives for every node in NODES the name of its parent node. For the root node 
it is required that PARENT[root] = "null"

There are a variety of specifiers that can be applied to model
components declared within a stochastic block:
<ul>
<li> Components can be declared in only some of the stages. This can be specified in two ways: 
- either by using the 'stages' keyword in the definition
  \code
    subject to Inventory{i in ASSETS} stages (STAGES diff {first(STAGES)}):
      ...
  \endcode
- or by specifying a stage block:
  \code
    set first := STAGES diff {first(STAGES)};
    stage first:
      subject to Inventory{i in ASSETS}:
        ...
    end stage
  \endcode
- also with the alternative syntax:
  \code
    set first := STAGES diff {first(STAGES)};
    stage first:{
      subject to Inventory{i in ASSETS}:
        ...
    }
  \endcode
<li> All model components declared in a stochastic block are by default
     indexed over the NODES set. Some parameters or variables only have one
     value for every stage not for every node within the stage. These can
     be defined by the 'deterministic' keyword:
    \code
      param Liability deterministic, >=0;
    \endcode

  @bug This is understood by the parser, but not implemented yet in the backend

<li> It is possible to explicitly refer to the node or stage of a component
(for example to refer to parameters declared outside the stochastic block):
\code
 param Liability{STAGES}; //an alternative deterministic parameter
 block alm stochastic using (NODES, PARENT, PROBS, STAGES): {
   subject to CashBalance:
      sum{i in ASSETS} xs[i] =  Liability[stage] + sum{i in ASSETS} xb[i];
   ...
 }
 \endcode

 Here the constraint CashBalance is repeated over all nodes in the
 stochastic block, whereas the keyword 'stage' in the constraint is replaced by
 the stage of the constraint.

 
 @bug Rather than implementing this as keywords 'node'/'stage' this
 should be done by declaring dummy variables for the NODES, STAGES
 set:
\code
 block alm stochastic using (nd in NODES, ANCESTORS, PROBS, st in STAGES):
  ..
      ...  Liability[st] +...
 \endcode

</ul>

\section exp Expectation Constraints

Add something on how to write Expectation type constraints

\section sec_problems Known problems

Read the \ref problems "Known Problems" page to see the current limitations
in SML.

*/

/**
\page problems Known Problems

@bug 
- Currently *all* parameters must be global (i.e. declared in the top
level block). This should be fixed once SML understands data files.

@bug 
- Variables that are defined over higher dimensional indexing sets
*must* have a dummy variable in their definition,i.e.
\code
  set NODES; 
  set ARCS within NODES cross NODES;
  var sparecap{(i,j) in ARCS};
\endcode
This is because SML needs to create a sum over all instances of the
variables (at least for the dummy objective). Without the dummy
variable (i,j) SML has no way of knowing that the set ARCS is
2-dimensional. This should be fixed once SML understands data files.

@bug 
- In a stochastic block all entities *must* have different names, even if
they are defined in different stages, i.e.
\code
  subject to CashBalance stages (TIME diff {first(TIME)}): ...
  subject to CashBalance1 stages ({first(TIME)}):...
\endcode
This can probably remedied by encoding the stages information in the
internally used global name somehow.


@bug
- At the moment all blocks have to contain at least a variable, otherwise
the following assertion is triggered:
\code
  smloops: MatrixSparseSimple.c:122: NewSparseMatrix: Assertion `element_sz == 0' failed.
\endcode
*/
/**
\page page_interface Solver Interface

The information about the problem and its structure as processed by
SML is stored in a tree of ExpandedModel objects. The solver can be called by

\code
SML_OOPS_driver(em)
\endcode
where 'em' is the ExpandedModel representing the root node.

@bug This behaviour should be changed. Rather we should do something like Amplsolver: 
- Either, compile SML into a library that provides a call of the form
\code
ExpandedModel *SML_process_model(char *model_file_name, char *data_file_name)
\endcode
- Or, pass the name of the solver to SML as a command line argument
\code
sml name-of-smlfile name-of-datafile name-of-solver 
\endcode
and SML would then finish with a final call of the form
\code
SML_<name-of-solver>_driver(em);
\endcode

Each ExpandedModel object represents one node of the Expanded Model
tree. It roughly corresponds to a child in a decomposition scheme
applied to solve the problem. It gives information about the variables
and constraints local to this node as well as links to any children.

It further provides an interface to AMPL that provides functions to
ask AMPL to evaluate constraint functions on this node.

*/
