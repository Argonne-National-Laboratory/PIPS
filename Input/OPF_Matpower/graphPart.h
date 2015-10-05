/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef graphPart_H
#define graphPart_H

class graphPart{

  public:
	int npartsReq;  // number of partition asked for 
	int actparts;   // actrual number of partitions
	int ncut_line;	// number of line cuts

	int *vtxdist;	// how to seperate the points in parmetis
	int *xadj;		// see Metis manual
	int *adjncy;	// see Metis manual
	int *id_part;   // part that each node should belong to
	
	int num_node;
	int num_line;
	int num_balCon; // number of balancing constraints. It should be at least 1. See Metis manual

	// weight of lines/nodes
	int *vwgt;
	int *adjwgt;
	
	int *vsize;
    double *ubvec;
	double *tpwgts;

  
	char *gname;
	char *dataname;

	int *line_status;	 // part that each line should belong to: -1, cut by metis; others: local line and part ID
    int *line_start;			 
    int *line_end;	
	char **line_name;	
	int *bus_name;	
	int *line_ID;		// corresponding line ID in oops, convert the line order in Metis graph file to the one used in mps file


	int *gen_status;	// part that each generator should belong to
	int *gen_bus;		// connected bus
	int num_gen;		// number of generators


	int DefineBlkFlag; 


    int local_gvnTxs;
	int local_numTxs;
	int *local_xadj;
	int *local_adjncy;

	int *local_adjwgt;

	int *local_vwgt;
	int *local_vsize;

	int* local_id_part;
	
		/*---------------------------------------------------------------------------
	  Methods
	---------------------------------------------------------------------------*/
	
	/* ------------------- constructors and destructors ---------------------- */
	graphPart();
	~graphPart();

	
	/* ---------------------------------- methods ------------------------------ */

  public:
	void computeGPfromGraphFile(const char gfilename[]);

  private:

	void _ReadGraph(const char gfilename[]); 
	void distributeLocalAndGlobalVar();



};


#endif
