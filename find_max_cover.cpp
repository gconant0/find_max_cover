#include <iostream>
#include <fstream>
#include <math.h>
#include "gen_dna_funcs.h"
#include "read_seq.h"
#include "linked_list.cpp"
#include "alignment.cpp"
#include "max_cover_graph.cpp"
#include <stdio.h>

#ifdef _MAKE_GENE_PLOT_
#include "plot.h"
#endif

//Timing functions do not exist under Windows
#ifdef _COMPILE_GCC
#include <sys/time.h>
#include <sys/resource.h>
#endif


//#define EVALUE_DEBUG


#ifdef _TIME_
unsigned long RunTime()
{
    struct rusage r;
    long sec, usec;

    getrusage(RUSAGE_SELF, &r);   /* constant defined in include files */ 
    sec = r.ru_utime.tv_sec;      /* whole seconds */
    usec = r.ru_utime.tv_usec;    /* microseconds */

    return((sec * 1000) + (usec / 1000));

} /* RunTime*/
#endif


#define ONE_MASK 1

void parse_args(int argc, char **argv, char filename[100], char plotfile[100], int &new_overlap, BOOL &calc_evalues, BOOL &use_score_graph);
void parse_inputfile(char *inputfile, double &gap_open, double &gap_extend, double &match, double &mismatch, int &mlen, int &nlen,
                                         double &kappa, double &lambda, double &H, DATATYPE &cdata, char *matrix_file, 
										 BOOL &blosum62, int &start_point);


void get_numstring(char *instring, char *num_string);
void get_alignments (int num_nodes, Node **the_nodes, DATATYPE cdata);
void correct_align_scores(int num_max_aligns, int overlap, Graph *the_graph, Node_list<Node> *max_combo,  double *&scores);
double calc_eval(double k, double lambda, double h, int aligns, int nlen, int mlen, double *scores);

#ifdef _MAKE_GENE_PLOT_
int read_lens(char *plotfile, int &num_genes, char **&names, int *&starts, int *&ends, int *&lens);
int make_plots(char *plotfile, int num_genes, char **names, int *starts, int *ends, int *lens, int combo_size, Node *best_node);
#endif


int main (int argc, char **argv) 
{
	int i, j, num_nodes, *stack1, *stack2, start, len, null=0,  new_overlap, 
		num_combo=0, *used_aligns, nlen, mlen, start_point=0, plotgenes, *plotstarts, *plotends, *plotlens;
	double lambda, k, h, *scores, best_Eval, score, match, mismatch, gap_open, gap_extend, best;
	long pre, post;
	char filename[100], alignment_file[100], plotfile[100], name[50],matrix_file[100], 
		inword[50], lastname[50], dump, readline[100], **plot_names;
	ifstream infile;
	BOOL already_in_list, calc_evalues, use_score_graph=FALSE, blosum62;
	DATATYPE cdata=NUCLEIC;
	Node **the_nodes, *a_node, *best_node, *save_node;
	Node_list<Node> *readnodes, *start_nodes, *search_pos;
	Graph *align_graph;
	
	if (argc>1) {

		parse_args(argc, argv, filename, plotfile, new_overlap, calc_evalues, use_score_graph);		

		
#ifdef _TIME_
		pre=RunTime();
#endif
		//If calculating E-values, get parameters first
	
		parse_inputfile(filename, gap_open, gap_extend, match, mismatch, mlen, nlen, 
				 k, lambda, h, cdata, matrix_file, blosum62, start_point);
	
			
			//Output current settings
			cout<<"Settings are:\n";
			if (calc_evalues == TRUE) {
				cout<<"Query size: "<<mlen<<" Reference sequence/DB size: "<<nlen<<endl;
				cout<<"kappa: "<<k<<" Lambda: "<<lambda<<" H: "<<h<<endl;
			}

			if(cdata == PROTEIN) {
				cout<<"Protein sequence. ";
				if (new_overlap != 0) {
					if (blosum62 ==TRUE) 
						cout<<"Using default BLOSUM62 matrix\n";
					else
						cout<<"Using matrix file "<<matrix_file<<endl;
				}
				else
					cout<<endl;
			}
			else {
				if (new_overlap != 0) 
					cout<<"Nucleotide sequence.  Match: "<<match<<" Mismatch: "<<mismatch<<endl;
				else
					cout<<"Nucleotide sequence.\n";
			}
			if (new_overlap != 0) 
				cout<<"Gap opening penalty: "<<gap_open<<" Gap extension penalty: "<<gap_extend<<endl;
		
			//Note that I assume that the gap opening penalty given is actually the
            //cost of a gap of length one: i.e. gap_open+gap_extend in my code
			gap_open-=gap_extend;
		
			
		infile.open(filename);
		strcpy(lastname, "\0");
		

		if (infile.fail()) {
			cerr<<"Error: Invalid input file "<<filename<<endl;
			return(-1);
		}

		
		for(i=0; i<start_point; i++)
			infile.getline(readline, 99);
		
		//Read the list of alignments
		score=-1e6;
		strcpy(alignment_file, "\0");
		num_nodes=1;	
		if ((calc_evalues == FALSE) && (use_score_graph == FALSE))
			infile>>name>>start>>len;	
		else {
			if (new_overlap == 0)	
				infile>>name>>start>>len>>score;
			else
				infile>>name>>start>>len>>score>>alignment_file;
		}
		
			
		a_node = new Node (len, start, name);
		a_node->set_score(score);
		a_node->set_filename(alignment_file);
		
		readnodes=new Node_list<Node>(0, a_node);
	
		start_nodes=readnodes;
		while (!infile.eof()) {
			if ((calc_evalues == FALSE) && (use_score_graph == FALSE))
				infile>>name>>start>>len;	
			else {
				if (new_overlap == 0)	
					infile>>name>>start>>len>>score;
				else
					infile>>name>>start>>len>>score>>alignment_file;
			}	
			
			//Check for duplicates in the alignment list
			
			if ((strcmp(lastname, name) != 0) && (!infile.eof())) {
				already_in_list=FALSE;
				search_pos=start_nodes;
				while ((search_pos != 0) && (already_in_list == FALSE)) {
					if ((search_pos->get_element()->get_start() == start)  && 
						(search_pos->get_element()->get_length() == len)) {
							already_in_list = TRUE;
							save_node=search_pos->get_element();
					}
					search_pos=search_pos->next;
				}
							
			
				if (already_in_list == FALSE) {
					a_node = new Node (len, start, name);
					a_node->set_score(score);
					a_node->set_filename(alignment_file);
		
					
					readnodes->next=new Node_list<Node>(readnodes, a_node);
					readnodes=readnodes->next;
					num_nodes++;
				}
				else {
				//If this version of the alignment has a better E-value, use that it
				//rather than the original alignment
					if (score>save_node->get_score())
						save_node->set_score(score);
						save_node->set_filename(alignment_file);
				}
			}	
			strcpy(lastname, name);
			if(!infile.eof())
				infile.get(dump);
		
		}
#ifdef _TIME_
		post=RunTime();
		cout<<"Runtime for reading "<<num_nodes<<" unique nodes: "<<post-pre<<endl;
#endif		
	
	//Create the graph nodes
	the_nodes=new Node*[num_nodes];
	stack1=new int[num_nodes];
	stack2=new int[num_nodes];
	

	readnodes=start_nodes;
	i=0;
	while(readnodes != 0) {
		the_nodes[i]=readnodes->get_element();
		readnodes=readnodes->next;
		i++;
	}


	if (use_score_graph == FALSE) 
		align_graph=new Length_Graph(num_nodes, the_nodes);
	else
		align_graph=new Score_Graph(num_nodes, the_nodes);

	if(((calc_evalues == TRUE) || (use_score_graph == TRUE)) && (new_overlap != 0)) {
		align_graph->set_aligner(blosum62, matrix_file, cdata, match, mismatch, gap_open, gap_extend);
		get_alignments (num_nodes, the_nodes, cdata);
	}
	


#ifdef _TIME_
	pre=RunTime();
#endif
	
	//Add an edge if two nodes overlap by no more than "overlap" residues
	for( i=0; i<num_nodes; i++) {
		for(j=0; j<num_nodes; j++)
		{
			if (j!=i) {
				if ((the_nodes[j]->get_start() > (the_nodes[i]->get_end()-new_overlap)) &&
					(the_nodes[j]->get_start() > the_nodes[i]->get_start())	)
					the_nodes[i]->add_edge(the_nodes[j]);
			}
		}

	}

#ifdef _TIME_	
	post=RunTime();
	cout<<"Runtime for initializing graph: "<<post-pre<<endl;
	pre=RunTime();
#endif	
	
	//Use depth-first search to identify the MCCA
	align_graph->count_in_edges();

	for(i=0; i<num_nodes; i++) {
		//If this is a starting node--go down it
		if (the_nodes[i]->get_num_edges() == 0) 
				align_graph->recurse_path(the_nodes[i]);
		
	}


	//Pass through the ending nodes to identify the one that ends the MCCA
	best=0;
	best_node=0;
	for(i=0; i<num_nodes; i++) {
		
		if (the_nodes[i]->get_num_edges() == 0) {
			if(the_nodes[i]->get_best() > best) {
				best=the_nodes[i]->get_best();
				best_node=the_nodes[i];
			}
		}
	}

#ifdef _TIME_
	post=RunTime();
	cout<<"Runtime for main computation: "<<post-pre<<endl;
#endif
	
	if(use_score_graph == FALSE) 
		cout<<"Best combination aligns "<<best<<" residues\n";
	else
		cout<<"Best combination has a score of "<<best<<endl;
	cout<<"Best alignment combination: ";
	readnodes=best_node->best_start;
	a_node=readnodes->get_element();
	while(readnodes != 0) {
		num_combo++;
		cout<<a_node->get_name()<<"\t";
		readnodes=readnodes->next;
		if (readnodes != 0)
			a_node=readnodes->get_element();
		
	}
	cout<<endl;
	
	cout<<"Combination contains "<<num_combo<<" alignments\n";
	used_aligns=new int [num_combo];
	readnodes=best_node->best_start;
	a_node=readnodes->get_element();
	i=0;
	while(readnodes != 0) {
		used_aligns[i++]=a_node->get_node_num();
		readnodes=readnodes->next;
		if (readnodes != 0)
			a_node=readnodes->get_element();
	}

	//Calculate E-values for MCCA if requested
	if (calc_evalues==TRUE) {
		correct_align_scores(num_combo, new_overlap, align_graph, best_node->best_start,  scores);
		if (num_combo>1) {
			best_Eval=calc_eval(k, lambda, h, num_combo, nlen, mlen, scores);
			if (best_Eval >=1e-5)
			best_Eval=-log(1.0-best_Eval);
			cout<<"E-value for combination is "<<best_Eval<<endl;
		}
		else {
				best_Eval=1e4;
				for( i=0; i<num_nodes; i++) {
					if (the_nodes[i]->get_best() ==best) {
						the_nodes[i]->calc_evalue(lambda, k, nlen, mlen);
						if (the_nodes[i]->get_evalue()<best_Eval)
							best_Eval=the_nodes[i]->get_evalue();
					}
				}
				cout<<"E-value is "<<best_Eval<<endl;
		}
			
	}

#ifdef _MAKE_GENE_PLOT_
	cout<<"Plotting : "<<plotfile<<endl;
	if (strcmp(plotfile, "\0") != 0 )
	{
		//A plot of the best alignment combination has been requested
	    if (read_lens(plotfile, plotgenes, plot_names, plotstarts, plotends, plotlens)==0) {
			make_plots(plotfile, plotgenes, plot_names, plotstarts, plotends, plotlens, num_combo, best_node);
			for(i=0; i<plotgenes; i++)
			    delete[] plot_names;
			delete[] plotstarts;
			delete[] plotends;
			delete[] plotlens;
	    }

	}
#endif
	return(0);
	}
	else {
		cerr<<"Usage: find_max_cover <infile> (-o:overlap) (-s) (-e) (-p:plotfile) \n";
		return(-1);
	}
}



void parse_args(int argc, char **argv, char filename[100], char plotfile[100], int &new_overlap, BOOL &calc_evalues, BOOL &use_score_graph)
{
	int i=2, j;
	char overstring[50];

	calc_evalues=FALSE;
	new_overlap=default_overlap;
	strcpy(filename, argv[1]);
	strcpy(plotfile, "\0");

	while (i<argc) {
		switch (argv[i][1]) {
			case 'o':
			case 'O':
				j=3;
				while(argv[i][j] != '\0')
					overstring[j-3]=argv[i][j++];
				overstring[j]='\0';
			
				new_overlap=string_to_int(overstring);
			break;
			
			case 'e':
			case 'E':
				cout<<"Calculating E-values for combinations\n";
				calc_evalues=TRUE;
			break;

			case 's':
			case 'S':
				use_score_graph=TRUE;
			break;
			case 'p':
			case 'P':
				j=3;
				while(argv[i][j] != '\0') {
					plotfile[j-3]=argv[i][j];
					j++;
				}	
				plotfile[j-3]='\0';
				cout<<"PLotfile is "<<plotfile<<endl;
				break;
		}
		i++;
	}

}



void parse_inputfile(char *inputfile, double &gap_open, double &gap_extend, 
					 double &match, double &mismatch, int &mlen, int &nlen,
					 double &kappa, double &lambda, double &H,
					 DATATYPE &cdata, char *matrix_file, BOOL &blosum62, int &start_point)

//Sets the options from the inputfile when either scores are used or E-values are required
//Also contains default values for match/mismatch, gap penalties, alignment matrices and kappa,
//lambda and H.  Note that the default values for these last three are almost certainly wrong for
//your application and will as a result yield non-meaningful E-values
{
	int i;
	char instring[100], nocasestring[100], numstring[100];
	BOOL pro_gap_set=FALSE;
	ifstream infile;

	start_point=1;
	cdata=NUCLEIC;
	gap_open=-12;
	gap_extend=-4;
	match=5;
	mismatch=-4;
	strcpy(matrix_file,"\0");
	kappa=0.138;
	lambda=0.6;
	H=0.449;
	blosum62=TRUE;

	infile.open(inputfile);
	
	infile.getline(instring, 99);
	if ( (strlen(instring) > 1) && (instring[strlen(instring)-1]+0 == 13))
	    instring[strlen(instring)-1]='\0';
	strcpy(nocasestring, instring);
	to_ucase(instring);

	if(!infile.fail()) {
		while(strcmp(instring, "BEGIN ALIGNMENTS") != 0) {
		
			if ((instring[0] > 64) && (instring[0] <91)) {
				switch (instring[0]) {
				case 'G':
					pro_gap_set=TRUE;
					get_numstring(instring, numstring);
					switch(instring[4]) {
					case 'O':
						gap_open=string_to_float(numstring);
						break;
					case 'E':
						gap_extend=string_to_float(numstring);
						break;
					}
					break;
				case 'K':
					get_numstring(instring, numstring);
					kappa=string_to_float(numstring);
					break;
				case 'L':
					get_numstring(instring, numstring);
					lambda=string_to_float(numstring);
					break;
				case 'H':
					get_numstring(instring, numstring);
					H=string_to_float(numstring);
					break;
				case 'M':
					switch(instring[1]) {
					case 'L':
						get_numstring(instring, numstring);
						mlen=string_to_int(numstring);
						break;
					case 'I':
						get_numstring(instring, numstring);
						mismatch=string_to_float(numstring);
						break;
					case 'A':
						switch(instring[3]) {
						case 'C':
							get_numstring(instring, numstring);
							match=string_to_float(numstring);
							break;
						case 'R':
							get_numstring(nocasestring, numstring);
							strcpy(matrix_file, numstring);
							blosum62=FALSE;
							break;
						}
						break;
					}
					break;
				case 'N':
					switch(instring[1]) {
					case 'L':
						get_numstring(instring, numstring);
						nlen=string_to_int(numstring);
						break;
					case 'U':  //Already set to nucleotide sequence
						break;
					}
					break;
				case 'P':
					cdata=PROTEIN;
					if (pro_gap_set==FALSE) {
						gap_open=-12;
						gap_extend=-2;
					}
					break;
			
				}
			}
		infile.getline(instring, 99);
		if ( (strlen(instring) > 1) && (instring[strlen(instring)-1]+0 == 13))
		    instring[strlen(instring)-1]='\0';
		
		strcpy(nocasestring, instring);
		to_ucase(instring);
		start_point++;
		}
		infile.close();
	}
	else 
		return;

}




void get_numstring(char *instring, char *num_string)
{
	int i=0, start, len;

	len=strlen(instring);
	while ((i<len) && (instring[i] != '=')) i++;

	if (i==len) {cerr<<"ERROR: Invalid argument in parameter file\n";}
	else {
		start=i+1;
		while ((instring[start] == ' ') || (instring[start] == '\t')) start++;
			
		for(i=start; i<strlen(instring); i++)
			num_string[i-start]=instring[i];
		num_string[strlen(instring)-start]='\0';
	}

	
}


void get_alignments (int num_nodes, Node **the_nodes, DATATYPE cdata)
{
	int i, ntaxa, nchars;
	BOOL pro_seq=FALSE;
	Read_Sequence *read_seq;
	Sequence_dataset *new_alignment;

	if (cdata==PROTEIN) 
		pro_seq=TRUE;

	for(i=0; i<num_nodes; i++) {
	  
		switch(guess_dataformat(the_nodes[i]->get_filename(), 
				strlen(the_nodes[i]->get_filename())))    
			{
				case NEXUS:
					read_seq=new Read_Nexus;
					 break;
				case PIR:
					read_seq=new Read_PIR;
					break;
				case PHYLIP:
					read_seq=new Read_Phylip_interleave;
					break;
				case FASTA:
					read_seq=new Read_FASTA;
					break;
	 
			}
	
		//Use local alignment code only to look up matrix entries
		
		new_alignment=read_seq->get_dataset(ntaxa, nchars, the_nodes[i]->get_filename(), pro_seq); 

		the_nodes[i]->set_my_align(new_alignment);
		delete read_seq;

	}


}


void correct_align_scores(int num_max_aligns, int overlap, Graph *the_graph, Node_list<Node> *max_combo,  double *&scores)
{
	int i, j, ntaxa, nchars, delta;
	BOOL pro_seq=FALSE;
	Node_list<Node> *runlist;
	

	scores=new double[num_max_aligns];


	if (overlap != 0) {
		//If overlap is allowed, we may need to correct the alignment scores so
		//that overlapping alignments do not count the same residue twice
		
		runlist=max_combo;
		i=1;
		scores[0]=max_combo->get_element()->get_score();
		runlist=runlist->next;
		while(runlist != 0) {
				scores[i]=runlist->get_element()->get_score()-
					the_graph->get_align_score_loss(runlist->last->get_element(), runlist->get_element());

			
#ifdef EVALUE_DEBUG
			cout<<"New score for alignment "<<i<<" is "<<scores[i]<<endl;
#endif
			runlist=runlist->next;
			i++;		
		}
	
		
	}
	else {
		scores=new double[num_max_aligns];
		runlist=max_combo;
		i=0;
		while(runlist != 0) {
			scores[i]=runlist->get_element()->get_score();
#ifdef EVALUE_DEBUG
			cout<<"Score for alignment "<<i<<" is "<<scores[i]<<endl;
#endif
			i++;
			runlist=runlist->next;
		}
	}
	

}


double calc_eval(double k, double lambda, double h, int aligns, int nlen, int mlen, double *scores)
//Uses Altschul and Gish integral to compute combined E-values (Methods in Enzym. 266, 460)
{
	int i, j, steps=4000, tsteps=2000;
	double Tval=0, tol=5e-2, pval, last_pval, max=10, tmax, step_size=0.07, tstep_size=0.015, 
		intval, y, t, rfact, r2fact;

	rfact=1;
	r2fact=1;
	for(i=1; i<aligns-2; i++) {
		rfact*=(double)i;
		r2fact*=(double)i;
	}
	rfact *=aligns*(aligns-1);

	for(i=0; i<aligns; i++) {
		Tval+=lambda*scores[i]-log(k*nlen*mlen);
#ifdef EVALUE_DEBUG
		cout<<"Score: "<<scores[i]<<" Sum: "<<lambda*scores[i]-log(k*nlen*mlen)<<" log: "<<log(k*nlen*mlen)<<endl;
#endif
	}

#ifdef EVALUE_DEBUG
	cout<<"T score is "<<Tval<<" facts are "<<rfact<<"\t"<<r2fact<<endl;
#endif

	pval=0;
	tmax=Tval+14;
	t=Tval;

	last_pval=10;


	while(fabs(pval-last_pval)/pval > tol) {
	//Repeat numerical integration until result converges
		last_pval=pval;
		pval=0;
	     
		if ((Tval >0) && (tmax >0))  
		  tsteps=(int)((tmax-Tval)/(double)tstep_size);		
		else if ((Tval < 0) && (tmax < 0))
		  tsteps=(int)(fabs(Tval-tmax)/(double)tstep_size);	
		else
		  tsteps=(int)((fabs(tmax)+fabs(Tval))/(double)tstep_size);	
		

		t=Tval+(tstep_size/2.0);
		for(i=0; i<tsteps; i++) {
			intval=0;
			y=(step_size/2.0);
			max=6+fabs(t);	
			steps = (int)(max/(double)step_size);
			for(j=0; j<steps; j++) {
				intval+=step_size*(pow(y, aligns-2)*exp(-exp((y-t)/aligns)));
				y+=step_size;		
			}
			pval+=tstep_size*((exp(-t)/(rfact*r2fact))*intval);
			t+=tstep_size;
		}
#ifdef EVALUE_DEBUG	
		cout<<"Inner steps: "<<steps<<" Outer: "<<tsteps<<endl;

		cout<<"Pval is "<<pval<<endl;
#endif

		step_size/=2;
		tstep_size/=2;
	}
	
	return(pval);

}


#ifdef _MAKE_GENE_PLOT_
//Uses the libplot drawing functions to make a diagram of the best alignment combination
int read_lens(char *plotfile, int &num_genes, char **&names, int *&starts, int *&ends, int *&lens)
{
	int i;
	char line[600];
	ifstream infile;

	num_genes=0;

	infile.open(plotfile);

	cout<<"Reading diagram information from file "<<plotfile<<endl;

	if (infile.fail()) 
	{
		cerr<<"Error: invalid plot file\n";
		return(-1);
	}
	
	while(!infile.eof()) {
		infile.getline(line, 599);
		num_genes++;
	}

	//Stupid c++ failing to enter eof till getline has failed
	num_genes-=1;

	
	infile.close();
	infile.clear();
	
	names=new char* [num_genes];
	starts=new int [num_genes];
	ends=new int [num_genes];
	lens=new int [num_genes];
	for(i=0; i<num_genes; i++)
		names[i]=new char [50];


	infile.open(plotfile);
	for(i=0; i<num_genes; i++) {
		infile>>names[i]>>starts[i]>>ends[i]>>lens[i];

	}
	infile.close();

	return(0);
}
	




int make_plots(char *plotfile, int num_genes, char **names, int *starts, 
			   int *ends, int *lens, int combo_size, Node *best_node)
{
	 int i, j, k, thandle, start_site, end_site, site_len, domain_start, domain_end,
	     *box_starts, *box_ends, *box_text, *box_y, q_start, 
	     q_end, *qlab1, *qlab2, *slab1, *slab2, qlaby[2], slaby[2], bigfont, smallfont;
  double spacex, spacey, box_height, box_dist, q_center, dist_per_site;
  char outputfile[100], loc_string[100];
  FILE *outfile;
  Node *a_node;
  Node_list<Node> *readnodes;

  strcpy(outputfile, plotfile);
  outputfile[strlen(plotfile)]='.';
  outputfile[strlen(plotfile)+1]='p';
  outputfile[strlen(plotfile)+2]='s';
  outputfile[strlen(plotfile)+3]='\0';
	
	cout<<"Writing plot "<<outputfile<<endl;
	outfile=fopen(outputfile, "w");
    
	spacex=1000;
	spacey=500.0;

	bigfont=spacey*0.05;
	smallfont=spacey*0.03;

	box_dist=spacey/(2*combo_size+1);

	box_height=box_dist;

	q_center=spacey/2;

	
	//Assume first entry in plotfile listing is the query sequence;

	start_site=0;
	end_site=lens[0];
	site_len=lens[0];
	domain_start=-1;
	domain_end=-1;
	box_starts=new int [combo_size];
	box_ends=new int [combo_size];
	box_text=new int [combo_size];
	box_y=new int [combo_size];
	qlab1=new int[combo_size];
	qlab2=new int[combo_size];
	slab1=new int[combo_size];
	slab2=new int[combo_size];

	readnodes=best_node->best_start;
	a_node=readnodes->get_element();
	i=0;
	while(readnodes != 0) {
		
		j=0;
		while ((j<num_genes) && (strcmp(names[j], a_node->get_name()) != 0)) {j++;};

		if (strcmp(names[j], a_node->get_name()) != 0)
		{
			cerr<<"Error: cannot find length of gene "<<a_node->get_name()<<endl;
			return(-1);
		}

		if ((starts[j] - a_node->get_start())>start_site) {
		//Does this aligned gene begin *before* the query?
		    
			start_site=starts[j]-a_node->get_start();
			cout<<"New start: "<<starts[j]<<" "<<start_site<<" "<<j<<endl;
			domain_start=j;
		}
		

		if (((lens[j]-ends[j])-(lens[0]- a_node->get_end())) > end_site) {
		//Does this aligned gene end *after* the query?
			end_site=(lens[j]-ends[j])-(lens[0]- a_node->get_end());
			domain_end=j;
		}

		readnodes=readnodes->next;
		if (readnodes != 0)
			a_node=readnodes->get_element();
	}

	end_site += start_site;
	dist_per_site=(double)(start_site+lens[0]+end_site)/spacex;

	q_start=(int)(start_site/dist_per_site);
	q_end=(int)(end_site/dist_per_site);


	qlaby[0]=(int)(q_center+0.375*box_height);
	qlaby[1]=(int)(q_center-0.375*box_height);

	slaby[0]=(int)(0.125*box_height);
	slaby[1]=(int)(-0.125*box_height);

	cout<<"Query box starts "<<q_start<<" ends "<<q_end<<endl;

	readnodes=best_node->best_start;
	a_node=readnodes->get_element();
	for(i=0; i<combo_size; i++)
	{
		j=0;
		while ((j<num_genes) && (strcmp(names[j], a_node->get_name()) != 0)) {j++;};

		if (strcmp(names[j], a_node->get_name()) != 0)
		{
			cerr<<"Error: cannot find length of gene "<<a_node->get_name()<<endl;
			return(-1);
		}

		box_starts[i]=(int) ((a_node->get_start()-starts[j]+ start_site)/dist_per_site);
		box_ends[i]=(int)(((a_node->get_start()-starts[j]+ start_site)+lens[j])/dist_per_site);
		box_text[i]=(int)((box_ends[i]-box_starts[i])/2+box_starts[i]);
		if (ONE_MASK & i == 0)
		    box_y[i]=(int) (q_center+0.5*box_height+ (i+1)*box_height);
		else
		    box_y[i]=(int) (q_center-0.5*box_height - (i)*box_height);

		cout<<"Box "<<i<<" starts "<<box_starts[i]<<" ends "<<box_ends[i]<<" text "<<box_text[i]<<endl;

		qlab1[i]=(int) ((a_node->get_start()+start_site)/dist_per_site);
		qlab2[i]=(int) ((a_node->get_end()+start_site)/dist_per_site);
		slab1[i]=qlab1[i];
		slab2[i]=qlab2[i];
		
		if (a_node->get_start()/dist_per_site < 10)
		    qlab1[i]+=10;
		if (a_node->get_end()/dist_per_site +10  > lens[0]/dist_per_site)
		    qlab2[i]-=10;


		if (starts[j]/dist_per_site < 10)
		    slab1[i]+=10;

		if (ends[j]/dist_per_site +10 > lens[j]/dist_per_site)
		    slab2[i] -=10;

		readnodes=readnodes->next;
		if (readnodes != 0)
		    a_node=readnodes->get_element();

	}


	/* set a Plotter parameter */
	/*pl_parampl ("PAGESIZE", "letter"); */ 
	
	/* create a Postscript Plotter that writes to standard output */
	if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0)
	{
	    fprintf (stderr, "Couldn't create Plotter\n");
	    return 1;
	}
	pl_selectpl (thandle);       /* select the Plotter for use */
	
	if (pl_openpl () < 0)       /* open Plotter */
	{
	    fprintf (stderr, "Couldn't open Plotter\n");
	    return 1;
	}
	pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
	pl_flinewidth (1.0);       /* line thickness in user coordinates */
	pl_pencolorname ("black");    /* path will be drawn in red */
	pl_erase ();                /* erase Plotter's graphics display */
	
	
	//Draw the Query Sequence
	cout<<"Start site: "<<start_site<<" End: "<<end_site<<endl;
	
	pl_fillcolor (65535, 65535, 65535); 
	pl_box(q_start, (int)(q_center+(box_height/2)), 
		q_end, (int)(q_center-(box_height/2)));

	pl_fontsize(bigfont);
	pl_fmove ((int)(((end_site-start_site)/2+start_site)/dist_per_site), q_center);
	pl_alabel ('c', 'c', names[0]); 


	readnodes=best_node->best_start;
	a_node=readnodes->get_element();

	for(i=0; i<combo_size; i++) {
		j=0;
		while ((j<num_genes) && (strcmp(names[j], a_node->get_name()) != 0)) {j++;};

		if (strcmp(names[j], a_node->get_name()) != 0)
		{
			cerr<<"Error: cannot find length of gene "<<a_node->get_name()<<endl;
			return(-1);
		}

		pl_fontsize(smallfont);
		int_to_string(loc_string, 99, a_node->get_start());
		pl_fmove(qlab1[i], qlaby[ONE_MASK & i]);
		pl_alabel('c', 'c', loc_string);

		int_to_string(loc_string, 99, a_node->get_end());
	       	pl_fmove(qlab2[i], qlaby[ONE_MASK & i]);
	       	pl_alabel('c', 'c', loc_string);

		int_to_string(loc_string, 99, starts[j]);
	       	pl_fmove(slab1[i], box_y[i]+slaby[ONE_MASK & i]);
	       	pl_alabel('c', 'c', loc_string);

		int_to_string(loc_string, 99, ends[j]);
	       	pl_fmove(slab2[i], box_y[i]+slaby[ONE_MASK & i]);
	       	pl_alabel('c', 'c', loc_string);
		
		
		pl_fontsize(bigfont);

		if (ONE_MASK & i == 0) {
		//Draw above the query
		    pl_box(box_starts[i], box_y[i], box_ends[i], box_y[i]+(int)box_height);
			   
		    pl_line ((int)((a_node->get_start()+start_site)/dist_per_site), box_y[i], 
			     (int)((a_node->get_start()+start_site)/dist_per_site), (int)(q_center+(box_height/2))); 
		    pl_line ((int)((a_node->get_end()+start_site)/dist_per_site), box_y[i], 
			     (int)((a_node->get_end()+start_site)/dist_per_site), (int)(q_center+(box_height/2))); 
		    
		    pl_fmove (box_text[i], box_y[i]+(int)(box_height/2));
		    pl_alabel ('c', 'c', names[j]); 
	
		}
		else {
		//Draw below the query
		    pl_box(box_starts[i], box_y[i], box_ends[i], box_y[i]-(int)box_height);
		    pl_line ((int)((a_node->get_start()+start_site)/dist_per_site), box_y[i], 
			     (int)((a_node->get_start()+start_site)/dist_per_site), (int)(q_center-(box_height/2))); 
		    pl_line ((int)((a_node->get_end()+start_site)/dist_per_site), box_y[i], 
			     (int)((a_node->get_end()+start_site)/dist_per_site), (int)(q_center-(box_height/2))); 

		    pl_fmove (box_text[i], box_y[i]-(int)(box_height/2));
		    pl_alabel ('c', 'c', names[j]); 
	
		}

		readnodes=readnodes->next;
		if (readnodes != 0)
			a_node=readnodes->get_element();

	}

  if (pl_closepl () < 0)      /* close Plotter */
    {
      fprintf (stderr, "Couldn't close Plotter\n");
      return 1;
    }

  delete[] box_starts;
  delete[] box_ends;
  delete[]box_text;
  delete[] box_y;
  delete[] qlab1;
  delete[] qlab2;
  delete[] slab1;
  delete[] slab2;


  pl_selectpl (0);            /* select default Plotter */
  if (pl_deletepl (thandle) < 0) /* delete Plotter we used */
   {
      fprintf (stderr, "Couldn't delete Plotter\n");
      return 1;
    }
  return 0;


}

#endif
