#include <string.h>
#include <iostream>
#include "gen_dna_funcs.h"
#include "score_matrix.h"
#include "alignment.h"


#ifndef ___MAX_COVER_GRAPH_H___
#define ___MAX_COVER_GRAPH_H___



//How many residues of overlap should be allowed
//between alignments by default?
#define default_overlap 0


//Linked list of nodes
template <class LIST_UNIT>
class Node_list {
	public:
		Node_list *next, *last;

		Node_list()   {};
		Node_list(Node_list *new_last, LIST_UNIT *new_element);
		LIST_UNIT* get_element()      {return(this_element);};
		~Node_list()   {};
	
	protected:
		LIST_UNIT *this_element;
};


class Node {
public:
	Node_list<Node> *best_path, *best_start;


	Node()                              {};
	Node(int len, int spos, char new_name[]);
	int get_length()                    {return(length);};
	int get_start()                     {return(start);};
	int get_end()						{return(start+length-1);};
	double get_evalue()					{return(evalue);};
	double get_score()					{return(score);};
	char* get_filename()				{return(file);};
	int get_num_edges()                 {return(num_out);};
	int get_num_in_edges()				{return(num_in);};
	char* get_name()					{return(name);};
	int get_node_num()                  {return(node_num);};
	void set_node_num(int num)          {node_num=num;};
	void set_my_align(Sequence_dataset *align)     {my_alignment=align;};
	int delete_edge(Node *old_edge_node);
	void add_edge(Node *new_edge_node);
	void add_in_edge(Node *new_edge_node);
	void reset_out_edges();
	Node* get_next_out();
	void reset_in_edges();
	Node* get_next_in();
	void calc_evalue(double lambda, double K, int nlen, int mlen);	
	void set_score(double s)	{score=s;};
	void set_filename(char* f);
	BOOL is_covered()               {return(node_covered);};
	void set_covered()			    {node_covered=TRUE;};
	void set_best(double newbest)      {best=newbest;};
	double get_best()              {return(best);};
	Sequence_dataset* get_alignment()    {return(my_alignment);};

protected:
	int length, start, num_out, num_in, node_num, *edges;
	char name[50], file[100];
	double score, evalue, best;
	BOOL node_covered;
	Node_list<Node> *out_edges, *trans_edges, *start_out,
			*in_edges, *start_in, *trans_in;

	Sequence_dataset *my_alignment;
};




class Graph {
	public:
		Graph();
		Graph(int n_nodes, Node **the_nodes); 
		int get_num_nodes()             {return(num_nodes);};
		Node* get_node(int n);
		void count_in_edges();
		int num_in_edges(int i);
		void recurse_path(Node *current);
		int get_best_path_len()         {return(max_len);};
		void reset_best_path()          {best_path=best_path_start;};
		double get_align_score_loss(Node *node1, Node *node2);
		void set_aligner(BOOL blosum62, char *matrix_file, DATATYPE cdata, double match, double mismatch, double gap_open, double gap_extend);
		~Graph();

	protected:
		int num_nodes, max_len;
		Node **nodes;
		Node_list<Node> *best_path, *best_path_start;
		Dynam_Prog *the_aligner;
		Score_matrix *the_matrix;

		virtual double get_edge_weight(Node *node1, Node *node2)=0;
		virtual double get_single_score(Node *node)=0;
};



class Length_Graph : public Graph {
public:
	Length_Graph() : Graph() {};
	Length_Graph(int n_nodes, Node **the_nodes) : Graph (n_nodes, the_nodes) {}; 
protected:
	double get_edge_weight(Node *node1, Node *node2);
	double get_single_score(Node *node);

};


class Score_Graph : public Graph{
public:
	Score_Graph() : Graph() {};
	Score_Graph(int n_nodes, Node **the_nodes) : Graph (n_nodes, the_nodes) {};
protected:
	double get_edge_weight(Node *node1, Node *node2);
	double get_single_score(Node *node);
};



#endif
