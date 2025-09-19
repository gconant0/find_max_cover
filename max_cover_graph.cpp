#include "max_cover_graph.h"


Node::Node(int len, int spos, char new_name[])
{
	evalue=1e4;
	length=len;
	start=spos;
	num_out=0;
	num_in=0;
	out_edges=start_out=trans_edges=0;
	in_edges=start_in=trans_in=0;
	best_path=best_start=0;
	strcpy(name, new_name);
	node_covered=FALSE;
}


void Node::set_filename(char* f) 
{
	strcpy(file, f);
}


void Node::calc_evalue(double lambda, double K, int nlen, int mlen)
//Calculates the E-value for this isolated node
{
	double b_score;

	b_score=(lambda*score-log(K))/log(2);
	evalue=(double)nlen*(double)mlen*pow(2, -b_score);
}

void Node::add_edge(Node *new_edge_node)
{
	if (out_edges == 0) {
		out_edges=new Node_list<Node>(0, new_edge_node);
		start_out=out_edges;
	}
	else {
		out_edges->next=new Node_list<Node>(out_edges, new_edge_node);
		out_edges=out_edges->next;
	}
	num_out++;
}



void Node::add_in_edge(Node *new_edge_node)
{
	if (in_edges == 0) {
		in_edges=new Node_list<Node>(0, new_edge_node);
		start_in=in_edges;
	}
	else {
		in_edges->next=new Node_list<Node>(in_edges, new_edge_node);
		in_edges=in_edges->next;
	}
	num_in++;
}



void Node::reset_out_edges()
//Resets the linked list that contains
//the edge list
{
	trans_edges=start_out;
}


Node* Node::get_next_out()
{
	Node* ret_val=0;

	if (trans_edges != 0) {
 		ret_val=trans_edges->get_element();
		trans_edges=trans_edges->next;
	}
	return(ret_val);
}



void Node::reset_in_edges()
//Resets the linked list that contains
//the edge list
{
	trans_in=start_in;
}



Node* Node::get_next_in()
{
	Node* ret_val=0;

	if (trans_in != 0) {
 		ret_val=trans_in->get_element();
		trans_in=trans_in->next;
	}
	return(ret_val);
}


template <class LIST_UNIT>
Node_list<LIST_UNIT>::Node_list(Node_list *new_last, LIST_UNIT *new_element)
{
	last=new_last;
	this_element=new_element;
	next=0;
}



Graph::Graph() 
//Default constructor for the Graph class:
//Should never be called                        
{
	num_nodes=0;
	cerr<<"Wrong Graph constructor\n";
}


Graph::Graph(int n_nodes, Node **the_nodes)
//Valid constructor for Graph class
{
	int i;
	
	num_nodes=n_nodes;
	nodes = new Node* [num_nodes];		

	for (i=0; i<num_nodes; i++) {
		nodes[i]=the_nodes[i];
		the_nodes[i]->set_node_num(i);	
	}

	max_len=0;
	best_path=best_path_start=0;
}




Node* Graph::get_node(int n)
//Function that prevents accessing an
//out-of-bounds node
{
	if (n<num_nodes)
		return(nodes[n]);
	else
		return(0);
}



void Graph::count_in_edges()
{
	int i;
	Node *current;
	
	for(i=0; i<num_nodes; i++) {
		if (nodes[i]->get_num_edges() != 0) {
			nodes[i]->reset_out_edges();
			current=nodes[i]->get_next_out();
		
			while (current!=0) {
				current->add_in_edge(nodes[i]);
				current=nodes[i]->get_next_out();
			}
		}

	}


}


int Graph::num_in_edges(int i)
{
	return(nodes[i]->get_num_in_edges());
}

void Graph::recurse_path(Node *current)
//Performs Depth-first search on the graph.
//If all of a node's in edges lead to nodes that
//have already been searched, the best combination 
//of alignments ending with the node "current" is 
//computed, and current is marked as complete.
//If the in-nodes are not complete, the program recurses
//to them and tries to compute their combination of
//alignments.
{
	double new_best, partial;
	Node *next_node, *best_par;
	Node_list<Node> *read_list;

	
	if (current->get_num_in_edges() != 0) {
		//Check if in-nodes are complete
		current->reset_in_edges();
		next_node=current->get_next_in();

		while (next_node != 0) {
			if (next_node->is_covered() == FALSE)
				recurse_path(next_node);
			next_node=current->get_next_in();
		}

		//Find best combination for this node
		new_best=0;
		best_par=0;

		current->reset_in_edges();
		next_node=current->get_next_in();
		while (next_node != 0) {
			partial=get_edge_weight(current, next_node);

			if (partial > new_best)
			{
				new_best=partial;
				best_par=next_node;
			}
			next_node=current->get_next_in();
		}

		if (new_best != 0) {
			current->set_best(new_best);
			current->set_covered();
			read_list=best_par->best_start;
			current->best_path=new Node_list<Node>(0, read_list->get_element());
			current->best_start=current->best_path;
			read_list=read_list->next;
			while(read_list != 0 )
			{
				current->best_path->next=new Node_list<Node>(current->best_path, read_list->get_element());
				current->best_path=current->best_path->next;
				read_list=read_list->next;
			}

			current->best_path->next=new Node_list<Node>(current->best_path, current
				);
			current->best_path=current->best_path->next;

		}
		else 
			cerr<<"Error in recurse_path: 0 length\n";
		
	}
	else {
		//Nodes without in edges are trivially computed
		current->set_covered();
		current->set_best(get_single_score(current));
		current->best_path=new Node_list<Node>(0, current);
		current->best_start=current->best_path;
	}


}



double Graph::get_align_score_loss(Node *node1, Node *node2)
{
	int i, delta;
	double score_loss1=0, score_loss2=0;

	delta=(node1->get_start()+ node1->get_length())-node2->get_start();	
		
	if (delta > 0) {
		score_loss1=
			the_aligner->get_matrix_entry( (*node2->get_alignment())[0][0], (*node2->get_alignment())[1][0]);

		for(i=1; i<delta; i++)
			score_loss1+=the_aligner->get_matrix_entry((*node2->get_alignment())[0][i], (*node2->get_alignment())[1][i], 
							(*node2->get_alignment())[0][i-1],(*node2->get_alignment())[1][i-1]);

		score_loss2=
			the_aligner->get_matrix_entry( (*node1->get_alignment())[0][node1->get_length()-1], 
				(*node1->get_alignment())[1][node1->get_length()-1]);

		for(i=1; i<delta; i++)
			score_loss2+=the_aligner->get_matrix_entry((*node1->get_alignment())[0][(*node1->get_alignment())[0].Sequence_size()-(i+1)],
							(*node1->get_alignment())[1][(*node1->get_alignment())[0].Sequence_size()-(i+1)], 
							(*node1->get_alignment())[0][(*node1->get_alignment())[0].Sequence_size()-i],
							(*node1->get_alignment())[1][(*node1->get_alignment())[0].Sequence_size()-i]);

	}

	if (score_loss2 < score_loss1)
		return(score_loss2);		

	else
		return(score_loss1);
}



void Graph::set_aligner(BOOL blosum62, char *matrix_file, DATATYPE cdata, double match, double mismatch, double gap_open, double gap_extend)
{
	if (cdata == PROTEIN) {
		if (blosum62 == FALSE) 
			the_matrix=new File_AA_matrix(matrix_file);			
		else
			the_matrix = new BLOSUM_62_matrix();
	}
	else
		the_matrix=new Nuc_matrix();

	the_aligner=new Local_Dynam_Prog(cdata, gap_extend, gap_open, the_matrix);
	}



Graph::~Graph()
{
	if (num_nodes > 0)
		delete[] nodes;
}


double Length_Graph::get_edge_weight(Node *node1, Node *node2) 
{
	double retval;

	retval=node2->get_best()+node1->get_length();

	if (node1->get_start() <= (node2->get_end())) 
			retval -= (node2->get_end())-node1->get_start()+1;

	return(retval);
}

double Length_Graph::get_single_score(Node *node)
{
	return(node->get_length());
}



double Score_Graph::get_edge_weight(Node *node1, Node *node2) 
{
	double retval;

	retval=node2->get_best()+node1->get_score();

	if (node1->get_start() <= (node2->get_end())) 
			retval -= get_align_score_loss(node2, node1);

	return(retval);
}


double Score_Graph::get_single_score(Node *node)
{
	return(node->get_score());
}

