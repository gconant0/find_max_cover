#include "gen_dna_funcs.h"
#include "score_matrix.h"
#include <fstream>
#include <iostream>

using namespace::std;

#ifndef ___ALIGNMENT_H___
#define ___ALIGNMENT_H___

enum directions {UP, LEFT, DIAG};
enum ALIGN_TYPE {GLOBAL, LOCAL};


class Matrix_entry {
public:
	Matrix_entry & operator= ( Matrix_entry & assign_from);

  int startx, starty;                    //Used for local alignments to ID start	
  int gap_extend_left, gap_extend_up;
  double score, left_score, up_score;
  BOOL tolast[3];
};



class Dynam_Prog {
 public:
  
  Dynam_Prog () {cerr<<"Call to default Dynam_Prog constructor\n";};  
  Dynam_Prog (DATATYPE data, double extend, double open, Score_matrix *matrix);

  Matrix_entry *path_line1, *path_line2, best_entry;

  void do_alignment(Molecule_Sequence *seq1, Molecule_Sequence *seq2);
  virtual double get_max_score()=0;
  int get_maxloc_x()			{return(maxloc1);};
  int get_maxloc_y()			{return(maxloc2);};
  int get_matrix_size()			{return(matrix_size);};
  double get_matrix_entry(int site1, int site2);
  double get_matrix_entry(int site1, int site2, int prevsite1, int prevsite2);
  double get_gap_open()				{return(gap_open);};
  double get_gap_extend()			{return(gap_extend);};
  DATATYPE get_datatype()			{return(curr_datatype);};
  void clear_path_arrays();
  

    virtual ~Dynam_Prog() {};
 protected:
  int matrix_size, size1, size2, maxloc1, maxloc2;
  int (*convert_char)(char);
  char code[24], matrix_file[100];
  //Sim matrix removed to used newly built class--GCC 6/3/04
  //double **sim_matrix, gap_open, gap_extend, blosum[24][24], align_max_score;
  double gap_open, gap_extend, blosum[24][24], align_max_score;  
  DATATYPE curr_datatype;
  Molecule_Sequence *sequence1, *sequence2; 
  Score_matrix *sim_matrix;
  
  
  void calc_similarity_matrix();
  void create_path_array();
  virtual void initialize_array()=0;  	
  virtual void get_direction(int loc_i, int loc_j)=0;
  virtual void get_first_element()=0;
  

};

class Global_Dynam_Prog : public Dynam_Prog {
 public:
  Global_Dynam_Prog () :Dynam_Prog() {};
  Global_Dynam_Prog (DATATYPE data, double extend, double open, Score_matrix *matrix) : Dynam_Prog(data, extend, open, matrix) {};
 
    ~Global_Dynam_Prog();

  Sequence_dataset * align_sequences(Sequence_dataset *seq1, Sequence_dataset *seq2);
  double get_max_score() {return(0);};


 protected:
  void initialize_array();
  void get_direction(int loc_i, int loc_j);
  void get_first_element();
};

class Local_Dynam_Prog : public Dynam_Prog {
public:
  Local_Dynam_Prog () :Dynam_Prog() {};
  Local_Dynam_Prog (DATATYPE data, double extend, double open, Score_matrix *matrix) : Dynam_Prog(data, extend, open, matrix) {};
    ~Local_Dynam_Prog();

  Sequence_dataset * align_sequences(Sequence_dataset *seq1, Sequence_dataset *seq2);
  double get_max_score() {return(align_max_score);};

protected:
  void initialize_array();
  void get_direction(int loc_i, int loc_j);
  void get_first_element();
};


class Align {
public:
	Sequence_dataset *aligned_sequences;

	Align ();
	Align (DATATYPE data, double extend, double open, ALIGN_TYPE type);
	Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, Score_matrix *mat, BOOL internal);
	Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, Score_matrix *mat);
	Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, double match, double mismatch);
	Align (DATATYPE data, char *matrix_loc, double extend, double open, ALIGN_TYPE type);
	void recurse_alignment(Molecule_Sequence *seq1, Molecule_Sequence *seq2, BOOL toplevel);
	double get_max_score()			{return(max_score);};
	int align_start(int sequence)	{return(starts[sequence]);};
	int align_end(int sequence)		{return(ends[sequence]);};
	~Align ();

protected:
	int seq1used[2], seq2used[2], starts[2], ends[2];
	double max_score;
	BOOL own_matrix;
	DATATYPE my_data;
	ALIGN_TYPE alignment_type;
	Molecule_Sequence *sequence1, *sequence2, *newseqs1[2], *newseqs2[2];
	Sequence_dataset *partial_align[2];
	Dynam_Prog *DP_Obj[2], *Local_DP_Obj;				
	Align *Align_obj[2];
	Score_matrix *the_matrix;
	int (*convert_char)(char);
	
	
	BOOL recover_alignment(int number, int pos);
	void assemble_alignment();
};
#endif
