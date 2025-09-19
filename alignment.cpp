#include <iostream>
#include <iomanip>
#include "alignment.h"

using namespace::std;

Matrix_entry & Matrix_entry::operator= (Matrix_entry & assign_from)
{
	int i;

	score=assign_from.score;
	
	for(i=0; i<3; i++)
		tolast[i]=assign_from.tolast[i];

	gap_extend_left=assign_from.gap_extend_left;
	gap_extend_up=assign_from.gap_extend_up;
	startx=assign_from.startx;
	starty=assign_from.starty;
	return(*this);
}




Dynam_Prog::Dynam_Prog (DATATYPE data, double extend, double open, Score_matrix *matrix)
{
  int i, j;

  sim_matrix=0;
  path_line1=path_line2=0;

  curr_datatype=data;
  gap_open=open;
  gap_extend=extend;
  align_max_score=0;
  maxloc1=maxloc2=0;
  sim_matrix=matrix;
  
  if (curr_datatype == NUCLEIC) { 
    convert_char=&readchar_to_base;
  }
  else { 
    convert_char=&readchar_to_aa;
  }
 
	
}




double Dynam_Prog::get_matrix_entry(int site1, int site2) 
{
	return((*sim_matrix)[site1][site2]);
}


double Dynam_Prog::get_matrix_entry(int site1, int site2, int prevsite1, int prevsite2)
{
	double val;
	if (curr_datatype == PROTEIN) {
		if (num_to_aa(site1) == '-') {
			if (num_to_aa(prevsite1) == '-')
				val=gap_extend;
			else
				val=gap_open+gap_extend;
		}
		else if (num_to_aa(site2) == '-') {
			if (num_to_aa(prevsite2) == '-')
				val=gap_extend;
			else
				val=gap_open+gap_extend;
		}
		else
			val=(*sim_matrix)[site1][site2];
	}
	else {
		if (num_to_base(site1) == '-') {
			if (num_to_base(prevsite1) == '-') 
				val=gap_extend;
			else
				val=gap_open+gap_extend;
		}
		else if (num_to_base(site2) == '-')  {
			if (num_to_base(prevsite2) == '-')
				val=gap_extend;
			else
				val=gap_open+gap_extend;
		}
		else
			val=(*sim_matrix)[site1][site2];

	}
	return(val);
}









void Dynam_Prog::create_path_array()
{
  int i;

  
  path_line1=new Matrix_entry [size1+1];
  path_line2=new Matrix_entry [size1+1];
 
 
  path_line1[0].score=0.0;
  path_line1[0].left_score=path_line1[0].up_score=0.0;
  path_line1[0].gap_extend_left=-1;
  path_line1[0].gap_extend_up=-1;
  
  for (i=0; i<3; i++)
      path_line1[0].tolast[i]=FALSE;
    
   initialize_array();
}


void Dynam_Prog::clear_path_arrays()
{
	if (path_line1 != 0)
		delete[] path_line1;

	if (path_line2 != 0)
		delete[] path_line2;

	path_line1=path_line2=0;

}




void Dynam_Prog::calc_similarity_matrix()
{
  

}




void Dynam_Prog::do_alignment(Molecule_Sequence *seq1, Molecule_Sequence *seq2)
{
	int i, j;
    Matrix_entry *dummy;	


	sequence1=seq1;
	sequence2=seq2;
	size1=sequence1->Sequence_size();
	size2=sequence2->Sequence_size();


	create_path_array();

	for(i=1; i<size2+1; i++) 
	{
		get_first_element();
		for(j=1; j<size1+1; j++)
			get_direction(j,i);
			
		
#if 0
		if (size1 < 15) {
		for(j=1; j<15; j++) {
			if (path_line2[j].tolast[DIAG] == TRUE)
				cout<<"  \\";
			else
				cout<<"   ";
			
			if (path_line2[j].tolast[UP] == TRUE)
				cout<<"|   ";
			else
				cout<<"    ";

		}
		cout<<endl;
		for(j=1; j<15; j++) {
			if (path_line2[j].tolast[LEFT] == TRUE)
				cout<<"<- ";
			else
				cout<<"  ";
					
			cout<<setw(3)<<path_line2[j].score<<"  ";
		}
		cout<<endl;
		}
#endif

		dummy=path_line1;
		path_line1=path_line2;
		path_line2=dummy;

    }
	//Does a final swap so the lines are in the expected order
	dummy=path_line1;
	path_line1=path_line2;
	path_line2=dummy;

	
}
  



void Global_Dynam_Prog::initialize_array() 
{
  int i, j;

  path_line1[1].score=path_line1[0].score+gap_extend+gap_open;
  
  for(j=0; j<3; j++)
    path_line1[1].tolast[j]=FALSE;
  
  path_line1[1].tolast[LEFT]=TRUE; 
  path_line1[1].gap_extend_left=1;
  path_line1[1].gap_extend_up=-1;
  path_line1[1].left_score=path_line1[1].score;			
	
  
  for(i=2; i<size1+1; i++)
    { 
		path_line1[i].score=path_line1[i-1].score+gap_extend;		
	    path_line1[i].gap_extend_left=path_line1[i-1].gap_extend_left;
		path_line1[i].left_score=path_line1[i].score;
		
		for(j=0; j<3; j++)
			path_line1[i].tolast[j]=FALSE;
      
		path_line1[i].tolast[LEFT]=TRUE;
		path_line1[i].gap_extend_up=-1;
    }
  
  
}


void Global_Dynam_Prog::get_direction(int loc_i, int loc_j)
{
	int i;
	double score_diag, score_left, score_up, max_score;	


	for (i=0; i<3; i++)
		path_line2[loc_i].tolast[i]=FALSE;
	path_line2[loc_i].gap_extend_left=-1;
	path_line2[loc_i].gap_extend_up=-1;
	


	score_diag=path_line1[loc_i-1].score+(*sim_matrix)[(*sequence1)[loc_i-1]][(*sequence2)[loc_j-1]];
	
	if (path_line2[loc_i-1].gap_extend_left!= -1) {
		//This isn't the second entry in the line
		score_left=path_line2[loc_i-1].score+gap_open+gap_extend;
		if (score_left > path_line2[loc_i-1].left_score+gap_extend) {
			path_line2[loc_i].left_score=score_left;
			path_line2[loc_i].gap_extend_left=loc_i;
		}
		else {
			score_left=path_line2[loc_i-1].left_score+gap_extend;
			path_line2[loc_i].left_score=score_left;
			path_line2[loc_i].gap_extend_left=path_line2[loc_i-1].gap_extend_left;
		}
 
	}
	else {
		//This is the first possible place to open a gap
		score_left=path_line2[loc_i-1].score+gap_open+gap_extend;
		path_line2[loc_i].left_score=score_left;
		path_line2[loc_i].gap_extend_left=loc_i;
	}

	if(path_line1[loc_i].gap_extend_up != -1) {
		//This isn't the second entry in the column
		score_up=path_line1[loc_i].score+gap_open+gap_extend;
		if (score_up > path_line1[loc_i].up_score+gap_extend) {
			path_line2[loc_i].up_score=score_up;
			path_line2[loc_i].gap_extend_up=loc_j;
		}
		else {
			score_up=path_line1[loc_i].up_score+gap_extend;
			path_line2[loc_i].up_score=score_up;
			path_line2[loc_i].gap_extend_up=path_line1[loc_i].gap_extend_up;
		}
	}
	else {
		//This is the first possible place to open a gap
		score_up=path_line1[loc_i].score+gap_open+gap_extend;
		path_line2[loc_i].up_score=score_up;
		path_line2[loc_i].gap_extend_up=loc_j;
	}

	max_score=score_diag;
	

	if(score_left > max_score) 
		max_score=score_left;

	if(score_up > max_score)
		max_score=score_up;

	path_line2[loc_i].score=max_score;

	if (score_diag == max_score)
		path_line2[loc_i].tolast[DIAG]=TRUE;

	if (score_left == max_score) {
		path_line2[loc_i].tolast[LEFT]=TRUE;
		if (path_line2[loc_i-1].gap_extend_left == -1) {
			path_line2[loc_i].gap_extend_left=loc_i;
			path_line2[loc_i].left_score=max_score;
		}
	}

	if (score_up == max_score) {
		path_line2[loc_i].tolast[UP]=TRUE;
		if(path_line1[loc_i].gap_extend_up == -1) {
			path_line2[loc_i].gap_extend_up=loc_j;
			path_line2[loc_i].up_score=max_score;
		}

	}
  
}


Global_Dynam_Prog::~Global_Dynam_Prog()
{
  int i;

  
  if (path_line1 != 0 )
      delete[] path_line1;
  if (path_line2 != 0 )
      delete[] path_line2;
    

}


void Global_Dynam_Prog::get_first_element()
{
	path_line2[0].tolast[DIAG]=FALSE;
	path_line2[0].tolast[LEFT]=FALSE;
	path_line2[0].tolast[UP]=TRUE;
	path_line2[0].gap_extend_left=-1;

	if(path_line1[0].tolast[UP] == TRUE) {
		//We're after postition 1 in the alignment
		path_line2[0].score=path_line1[0].score+gap_extend;
		path_line2[0].gap_extend_up=path_line1[0].gap_extend_up;
		path_line2[0].up_score=path_line2[0].score;
	}
	else {
		path_line2[0].score=path_line1[0].score+gap_open+gap_extend;
		path_line2[0].gap_extend_up=1;
		path_line2[0].up_score=path_line2[0].score;
	}


}







void Local_Dynam_Prog::initialize_array() 
{
  int i, j;

  path_line1[0].score=0;

  //Our best alignment so far is null
  path_line1[0].startx=path_line1[0].starty=-1;

  for (j=0; j<3; j++)
	path_line1[0].tolast[j]=FALSE;
  path_line1[0].gap_extend_left=path_line1[0].gap_extend_up=-1;

  for(i=1; i<size1+1; i++)
    { 
      path_line1[i].score=0;
	  
	  //Still a null best alignment
	  path_line1[i].startx=-1;
	  path_line1[i].starty=-1;

      
      for (j=0; j<3; j++)
		path_line1[i].tolast[j]=FALSE;
      path_line1[i].gap_extend_left=path_line1[i].gap_extend_up=-1;

    }
  
}



void Local_Dynam_Prog::get_direction(int loc_i, int loc_j)
{
	int i;
	double score_diag, score_left, score_up, max_score;	


	for (i=0; i<3; i++)
		path_line2[loc_i].tolast[i]=FALSE;

	path_line2[loc_i].gap_extend_left=path_line2[loc_i].gap_extend_up=-1;

	score_diag=path_line1[loc_i-1].score+(*sim_matrix)[(*sequence1)[loc_i-1]][(*sequence2)[loc_j-1]];
	
	if (path_line2[loc_i-1].gap_extend_left!= -1) {
		//This isn't the second entry in the line
		score_left=path_line2[loc_i-1].score+gap_open+gap_extend;
		if (score_left > path_line2[loc_i-1].left_score+gap_extend) {
			path_line2[loc_i].left_score=score_left;
			path_line2[loc_i].gap_extend_left=loc_i;
		}
		else {
			score_left=path_line2[loc_i-1].left_score+gap_extend;
			path_line2[loc_i].left_score=score_left;
			path_line2[loc_i].gap_extend_left=path_line2[loc_i-1].gap_extend_left;
		}
 
	}
	else {
		//This is the first possible place to open a gap
		score_left=path_line2[loc_i-1].score+gap_open+gap_extend;
		path_line2[loc_i].left_score=score_left;
		path_line2[loc_i].gap_extend_left=loc_i;
	}

	if(path_line1[loc_i].gap_extend_up != -1) {
		//This isn't the second entry in the column
		score_up=path_line1[loc_i].score+gap_open+gap_extend;
		if (score_up > path_line1[loc_i].up_score+gap_extend) {
			path_line2[loc_i].up_score=score_up;
			path_line2[loc_i].gap_extend_up=loc_j;
		}
		else {
			score_up=path_line1[loc_i].up_score+gap_extend;
			path_line2[loc_i].up_score=score_up;
			path_line2[loc_i].gap_extend_up=path_line1[loc_i].gap_extend_up;
		}
	}
	else {
		//This is the first possible place to open a gap
		score_up=path_line1[loc_i].score+gap_open+gap_extend;
		path_line2[loc_i].up_score=score_up;
		path_line2[loc_i].gap_extend_up=loc_j;
	}

	
	max_score=score_diag;

	if(score_left > max_score) 
		max_score=score_left;

	if(score_up > max_score)
		max_score=score_up;


	if (max_score > 0) {
		path_line2[loc_i].score=max_score;

		if (score_left == max_score) {
			path_line2[loc_i].tolast[LEFT]=TRUE;
			if(path_line2[loc_i-1].gap_extend_left == -1) {
				path_line2[loc_i].gap_extend_left=loc_i;
				path_line2[loc_i].left_score=max_score;
			}
			//We can't *start* a local alignment with a gap--this must be a continuation
			path_line2[loc_i].startx=path_line2[loc_i-1].startx;
			path_line2[loc_i].starty=path_line2[loc_i-1].starty;
		}

		if (score_diag == max_score) {
			path_line2[loc_i].tolast[DIAG]=TRUE;

			//Are we starting an alignment or continuing one?
			if (path_line1[loc_i-1].startx != -1) {
				path_line2[loc_i].startx=path_line1[loc_i-1].startx;
				path_line2[loc_i].starty=path_line1[loc_i-1].starty;
			}
			else {
				path_line2[loc_i].startx=loc_i;
				path_line2[loc_i].starty=loc_j;
			}
		}

		if (score_up == max_score) {
			path_line2[loc_i].tolast[UP]=TRUE;
			if(path_line1[loc_i].gap_extend_up == -1) {
				path_line2[loc_i].gap_extend_up=loc_j;
				path_line2[loc_i].up_score=max_score;
			}
			//We can't *start* a local alignment with a gap--this must be a continuation
			path_line2[loc_i].startx=path_line1[loc_i].startx;
			path_line2[loc_i].starty=path_line1[loc_i].starty;
		}
  
	}
	else {
		path_line2[loc_i].score=0;
		path_line2[loc_i].startx=path_line2[loc_i].starty=-1;
	}

	if (path_line2[loc_i].score >= align_max_score) {
		align_max_score=path_line2[loc_i].score;
		maxloc1=loc_i;
		maxloc2=loc_j;
		best_entry=path_line2[loc_i];
	}
  
}






void Local_Dynam_Prog::get_first_element()
{
	path_line2[0].tolast[UP]=path_line2[0].tolast[LEFT]=path_line2[0].tolast[DIAG]=FALSE;
	path_line2[0].score=0;
	path_line2[0].gap_extend_left=path_line2[0].gap_extend_up=-1;
	path_line2[0].startx=-1;
	path_line2[0].starty=-1;


}


Local_Dynam_Prog::~Local_Dynam_Prog()
{
  int i;

  
  if (path_line1 != 0 )
      delete[] path_line1;
  if (path_line2 != 0 )
      delete[] path_line2;
    

}



Align::Align ()
{
	cerr<<"ERROR: call to default Align constructor\n";
}


Align::Align (DATATYPE data, double extend, double open, ALIGN_TYPE type)
{
	//Top Level constructor--only call to create inital alignment sim matrix
	//(i.e. BLOSUM 62)
	int i,j;
	
	my_data=data;
	alignment_type=type;
	Align_obj[0]=Align_obj[1]=0;
	Local_DP_Obj=0;
	partial_align[0]=partial_align[1]=0;
	newseqs1[0]=newseqs1[1]=0;
	newseqs2[0]=newseqs2[1]=0;
	aligned_sequences=0;
	own_matrix=TRUE;
	
	if (data == NUCLEIC) { 
		convert_char=&readchar_to_base;  
		the_matrix = new Nuc_matrix();
	}
	else { 
		convert_char=&readchar_to_aa;
		the_matrix= new BLOSUM_62_matrix(); 
	}

	
	//Our alignment reconstruction is always global--local alignments are
	//performed with the Local_DP_Obj and then the sequences trimmed and
	//globally aligned
	DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, the_matrix);
	DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, the_matrix);
	
	if (alignment_type == LOCAL) 
		Local_DP_Obj=new Local_Dynam_Prog(data, extend, open, the_matrix);
}
	


Align::Align (DATATYPE data, char *matrix_loc, double extend, double open, ALIGN_TYPE type)
{
	//Top Level constructor--only call to create inital alignment sim matrix
	//(i.e. BLOSUM 62)
	int i,j;
	
	my_data=data;
	alignment_type=type;
	Align_obj[0]=Align_obj[1]=0;
	Local_DP_Obj=0;
	partial_align[0]=partial_align[1]=0;
	newseqs1[0]=newseqs1[1]=0;
	newseqs2[0]=newseqs2[1]=0;
	aligned_sequences=0;
	own_matrix=TRUE;

	if (data == NUCLEIC) 
		cerr<<"Error: File-based matrices only valild with AA matrices\n";
		//    convert_char=&readchar_to_base;  
	else {
		convert_char=&readchar_to_aa;
		the_matrix= new File_AA_matrix(matrix_loc);
    }

	if (data == NUCLEIC) 
		convert_char=&readchar_to_base;  
	else 
		convert_char=&readchar_to_aa;

	//Our alignment reconstruction is always global--local alignments are
	//performed with the Local_DP_Obj and then the sequences trimmed and
	//globally aligned
	DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, the_matrix);
	DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, the_matrix);

	if (alignment_type == LOCAL) 
		Local_DP_Obj=new Local_Dynam_Prog(data, extend, open, the_matrix);
	
}


Align::Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, double match, double mismatch)
{
	//Top Level constructor--only call to create inital alignment sim matrix
	//with match and mismatch penalties
	int i,j;
	
	my_data=data;
	alignment_type=type;
	Align_obj[0]=Align_obj[1]=0;
	Local_DP_Obj=0;
	partial_align[0]=partial_align[1]=0;
	newseqs1[0]=newseqs1[1]=0;
	newseqs2[0]=newseqs2[1]=0;
	aligned_sequences=0;
	own_matrix=TRUE;

	if (data == NUCLEIC) { 
		convert_char=&readchar_to_base; 
		the_matrix = new Nuc_matrix(match, mismatch);
	}
	else { 
		convert_char=&readchar_to_aa;
		the_matrix= new Simple_AA_matrix(match, mismatch);
	}
 
	

	//Our alignment reconstruction is always global--local alignments are
	//performed with the Local_DP_Obj and then the sequences trimmed and
	//globally aligned
	DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, the_matrix);
	DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, the_matrix);

	if (alignment_type == LOCAL) 
		Local_DP_Obj=new Local_Dynam_Prog(data, extend, open, the_matrix);


}

Align::Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, Score_matrix *mat)
{
	//Internal constructor--use to avoid recreating sim_matrix from file

	int i,j, mat_size;
	
	my_data=data;
	alignment_type=type;
	Align_obj[0]=Align_obj[1]=0;
	Local_DP_Obj=0;
	partial_align[0]=partial_align[1]=0;
	newseqs1[0]=newseqs1[1]=0;
	newseqs2[0]=newseqs2[1]=0;
	aligned_sequences=0;
	the_matrix=mat;
	own_matrix=FALSE;


	if (data == NUCLEIC) {
		convert_char=&readchar_to_base; 
	}
	else {
		convert_char=&readchar_to_aa;
	}


	//Our alignment reconstruction is always global--local alignments are
	//performed with the Local_DP_Obj and then the sequences trimmed and
	//globally aligned.  Since this is an internal constructor, we don't
	//need to worry about the possiblity of local alignments.
	DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, mat);
	DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, mat);
		

}



Align::Align (DATATYPE data, double extend, double open, ALIGN_TYPE type, Score_matrix *mat, BOOL internal)
{
	//Internal constructor--use to avoid recreating sim_matrix from file
	
	int i,j, mat_size;
	

	my_data=data;
	alignment_type=type;
	Align_obj[0]=Align_obj[1]=0;
	Local_DP_Obj=0;
	partial_align[0]=partial_align[1]=0;
	newseqs1[0]=newseqs1[1]=0;
	newseqs2[0]=newseqs2[1]=0;
	aligned_sequences=0;
	if (internal == TRUE) {
		the_matrix=mat;
		own_matrix=FALSE;
	
	
		if (data == NUCLEIC) {
			convert_char=&readchar_to_base; 
		}
		else {
			convert_char=&readchar_to_aa;
		}
	
	
		//Our alignment reconstruction is always global--local alignments are
		//performed with the Local_DP_Obj and then the sequences trimmed and
		//globally aligned.  Since this is an internal constructor, we don't
		//need to worry about the possiblity of local alignments.
		DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, mat);
		DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, mat);
	
	}
	else {
		own_matrix=TRUE;
		the_matrix=mat;
		
		if (data == NUCLEIC) 
			convert_char=&readchar_to_base; 		
		else  
			convert_char=&readchar_to_aa;		
		
		//Our alignment reconstruction is always global--local alignments are
		//performed with the Local_DP_Obj and then the sequences trimmed and
		//globally aligned
		DP_Obj[0]=new Global_Dynam_Prog(data, extend, open, the_matrix);
		DP_Obj[1]=new Global_Dynam_Prog(data, extend, open, the_matrix);
		
		if (alignment_type == LOCAL) 
			Local_DP_Obj=new Local_Dynam_Prog(data, extend, open, the_matrix);
		
	}
}





void Align::recurse_alignment(Molecule_Sequence *seq1, Molecule_Sequence *seq2, BOOL toplevel)
{
	int i, half_len, other_half_len, one_len, best_score_pos,
		new_size1, new_size2;
	double best_score, new_score, left_partial_score, up_partial_score;
	BOOL done=FALSE;
	Molecule_Sequence *oldseq1, *oldseq2, *pass_seq1, *pass_seq2;


	sequence1=seq1;
	sequence2=seq2;

    //if (sequence1->Sequence_size() <1) {
    //    cerr<<"Call to recurse_alignment with 0 sequence: "<<sequence1->Sequence_size()<<endl;
    //}

	if ((toplevel == TRUE) && (alignment_type == LOCAL)) {
		//We do a primary local alignment to determine where we
		//need to align
	
		Local_DP_Obj->do_alignment(sequence1, sequence2);
		if (Local_DP_Obj->get_max_score() >0 ) {
			starts[0]=Local_DP_Obj->best_entry.startx;
			starts[1]=Local_DP_Obj->best_entry.starty;

			ends[0]=Local_DP_Obj->get_maxloc_x();
			ends[1]=Local_DP_Obj->get_maxloc_y();
		
			oldseq1=sequence1;
			oldseq2=sequence2;

			new_size1=ends[0]-starts[0]+1;
			new_size2=ends[1]-starts[1]+1;


			sequence1= new Molecule_Sequence(new_size1);
			sequence2= new Molecule_Sequence(new_size2);

			for(i=0; i<new_size1; i++)
				sequence1->Assign_site(i, (*oldseq1)[i+starts[0]-1]);


			for(i=0; i<new_size2; i++)
				sequence2->Assign_site(i, (*oldseq2)[i+starts[1]-1]);
		
			sequence1->Assign_name(oldseq1->Sequence_name());
			sequence2->Assign_name(oldseq2->Sequence_name());
#if 0
			cout<<"Original size 1: "<<oldseq1->Sequence_size()<<" And 2: "<<oldseq2->Sequence_size()<<endl;
			cout<<"Local alignment ends at "<<ends[0]<<", "<<ends[1]<<endl;
			cout<<"First: "<<num_to_aa(*sequence1)[0])<<" 2: "<<num_to_aa((*sequence2)[0])<<endl;
			cout<<"Last 1-1: "<<num_to_aa((*sequence1)[sequence1->Sequence_size()-2])
				<<" Last 2-1: "<<num_to_aa((*sequence2)[sequence2)->Sequence_size()-2])<<endl;
			cout<<"Last 1: "<<num_to_aa((*sequence1)[sequence1->Sequence_size()-1])
				<<" Last 2: "<<num_to_aa((*sequence2)[sequence2->Sequence_size()-1])<<endl;
#endif
				Local_DP_Obj->clear_path_arrays();
		}
		else {
			cerr<<"ERROR: no local alignments\n";
			return;
		}

	}



	if (sequence1->Sequence_size() < 1)  {
		aligned_sequences=new Sequence_dataset(2, sequence2->Sequence_size());
	
		for(i=0; i<sequence2->Sequence_size(); i++) {
			(*aligned_sequences)[0].Assign_site(i, convert_char('-'));
			(*aligned_sequences)[1].Assign_site(i, (*sequence2)[i]);
		}
	
	}
	else if (sequence2->Sequence_size() < 1) {
		aligned_sequences=new Sequence_dataset(2, sequence1->Sequence_size());

		for(i=0; i<sequence1->Sequence_size(); i++) {
			(*aligned_sequences)[1].Assign_site(i, convert_char('-'));
			(*aligned_sequences)[0].Assign_site(i, (*sequence1)[i]);
		}
		
	}
	else if (sequence2->Sequence_size() < 2) {
		//We can finish this alignment without recursing
		newseqs1[0]=sequence1;
		newseqs2[0]=sequence2;
		DP_Obj[0]->do_alignment(newseqs1[0], newseqs2[0]);
		best_score=DP_Obj[0]->path_line2[newseqs1[0]->Sequence_size()].score;
		best_score_pos=newseqs1[0]->Sequence_size();

		done=recover_alignment(0,best_score_pos);
		
		DP_Obj[0]->clear_path_arrays();
		assemble_alignment();
		newseqs1[0]=newseqs2[0]=0;
		
	}
	else {

		half_len=sequence2->Sequence_size()/2;
		other_half_len=sequence2->Sequence_size()-half_len;
		one_len=sequence1->Sequence_size();


		newseqs2[0]=new Molecule_Sequence (half_len, my_data);
		newseqs2[1]=new Molecule_Sequence (other_half_len, my_data);
		newseqs1[1]=new Molecule_Sequence (one_len, my_data);
		newseqs1[0]=sequence1;

		for(i=0; i<one_len; i++) 
			newseqs1[1]->Assign_site(i, (*sequence1)[one_len-i-1]);

		for(i=0; i<half_len; i++) 
			newseqs2[0]->Assign_site(i, (*sequence2)[i]);

		for(i=0; i<other_half_len; i++) 
			newseqs2[1]->Assign_site(i, (*sequence2)[sequence2->Sequence_size()-i-1]);
		
		DP_Obj[0]->do_alignment(newseqs1[0], newseqs2[0]);
		DP_Obj[1]->do_alignment(newseqs1[1], newseqs2[1]);

		

		best_score=DP_Obj[0]->path_line2[0].score+DP_Obj[1]->path_line2[sequence1->Sequence_size()].score;

		if ((DP_Obj[0]->path_line2[0].tolast[UP]==TRUE) && 
				(DP_Obj[1]->path_line2[sequence1->Sequence_size()].tolast[UP]==TRUE)) 
				best_score -= DP_Obj[0]->get_gap_open();

		

		up_partial_score=DP_Obj[0]->path_line2[0].up_score+
			DP_Obj[1]->path_line2[sequence1->Sequence_size()].up_score-DP_Obj[1]->get_gap_open();

		if (up_partial_score > best_score) {
			best_score=up_partial_score;
			DP_Obj[1]->path_line2[sequence1->Sequence_size()].up_score-=DP_Obj[1]->get_gap_open();
			if (DP_Obj[1]->path_line2[sequence1->Sequence_size()].up_score >=
				DP_Obj[1]->path_line2[sequence1->Sequence_size()].score) {
				DP_Obj[1]->path_line2[sequence1->Sequence_size()].tolast[UP]=TRUE;
				DP_Obj[0]->path_line2[0].tolast[UP]=TRUE;
				if (DP_Obj[1]->path_line2[sequence1->Sequence_size()].up_score >
				DP_Obj[1]->path_line2[sequence1->Sequence_size()].score) {
 					DP_Obj[1]->path_line2[sequence1->Sequence_size()].score=
						DP_Obj[1]->path_line2[sequence1->Sequence_size()].up_score;
				
					DP_Obj[1]->path_line2[sequence1->Sequence_size()].tolast[DIAG]=
						DP_Obj[1]->path_line2[sequence1->Sequence_size()].tolast[LEFT]=FALSE;
					DP_Obj[0]->path_line2[0].tolast[DIAG]=DP_Obj[0]->path_line2[0].tolast[LEFT]=FALSE;
				}
				
				
			}
			
		}	
	

		best_score_pos=0;


		for(i=1; i<sequence1->Sequence_size()+1; i++) {
			new_score=DP_Obj[0]->path_line2[i].score+
					DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score;

			if ((DP_Obj[0]->path_line2[i].tolast[UP]==TRUE) && 
				(DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[UP]==TRUE)) 
				new_score -= DP_Obj[0]->get_gap_open();

			if ((DP_Obj[0]->path_line2[i].tolast[LEFT]==TRUE) && 
				(DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[LEFT]==TRUE)) 
				new_score -= DP_Obj[0]->get_gap_open();

			left_partial_score=DP_Obj[0]->path_line2[i].left_score+
			DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].left_score-DP_Obj[1]->get_gap_open();

			up_partial_score=DP_Obj[0]->path_line2[i].up_score+
			DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].up_score-DP_Obj[1]->get_gap_open();

			if (up_partial_score > new_score) {
				new_score=up_partial_score;
				


				DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].up_score-=DP_Obj[1]->get_gap_open();
				if (DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].up_score >=
					DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score) {
					DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[UP]=TRUE;
					DP_Obj[0]->path_line2[i].tolast[UP]=TRUE;
					if (DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].up_score >
						DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score) {
 							DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score=
								DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].up_score;
				
							DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[DIAG]=
								DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[LEFT]=FALSE;
							DP_Obj[0]->path_line2[i].tolast[DIAG]=DP_Obj[0]->path_line2[i].tolast[LEFT]=FALSE;
					}
						
				}
			
			}	
			else if (left_partial_score > new_score) {
				if (sequence1->Sequence_size()-i != 0) {
					new_score=left_partial_score;
					DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].left_score-=DP_Obj[1]->get_gap_open();

					if (DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].left_score >=
						DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score) {
						DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[LEFT]=TRUE;
						DP_Obj[0]->path_line2[i].tolast[LEFT]=TRUE;
						if (DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].left_score >
							DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score) {
 								DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].score=
									DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].left_score;
				
								DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[DIAG]=
									DP_Obj[1]->path_line2[sequence1->Sequence_size()-i].tolast[UP]=FALSE;
								DP_Obj[0]->path_line2[i].tolast[DIAG]=DP_Obj[0]->path_line2[i].tolast[UP]=FALSE;
						}
						
					}
			
				}
			}	




			if (new_score >=best_score) {
				best_score=new_score;
				best_score_pos=i;

			}

		}

		if (toplevel == TRUE) {
			max_score=best_score;
#if 0
			cout<<"Max score found at "<<best_score_pos<<" : "
				<<max_score<<" = "<<DP_Obj[0]->path_line2[best_score_pos].score<<" + "
				<<DP_Obj[1]->path_line2[sequence1->Sequence_size()-best_score_pos].score<<endl;
#endif
		}

		done=recover_alignment(0,best_score_pos);
		done=recover_alignment(1,newseqs1[1]->Sequence_size()-best_score_pos);
		
		
		if (done == FALSE)
		{
			//Recursively align top
			if ((newseqs2[0]->Sequence_size()-seq2used[0] > 0) ||
				(best_score_pos-seq1used[0] > 0)) {
				new_size1=best_score_pos-seq1used[0];

			
				pass_seq1=new Molecule_Sequence (new_size1);
				for(i=0; i<new_size1; i++)
					pass_seq1->Assign_site(i, (*sequence1)[i]);
			
			
				new_size2=newseqs2[0]->Sequence_size()-seq2used[0];

				pass_seq2=new Molecule_Sequence (new_size2);
				for(i=0; i<new_size2; i++)
					pass_seq2->Assign_site(i, (*newseqs2[0])[i]);
	
			
				
				Align_obj[0]= new Align(DP_Obj[0]->get_datatype(), DP_Obj[0]->get_gap_extend(), 
					DP_Obj[0]->get_gap_open(), alignment_type, the_matrix);

				delete DP_Obj[0];
				DP_Obj[0]=0;

				Align_obj[0]->recurse_alignment(pass_seq1, pass_seq2, FALSE);
				delete pass_seq1;
				delete pass_seq2;
				pass_seq1=pass_seq2=0;
				
			}

			//Recursively align bottom
			if ((newseqs2[1]->Sequence_size()-seq2used[1] > 0) ||
				((newseqs1[1]->Sequence_size()-best_score_pos)-seq1used[1] > 0)) {
				new_size1=(newseqs1[1]->Sequence_size()-best_score_pos)-seq1used[1];

				pass_seq1=new Molecule_Sequence (new_size1);
				for(i=0; i<new_size1; i++)
					pass_seq1->Assign_site(i, (*newseqs1[1])[i]);

		
				new_size2=newseqs2[1]->Sequence_size()-seq2used[1];

				pass_seq2=new Molecule_Sequence (new_size2);
				for(i=0; i<new_size2; i++)
					pass_seq2->Assign_site(i, (*newseqs2[1])[i]);
	
			


			
				Align_obj[1]= new Align(DP_Obj[1]->get_datatype(), DP_Obj[1]->get_gap_extend(), 
					DP_Obj[1]->get_gap_open(), alignment_type, the_matrix);

				
				delete DP_Obj[1];
				DP_Obj[1]=0;

				Align_obj[1]->recurse_alignment(pass_seq1, pass_seq2, FALSE);
				delete pass_seq1;
				delete pass_seq2;
				pass_seq1=pass_seq2=0;
			}
			

		}
		delete newseqs2[0];
		newseqs2[0]=0;
		delete newseqs1[1];
		delete newseqs2[1];
		newseqs1[1]=newseqs2[1]=0;
		

		

		assemble_alignment();

		if (toplevel == TRUE) {
			(*aligned_sequences)[0].Assign_name(sequence1->Sequence_name());
			(*aligned_sequences)[1].Assign_name(sequence2->Sequence_name());
			if (alignment_type == LOCAL) {
				delete sequence1;
				delete sequence2;
			}

		}
		
		//if (Align_obj[0] != 0) delete Align_obj[0];
		//if (Align_obj[1] != 0) delete Align_obj[1];
		//Align_obj[0]=Align_obj[1]=0;
	}
	


}



BOOL Align::recover_alignment(int number, int pos)
{
	int i, best_score_pos, ypos, val;
	BOOL done = FALSE;
	List<int> *list1, *list2;


	list1 = new List<int>;
	list2 = new List<int>;

	seq1used[number]=0;
	seq2used[number]=0;

	best_score_pos=pos;
	ypos=newseqs2[number]->Sequence_size()-1;

	if (DP_Obj[number]->path_line2[best_score_pos].tolast[UP] == TRUE) {
	
		while (ypos+1 >= DP_Obj[number]->path_line2[best_score_pos].gap_extend_up) {
			val=convert_char('-');
			list1->add_to_list(val);
			val=(*newseqs2[number])[ypos--];
			list2->add_to_list(val);
			seq2used[number]++;
		}
	}
	else if (DP_Obj[number]->path_line2[best_score_pos].tolast[DIAG] == TRUE) {	
		val=(*newseqs1[number])[best_score_pos-1];
	    list1->add_to_list(val);
		seq1used[number]++;
		val=(*newseqs2[number])[ypos--];
		list2->add_to_list(val);
		seq2used[number]++;
		best_score_pos--;
	}
	else if (DP_Obj[number]->path_line2[best_score_pos].tolast[LEFT] == TRUE) {
			i=best_score_pos;
			while(i>=DP_Obj[number]->path_line2[best_score_pos].gap_extend_left) {
				val=(*newseqs1[number])[i-1];
				list1->add_to_list(val);
				seq1used[number]++;
				val=convert_char('-');
				list2->add_to_list(val);
				i--;
			}
			if (DP_Obj[number]->path_line2[i].tolast[UP] ==TRUE) {	
				while (ypos+1 >=DP_Obj[number]->path_line2[best_score_pos].gap_extend_up) {
					val=convert_char('-');
					list1->add_to_list(val);
					val=(*newseqs2[number])[ypos--];
					list2->add_to_list(val);
					seq2used[number]++;
				}
			}
			 else if (DP_Obj[number]->path_line2[i].tolast[DIAG] ==TRUE) {
				 val=(*newseqs1[number])[i-1];
				list1->add_to_list(val);
				 val=(*newseqs2[number])[ypos--];
				list2->add_to_list(val);
				seq1used[number]++;
				seq2used[number]++;
			}
			else
				//We've finished this alignment
				done = TRUE;
			best_score_pos=i-1;
				
	}
	
	
	if ((done == FALSE) && (seq2used[number] < 2))  {
		//if we've not finished and still have data for the second sequence
			if (DP_Obj[number]->path_line1[best_score_pos].tolast[UP] == TRUE) {
				while (ypos+1 >=DP_Obj[number]->path_line1[best_score_pos].gap_extend_up) {
					val=convert_char('-');
					list1->add_to_list(val);
					val=(*newseqs2[number])[ypos--];
					list2->add_to_list(val);
					seq2used[number]++;
					
				}
			}

			else if (DP_Obj[number]->path_line1[best_score_pos].tolast[DIAG] == TRUE) {	
				val=(*newseqs1[number])[best_score_pos-1];
			    list1->add_to_list(val);
				seq1used[number]++;
				val=(*newseqs2[number])[ypos--];
				list2->add_to_list(val);
				seq2used[number]++;
				best_score_pos--;
			}
			else if (DP_Obj[number]->path_line1[best_score_pos].tolast[LEFT] == TRUE) {
					i=best_score_pos;
					while(i>=DP_Obj[number]->path_line1[best_score_pos].gap_extend_left) {
						val=(*newseqs1[number])[i-1];
						list1->add_to_list(val);
						seq1used[number]++;
						val=convert_char('-');
						list2->add_to_list(val);
						i--;
					}
					if (DP_Obj[number]->path_line1[i].tolast[UP] ==TRUE) {								
						while (ypos+1 >=DP_Obj[number]->path_line2[best_score_pos].gap_extend_up) {
							val=convert_char('-');
							list1->add_to_list(val);
							val=(*newseqs2[number])[ypos--];
							list2->add_to_list(val);
							seq2used[number]++;
						}
					}
					else if (DP_Obj[number]->path_line1[i].tolast[DIAG] ==TRUE) {
						val=(*newseqs1[number])[i-1];
						list1->add_to_list(val);
						val=(*newseqs2[number])[ypos--];
						list2->add_to_list(val);
						seq1used[number]++;
						seq2used[number]++;
					}
					else 
						//We've finished this alignment
						done = TRUE;
					best_score_pos=i-1;
				
			}
			else {
				//For local alignments which may end at this point
				done=TRUE;	
			}


	}
		
	
	partial_align[number] = new Sequence_dataset(2, list1->get_list_length());

	list1->return_to_start();
	list2->return_to_start();
     
 
	for(i=list1->get_list_length()-1; i>=0; i--)
	{
		(*partial_align[number])[0].Assign_site(i, *list1->get_current_item());
		(*partial_align[number])[1].Assign_site(i, *list2->get_current_item());
		list1->get_next();
		list2->get_next();
	}

	delete list1;
	delete list2;
#if 0
		cout<<"Partial Alignment "<<number<<endl;
		for(i=0; i<(*partial_align[number])[0].Sequence_size(); i++)
			cout<<num_to_aa((*partial_align[number])[0][i]);
		cout<<endl;

		for(i=0; i<(*partial_align[number])[0].Sequence_size(); i++)
			cout<<num_to_aa((*partial_align[number])[1][i]);
		cout<<endl;
#endif
	return(done);
}




void Align::assemble_alignment()
{
	int i, j, size, offset, len;


	size=(*partial_align[0])[0].Sequence_size();
	if (partial_align[1] != 0)
		size+=(*partial_align[1])[0].Sequence_size();


	if (Align_obj[0] != 0 )
		size += (*Align_obj[0]->aligned_sequences)[0].Sequence_size();

	if (Align_obj[1] != 0 )
		size += (*Align_obj[1]->aligned_sequences)[0].Sequence_size();


	aligned_sequences= new Sequence_dataset(2, size, my_data);

	offset=0;

	if (Align_obj[0] != 0 ) {
		for(i=0; i<(*Align_obj[0]->aligned_sequences)[0].Sequence_size(); i++) {
			(*aligned_sequences)[0].Assign_site(i, (*Align_obj[0]->aligned_sequences)[0][i]);
			(*aligned_sequences)[1].Assign_site(i, (*Align_obj[0]->aligned_sequences)[1][i]);

		}
		offset +=(*Align_obj[0]->aligned_sequences)[0].Sequence_size();
	}
	
	
 
	for(i=0; i<(*partial_align[0])[0].Sequence_size(); i++)
       {
		(*aligned_sequences)[0].Assign_site(i+offset, (*partial_align[0])[0][i]);
		(*aligned_sequences)[1].Assign_site(i+offset, (*partial_align[0])[1][i]);
       }
	offset+=(*partial_align[0])[0].Sequence_size();
	delete partial_align[0];
	partial_align[0]=0;


	if (partial_align[1] != 0) {
		len=(*partial_align[1])[0].Sequence_size();
		for(i=0; i<len; i++)
		{	
			(*aligned_sequences)[0].Assign_site(i+offset, (*partial_align[1])[0][len-i-1]);
			(*aligned_sequences)[1].Assign_site(i+offset, (*partial_align[1])[1][len-i-1]);
		}
		offset+=(*partial_align[1])[0].Sequence_size();
		delete partial_align[1];
		partial_align[1]=0;
	}



    if (Align_obj[1] != 0 ) {
		len=(*Align_obj[1]->aligned_sequences)[0].Sequence_size();
		for(i=0; i<len; i++) {
			(*aligned_sequences)[0].Assign_site(i+offset, (*Align_obj[1]->aligned_sequences)[0][len-i-1]);
			(*aligned_sequences)[1].Assign_site(i+offset, (*Align_obj[1]->aligned_sequences)[1][len-i-1]);

		}
	}

}



Align::~Align()
{
	if (aligned_sequences != 0)
		delete aligned_sequences;
	
#if 1
	if(	newseqs1[1] !=0) delete newseqs1[1];
	if(newseqs2[0] != 0) delete newseqs2[0];
	if(newseqs2[1] != 0) delete newseqs2[1];
	if(partial_align[0] != 0) delete partial_align[0];
	if(partial_align[1] != 0) delete partial_align[1];
#endif
	
	
	if (Align_obj[0] !=0)
		delete Align_obj[0];

	if (Align_obj[1] !=0)
		delete Align_obj[1];

	if (DP_Obj[0] !=0)
		delete DP_Obj[0];

	if (DP_Obj[1] !=0)
		delete DP_Obj[1];
	
	if(Local_DP_Obj != 0)
		delete Local_DP_Obj;

	if ((own_matrix == TRUE) && (the_matrix != 0)) delete the_matrix;
	
}




