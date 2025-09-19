# Master Makefile for Sequence Utilities for Genome analysis
# Type make to build

O = o
SRC_DIR = ./

cc = gcc
CC = g++

INCLUDE = -I. 
CFLAGS = -g                 		                     # compiler switches to be applied to every module
OPTIM_SPEED = -O3             	  	                     # switches that give speed priority over size
OPTIM_SIZE = -O1              	  	                     # switches that give size priority over speed

OPTIONS = $(CFLAGS) $(INCLUDE)


all:  find_max_cover 



FIND_MAX_COVER_OBJS = find_max_cover.$(O)	\
		      read_seq.$(O)    \
			score_matrix.$(O)   \
	              gen_dna_funcs.$(O)    


find_max_cover : $(FIND_MAX_COVER_OBJS) 
	$(CC) $(LINUX_BUILD)  $(LIBRARY_DIR) \
	 -o find_max_cover $(OPTIONS) $(FIND_MAX_COVER_OBJS) 


%.o: %.cpp
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<







