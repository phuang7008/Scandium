# (note: the <tab> in the command line is necessary for make to work)
# variable definitions
CC=g++
CXXFLAGS=-g -Wall -std=c++0x
CODEDIR=/stornext/snfs5/next-gen/scratch/phuang/git_repo/dev/batch_analysis/src
MYSQLDIR=-I/usr/include/mysql
MYSQLLIB=/usr/lib64

# define any libraries to be link into executables:
LIBS = -lmysqlclient

# define dependent header files
_DEPS = hts_data.hpp exon.hpp annotation.hpp user_inputs.hpp utils.hpp
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C++ source files
SRCS = annotation.cpp main.cpp user_inputs.cpp utils.cpp

# define the C++ object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
OBJS = $(SRCS:.c=.o)

# the build target executable:
MAIN = batch_analysis

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
.PHONY:depend clean

all: $(MAIN)
	@echo  batch_analysis executable has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CXXFLAGS) $(MYSQLDIR) -o $(MAIN) $(OBJS) $(LIBS) -L$(MYSQLLIB)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(MYSQLDIR) -c $<  -o $@

# clean up everything
clean:
	$(RM) *.o *~ $(MAIN)

# DO NOT DELETE THIS LINE -- make depend needs it
