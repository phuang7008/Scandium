# (note: the <tab> in the command line is necessary for make to work)
# variable definitions
CC=g++
CXXFLAGS=-g -Wall -std=c++0x
CODEDIR=/stornext/snfs5/next-gen/scratch/phuang/dev/c_graphing/src

# define dependent header files
_DEPS = user_inputs.hpp utils.hpp
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C++ source files
SRCS = graphing.cpp user_inputs.cpp utils.cpp

# define the C++ object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
OBJS = $(SRCS:.c=.o)

# the build target executable:
MAIN = draw_graph

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
.PHONY:depend clean

all: $(MAIN)
	@echo  draw_graph executable has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CXXFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) -c $<  -o $@

# clean up everything
clean:
	$(RM) *.o *~ $(MAIN)

# DO NOT DELETE THIS LINE -- make depend needs it
