IDIR=/usr/lib/gcc/x86_64-redhat-linux/4.4.7/include
CODEDIR=/stornext/snfs5/next-gen/scratch/phuang/dev/sequenceStats
HTSDIR=$(CODEDIR)/htslib
MYSQLDIR=/usr/include/mysql
MYSQLLIB=/usr/lib64
CC=gcc
#CFLAGS= -g -Wall -O2 -fopenmp
CFLAGS= -g -Wall -O0 -fopenmp
#CFLAGS= -g3 -gdwarf-2 -fvar-tracking -fvar-tracking-assignments -Wall -O0 -fopenmp

# define any directories containing header files other than /usr/bin
INCLUDES= -I$(IDIR) -I$(HTSDIR) -I$(MYSQLDIR)

# define library paths in addition to /usr/lib
LFLAGS = -L$(HTSDIR) -L$(MYSQLLIB)

# define any libraries to be link into executables:
#LIBS =  $(HTSDIR)/libhts.a -lm -lz -lpthread -lmysqlclient -lrt -lssl -lcrypto -ldl
LIBS =  $(HTSDIR)/libhts.a -lm -lz -lpthread -lmysqlclient_r -lrt -lssl -lcrypto -ldl

# define dependent header files
_DEPS = stats.h  targets.h  terms.h  utils.h for_mysql.h reports.h
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C ssource files
SRCS = main.c stats.c targets.c  utils.c for_mysql.c reports.c

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = coverage

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:	$(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
#$(CODEDIR)/%.o: %.c $(DEPS)
#.c.o:
.c.o: $(DEPS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
