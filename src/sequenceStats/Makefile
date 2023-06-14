IDIR=/opt/rh/devtoolset-8/root/usr/include
CODEDIR=/stornext/snfs130/NGIRD/scratch/phuang/repo/scandium/src/sequenceStats
HTSDIR=/stornext/snfs130/NGIRD/scratch/phuang/
MYSQLDIR=/usr/include/mysql
MYSQLLIB=/usr/lib64
MYLIB=/stornext/snfs130/NGIRD/scratch/phuang/repo/common_libs
CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
#CC=/hgsc_software/gcc/gcc-7.3.0/bin/gcc
#CFLAGS= -g -Wall -O2 -fopenmp
CFLAGS= -g -Wall -Wextra -O0 -fopenmp
#CFLAGS= -g3 -gdwarf-2 -fvar-tracking -fvar-tracking-assignments -Wall -O0 -fopenmp

ifeq ($(shell hostname), analysis-app1.sug.hgsc.bcm.edu)
	CC=/usr/bin/gcc
	MYSQLLIB=/usr/lib64/mysql
else
	CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
endif

# define any directories containing header files other than /usr/bin
include_dirs= -I$(IDIR) -I$(HTSDIR) -I$(MYSQLDIR) -I$(MYLIB)

# define library paths in addition to /usr/lib
LFLAGS = -L$(HTSDIR) -L$(MYSQLLIB) -L$(MYLIB)

# define any libraries to be link into executables:
#LIBS =  $(HTSDIR)/libhts.a -lm -lz -lpthread -lmysqlclient -lrt -lssl -lcrypto -ldl
#LIBS =  $(HTSDIR)/libhts.a $(MYLIB)/libcov_tracking.a -lm -lz -lpthread -lmysqlclient_r -lrt -lssl -lcrypto -ldl -lcurl -llzma -lbz2
ifeq ($(shell hostname), analysis-app1.sug.hgsc.bcm.edu)
	LIBS =  $(HTSDIR)/htslib/lib/libhts.a $(MYLIB)/libcov_tracking.a -lm -lz -lpthread -lmysqlclient -lrt -lssl -lcrypto -ldl -lcurl -llzma -lbz2
else
	LIBS =  $(HTSDIR)/libhts.a -lm -lz -lpthread -lmysqlclient_r -lrt -lssl -lcrypto -ldl -lcurl -llzma -lbz2
endif

# define dependent header files
#_DEPS = stats.h  targets.h  terms.h  utils.h annotation.h reports.h user_defined_annotation.h user_inputs.h
_DEPS = $(wildcard *.h)
DEPS=$(patsubst %,$(CODEDIR)/%,$(_DEPS))

# define the C ssource files
#SRCS = main.c stats.c targets.c utils.c annotation.c reports.c user_defined_annotation.c user_inputs.c
SRCS = $(wildcard *.c)

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
MAIN = scandium

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:	$(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(include_dirs) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
#$(CODEDIR)/%.o: %.c $(DEPS)
#.c.o:
#.c.o: $(DEPS)
%.o: %.c
	$(CC) $(CFLAGS) $(include_dirs) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(include_dirs) $^

# DO NOT DELETE THIS LINE -- make depend needs it
