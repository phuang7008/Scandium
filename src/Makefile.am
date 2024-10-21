bin_PROGRAMS=scandium
scandium_SOURCES=annotation.c main.c stats.c user_defined_annotation.c user_inputs.c utils.c reports.c targets.c annotation.h stats.h user_defined_annotation.h user_inputs.h utils.h reports.h targets.h terms.h

scandium_LDFLAGS = -lm -lz -lpthread -lmysqlclient -lcrypto -ldl -llzma -lbz2 -lcurl

scandium_CFLAGS = -Wall -g -fopenmp
scandium_CPPFLAGS=-I$(includedir) -I$(abs_top_builddir)/htslib -I$(abs_top_builddir)/lib
scandium_LDADD = $(abs_top_builddir)/htslib/lib/libhts.a $(abs_top_builddir)/lib/libcov_tracking.a
