libsaftdir = $(top_srcdir)/src/libsaft

INCLUDES =				\
	-I$(libsaftdir)

noinst_PROGRAMS =			\
	test_fasta			\
	test_search2seqs		\
	test_mean_var			\
	test_pgamma			\
	test_BH

test_fasta_LDADD = $(libsaftdir)/libsaft.la
test_fasta_SOURCES = test_fasta.c

test_search2seqs_LDADD = $(libsaftdir)/libsaft.la
test_search2seqs_SOURCES = test_search2seqs.c

test_mean_var_LDADD = $(libsaftdir)/libsaft.la
test_mean_var_SOURCES = test_mean_var.c

test_pgamma_LDADD = $(libsaftdir)/libsaft.la
test_pgamma_SOURCES = test_pgamma.c

test_BH_LDADD = $(libsaftdir)/libsaft.la
test_BH_SOURCES = test_BH.c

MAINTAINERCLEANFILES =			\
	Makefile.in
