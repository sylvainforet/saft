libsaftdir = $(top_srcdir)/src/libsaft

INCLUDES =				\
	-I$(libsaftdir)

noinst_PROGRAMS =			\
	test_fasta			\
	test_search2seqs		\
	test_stats

test_fasta_LDADD =			\
	$(libsaftdir)/libsaft.la
test_fasta_SOURCES =			\
	test_fasta.c

test_search2seqs_LDADD =		\
	$(libsaftdir)/libsaft.la
test_search2seqs_SOURCES =		\
	test_search2seqs.c

test_stats_LDADD =			\
	$(libsaftdir)/libsaft.la
test_stats_SOURCES =			\
	test_stats.c

MAINTAINERCLEANFILES =			\
	Makefile.in
