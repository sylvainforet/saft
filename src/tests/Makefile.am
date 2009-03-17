libwaftdir = $(top_srcdir)/src/libwaft

INCLUDES =				\
	-I$(libwaftdir)

noinst_PROGRAMS =			\
	test_fasta			\
	test_search2seqs

test_fasta_LDADD =			\
	$(libwaftdir)/libwaft.la
test_fasta_SOURCES =			\
	test_fasta.c

test_search2seqs_LDADD =		\
	$(libwaftdir)/libwaft.la
test_search2seqs_SOURCES =		\
	test_search2seqs.c

MAINTAINERCLEANFILES =			\
	Makefile.in
