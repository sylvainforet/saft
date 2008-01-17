libbfastdir = $(top_srcdir)/src/libbfast

INCLUDES =				\
	-I$(libbfastdir)

noinst_PROGRAMS =			\
	test_fasta			\
	test_search2seqs

test_fasta_LDADD =			\
	$(libbfastdir)/libbfast.la
test_fasta_SOURCES =			\
	test_fasta.c

test_search2seqs_LDADD =		\
	$(libbfastdir)/libbfast.la
test_search2seqs_SOURCES =		\
	test_search2seqs.c

MAINTAINERCLEANFILES =			\
	Makefile.in
