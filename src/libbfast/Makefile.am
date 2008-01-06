lib_LTLIBRARIES =		\
	libbfast.la

libbfast_sources =		\
	bfasterror.c		\
	bfastfasta.c		\
	bfastsequence.c

libbfast_headers =		\
	bfasterror.h		\
	bfastfasta.h		\
	bfastsequence.h

# TODO add library version info like that:
#libbfast_la_LDFLAGS =	\
#	-version-info $(LT_VERSION_INFO)

libbfast_la_SOURCES =		\
	$(libbfast_sources)	\
	$(libbfast_headers)

libbfastincludedir =		\
	$(includedir)/bfast

libbfastinclude_HEADERS =	\
	$(libbfast_headers)

MAINTAINERCLEANFILES =		\
	Makefile.in
