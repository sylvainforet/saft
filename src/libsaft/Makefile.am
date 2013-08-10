ACLOCAL_AMFLAGS =		\
	-I m4

AM_CPPFLAGS =			\
	$(GSL_CFLAGS)

lib_LTLIBRARIES =		\
	libsaft.la

libsaft_sources =		\
	safterror.c		\
	saftfasta.c		\
	safthash.c		\
	saftsearch.c		\
	saftsequence.c		\
	saftstats.c

libsaft_headers =		\
	safterror.h		\
	saftfasta.h		\
	safthash.h		\
	saftsearch.h		\
	saftsequence.h		\
	saftstats.h

# FIXME add library version info like that:
#libsaft_la_LDFLAGS =	\
#	-version-info $(LT_VERSION_INFO)
#
libsaft_la_LIBADD =		\
	$(GSL_LIBS)

libsaft_la_SOURCES =		\
	$(libsaft_sources)	\
	$(libsaft_headers)

libsaftincludedir =		\
	$(includedir)/saft

libsaftinclude_HEADERS =	\
	$(libsaft_headers)

MAINTAINERCLEANFILES =		\
	Makefile.in
