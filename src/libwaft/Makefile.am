lib_LTLIBRARIES =		\
	libwaft.la

libwaft_sources =		\
	wafterror.c		\
	waftfasta.c		\
	wafthash.c		\
	waftsearch.c		\
	waftsequence.c

libwaft_headers =		\
	wafterror.h		\
	waftfasta.h		\
	wafthash.h		\
	waftsearch.h		\
	waftsequence.h

# TODO add library version info like that:
#libwaft_la_LDFLAGS =	\
#	-version-info $(LT_VERSION_INFO)

libwaft_la_SOURCES =		\
	$(libwaft_sources)	\
	$(libwaft_headers)

libwaftincludedir =		\
	$(includedir)/waft

libwaftinclude_HEADERS =	\
	$(libwaft_headers)

MAINTAINERCLEANFILES =		\
	Makefile.in
