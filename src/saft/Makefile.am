ACLOCAL_AMFLAGS =		\
	-I m4
libsaftdir = $(top_srcdir)/src/libsaft

AM_CPPFLAGS =			\
	-I$(libsaftdir)		\
	-DSAFT_VERSION=\"$(PACKAGE_VERSION)\"

bin_PROGRAMS =	saft

saft_LDADD =			\
	$(libsaftdir)/libsaft.la

saft_SOURCES =			\
	saftmain.c

MAINTAINERCLEANFILES =		\
	Makefile.in
