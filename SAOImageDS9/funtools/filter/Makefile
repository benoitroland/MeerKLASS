#
# This file is a Makefile for Filter.  If it has the name "Makefile.in"
# then it is a template for a Makefile;  to generate the actual Makefile,
# run "./configure", which is a configuration script generated by the
# "autoconf" program (constructs like "@foo@" will get replaced in the
# actual Makefile.
#

# Current FILT version;  used in various names.
PACKAGE = filter
VERSION = 1.1

DISTNAME = filter-${VERSION}
DISTDIR = ../export/${DISTNAME}
FTPDIR	= ../ftp

#----------------------------------------------------------------
# Things you can change to personalize the Makefile for your own
# site (you can make these changes in either Makefile.in or
# Makefile, but changes to Makefile will get lost if you re-run
# the configuration script).
#----------------------------------------------------------------

# Default top-level directories in which to install architecture-
# specific files (exec_prefix) and machine-independent files such
# as scripts (prefix).  The values specified here may be overridden
# at configure-time with the --exec-prefix and --prefix options
# to the "configure" script.

prefix =		/usr/local/src/SAOImageDS9
exec_prefix =		/usr/local/src/SAOImageDS9

# The following definition can be set to non-null for special systems
# like AFS with replication.  It allows the pathnames used for installation
# to be different than those used for actually reference files at
# run-time.  INSTALL_ROOT is prepended to $prefix and $exec_prefix
# when installing files.
INSTALL_ROOT =

# Directory in which to install the .a, .so, and .o files:
LIB_INSTALL_DIR =	$(INSTALL_ROOT)$(exec_prefix)/lib

# Directory in which to install the programs:
BIN_INSTALL_DIR =	$(INSTALL_ROOT)$(exec_prefix)/bin

# Directory in which to install the include files:
INCLUDE_INSTALL_DIR =	$(INSTALL_ROOT)$(prefix)/include

# util files are in the util directory at same level
UTIL_INC =		-I../util
#UTIL_LIBS =		-L../util -lutil

# fitsy files are in the fitsy directory at same level
FITSY_INC =		-I../fitsy
#FITSY_LIBS =		-L../fitsy

# wcs files are in the wcs subdirectory
WCS_INC =		-I../wcs
#WCS_LIBS =		-L../wcs -lwcs

# extra includes for compiling
INCLUDES =		$(UTIL_INC) $(FITSY_INC) $(WCS_INC) 

# extra libs
EXTRA_LIBS =		-L.. -lfuntools 

# the full set of libraries for linking
LIBS =			$(WCS_LIBS) $(FITSY_LIBS) $(EXTRA_LIBS) -lm 

# To change the compiler switches, for example to change from -O
# to -g, change the following line:
CFLAGS = 		-g -O2

# To add ld switches, change the following line:
LDFLAGS =		

# whether to set the FPU in double round mode
FPU_DOUBLE =		0

# Some versions of make, like SGI's, use the following variable to
# determine which shell to use for executing commands:
SHELL =			/bin/sh

# extra modules added by configure.in to fix OS bugs
EXTRA_OBJS = 		

# special definitions for filtering
OBJPATH	=		\"$(LIB_INSTALL_DIR):.\"
FILTER_CC =		NULL
FILTER_CFLAGS =		NULL

# There are just too many different versions of "install" around;
# better to use the install-sh script that comes with the distribution,
# which is slower but guaranteed to work.

INSTALL = 		./install-sh -c
INSTALL_PROGRAM =	${INSTALL}
INSTALL_DATA =		${INSTALL} -m 644

#----------------------------------------------------------------
# The information below is modified by the configure script when
# Makefile is generated from Makefile.in.  You shouldn't normally
# modify any of this stuff by hand.
#----------------------------------------------------------------

AC_FLAGS =		-DHAVE_CONFIG_H
RANLIB =		ranlib

#----------------------------------------------------------------
# The information below should be usable as is.  The configure
# script won't modify it and you shouldn't need to modify it
# either.
#----------------------------------------------------------------

CC =		gcc

CC_SWITCHES =	-I. ${CFLAGS} ${CPPFLAGS} -DFILT_VERSION="\"$(VERSION)\""  \
		${INCLUDES} ${AC_FLAGS}

DEPEND_SWITCHES = -I. ${CFLAGS} -I. ${INCLUDES} ${AC_FLAGS}

SRCS =		evregions.c filter.c filtprog.c filtprog_c.c \
		imregions.c symbols.c lex.filt.c column.c dl.c tl.c \
		lex.idx.c idx.tab.c idxacts.c

OBJS =		evregions.o filter.o filtprog.o filtprog_c.o \
		imregions.o symbols.o lex.filt.o column.o dl.o tl.o \
		lex.idx.o idx.tab.o idxacts.o

# these are all the modules going into the "normal" filt library
LIBOBJS =	${OBJS}

# the default library for this package
DEFLIB =	libfilter.a
# the actual library we are building (if this is a subpackage)
LIB =		/usr/local/src/SAOImageDS9/funtools/libfuntools.a

TESTPROGS =	tfilt

all:		$(LIB) $(PROGS)

testall:	$(TESTPROGS)

All:		all testall

install:	install-binaries

$(LIB):		$(LIBOBJS)
		ar cruv $(LIB) $(LIBOBJS)
		$(RANLIB) $(LIB)

tfilt:		$(LIB) tfilt.o
		$(CC) $(LDFLAGS) tfilt.o -o tfilt $(LIB) ${LIBS}

# Note: before running ranlib below, must cd to target directory because
# some ranlibs write to current directory, and this might not always be
# possible (e.g. if installing as root).

install-binaries: $(LIB) $(PROGS)
	@for i in $(LIB_INSTALL_DIR) $(INCLUDE_INSTALL_DIR) $(BIN_INSTALL_DIR) ; \
	    do \
	    if [ ! -d $$i ] ; then \
		echo "Making directory $$i"; \
		mkdir $$i; \
		chmod 755 $$i; \
		else true; \
		fi; \
	    done;
	@if [ "$(LIB)" = "$(DEFLIB)" ] ; then \
	    XLIB=`basename $(LIB)`; \
	    echo "Installing $$XLIB"; \
	    $(INSTALL_DATA) $$XLIB $(LIB_INSTALL_DIR)/$$XLIB; \
	    (cd $(LIB_INSTALL_DIR); $(RANLIB) $$XLIB); \
	    chmod 555 $(LIB_INSTALL_DIR)/$$XLIB; \
        fi;
#	@for i in imregions.o evregions.o ; \
#	    do \
#		if [ -f $$i ] ; then \
#		  echo "Installing $$i" ; \
#		  $(INSTALL_DATA) $$i $(LIB_INSTALL_DIR)/$$i ; \
#		fi; \
#	    done;
	@for i in filter.h column.h dl.h tl.h idx.h; \
	    do \
		if [ -f $$i ] ; then \
		  echo "Installing $$i" ; \
		  $(INSTALL_DATA) $$i $(INCLUDE_INSTALL_DIR)/$$i ; \
		fi; \
	    done;

Makefile:	Makefile.in
		$(SHELL) config.status

RM = 		rm -f

clean:
		$(RM) *.a *.so *.o *.exe core errs *pure* .nfs* \
			foo* *~ *.log \#* TAGS *.E a.out errors \
			$(PROGS) $(TESTPROGS) \
			gmon.out *.pg *.bak \
			config.info config.log \
			doc/*~
		$(RM) -r autom4te.cache

distclean:	clean
		$(RM) Makefile config.status config.cache config.log

depend:
		makedepend -- $(DEPEND_SWITCHES) -- $(SRCS)

.c.o:
		$(CC) -c $(CC_SWITCHES) $<

filtprog_c.o:	filtprog_c.c regions_h.h xalloc_c.h swap_c.h \
		events_c.h image_c.h evregions_c.h imregions_c.h filter.h
		$(CC) -c $(CC_SWITCHES) -DOBJPATH="$(OBJPATH)" \
		-DFILTER_CC="$(FILTER_CC)" -DFILTER_CFLAGS="$(FILTER_CFLAGS)" \
		filtprog_c.c

# header and source include files for program generation
regions_h.h:	regions.h inc.sed
		$(RM) regions_h.h
		./inc.sed REGIONS_H < regions.h > regions_h.h

xalloc_c.h:	../util/xalloc.c inc.sed
		$(RM) xalloc_c.h
		./inc.sed XALLOC_C < ../util/xalloc.c > xalloc_c.h

swap_c.h:	swap.c inc.sed
		$(RM) swap_c.h
		./inc.sed SWAP_C < swap.c > swap_c.h

events_c.h:	evfilter.c inc.sed
		$(RM) events_c.h
		./inc.sed EVENTS_C < evfilter.c > events_c.h

image_c.h:	imfilter.c inc.sed
		$(RM) image_c.h
		./inc.sed IMAGE_C < imfilter.c > image_c.h

evregions.c:	evregions_c.tmpl Makefile
		$(RM) evregions.c
		sed 's/_FPUDBL_/$(FPU_DOUBLE)/' < evregions_c.tmpl > evregions.c

imregions.c:	imregions_c.tmpl Makefile
		$(RM) imregions.c
		sed 's/_FPUDBL_/$(FPU_DOUBLE)/' < imregions_c.tmpl > imregions.c

evregions_c.h:	evregions.c inc.sed
		$(RM) evregions_c.h
		./inc.sed EVREGIONS_C < evregions.c > evregions_c.h

imregions_c.h:	imregions.c inc.sed
		$(RM) imregions_c.h
		./inc.sed IMREGIONS_C < imregions.c > imregions_c.h

headers:	regions_h.h xalloc_c.h swap_c.h \
		events_c.h image_c.h evregions_c.h imregions_c.h

# remake the parser
LEX =		flex -i
parser:		filt.l
		@($(LEX) -Pfilt filt.l; \
		sed "s/yytext_ptr/filttext_ptr/g" < lex.filt.c > nlex.filt.c; \
		mv nlex.filt.c lex.filt.c)

iparser:	idx.l
		@($(LEX) -Pidx idx.l; \
		sed "s/yytext_ptr/idxtext_ptr/g" < lex.idx.c > nlex.idx.c; \
		mv nlex.idx.c lex.idx.c)

YACC =		bison -d -v
igrammar:	idx.y
		$(YACC) -p idx idx.y

pure:		tfilt.pure

tfilt.pure:	$(FILDLIB) tfilt.o
		purify $(CC) $(LDFLAGS) tfilt.o -o tfilt.pure \
		$(LIB) $(LIBS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
