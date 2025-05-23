# Makefile for library files used by GNU text utilities.
# Copyright (C) 1991, 1992, 1993 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

SHELL = /bin/sh

srcdir = .


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

# There are just too many different versions of "install" around;
# better to use the install-sh script that comes with the distribution,
# which is slower but guaranteed to work.

INSTALL = 		./install-sh -c
INSTALL_PROGRAM =	${INSTALL}
INSTALL_DATA =		${INSTALL} -m 644

CC = gcc
#CC = gcc -fPIC -Wall -Wno-comment
#CC = gcc -fbounds-checking -DBOUNDS_CHECKING
AR = ar
RANLIB = ranlib
#DEFS = -DCONFIG_BROKETS -DHAVE_CONFIG_H
DEFS = -DHAVE_CONFIG_H
#CFLAGS = -g
# CFLAGS = -g -O2
#CFLAGS = -ansi -pedantic -Wimplicit -Wreturn-type -Wunused -Wformat -Wuninitialized -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Werror
CFLAGS = -g -O2
LDFLAGS = 
YACC = @YACC@

# extra libs
EXTRA_LIBS	= 

# extra objs (usually set to $ZOBJS)
EXTRA_OBJS	= 

# the full set of libraries for linking
LIBS	= $(EXTRA_LIBS)

OBJECTS=cardpar.o	\
	cardfmt.o	\
	cardins.o	\
	cardfind.o	\
	headget.o	\
	headset.o	\
	headins.o	\
	headinit.o	\
	headdata.o	\
	headfind.o	\
	headfile.o	\
	headimage.o	\
	headtable.o	\
	headsimple.o	\
	headwcs.o	\
	inverse.o	\
	dataimage.o	\
	fitshead.o	\
	fitsimage.o	\
	ftacht.o	\
	ftacht2.o	\
	nd.o		\
	tableheader.o	\
	$(EXTRA_OBJS)

# comment this out to use default defintions for stand-alone compile
XFILE =	-DUSE_XFILEIO=1  -I../util

INCL =	$(XFILE) -I$(srcdir)

# the default library for this package
DEFLIB =	libfitsy.a
# the actual library we are building (if this is a subpackage)
LIB =	/usr/local/src/SAOImageDS9/funtools/libfuntools.a

all:	$(LIB)

progs:	$(PROGS)

All:	$(LIB) $(PROGS)

.c.o:
	$(CC) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) $<

$(LIB):	$(OBJECTS)
	$(AR) cruv $@ $(OBJECTS)
	$(RANLIB) $@

objects:$(OBJECTS)

fitsy:	Force
	$(MAKE) clean
	$(MAKE) objects XFILE=""
	$(MAKE) libfitsy.a

libfitsy.a : Force
	$(AR) cruv $@ $(OBJECTS)
	$(RANLIB) $@


# these are the zlib routines:
adler32.o:	$(ZDIR)/adler32.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/adler32.c
compress.o:	$(ZDIR)/compress.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/compress.c
crc32.o:	$(ZDIR)/crc32.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/crc32.c
gzio.o:		$(ZDIR)/gzio.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/gzio.c
uncompr.o:	$(ZDIR)/uncompr.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/uncompr.c
deflate.o:	$(ZDIR)/deflate.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/deflate.c
trees.o:	$(ZDIR)/trees.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/trees.c
zutil.o:	$(ZDIR)/zutil.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/zutil.c
inflate.o:	$(ZDIR)/inflate.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/inflate.c
infblock.o:	$(ZDIR)/infblock.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/infblock.c
inftrees.o:	$(ZDIR)/inftrees.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/inftrees.c
infcodes.o:	$(ZDIR)/infcodes.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/infcodes.c
infutil.o:	$(ZDIR)/infutil.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/infutil.c
inffast.o:	$(ZDIR)/inffast.c
		$(CC) -I$(ZDIR) -c $(CPPFLAGS) $(DEFS) $(INCL) $(CFLAGS) \
		$(ZDIR)/inffast.c

install: install-binaries

install-binaries: $(LIB)
	@for i in $(LIB_INSTALL_DIR) $(INCLUDE_INSTALL_DIR) ; \
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
	@echo "Installing fitsy.h"
	@$(INSTALL_DATA) fitsy.h $(INCLUDE_INSTALL_DIR)/fitsy.h
	@echo "Installing xfile.h"
	@$(INSTALL_DATA) xfile.h $(INCLUDE_INSTALL_DIR)/xfile.h
#	@for i in $(PROGS) ; \
#	    do \
#		if [ -f $$i ] ; then \
#		    echo "Installing $$i" ; \
#		    $(INSTALL_PROGRAM) $$i $(BIN_INSTALL_DIR)/$$i ; \
#		fi; \
#	    done;

uninstall:

TAGS: $(SOURCES)
	etags $(SOURCES)

check:

RM =	rm -f

clean:
	$(RM) *.a *.so *.o *.exe core errs *pure* .nfs* \
	foo* *~ *.log \#* TAGS *.E a.out errors \
	$(ZDIR)/*.o \
	$(PROGS) $(TESTPROGS) \
	gmon.out *.pg *.bak \
	config.info config.log \
	doc/*~
	$(RM) -r autom4te.cache

mostlyclean: clean

distclean: clean
	$(RM) Makefile config.cache config.status config.log

realclean: distclean
	$(RM) TAGS

