prefix		= /usr/local/src/SAOImageDS9
exec_prefix	= /usr/local/src/SAOImageDS9

bindir		= ${exec_prefix}/bin
libdir		= $(exec_prefix)/lib
#libdir		= ${exec_prefix}/lib
includedir	= ${prefix}/include
datarootdir	= ${prefix}/share
datadir		= ${datarootdir}
mandir		= ${datarootdir}/man

OS		= unix
ARCH		= unknown
TKFLAGS		= 
SSLCONFIG	= ./config
JOBS		= 96
XML2CONFIG	= /usr/bin/xml2-config
XSLTCONFIG	= /usr/bin/xslt-config

include ./make.include

WITHTCL 	= --with-tcl=$(prefix)/tcl$(TCL_VERSION)/unix
WITHTK  	= --with-tk=$(prefix)/tk$(TCL_VERSION)/unix $(WITHTCL)
CACHE		= --config-cache --cache-file=$(prefix)/config.cache

#--------------------------build

.PHONY	: saods9 tcl tk

saods9	: dirs tcl tk \
	openssl \
	xpa funtools ast \
	util vector fitsy \
	tcliis tclsignal \
	tclzipfs tls tclxml \
	tkblt tclfitsy tktable tkimg \
	tkmpeg tksvg tkagif tkhtml1 tkcon \
	tksao ds9

tcl	: tcl$(TCL_VERSION)/unix/Makefile
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tcl$(TCL_VERSION)/unix -j $(JOBS) install

tcl$(TCL_VERSION)/unix/Makefile :
	@echo ""
	@echo "*** $@ ***"
	cd tcl$(TCL_VERSION)/unix; \
	./configure \
	$(PREFIX) \
	--disable-shared --enable-symbols \
	$(CACHE)

tk	: tk$(TCL_VERSION)/unix/Makefile
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tk$(TCL_VERSION)/unix -j $(JOBS) install

tk$(TCL_VERSION)/unix/Makefile :
	@echo ""
	@echo "*** $@ ***"
	cd tk$(TCL_VERSION)/unix; \
	./configure \
	$(TKFLAGS) \
	$(PREFIX) \
	--disable-shared --enable-symbols \
	$(CACHE)

#-------------------------- tclzipfs

tclzipfs: tclzipfs/Makefile
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclzipfs install

tclzipfs/Makefile :
	@echo ""
	@echo "*** $@ ***"
	cd tclzipfs; \
	./configure \
	$(WITHTCL) \
	$(PREFIX) \
	$(STD) \
	$(CACHE)

#-------------------------- tclsignal

tclsignal: tclsignal/Makefile
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclsignal install

tclsignal/Makefile :
	@echo ""
	@echo "*** $@ ***"
	cd tclsignal; \
	./configure \
	$(WITHTCL) \
	$(PREFIX) \
	$(STD) \
	$(CACHE)

#--------------------------clean

.PHONY	: tclclean tkclean osclean tclzipfsclean tclsignalclean

tclclean :
	@echo ""
	@echo "*** tclclean ***"
	$(MAKE) -C tcl$(TCL_VERSION)/unix clean

tkclean :
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tk$(TCL_VERSION)/unix clean

osclean : tclzipfsclean tclsignalclean

tclzipfsclean:
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclzipfs clean

tclsignalclean :
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclsignal clean

#--------------------------distclean

.PHONY	: tcldistclean tkdistclean
.PHONY	: osdistclean tclzipfsdistclean tclsignaldistclean

tcldistclean :
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tcl$(TCL_VERSION)/unix distclean

tkdistclean :
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tk$(TCL_VERSION)/unix distclean

osdistclean : tclzipfsdistclean tclsignaldistclean

tclzipfsdistclean:
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclzipfs distclean
	rm -f tclzipfs/tclzipfsConfig.sh

tclsignaldistclean :
	@echo ""
	@echo "*** $@ ***"
	$(MAKE) -C tclsignal distclean
	rm -f tclsignal/tclsignalConfig.sh

#--------------------------distribution

.PHONY	: dist

dist	:
	cd bin; tar cvf - ds9* | gzip > $(prefix)/dist/ds9.$(ARCH).$(DS9_VERSION).tar.gz
	cd bin; tar cvf - xpa* | gzip > $(prefix)/dist/xpa.$(ARCH).$(XPA_VERSION).tar.gz
