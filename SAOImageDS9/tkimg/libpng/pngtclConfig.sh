# pngtclConfig.sh --
#
# This shell script (for sh) is generated automatically by pngtcl's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for pngtcl extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# pngtcl's version number.
pngtcl_VERSION='1.6.37'
pngtcl_MAJOR_VERSION=''
pngtcl_MINOR_VERSION=''
pngtcl_RELEASE_LEVEL=''

# The name of the pngtcl library (may be either a .a file or a shared library):
pngtcl_LIB_FILE=libpngtcl1.6.37.a

# String to pass to linker to pick up the pngtcl library from its
# build directory.
pngtcl_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/libpng -lpngtcl1.6.37'

# String to pass to linker to pick up the pngtcl library from its
# installed directory.
pngtcl_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/pngtcl1.6.37 -lpngtcl1.6.37'

# The name of the pngtcl stub library (a .a file):
pngtcl_STUB_LIB_FILE=libpngtclstub1.6.37.a

# String to pass to linker to pick up the pngtcl stub library from its
# build directory.
pngtcl_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/libpng -lpngtclstub1.6.37'

# String to pass to linker to pick up the pngtcl stub library from its
# installed directory.
pngtcl_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/pngtcl1.6.37 -lpngtclstub1.6.37'

# String to pass to linker to pick up the pngtcl stub library from its
# build directory.
pngtcl_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tkimg/libpng/libpngtclstub1.6.37.a'

# String to pass to linker to pick up the pngtcl stub library from its
# installed directory.
pngtcl_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/pngtcl1.6.37/libpngtclstub1.6.37.a'

# Location of the top-level source directories from which [incr Tcl]
# was built.  This is the directory that contains generic, unix, etc.
# If [incr Tcl] was compiled in a different place than the directory
# containing the source files, this points to the location of the sources,
# not the location where [incr Tcl] was compiled.
pngtcl_SRC_DIR='.'
