# tifftclConfig.sh --
#
# This shell script (for sh) is generated automatically by tifftcl's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tifftcl extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tifftcl's version number.
tifftcl_VERSION='4.1.0'
tifftcl_MAJOR_VERSION=''
tifftcl_MINOR_VERSION=''
tifftcl_RELEASE_LEVEL=''

# The name of the tifftcl library (may be either a .a file or a shared library):
tifftcl_LIB_FILE=libtifftcl4.1.0.a

# String to pass to linker to pick up the tifftcl library from its
# build directory.
tifftcl_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/libtiff -ltifftcl4.1.0'

# String to pass to linker to pick up the tifftcl library from its
# installed directory.
tifftcl_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tifftcl4.1.0 -ltifftcl4.1.0'

# The name of the tifftcl stub library (a .a file):
tifftcl_STUB_LIB_FILE=libtifftclstub4.1.0.a

# String to pass to linker to pick up the tifftcl stub library from its
# build directory.
tifftcl_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/libtiff -ltifftclstub4.1.0'

# String to pass to linker to pick up the tifftcl stub library from its
# installed directory.
tifftcl_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tifftcl4.1.0 -ltifftclstub4.1.0'

# String to pass to linker to pick up the tifftcl stub library from its
# build directory.
tifftcl_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tkimg/libtiff/libtifftclstub4.1.0.a'

# String to pass to linker to pick up the tifftcl stub library from its
# installed directory.
tifftcl_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tifftcl4.1.0/libtifftclstub4.1.0.a'

# Location of the top-level source directories from which tifftcl
# was built.  This is the directory that contains generic, unix, etc.
# If tifftcl was compiled in a different place than the directory
# containing the source files, this points to the location of the
# sources, not the location where tifftcl was compiled. This can
# be relative to the build directory.

tifftcl_SRC_DIR='.'
