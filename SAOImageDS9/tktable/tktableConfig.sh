# tktableConfig.sh --
#
# This shell script (for sh) is generated automatically by tktable's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tktable extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tktable's version number.
tktable_VERSION='2.10'

# The name of the tktable library (may be either a .a file or a shared library):
tktable_LIB_FILE=libtktable2.10.a

# String to pass to linker to pick up the tktable library from its
# build directory.
tktable_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tktable -ltktable2.10'

# String to pass to linker to pick up the tktable library from its
# installed directory.
tktable_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tktable2.10 -ltktable2.10'

# The name of the tktable stub library (a .a file):
#tktable_STUB_LIB_FILE=libtktablestub2.10.a

# String to pass to linker to pick up the tktable stub library from its
# build directory.
#tktable_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tktable -ltktablestub2.10'

# String to pass to linker to pick up the tktable stub library from its
# installed directory.
#tktable_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tktable2.10 -ltktablestub2.10'

# String to pass to linker to pick up the tktable stub library from its
# build directory.
#tktable_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tktable/libtktablestub2.10.a'

# String to pass to linker to pick up the tktable stub library from its
# installed directory.
#tktable_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tktable2.10/libtktablestub2.10.a'
