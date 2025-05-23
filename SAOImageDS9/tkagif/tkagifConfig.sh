# tkagifConfig.sh --
#
# This shell script (for sh) is generated automatically by tkagif's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tkagif extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tkagif's version number.
tkagif_VERSION='1.0'

# The name of the tkagif library (may be either a .a file or a shared library):
tkagif_LIB_FILE=libtkagif1.0.a

# String to pass to linker to pick up the tkagif library from its
# build directory.
tkagif_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkagif -ltkagif1.0'

# String to pass to linker to pick up the tkagif library from its
# installed directory.
tkagif_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkagif1.0 -ltkagif1.0'

# The name of the tkagif stub library (a .a file):
#tkagif_STUB_LIB_FILE=libtkagifstub1.0.a

# String to pass to linker to pick up the tkagif stub library from its
# build directory.
#tkagif_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkagif -ltkagifstub1.0'

# String to pass to linker to pick up the tkagif stub library from its
# installed directory.
#tkagif_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkagif1.0 -ltkagifstub1.0'

# String to pass to linker to pick up the tkagif stub library from its
# build directory.
#tkagif_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tkagif/libtkagifstub1.0.a'

# String to pass to linker to pick up the tkagif stub library from its
# installed directory.
#tkagif_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tkagif1.0/libtkagifstub1.0.a'
